"""
The :py:mod:`h2_fuel` module contains a class to read the required data and
a class to evaluate the power-to-fuel system.
"""

import os
import pandas as pd
import numpy as np
import pvlib


class ReadData:
    """

    This class enables to read data from the data files.

    Parameters
    ----------
    filename_climate : str
        The directory of the file with information on the
        climate data.

    """

    def __init__(self, filename_climate):
        self.filename_climate = filename_climate
        self.path = os.path.dirname(os.path.abspath(__file__))

    def load_climate(self):
        """

        This method loads the hourly solar irradiance data
        and ambient temperature data,
        situated in the 'sol_irr' and 'T_amb' columns of the
        climate data file.

        Returns
        -------
        sol_irr : ndarray
            The hourly solar irradiance data for a Typical
            Meteorological Year. (8760 elements)
        t_amb : ndarray
            The hourly ambient temperature data for a Typical
            Meteorological Year. (8760 elements)

        """
        data = pd.read_csv(self.filename_climate)
        sol_irr = data['sol_irr'].to_numpy()
        t_amb = data['T_amb'].to_numpy()

        return sol_irr, t_amb

    def load_parameters(self):
        """

        This method loads the deterministic values of the model
        parameters, defined in the design_space file. This is
        useful when the deterministic performance of a specific
        design needs to be evaluated.

        Returns
        -------
        param_dict : dict
            Dictionary with the names of the model parameters
            and the corresponding deterministic values.

        """
        param_dict = {}
        design_space = os.path.join(self.path, 'design_space')

        # read the deterministic values for the parameters in `design_space`
        with open(design_space, 'r') as file:
            for line in file:
                tmp = line.split()
                if tmp[1] == 'par':
                    param_dict[tmp[0]] = float(tmp[2])

        return param_dict


class Evaluation:
    """

    This class evaluates the photovoltaic-electrolyzer system.
    For a given design, the solar irradiance, ambient temperature
    and the characterization of the model parameters,
    the levelized cost of hydrogen and the annual hydrogen production
    are quantified.

    Parameters
    ----------
    sol_irr : ndarray
        The hourly solar irradiance for the evaluated year.
    t_amb : ndarray
        The hourly ambient temperature for the evaluated year.
    parameters : dict
        Dictionary with the model parameters and design variables values.

    """

    def __init__(self, sol_irr, t_amb, par):
        self.par = par

        # the solar irradiance and ambient temperature are scaled with the
        # corresponding uncertainty
        self.sol_irr = sol_irr * self.par['u_sol_irr']
        self.t_amb = t_amb + self.par['u_t_amb']

        # the result dictionary
        self.res = {}

        # the system lifetime
        self.par['life_sys'] = 20.

        # the amount of hydrogen produced
        self.res['m_h2'] = 0.

        # initialize the pv power that is consumed by the electrolyzer array
        self.res['p_pv_consumed'] = 0.

        # initialize the operating hours of the electrolyzer array
        self.res['running_hours_pemel'] = 0.

        # the number of PEM electrolyzer cells, corresponding to the
        # nominal capacity of the considered PEM cell and the provided
        # PEM capacity
        self.n_pemel_array = self.par['n_pemel'] / 0.4

        # generate the fitted polynomial on the electrolyzer array
        # power - current relation
        self.polyfit_pemel()

    #############################
    # photovoltaic array module #
    #############################

    def quantify_mpp(self, sol_irr, t_amb, pv_system):
        """

        Quantify the maximum power of the photovoltaic array
        for a given solar irradiance and ambient temperature.

        Parameters
        ----------
        sol_irr : float
            The solar irradiance [W/m2].
        t_amb : float
            The ambient temperature [C].
        pv_system : pandas.core.series.Series
            The pv system characteristics
        Returns
        -------
        pmp : float
            The maximum power.

        """

        # quantify the parameters for the pv system using De Soto method
        pv_inputs = pvlib.pvsystem.calcparams_desoto(sol_irr,
                                                     t_amb,
                                                     pv_system['alpha_sc'],
                                                     pv_system['a_ref'],
                                                     pv_system['I_L_ref'],
                                                     pv_system['I_o_ref'],
                                                     pv_system['R_sh_ref'],
                                                     pv_system['R_s'],
                                                     EgRef=1.121,
                                                     dEgdT=-0.0002677,
                                                     irrad_ref=1000.,
                                                     temp_ref=25.)

        # determine the maximum power for the given pv system
        pmp = pvlib.pvsystem.max_power_point(pv_inputs[0],
                                             pv_inputs[1],
                                             pv_inputs[2],
                                             pv_inputs[3],
                                             pv_inputs[4],
                                             method='newton')['p_mp']

        return pmp

    def photovoltaic(self):
        """

        The hourly photovoltaic power is quantified via the PVlib package.
        Using this package, first the characteristics for a typical
        photovoltaic panel are defined. Based on these characteristics,
        the maximum power point is quantified for each hour, based on the
        corresponding solar irradiance and ambient temperature. Finally, the
        hourly power production is scaled by the considered photovoltaic array
        capacity.

        """

        p_pv = np.zeros(len(self.sol_irr))

        # get the specific photovoltaic panel characteristics
        pv_database = pvlib.pvsystem.retrieve_sam('CECmod')
        pv_system = pv_database.SunPower_SPR_X19_240_BLK

        p_mpp_ref = self.quantify_mpp(1000., 25., pv_system)  # W

        # maximum power point determination for each hour in the timeframe
        for i, irr in enumerate(self.sol_irr):
            if irr > 0.:
                p_mpp = self.quantify_mpp(irr, self.t_amb[i], pv_system)
                p_mpp_array = p_mpp / p_mpp_ref * self.par['n_pv']
                p_pv[i] = min(p_mpp_array, self.par['n_dcdc_pv']) * 1e3  # W
            else:
                p_pv[i] = 0.

        # store the hourly pv power in the result dictionary
        self.res['p_pv'] = p_pv

    #############################
    # electrolyzer array module #
    #############################

    def pemel(self, i_pemel):
        """
        The electrolyzer model, based on the work of Saeed et al. [1]. For a
        given current, the model determines the operating voltage by
        considering the activation, concentration and ohmic overpotentials.
        The model quantifies the operating voltage, power, efficiency and
        hydrogen production.

        [1] Saeed, E. W., & Warkozek, E. G. (2015). Modeling and Analysis of
            Renewable PEM Fuel Cell System. Energy Procedia, 74, 87–101.
            https://doi.org/10.1016/j.egypro.2015.07.527

        Parameters
        ----------
        i_pemel : float
            The electrolyzer input current [A].

        Returns
        -------
        res : dict
            Dictionary with the operating conditions of the electrolyzer for a
            given current. It contains items on the operating voltage, power,
            efficiency and hydrogen mass flow rate.

        """
        par_pemel = {'T': 353.,
                     'a': 1.,
                     'p_o2': 1.,
                     'p_h2': 1.,
                     'p_h2o': 1.,
                     'i_L': 2.,
                     'A': 100.,
                     'i_0': 1e-4,
                     'n': 2.,
                     't_mem': 50e-4,
                     'alpha': 0.3,
                     'R': 8.3143,
                     'F': 96485.,
                     'HHV': 141.7e6,
                     }

        res = {}
        i = i_pemel / par_pemel['A']

        # minimum operating voltage of electrolyzer
        e_0 = (1.48 - 0.85e-3 * (par_pemel['T'] - 298.15) + 4.3085e-5 *
               par_pemel['T'] * np.log(par_pemel['p_h2'] *
                                       np.sqrt(par_pemel['p_o2']) /
                                       par_pemel['p_h2o']))

        # activation overpotential
        v_act = (np.log(i / par_pemel['i_0']) /
                 (par_pemel['alpha'] * par_pemel['n'] * par_pemel['F']) *
                 par_pemel['R'] * par_pemel['T'])

        # ohmic overpotential
        lambda_mem = (0.043 + 17.81 * par_pemel['a'] -
                      39.85 * par_pemel['a']**2. +
                      36. * par_pemel['a']**3.)
        sigma_mem = ((0.005139 * lambda_mem - 0.00326) *
                     np.exp(1268 * (1. / 303. - 1. / par_pemel['T'])))
        v_ohm = i * par_pemel['t_mem'] / sigma_mem

        # the concentration overpotential
        v_con = - (par_pemel['R'] * par_pemel['T'] /
                   (par_pemel['n'] * par_pemel['F']) *
                   np.log(1. - i / par_pemel['i_L']))

        # model outputs
        res['v_pemel'] = (e_0 + v_act + v_ohm + v_con) * self.n_pemel_array
        res['m_pemel'] = self.current_to_mh2(i_pemel) * self.n_pemel_array
        res['p_pemel'] = i_pemel * res['v_pemel']
        res['eff_pemel'] = (res['m_pemel'] * par_pemel['HHV'] /
                          (res['p_pemel'] * 3600.))
        return res

    def current_to_mh2(self, current):
        """
        When current is provided, this function determines the
        corresponding hydrogen mass flow rate per hour.

        Parameters
        ----------
        current : float
            The electrolyzer input current [A].

        Returns
        -------
        m_h2 : float
            The produced hydrogen mass flow rate [kg/h].

        """
        far_cons = 96485.
        m_h2 = current / (2. * far_cons) * 2.02e-3 * 3600.

        return m_h2

    def polyfit_pemel(self):
        """
        The electrolyzer stack is evaluated over a range of input currents.
        Following these evaluations, a polynomial is fitted on the
        power - current relation of the electrolyzer. This polynomial enables
        to rapidly determine the input current when a certain amount of power
        is available. Since this relation is fairly linear, the polynomial
        should reach good agreement with the actual power - current relation,
        while maintaining the level of fidelity of the actual model.

        """

        # evaluate the electrolyzer stack for a set of currents
        i_list = np.arange(start=3, stop=200, step=4)
        p_pemel = np.zeros(len(i_list))
        for index, i in enumerate(i_list):
            res = self.pemel(i)
            p_pemel[index] = res['p_pemel']

        # generate a polynomial fitted on the power - current points
        self.p_to_i_pemel = polyfit_func(p_pemel, i_list)

    def charge_pemel(self, p_pemel):
        """
        For a given power supplied to the electrolyzer, this function
        determines the actual hydrogen produced. First, the method evaluates
        if the power supplied lies within the operating bounds of the
        electrolyzer stack. If the power is situated below the lower limit,
        the electrolyzer does not run. Instead, when the power is situated
        above the upper limit, the electrolyzer operates at nominal conditions.
        At nominal conditions, the current is known and the hydrogen mass flow
        rate is quantified. Otherwise, the input current is determined through
        the fitted polynomial on the power - current relation of the
        electrolyzer stack. As the power is an output from the electrolyzer
        model, considering this polynomial avoids the use of root finding
        methods and is therefore more robust in optimization and uncertainty
        quantification approaches. Finally, when hydrogen is produced, the
        running hours of the electrolyzer stack is increased by 1 and the
        hydrogen mass flow rate is returned.

        Parameters
        ----------
        p_pemel : float
            The power available for electrolysis [W].

        Returns
        -------
        m_h2 : float
            The produced hydrogen mass flow rate [kg/h].

        """
        # the operating bounds
        op_lower_lim = self.par['n_pemel'] * 10.
        op_upper_lim = self.par['n_pemel'] * 1e3

        # check if power is higher than the lowest operating point
        if p_pemel > op_lower_lim:

            # if the power exceeds the upper bound, operate at upper bound
            p_pemel_applied = min(p_pemel, op_upper_lim)

            # current at this applied power
            i_pemel = self.p_to_i_pemel(p_pemel_applied)
            m_h2 = self.current_to_mh2(i_pemel) * self.n_pemel_array
            self.res['p_pv_consumed'] += p_pemel_applied

            # increase the operating hours by 1
            self.res['running_hours_pemel'] += 1.

        # no hydrogen production when the power falls outside the operating
        # bounds
        else:
            m_h2 = 0.

        return m_h2

    #####################
    # evaluation module #
    #####################

    def evaluation(self):
        """

        This is the main method of the Evaluation class.
        In this method, the hourly photovoltaic power is
        quantified first. Then, for each hour, the hydrogen
        is determined. Finally, the electrolyzer lifetime and
        the system cost are determined.

        """

        # get the hourly photovoltaic array power
        self.photovoltaic()

        # evaluate hourly hydrogen production
        for p_pv in self.res['p_pv']:
            self.res['m_h2'] += self.charge_pemel(p_pv)

        # determine the electrolyzer lifetime
        self.lifetime()

        # determine the system cost and levelized cost of hydrogen
        self.cost()

    def lifetime(self):
        """

        The lifetime method determines the lifetime of
        the electrolyzer array, based on the number of
        operating hours during the evaluated year.

        """

        # lifetime of the electrolyzer array
        if self.res['running_hours_pemel'] == 0.:
            self.res['life_pemel'] = 1e8
        else:
            self.res['life_pemel'] = (self.par['life_pemel'] /
                                      self.res['running_hours_pemel'])

    def cost(self):
        """

        Based on the capital recovery factor, the CAPEX,
        OPEX and replacement cost of the system components,
        the levelized cost of hydrogen is determined. The formula
        for the annualized system cost is adopted from Zakeri et al. [2].

        [2] Zakeri, B., & Syri, S. (2015). Electrical energy storage systems:
            A comparative life cycle cost analysis. Renewable and Sustainable
            Energy Reviews, 42, 569–596.
            https://doi.org/10.1016/j.rser.2014.10.011
        """

        # the capital recovery factor
        inv_rate = ((self.par['int_rate'] - self.par['infl_rate']) /
                    (1. + self.par['infl_rate']))
        crf = (((1. + inv_rate)**self.par['life_sys'] - 1.) /
               (inv_rate * (1. + inv_rate)**self.par['life_sys']))**(-1)

        # annual cost of photovoltaic array and DC-DC converter
        pv_cost = self.par['n_pv'] * (crf * self.par['capex_pv'] +
                                      self.par['opex_pv'])
        pv_dcdc_cost = self.par['n_dcdc_pv'] * (self.par['capex_dcdc'] *
                                                (crf + self.par['opex_dcdc']))
        components_cost = pv_cost + pv_dcdc_cost

        # annual cost of electrolyzer array
        pemel_cost = self.par['n_pemel'] * (self.par['capex_pemel'] *
                                            (crf + self.par['opex_pemel']))
        components_cost += pemel_cost

        # annual replacement cost of the electrolyzer array
        arc = crf * sum([(1. + inv_rate)**(-(i + 1.) *
                                           self.res['life_pemel']) *
                         self.par['n_pemel'] *
                         self.par['repl_pemel'] *
                         self.par['capex_pemel'] for
                         i in range(int(self.par['life_sys'] /
                                        self.res['life_pemel']))])

        # the levelized cost of hydrogen
        cost = arc + components_cost
        if self.res['m_h2'] < 1e-5:
            self.res['lcoh'] = 1e8
        else:
            self.res['lcoh'] = cost / self.res['m_h2']

    def print_results(self):
        """

        This method prints the levelized cost of hydrogen,
        the hydrogen production, the annual energy produced
        by the photovoltaic array and the energy consumed by
        the electrolyzer array.

        """

        print('outputs:')
        print('LCOH:'.ljust(30) + '%.2f euro/kg' % self.res['lcoh'])
        print('m_h2:'.ljust(30) + '%.2f kg' % self.res['m_h2'])
        print('PV electricity generated:'.ljust(30) +
              '%.2f MWh' % (sum(self.res['p_pv']) / 1e6))
        print(
            'PV electricity consumed:'.ljust(30) + '%.2f MWh' %
            (self.res['p_pv_consumed'] / 1e6))
        print('self-consumption ratio:'.ljust(30) + '%.2f %%' %
              (1e2 * self.res['p_pv_consumed'] /
               sum(self.res['p_pv'])))
        print('lifetime electrolyzer:'.ljust(30) + '%.2f year' %
              self.res['life_pemel'])


def polyfit_func(x_in, y_in, threshold=0.99999999):
    """
    The function fits a polynomial to the points of x_in and y_in. The
    polynomial starts with order 1. To evaluate its performance, the R-squared
    performance indicator is quantified. If the value for R-squared does
    not reach the defined threshold, the polynomial order is increased and
    the polynomial is fitted again on the points, until the threshold is
    satisfied. Once satisfied, the function returns the polynomial.

    Parameters
    ----------
    x_in : ndarray
        The x-coordinates for the sample points.
    y_in : ndarray
        The y-coordinates for the sample points.
    threshold : float, optional
        The threshold for the R-squared parameter. The default is 0.99999.

    Returns
    -------
    poly_func : numpy.poly1d
        A one-dimensional polynomial.

    """
    order = 0
    r_squared = 0.
    while r_squared < threshold:
        order += 1

        # the polynomial
        poly_coeff = np.polyfit(x_in, y_in, order)
        poly_func = np.poly1d(poly_coeff)

        # r-squared
        yhat = poly_func(x_in)
        ybar = np.sum(y_in) / len(y_in)
        ssreg = np.sum((yhat - ybar)**2.)
        sstot = np.sum((y_in - ybar)**2.)
        r_squared = ssreg / sstot

    return poly_func
