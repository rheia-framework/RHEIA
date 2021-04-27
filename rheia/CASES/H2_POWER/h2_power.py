"""
The :py:mod:`h2_power` module contains a class to read the required data and
a class to evaluate the power-to-power system.
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pvlib


class ReadData:
    """

    This class enables to read data from the data files.

    Parameters
    ----------
    filename_climate : str
        The directory of the file with information on the
        solar irradiance.
    filename_climate : str
        The directory of the file with information on the
        demand.

    """

    def __init__(self, filename_climate, filename_demand):
        self.filename_climate = filename_climate
        self.filename_demand = filename_demand
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

    def load_demand(self):
        """

        This method loads the hourly electricity demand data,
        situated in the 'total electricity' column of the demand data file.

        Returns
        -------
        load_elec : ndarray
            The hourly electricity demand.

        """
        data = pd.read_csv(self.filename_demand)
        load_elec = data['total electricity'].to_numpy() * 1e3

        return load_elec

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

    This class evaluates the photovoltaic-hydrogen system.
    For a given design, the solar irradiance, electricity demand
    and the characterization of the model parameters,
    the levelized cost of electricity and the self-sufficiency ratio
    are quantified.

    Parameters
    ----------
    sol_irr : ndarray
        The hourly solar irradiance for the evaluated year.
    t_amb : ndarray
        The hourly ambient temperature for the evaluated year.
    L_elec : ndarray
        The hourly electricity demand for the evaluated year.
    parameters : dict
        Dictionary with the model parameters and design variables values.

    """

    def __init__(self, sol_irr, t_amb, L_elec, par):
        self.par = par

        # the solar irradiance and ambient temperature are scaled with the
        # corresponding uncertainty
        self.sol_irr = sol_irr * self.par['u_sol_irr']
        self.t_amb = t_amb + self.par['u_t_amb']

        # the electric load, scaled with the uncertainty on the load
        self.load_elec = L_elec * self.par['u_load_elec']
        # the result dictionary
        self.res = {}

        # the system lifetime
        self.par['life_sys'] = 20.

        # initialize the operating hours of the electrolyzer and fuel cell
        self.res['running_hours_pemel'] = 0.
        self.res['running_hours_pemfc'] = 0.

        # initialize the storage tank size and its starting status
        self.m_h2_max = self.tank()
        self.m_h2_min = 0.05 * self.m_h2_max
        self.m_h2 = self.m_h2_min

        # the number of PEM electrolyzer cells and PEM fuel cells,
        # corresponding to the nominal capacity of the considered
        # PEM electrolyzer and PEM fuel cell and the provided
        # system capacities
        self.n_pemel_array = self.par['n_pemel'] / 0.4
        self.n_pemfc_array = self.par['n_pemfc'] / (5. / 35.)

        # generate the fitted polynomial on the electrolyzer array and
        # fuel cell array power - current relation
        self.polyfit_pemel()
        self.polyfit_pemfc()

        # instantiate the profiles for grid electricity price
        self.elec_profiles()

    def elec_profiles(self):
        """
        Set the grid electricity price for buying and selling electricity.
        A contract with fixed electricity price is considered, for which the
        price for buying electricity consists of three segments: the energy
        price itself (i.e. 'elec cost'), the profit made on this price by the
        electricity provider (i.e. 'elec_cost_profit') and the fraction of the
        energy price to the final retail price (i.e. 'elec_cost_ratio', e.g.
        when this value equal 0.3, then the energy price corresponds to 30% of
        the final bill, while 70% corresponds to transportation cost,
        taxes,...). The price for selling electricity back to the grid
        corresponds to the energy price.

        """

        self.elec_profile = np.ones(len(self.sol_irr)) * (
            (self.par['elec_cost'] + self.par['elec_cost_profit']) /
            self.par['elec_cost_ratio']) / 1e6

        self.elec_profile_sale = np.ones(len(self.sol_irr)) * (
            self.par['elec_cost'] / 1e6)

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
                p_pv[i] = p_mpp / p_mpp_ref * self.par['n_pv'] * 1e3  # W
            else:
                p_pv[i] = 0.

        # store the hourly pv power in the result dictionary
        self.res['p_pv'] = p_pv

        # the dc-dc converter capacity in kW
        self.res['n_dcdc_pv'] = max(p_pv) / 1e3

    def net_power(self):
        """

        Determine the hourly net power.
        This corresponds to the net power available (or still required) on the
        DC bus bar after extracting the load from the photovoltaic power.

        Returns
        -------
        p_net : ndarray
            The hourly net power. [W]

        """

        p_net = np.zeros(len(self.sol_irr))

        # determine the hourly photovoltaic power
        self.photovoltaic()

        for i, p_pv in enumerate(self.res['p_pv']):

            # the excess pv power after extracting the electric load
            p_net[i] = p_pv - self.load_elec[i]

        return p_net

    ###############
    # tank module #
    ###############

    def tank(self):
        """

        The maximum storage capacity of the hydrogen tank.

        Returns
        -------
        m_max : float
            The hydrogen storage capacity. [kg]

        """

        # conversion from energy (kWh) into mass (kg)
        m_max = self.par['n_tank'] / 33.33

        return m_max

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
        The input current is determined through the fitted polynomial on the
        power - current relation of the electrolyzer stack. As the power is an
        output from the electrolyzer model, considering this polynomial avoids
        the use of root finding methods and is therefore more robust in
        optimization and uncertainty quantification approaches.

        When the hydrogen mass is determined that can be generated, the current
        hydrogen tank capacity is addressed. If the capacity exceeds the
        maximum storage tank capacity, the current applied to the electrolyzer
        is re-evaluated, such that this current matches the hydrogen production
        that leads to a full hydrogen tank.

        Finally, when hydrogen is produced, the
        running hours of the electrolyzer stack is increased by 1 and the
        power consumed by the electrolyzer array is returned.

        Parameters
        ----------
        p_pemel : float
            The power available for electrolysis [W].

        Returns
        -------
        p_consumed : float
            The power consumed by the electrolyzer array [W].

        """

        # the operating bounds
        op_lower_lim = self.par['n_pemel'] * 10.
        op_upper_lim = self.par['n_pemel'] * 1e3

        # check if power is higher than the lowest operating point
        if (self.m_h2 < self.m_h2_max and
                p_pemel > op_lower_lim and
                self.par['n_pemel'] > 1e-1):

            # if the power exceeds the upper bound, operate at upper bound
            p_pemel_applied = min(p_pemel, op_upper_lim)

            # current at this applied power
            i_pemel = self.p_to_i_pemel(p_pemel_applied)
            m_h2 = self.current_to_mh2(i_pemel) * self.n_pemel_array
            p_consumed = p_pemel_applied

            # produced hydrogen is added to the tank
            self.m_h2 += m_h2

            # if the new hydrogen capacity exceeds the maximum storage capacity
            if self.m_h2 > self.m_h2_max:

                # define current that results in a full storage tank
                excess = self.m_h2 - self.m_h2_max
                addition = m_h2
                allowed_current = (1. - excess / addition) * i_pemel

                # check if this current is still higher than lower limit
                if allowed_current > self.p_to_i_pemel(op_lower_lim):
                    p_consumed = self.pemel(allowed_current)['p_pemel']
                    self.m_h2 = self.m_h2_max
                    self.res['running_hours_pemel'] += 1.
                else:
                    self.m_h2 -= m_h2
                    p_consumed = 0.
            else:
                self.res['running_hours_pemel'] += 1.
        else:
            p_consumed = 0.

        return p_consumed

    ##########################
    # fuel cell array module #
    ##########################

    def polyfit_pemfc(self):
        """
        The fuel cell stack is evaluated over a range of input currents.
        Following these evaluations, a polynomial is fitted on the
        power - current relation of the fuel cell. This polynomial enables
        to rapidly determine the input current when a certain amount of power
        needs to be provided. Since this relation is fairly linear, the
        polynomial should reach good agreement with the actual power - current
        relation, while maintaining the level of fidelity of the actual model.

        """

        # evaluate the electrolyzer stack for a set of currents
        i_list = np.arange(start=1, stop=325, step=4)
        p_pemfc = np.zeros(len(i_list))
        for index, i in enumerate(i_list):
            res = self.pemfc(i)
            p_pemfc[index] = res['p_pemfc']

        # generate a polynomial fitted on the power - current points
        self.p_to_i_pemfc = polyfit_func(p_pemfc, i_list)

    def pemfc(self, i_pemfc):
        """
        The PEM fuel cell model, based on the work of Murugesan et al. [2].
        The model determines the current-voltage characteristic and provides
        the voltage, power, efficiency and hydrogen consumption.

        [2] Murugesan, K., & Senniappan, V. (2013). Investigation of water
        management dynamics on the performance of a Ballard-Mark-V proton
        exchange membrane fuel cell stack system. International Journal of
        Electrochemical Science, 8(6), 7885–7904.

        Parameters
        ----------
        i_pemfc : float
            The fuel cell output current.

        Returns
        -------
        res : dict
            Dictionary with the operating conditions of the fuel cell for a
            given current. It contains items on the operating voltage, power,
            efficiency and hydrogen mass flow rate.

        """

        par_pemfc = {'p_o2': 3.,
                     'p_h2': 3.,
                     'i_L': 1.5,
                     'T': 353.,
                     'A': 232.,
                     't_mem': 178e-4,
                     'F': 96485.,
                     'B': 0.016,
                     'HHV': 141.7e6,
                     'lambda_mem': 23.,
                     'eps_1': -0.948,
                     'eps_2': 0.00354,
                     'eps_3': 7.6e-5,
                     'eps_4': -1.93e-4,
                     'b2': 1268.,
                     'b11': 0.005139,
                     'b12': 0.00326,
                     }

        res = {}
        i = i_pemfc / par_pemfc['A']

        # minimum operating voltage of fuel cell
        e_0 = (1.229 - 0.85e-3 * (par_pemfc['T'] - 298.15) + 4.31e-5 *
               par_pemfc['T'] * np.log(par_pemfc['p_h2'] *
                                       np.sqrt(par_pemfc['p_o2'])))

        # activation overpotential
        c_o2 = par_pemfc['p_o2'] / (5.08e6 * np.exp(498. / par_pemfc['T']))
        v_act = - (par_pemfc['eps_1'] + par_pemfc['eps_2'] * par_pemfc['T'] +
                   par_pemfc['eps_3'] * par_pemfc['T'] * np.log(c_o2) +
                   par_pemfc['eps_4'] * np.log(i_pemfc) * par_pemfc['T'])

        # ohmic overpotential
        b_1 = par_pemfc['b11'] * par_pemfc['lambda_mem'] - par_pemfc['b12']
        sigma_mem = b_1 * np.exp(par_pemfc['b2'] *
                                 (1 / 303. - 1. / par_pemfc['T']))
        r_m = par_pemfc['t_mem'] / sigma_mem
        v_ohm = i_pemfc * r_m / par_pemfc['A']

        # concentration overpotential
        v_con = - par_pemfc['B'] * np.log(1. - i / par_pemfc['i_L'])

        # model outputs
        res['v_pemfc'] = (e_0 - v_act - v_ohm - v_con) * self.n_pemfc_array
        res['m_pemfc'] = self.current_to_mh2(i_pemfc) * self.n_pemfc_array
        res['p_pemfc'] = i_pemfc * res['v_pemfc']
        res['eff_pemfc'] = (res['m_pemfc'] * par_pemfc['HHV'] /
                            (res['p_pemfc'] * 3600.))
        return res

    def charge_pemfc(self, p_pemfc):
        """

        This method evaluates if the power required from the fuel cell
        lies within the operating range of the DC-DC converter and of the
        fuel cell array. If yes, the power is supplied and the
        consumed hydrogen is quantified. If the required power is larger than
        the fuel cell capacity, than the nominal power is supplied.
        If the consumed hydrogen is larger than the available hydrogen in
        the hydrogen storage tank, the supplied power is recalculated,
        such that the consumed hydrogen matches the initial available hydrogen
        left in the hydrogen storage tank.
        Finally, the operating hours of the fuel cell array is
        increased by 1.

        Parameters
        ----------
        p_pemfc : float
            The power demanded from the fuel cell array.

        Returns
        -------
        p_produced : float
            The actual power produced by the fuel cell array.
        """

        # the operating bounds
        op_lower_lim = self.par['n_pemfc'] * 10.
        op_upper_lim = self.par['n_pemfc'] * 1e3

        # check if power is higher than the lowest operating point
        if (self.m_h2 > self.m_h2_min and
                p_pemfc > op_lower_lim and
                self.par['n_pemfc'] > 1e-1):

            # if the power exceeds the upper bound, operate at upper bound
            p_pemfc_applied = min(p_pemfc, op_upper_lim)

            # current at this required power
            i_pemfc = self.p_to_i_pemfc(p_pemfc_applied)
            m_h2 = self.current_to_mh2(i_pemfc) * self.n_pemfc_array
            p_produced = p_pemfc_applied

            # produced hydrogen is extracted from the tank
            self.m_h2 -= m_h2

            # if the new hydrogen capacity is below the minimum storage
            # capacity
            if self.m_h2 < self.m_h2_min:

                # define current that results in an empty storage tank
                lack = self.m_h2_min - self.m_h2
                extracted = m_h2
                allowed_current = (1. - lack / extracted) * i_pemfc

                # check if this current is still higher than lower limit
                if allowed_current > self.p_to_i_pemfc(op_lower_lim):
                    p_produced = self.pemfc(allowed_current)['p_pemfc']
                    self.m_h2 = self.m_h2_min
                    self.res['running_hours_pemfc'] += 1.
                else:
                    self.m_h2 += m_h2
                    p_produced = 0.
            else:
                self.res['running_hours_pemfc'] += 1.
        else:
            p_produced = 0.

        return p_produced

    #####################
    # evaluation module #
    #####################

    def evaluation(self):
        """

        This is the main method of the Evaluation class.
        For each hour, the power management strategy is applied.
        If the net power is positive, the electrolyzer is charged. If, after
        the charging, there is still power available, this power is sold to
        the grid. Instead, when the net power is negative, the fuel cell is
        charged. If the power generated by the fuel cell is insufficient, the
        remaining power is extracted from the grid.
        Finally, the self-sufficiency ratio and the system cost are determined.

        """

        n_dcdc_pemel = np.zeros(len(self.sol_irr))
        n_dcdc_pemfc = np.zeros(len(self.sol_irr))
        n_dcac = np.zeros(len(self.sol_irr))

        self.res['m_h2_array'] = np.ones(len(self.sol_irr))
        self.res['grid_e_buy'] = np.ones(len(self.sol_irr))
        self.res['grid_e_sold'] = np.ones(len(self.sol_irr))

        # the net power, i.e. after using the pv power to cover the load
        p_net = self.net_power()

        for t, p_net in enumerate(p_net):
            e_grid_buy = 0.
            e_grid_sold = 0.

            # if excess power available, charge electrolyzer
            if p_net > 0.:
                p_consumed = self.charge_pemel(p_net)
                n_dcdc_pemel[t] = p_consumed
                p_rem = p_net - p_consumed

                # sell remaining electricity to the grid
                e_grid_sold += p_rem

                # power that goes through the DC-AC converter,
                # i.e. the load is covered by the pv array +
                # additional energy is sold to the grid
                n_dcac[t] = self.load_elec[t] + e_grid_sold

            # if load remaining, charge fuel cell
            else:
                p_produced = self.charge_pemfc(abs(p_net))
                n_dcdc_pemfc[t] = p_produced
                p_req = abs(p_net) - p_produced

                # buy remaining electricity from the grid
                e_grid_buy += p_req

                # power that goes through the DC-AC converter,
                # i.e. the load that was covered by the pv array
                n_dcac[t] = (self.load_elec[t] - e_grid_buy)

            # store evolution of storage tank and electricity bought/sold
            self.res['m_h2_array'][t] = ((self.m_h2 - self.m_h2_min) /
                                         (self.m_h2_max - self.m_h2_min))
            self.res['grid_e_buy'][t] = e_grid_buy
            self.res['grid_e_sold'][t] = e_grid_sold

        # define the capacity of the converters and inverter
        self.res['n_dcac'] = max(n_dcac) / 1e3
        self.res['n_dcdc_pemel'] = max(n_dcdc_pemel) / 1e3
        self.res['n_dcdc_pemfc'] = max(n_dcdc_pemfc) / 1e3

        # determine the electrolyzer and fuel cell lifetime
        self.lifetime()

        # determine the self-sufficiency ratio
        self.self_sufficiency_ratio()

        # determine the system cost
        self.cost()

    def lifetime(self):
        """

        The lifetime method determines the lifetime of
        the electrolyzer array and fuel cell array, based on the number of
        operating hours for each component during the evaluated period.

        """
        # lifetime of the electrolyzer array
        if self.res['running_hours_pemel'] == 0.:
            self.res['life_pemel'] = 1e8
        else:
            self.res['life_pemel'] = (self.par['life_pemel'] /
                                      self.res['running_hours_pemel'])

        # lifetime of the fuel cell array
        if self.res['running_hours_pemfc'] == 0.:
            self.res['life_pemfc'] = 1e8
        else:
            self.res['life_pemfc'] = (self.par['life_pemfc'] /
                                      self.res['running_hours_pemfc'])

    def self_sufficiency_ratio(self):
        """

        The self-sufficiency ratio is quantified. The self-sufficiency ratio
        corresponds to the fraction of the electric load that is covered by
        the photovoltaic-hydrogen system.

        """
        self.res['ssr'] = (1. - sum(self.res['grid_e_buy']) /
                           sum(self.load_elec))

    def cost(self):
        """

        Based on the capital recovery factor, the CAPEX,
        OPEX and replacement cost of the system components,
        the levelized cost of electricity is determined. The formula
        for the annualized system cost is adopted from Coppitters et al. [3].

        [3] Coppitters, D., De Paepe, W., & Contino, F. (2020). Robust design
            optimization and stochastic performance analysis of a
            grid-connected photovoltaic system with battery storage and
            hydrogen storage. Energy, 213, 118798.
            https://doi.org/10.1016/j.energy.2020.118798
        """

        # the capital recovery factor
        inv_rate = ((self.par['int_rate'] - self.par['infl_rate']) /
                    (1. + self.par['infl_rate']))
        crf = (((1. + inv_rate)**self.par['life_sys'] - 1.) /
               (inv_rate * (1. + inv_rate)**self.par['life_sys']))**(-1)

        # annual cost of photovoltaic array and DC-DC converter
        pv_cost = self.par['n_pv'] * (crf * self.par['capex_pv'] +
                                      self.par['opex_pv'])
        pv_dcdc_cost = self.res['n_dcdc_pv'] * (self.par['capex_dcdc'] *
                                                (crf + self.par['opex_dcdc']))
        components_cost = pv_cost + pv_dcdc_cost

        # annual cost of electrolyzer array and DC-DC converter
        pemel_cost = self.par['n_pemel'] * (self.par['capex_pemel'] *
                                            (crf + self.par['opex_pemel']))
        pemel_dcdc_cost = (self.res['n_dcdc_pemel'] *
                           (self.par['capex_dcdc'] *
                            (crf + self.par['opex_dcdc'])))
        components_cost += pemel_cost + pemel_dcdc_cost

        # annual cost of fuel cell array and DC-DC converter
        pemfc_cost = self.par['n_pemfc'] * (self.par['capex_pemfc'] * crf +
                                            self.par['opex_pemfc'] *
                                            self.res['running_hours_pemfc'])
        pemfc_dcdc_cost = self.res['n_dcdc_pemfc'] * (self.par['capex_dcdc'] *
                                                      (crf +
                                                       self.par['opex_dcdc']))
        components_cost += pemfc_cost + pemfc_dcdc_cost

        # annual cost of hydrogen storage tank
        tank_cost = self.par['n_tank'] * (self.par['capex_tank'] *
                                          (crf + self.par['opex_tank']))
        components_cost += tank_cost

        # annual cost of DC-AC inverter
        dcac_cost = self.res['n_dcac'] * (self.par['capex_dcac'] *
                                          (crf + self.par['opex_dcac']))
        components_cost += dcac_cost

        # annual replacement cost of the electrolyzer and fuel cell
        arc = crf * sum([(1. + inv_rate)**(-(i + 1.) *
                                           self.res['life_pemel']) *
                         self.par['n_pemel'] *
                         self.par['repl_pemel'] *
                         self.par['capex_pemel'] for i in
                         range(int(self.par['life_sys'] /
                                   self.res['life_pemel']))])
        arc += crf * sum([(1. + inv_rate)**(-(i + 1.) *
                                            self.res['life_pemfc']) *
                          self.par['n_pemfc'] *
                          self.par['repl_pemfc'] *
                          self.par['capex_pemfc'] for i in
                          range(int(self.par['life_sys'] /
                                    self.res['life_pemfc']))])

        grid_e_cost = sum(self.res['grid_e_buy'] * self.elec_profile)
        grid_e_gain = sum(self.res['grid_e_sold'] * self.elec_profile_sale)

        # annual cost
        cost = arc + components_cost + grid_e_cost - grid_e_gain

        # levelized cost of electricity in euro/MWh
        self.res['lcoe'] = cost / (sum(self.load_elec)) * 1e6

    def print_results(self):
        """

        This method prints the levelized cost of electricity,
        the self-sufficiency ratio and the annual energy produced
        by the photovoltaic array.

        """
        print('outputs:')
        print('LCOE:'.ljust(30) + '%.2f euro/MWh' % self.res['lcoe'])
        print('SSR:'.ljust(30) + '%.2f %%' % (self.res['ssr'] * 100.))
        print('PV electricity generated:'.ljust(30) +
              '%.2f MWh' % (sum(self.res['p_pv']) / 1e6))
        print('grid electricity bought:'.ljust(30) +
              '%.2f MWh' % (sum(self.res['grid_e_buy']) / 1e6))
        print('grid electricity sold:'.ljust(30) +
              '%.2f MWh' % (sum(self.res['grid_e_sold']) / 1e6))
        print(
            'lifetime electrolyzer:'.ljust(30) +
            '%.2f year' %
            self.res['life_pemel'])
        print(
            'lifetime fuel cell:'.ljust(30) +
            '%.2f year' %
            self.res['life_pemfc'])

        plt.plot(self.res['m_h2_array'])
        plt.show(block=False)


def polyfit_func(x_in, y_in, threshold=0.9999999):
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
