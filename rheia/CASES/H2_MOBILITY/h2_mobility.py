"""
The :py:mod:`h2_mobility` module contains a class to read the required data and
a class to evaluate the power-to-mobility system.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pvlib


class ReadData:
    """

    This class enables to read data from the data files.

    Parameters
    ----------
    filename_climate : str
        The directory of the file with information on the
        solar irradiance.

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
            Meteorological Year.
        t_amb : ndarray
            The hourly ambient temperature data for a Typical
            Meteorological Year.

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

    This class evaluates the power-to-mobility system.
    For a given design, the solar irradiance, ambient temperature
    and the characterization of the model parameters,
    the levelized cost of driving, carbon intensity and the annual
    grid consumption are quantified.

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

        self.length = len(self.sol_irr)

        # the result dictionary
        self.res = {}

        # the system lifetime
        self.par['life_sys'] = 20.

        # initialize the operating hours of the electrolyzer array
        self.res['running_hours_pemel'] = 0.

        # initialize the storage tank size and its starting status
        self.m_h2_max = self.tank()
        self.m_h2_min = 0.05 * self.m_h2_max
        self.m_h2 = self.m_h2_min

        # instantiate the profiles for grid electricity price
        self.demand_profiles()

        # the number of PEM electrolyzer cells, corresponding to the
        # nominal capacity of the considered PEM cell and the provided
        # PEM capacity
        self.n_pemel_array = self.par['n_pemel'] / 0.4

        # the fitted polynomials on the electrolyzer and compressor
        self.polyfit_pemel()
        self.polyfit_pemel_compr()

    def demand_profiles(self):
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

        In addition, the demand profiles from the hydrogen buses and diesel
        buses is determined, based on the European daily refueling profile [1].

        [1] T. A. Gunawan, I.Williamson, D. Raine, and R. F. Monaghan,
        “Decarbonising city bus networks in ireland with renewable hydrogen,”
        International Journal of Hydrogen Energy, 2020.

        """

        # electricity cost profile [euro/Wh]
        self.elec_profile = np.ones(self.length) * (
            (self.par['elec_cost'] +
             self.par['elec_cost_profit']) /
            self.par['elec_cost_ratio']) / 1e6

        # electricity selling profile [euro/Wh]
        self.elec_profile_sale = np.ones(
            self.length) * self.par['elec_cost'] / 1e6

        self.diesel_profile = np.ones(self.length) * self.par['diesel_cost']

        # number of km driven per day per bus
        self.par['n_km_bus'] = 250.

        # number of buses in the fleet
        self.par['n_bus'] = 40.

        # energy consumed by the diesel buses and hydrogen buses per day
        # [kWh/day]
        energy_h2 = (self.par['cons_h2_bus'] * self.par['n_km_bus'] *
                     self.par['n_h2_bus'])
        energy_diesel = (self.par['cons_diesel_bus'] * self.par['n_km_bus'] *
                         (self.par['n_bus'] - self.par['n_h2_bus']))

        h2_required = energy_h2 / 33.33  # kg
        diesel_required = energy_diesel / 10.  # litre

        # the daily refueling profile of the buses.
        fill_profile = np.array([0.09, 0.015, 0.005, 0.04, 0.04, 0., 0.01,
                                 0.01, 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                 0.08, 0.08, 0.13, 0.13, 0.13, 0.13, 0.11])

        # daily refueling profile for the hydrogen buses and diesel buses
        day_h2 = np.ones(24) * fill_profile * h2_required
        day_diesel = np.ones(24) * fill_profile * diesel_required

        # annual refueling profiles
        self.load_h2 = list(day_h2) * int(365 * self.length / 8760)
        self.load_diesel = list(day_diesel) * int(365 * self.length / 8760)

        # dispenser capacity such that the hourly hydrogen demand can be
        # complied with
        dispenser_mass_flow_rate = 33.33
        self.par['n_disp'] = max(self.load_h2) / dispenser_mass_flow_rate

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

        # determine the maximum power point at reference conditions
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

    #############################
    # electrolyzer array module #
    #############################

    def pemel(self, i_pemel):
        """
        The electrolyzer model, based on the work of Saeed et al. [2]. For a
        given current, the model determines the operating voltage by
        considering the activation, concentration and ohmic overpotentials.
        The model quantifies the operating voltage, power, efficiency and
        hydrogen production.

        [2] Saeed, E. W., & Warkozek, E. G. (2015). Modeling and Analysis of
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

    def mh2_to_power(self, m_h2):
        """
        When the hydrogen mass flow rate is provided, this function determines
        the corresponding required power per hour.

        Parameters
        ----------
        m_h2 : float
            The produced hydrogen mass flow rate [kg/h].

        Returns
        -------
        power : float
            The required power to produce the hydrogen [W].

        """

        far_cons = 96485.
        current = m_h2 * (2. * far_cons) / (2.02e-3 * 3600. *
                                            self.n_pemel_array)
        power = self.pemel(current)['p_pemel']

        return power

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

    #####################
    # compressor module #
    #####################

    def compressor(self, m_h2):
        """
        The compressor module defined the required compression power to
        compress the hydrogen mass flow rate [3].

        [3] Zhao, L., Brouwer, J., & Samuelsen, S. (2014). Dynamic analysis of
        a self-sustainable renewable hydrogen fueling station. ASME 2014 12th
        International Conference on Fuel Cell Science, Engineering and
        Technology, FUELCELL 2014 Collocated with the ASME 2014 8th
        International Conference on Energy Sustainability.
        https://doi.org/10.1115/FuelCell2014-6330

        Parameters
        ----------
        m_h2 : float
            Hydrogen mass flow rate [kg/h].

        Returns
        -------
        power : float
            The required compression power [W].

        """

        # convert the flow rate into kg/s
        m_h2 *= 1. / 3600.

        par_c = {
            'T_in': 353.,
            'p_in': 20.,
            'p_out': 440.,
            'eta_c': 0.85,
            'R': 4.124,
            'n': 1.609,
        }

        power = (m_h2 *
                 par_c['n'] *
                 par_c['R'] *
                 par_c['T_in'] *
                 ((par_c['p_out'] /
                   par_c['p_in'])**((par_c['n'] -
                                     1.) /
                                    par_c['n']) -
                     1.) *
                 1000. /
                 (par_c['eta_c'] *
                     (par_c['n'] -
                      1.)))

        return power

    def polyfit_pemel_compr(self):
        """
        The power consumption by the electrolyzer stack and compressor are
        evaluated over a range of hydrogen mass flow rates. Following these
        evaluations, a polynomial is fitted on the mass flow rate - power
        relation. This polynomial enables to rapidly determine the input
        mass flow rate when a certain amount of power is available.
        Since this relation is fairly linear, the polynomial should reach good
        agreement with the actual mass flow rate - power relation, while
        maintaining the level of fidelity of the actual model.

        """

        # the electrolyzer array operating limits
        pemel_lower_lim = self.par['n_pemel'] * 10.
        pemel_upper_lim = self.par['n_pemel'] * 1e3

        # the operating current at these limits
        current_ub = self.p_to_i_pemel(pemel_upper_lim)
        current_lb = self.p_to_i_pemel(pemel_lower_lim)

        # the characteristics of the electrolyzer at the limits
        pemel_ub = self.pemel(current_ub)
        pemel_lb = self.pemel(current_lb)

        # the compression power at the hydrogen production limits
        p_compr_ub = self.compressor(pemel_ub['m_pemel'])
        p_compr_lb = self.compressor(pemel_lb['m_pemel'])

        # the compression and electrolyzer power at the limits
        pemel_compr_ub = pemel_ub['p_pemel'] + p_compr_ub
        pemel_compr_lb = pemel_lb['p_pemel'] + p_compr_lb

        # the operating bounds for the compressor and electrolyzer
        self.bounds = {'pemel_lb': pemel_lb,
                       'pemel_ub': pemel_ub,
                       'pemel_compr_lb': pemel_compr_lb,
                       'pemel_compr_ub': pemel_compr_ub
                       }

        # evaluate the electrolyzer and compressor for a set of mass flow rates
        step = pemel_ub['m_pemel'] / 50.
        m_h2_list = np.arange(start=step, stop=pemel_ub['m_pemel'] - step,
                              step=step)
        p_pemel_compr = np.zeros(len(m_h2_list))
        for index, m_h2 in enumerate(m_h2_list):
            p_pemel = self.mh2_to_power(m_h2)
            p_compr = self.compressor(m_h2)
            p_pemel_compr[index] = p_pemel + p_compr

        # generate a polynomial fitted on the power - mass flow rate points
        self.p_to_m_pemel_comp = polyfit_func(p_pemel_compr, m_h2_list)

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

    ##############################
    # management strategy module #
    ##############################

    def p_for_inst_demand(self, m_h2):
        """
        When a hydrogen demand is provided, this method determines the power
        to generate the hydrogen in the electrolyzer array and the power to
        compress the hydrogen in the compressor. The method indicates when the
        electrolyzer array capacity is not sufficient to produce the desired
        hydrogen mass flow rate.

        Parameters
        ----------
        m_h2 : float
            Hydrogen mass flow rate [kg/h].

        Returns
        -------
        p_pemel : float
            The required electrolyzer array power [W].
        p_compr : float
            The required compressor power [W].
        bool
            True when the electrolyzer capacity is insufficient to produce the
            hydrogen.

        """

        # check if the desired hydrogen can be produced by the PEM
        if m_h2 > self.bounds['pemel_ub']['m_pemel']:
            return 0., 0., True

        # determine the power needed from electrolyzer and compressor to
        # deliver the required hydrogen mass flow rate
        p_pemel = self.mh2_to_power(m_h2)
        p_compr = self.compressor(m_h2)
        self.res['running_hours_pemel'] += 1.
        return p_pemel, p_compr, False

    def prod_mh2(self, p_in):
        """
        When there is power available to produce hydrogen, this method
        distributes this power over the electrolyzer array and compressor,
        such that the hydrogen is produced and compressed with this given
        power. If the given power will lead to an excess of hydrogen (i.e
        the storage tank is full before the energy is fully consumed), the
        power supplied to the electrolyzer array and compressor is recalculated
        such that no excess hydrogen is produced.

        Parameters
        ----------
        p_in : float
            The power available to produce hydrogen [W].

        Returns
        -------
        p_pemel : float
            The power consumed by the electrolyzer array [W].
        p_compr : float
            The power consumed by the compressor [W].

        """

        no_run = False

        # produce H2 only if tank is not yet full
        # nothing happens when the desired power is under the lower limit
        if self.m_h2 >= self.m_h2_max or p_in < self.bounds['pemel_compr_lb']:
            p_pemel = 0.
            p_compr = 0.
            no_run = True  # the electrolyzer did not run

        # if the power is higher than the upper limit, work at the upper limit
        elif p_in > self.bounds['pemel_compr_ub']:

            # this is what can be produced at the upper limit for H2
            m_h2 = self.bounds['pemel_ub']['m_pemel']

            # produced hydrogen is added to the tank
            self.m_h2 += m_h2

        else:

            # quantify the hydrogen created with the available power
            m_h2 = self.p_to_m_pemel_comp(p_in)

            # produced hydrogen is added to the tank
            self.m_h2 += m_h2

        # if the new hydrogen capacity exceeds the maximum storage capacity
        if self.m_h2 > self.m_h2_max:

            # the hydrogen that is still allowed in the tank
            m_h2 -= (self.m_h2 - self.m_h2_max)

            # the power for the PEM to generate this hydrogen
            p_rev = self.mh2_to_power(m_h2)

            # check if this power is higher than the lower limit for the PEM
            if p_rev > self.bounds['pemel_lb']['p_pemel']:
                self.m_h2 = self.m_h2_max

            # the space left in the tank is too small
            else:
                self.m_h2 -= (m_h2 + self.m_h2 - self.m_h2_max)
                p_pemel = 0.
                p_compr = 0.
                no_run = True

        # if power is applied to the PEM and compessor
        if not no_run:
            p_compr = self.compressor(m_h2)
            p_pemel = self.mh2_to_power(m_h2)
            self.res['running_hours_pemel'] += 1.

        return p_pemel, p_compr

    def extract_h2_from_tank(self, demand):
        """
        Extract the hydrogen demand from the storage tank. If more hydrogen is
        demanded than available in the tank, extract only the available amount
        of hydrogen.

        Parameters
        ----------
        demand : float
            The hydrogen demand. [kg]

        Returns
        -------
        demand_left : float
            The amount of hydrogen demand that is not covered by the tank. [kg]

        """

        # extract the hydrogen from the tank
        self.m_h2 -= demand

        if self.m_h2 < self.m_h2_min:

            # if the demand was too high, set the tank to its minimum level and
            # define the demand that was not covered by the tank
            demand_left = self.m_h2_min - self.m_h2
            self.m_h2 = self.m_h2_min
        else:

            # if the tank covers the demand, there is no hydrogen demand left
            demand_left = 0.

        return demand_left

    ##########################################
    # model evaluation
    ##########################################

    def evaluation(self):
        """

        This is the main method of the Evaluation class.
        For each hour, the power management strategy is applied.
        The hydrogen demand is extracted from the hydrogen tank, when
        sufficient hydrogen is available in the tank. When the hydrogen in the
        tank does not comply with the hydrogen demand, the power to run the
        electrolyzer array and compressor is calculated to generate and
        compress the remaining hydrogen. To generate this power,
        the photovoltaic power is called upon first. If
        necessary, the remaining power is covered by the grid.
        When excess photovoltaic power is present, this power is used to
        generate and compress hydrogen in the electrolyzer array and
        compressor, respectively. Finally, the lifetime, cost and
        CO2-emission of the system are quantified.

        Returns
        -------
        bool
            True when the capacity of the electrolyzer array is sufficient to
            cover the instantaneous hydrogen demand during the year.

        """

        n_compr = np.zeros(self.length)
        n_dcdc_pem = np.zeros(self.length)
        n_dcac = np.zeros(self.length)

        self.res['m_h2_array'] = np.ones(self.length)
        self.res['grid_e_buy'] = np.ones(self.length)
        self.res['grid_e_sold'] = np.ones(self.length)
        self.res['grid_co2'] = 0.

        # get the hourly photovoltaic array power
        self.photovoltaic()

        for t in range(self.length):
            e_grid_buy = 0.
            e_grid_sold = 0.

            # define if there is any H2 demand left after assessing the tank
            demand_left = self.extract_h2_from_tank(self.load_h2[t])

            if demand_left > 0.:

                # power needed by the PEM and compressor to generate the demand
                p_pemel, p_compr, check = self.p_for_inst_demand(demand_left)
                n_dcdc_pem[t] = p_pemel
                n_compr[t] = p_compr

                # check if PEM capacity can cover the demand
                if check:
                    return False

                # positive if the PV energy can comply with the PEM power
                net_p_1 = self.res['p_pv'][t] - p_pemel

                if net_p_1 > 0.:
                    # if yes, check if remaining PV power can comply with the
                    # compressor
                    net_p_2 = net_p_1 - p_compr

                    if net_p_2 > 0.:
                        # when still excess power available, sell to the grid
                        e_grid_sold = net_p_2
                        n_dcac[t] = e_grid_sold

                    else:
                        # if the compressor demand cannot be
                        # covered, buy remaining demand from the grid
                        e_grid_buy = abs(net_p_2)

                else:
                    # if PEM and compressor cannot be covered by
                    # PV energy, buy required electricity from the grid
                    p_to_buy = abs(net_p_1) + p_compr
                    e_grid_buy = p_to_buy
                    n_dcac[t] = abs(net_p_1)

            else:  # excess PV energy available to generate hydrogen

                if self.res['p_pv'][t] > 0.:
                    # quantify how much of the remaining PV energy can be used
                    # by the PEM and compressor to generate hydrogen
                    p_pemel, p_compr = self.prod_mh2(self.res['p_pv'][t])
                    n_compr[t] = p_compr

                    # the remaining excess PV energy can be sold
                    e_grid_sold = self.res['p_pv'][t] - p_pemel - p_compr
                    n_dcac[t] = e_grid_sold + p_compr
                    n_dcdc_pem[t] = p_pemel

            # store evolution of storage tank, electricity bought/sold
            # and the CO2 amount from buying grid electricity
            self.res['m_h2_array'][t] = ((self.m_h2 - self.m_h2_min) /
                                         (self.m_h2_max - self.m_h2_min))
            self.res['grid_e_buy'][t] = e_grid_buy
            self.res['grid_e_sold'][t] = e_grid_sold
            self.res['grid_co2'] += e_grid_buy * self.par['co2_elec']

        # define the capacity of the converters and compression
        self.res['n_compr'] = max(n_compr) / 1e3
        self.res['n_dcdc_pemel'] = max(n_dcdc_pem) / 1e3
        self.res['n_dcac'] = max(n_dcac) / 1e3

        # the cost of annual diesel consumption
        self.res['diesel_cost'] = self.load_diesel * self.diesel_profile

        # determine the lifetime of the electrolyzer
        self.lifetime()

        # determine the system cost
        self.cost()

        # determine the system annual CO2 emission
        self.lca()

        return True

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
        the levelized cost of mobility [euro/km] is determined. The formula
        for the annualized system cost is adopted from Coppitters et al. [5].

        [5] Coppitters, D., De Paepe, W., & Contino, F. (2020). Robust design
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

        # annual cost of hydrogen storage tank
        tank_cost = self.par['n_tank'] * (self.par['capex_tank'] *
                                          (crf + self.par['opex_tank']))
        components_cost += tank_cost

        # annual cost of compressor
        compressor_cost = (self.par['capex_compr'] *
                           self.res['n_compr'] *
                           (crf + self.par['opex_compr']))
        components_cost += compressor_cost

        # annual cost of dispenser
        dispenser_cost = self.par['n_disp'] * (self.par['capex_disp'] *
                                               (crf + self.par['opex_disp']))
        components_cost += dispenser_cost

        # annual cost of DC-AC inverter
        dcac_cost = self.res['n_dcac'] * (self.par['capex_dcac'] *
                                          (crf + self.par['opex_dcac']))
        components_cost += dcac_cost

        # annual cost of buses
        h2_bus_cost = (self.par['capex_h2_bus'] * crf + self.par['n_km_bus'] *
                       365. * self.par['opex_h2_bus']) * self.par['n_h2_bus']
        diesel_bus_cost = (self.par['capex_diesel_bus'] * crf +
                           self.par['opex_diesel_bus'] * 365. *
                           self.par['n_km_bus']) * (self.par['n_bus'] -
                                                    self.par['n_h2_bus'])
        components_cost += h2_bus_cost + diesel_bus_cost

        # annual replacement cost of the electrolyzer and buses
        arc = crf * sum([(1. + inv_rate)**(-(i + 1.) *
                                           self.res['life_pemel']) *
                         self.par['n_pemel'] *
                         self.par['repl_pemel'] *
                         self.par['capex_pemel'] for
                         i in range(int(self.par['life_sys'] /
                                        self.res['life_pemel']))])
        arc += crf * sum([(1. + inv_rate)**(-(i + 1.) * 10.) *
                          self.par['n_h2_bus'] * self.par['repl_h2_bus']
                          for i in range(int(self.par['life_sys'] / 10.))])
        arc += crf * sum([(1. + inv_rate)**(-(i + 1.) * 10.) *
                          (self.par['n_bus'] - self.par['n_h2_bus']) *
                          self.par['repl_diesel_bus'] for i in
                          range(int(self.par['life_sys'] / 10.))])

        grid_e_cost = sum(self.res['grid_e_buy'] * self.elec_profile)
        grid_e_gain = sum(self.res['grid_e_sold'] * self.elec_profile_sale)

        # total annual cost
        cost = (arc + components_cost + grid_e_cost - grid_e_gain +
                sum(self.res['diesel_cost']))

        # levelized cost of mobility in euro/km
        self.res['lcom'] = cost / (self.par['n_km_bus'] * self.par['n_bus'] *
                                   365.)

    def lca(self):
        """
        The life cycle assessment is performed based on the CO2 emissions from
        constructing the system components and the emissions related to
        consuming grid electricity and diesel. The annual CO2-equivalent
        emissions of the system is divided by the annual distance covered by
        the bus fleet, resulting in the levelized cost of mobility.

        """

        # annual CO2 emission of photovoltaic array production
        pv_co2 = self.par['n_pv'] * self.par['co2_pv']
        pv_dcdc_co2 = self.res['n_dcdc_pv'] * self.par['co2_dcdc']
        comp_co2 = pv_co2 + pv_dcdc_co2

        # annual CO2 emission of electrolyzer array production
        pemel_co2 = (self.par['n_pemel'] * self.par['co2_pemel'] *
                     (1 + int(self.par['life_sys'] / self.res['life_pemel'])))
        pemel_dcdc_co2 = self.res['n_dcdc_pemel'] * self.par['co2_dcdc']
        comp_co2 += pemel_co2 + pemel_dcdc_co2

        # annual CO2 emission of hydrogen storage tank production
        tank_co2 = self.par['n_tank'] * self.par['co2_tank']
        comp_co2 += tank_co2

        # annual CO2 emission of compressor production
        compressor_co2 = self.res['n_compr'] * self.par['co2_compr']
        comp_co2 += compressor_co2

        # annual CO2 emission of dispenser production
        dispenser_co2 = self.par['n_disp'] * self.par['co2_disp']
        comp_co2 += dispenser_co2

        # annual CO2 emission of DC-AC inverter production
        dcac_co2 = self.res['n_dcac'] * self.par['co2_dcac']
        comp_co2 += dcac_co2

        # annual CO2 emission of diesel and hydrogen engine production
        diesel_engine_co2 = (self.par['co2_diesel_engine'] * 200. *  # 200 kW
                             (self.par['n_bus'] - self.par['n_h2_bus']))
        h2_engine_co2 = self.par['co2_fc_engine'] * 200. * self.par['n_h2_bus']
        comp_co2 += ((diesel_engine_co2 + h2_engine_co2) *
                     (1 + int(self.par['life_sys'] / 11.)))  # 10y lifetime

        # annual CO2 emission of diesel consumption
        diesel_co2 = sum(self.load_diesel) * self.par['co2_diesel']

        # annual CO2 emission of system
        co2 = (comp_co2 / self.par['life_sys'] + self.res['grid_co2'] +
               diesel_co2)

        # CO2 emitted per km driven by the fleet
        self.res['ci'] = co2 / (self.par['n_km_bus'] * self.par['n_bus'] *
                                  365. * self.length / 8760.)

    def print_results(self, succes=True):
        """

        This method prints the levelized cost of electricity,
        the self-sufficiency ratio and the annual energy produced
        by the photovoltaic array.

        """
        if not succes:
            print("""Evaluation failed: the electrolyzer array capacity
                     of %f kW was not sufficient to cover the instantaneous
                     hydrogen demand.""" % self.par['n_pemel'])
        else:
            print('outputs:')
            print('LCOE:'.ljust(30) + '%.5f euro/km' % self.res['lcom'])
            print('CI:'.ljust(30) + '%.5f kg co2-eq/km' % self.res['ci'])
            print(
                'PV electricity generated:'.ljust(30) +
                '%.5f MWh' %
                (sum(self.res['p_pv']) / 1e6))
            print('grid energy bought:'.ljust(30) + '%.5f MWh' %
                  (sum(self.res['grid_e_buy']) / 1e6))
            print('grid energy sold:'.ljust(30) + '%.5f MWh' %
                  (sum(self.res['grid_e_sold']) / 1e6))
            print(
                'compressor capacity:'.ljust(30) +
                '%.5f kW' %
                self.res['n_compr'])
            print('life electrolyzer:'.ljust(30) + '%.5f year' %
                  self.res['life_pemel'])

            plt.plot(self.res['m_h2_array'])
            plt.show(block=False)


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
