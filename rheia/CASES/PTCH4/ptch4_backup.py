"""
The :py:mod:`h2_fuel` module contains a class to read the required data and
a class to evaluate the power-to-fuel system.
"""

import os
import pandas as pd
import numpy as np
import pvlib
import matplotlib.pyplot as plt
try:
    import brightway2 as bw
    from brightway2 import *
except ImportError:
    pass


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
        and ambient temperature data and wind data,
        situated in the 'sol_irr' and 'T_amb' and 'wind_pow'columns of the
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
        wind_pow = data['wind_pow'].to_numpy()

        return sol_irr, t_amb, wind_pow

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
        design_space = os.path.join(self.path, 'design_space.csv')

        # read the deterministic values for the parameters in `design_space`
        with open(design_space, 'r') as file:
            for line in file:
                tmp = line.split(",")
                if tmp[1] == 'par':
                    param_dict[tmp[0]] = float(tmp[2])

        return param_dict


class Evaluation:
    """

    This class evaluates the system.
    For a given design, the solar irradiance, ambient temperature
    and the characterization of the model parameters,
    the levelized cost of hydrogen and the PB impacts are evaluated

    Parameters
    ----------
    sol_irr : ndarray
        The hourly solar irradiance for the evaluated year.
    t_amb : ndarray
        The hourly ambient temperature for the evaluated year.
    parameters : dict
        Dictionary with the model parameters and design variables values.

    """

    def __init__(self, sol_irr, t_amb, wind, bog, full_lca, dict_lcia, dict_lcia_pb, all_dbs, scenario, par):
        self.par = par
        self.full_lca = full_lca
        self.dict_lcia = dict_lcia
        self.dict_lcia_pb = dict_lcia_pb
        self.all_dbs = all_dbs
        
        self.bog = bog
        
        self.smr = True
        
        # the system lifetime
        self.par['life_sys'] = 30.


        self.pb_names = [
                'climate change CO2 concentration',
                'climate change energy imbalance',
                'stratospheric ozone depletion',
                'ocean acidification',
                'biogeochemical flows P',
                'biogeochemical flows N',
                'land system change global',
                'freshwater use global',
                'atmospheric aerosol loading',
                'biosphere integrity',]

                       
        self.pb_limits = np.array([ 72., 1., 15., 0.688, 6.2, 62., 25., 4000., 0.11, 195000.])
        
        self.pb_limits_uncertain = np.array([ 172., 1.5, 29., 2.75*0.3, 11.2, 82., 46., 6000., 0.36, 7.*195000.])
        
        
        self.dbs = [
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg1150_2030',
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg1150_2040',
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg1150_2050',
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg500_2030',
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg500_2040',
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg500_2050',
               'ecoinvent_cutoff_3.9_remind_SSP2-Base_2030',
               'ecoinvent_cutoff_3.9_remind_SSP2-Base_2040',
               'ecoinvent_cutoff_3.9_remind_SSP2-Base_2050',
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg1150_2030',
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg1150_2040',
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg1150_2050',
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg500_2030',
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg500_2040',
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg500_2050',
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2030',
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2040',
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2050',
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg1150_2030',
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg1150_2040',
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg1150_2050',
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg500_2030',
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg500_2040',
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg500_2050',
               'ecoinvent_cutoff_3.9_remind_SSP5-Base_2030',
               'ecoinvent_cutoff_3.9_remind_SSP5-Base_2040',
               'ecoinvent_cutoff_3.9_remind_SSP5-Base_2050',
               ]

        self.pb = {}
        self.pb_limits_pros = {}
        
        self.pb_gdp = {
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg1150_2030': 3.3159E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg1150_2040': 3.1828E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg1150_2050': 3.13843E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg500_2030': 3.741E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg500_2040': 3.90351E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg500_2050': 3.84133E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP2-Base_2030': 3.03636E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP2-Base_2040': 2.69498E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP2-Base_2050': 2.47912E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg1150_2030': 3.78721E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg1150_2040': 3.63542E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg1150_2050': 3.60115E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg500_2030': 5.30869E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg500_2040': 5.07707E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg500_2050': 4.5867E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2030': 3.52249E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2040': 3.15423E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2050': 2.96551E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg1150_2030': 2.65014E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg1150_2040': 2.28696E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg1150_2050': 2.04255E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg500_2030': 2.83382E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg500_2040': 2.6826E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg500_2050': 2.48155E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP5-Base_2030': 2.52959E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP5-Base_2040': 1.94297E-06 * 1.874932065,
               'ecoinvent_cutoff_3.9_remind_SSP5-Base_2050': 1.57375E-06 * 1.874932065,}

        self.pb_pop_all_h2 = {
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2010': 1.82E-04,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg1150_2030': 2.57E-03,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg1150_2040': 9.69E-04,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg1150_2050': 7.45E-04,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg500_2030': 2.71E-03,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg500_2040': 1.72E-03,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg500_2050': 8.73E-04,
               'ecoinvent_cutoff_3.9_remind_SSP2-Base_2030': 3.29E-03, #eigenlijk geen H2 in dit scenario!
               'ecoinvent_cutoff_3.9_remind_SSP2-Base_2040': 3.29E-03,
               'ecoinvent_cutoff_3.9_remind_SSP2-Base_2050': 8.17E-04,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg1150_2030': 2.12E-03,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg1150_2040': 6.44E-04,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg1150_2050': 5.22E-04,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg500_2030': 3.80E-04,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg500_2040': 9.12E-05,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg500_2050': 7.30E-05,
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2030': 2.16E-03,
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2040': 6.64E-04,
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2050': 6.35E-04,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg1150_2030': 1.27E-02,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg1150_2040': 3.79E-03,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg1150_2050': 9.87E-04,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg500_2030': 3.37E-03,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg500_2040': 3.58E-04,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg500_2050': 1.50E-04,
               'ecoinvent_cutoff_3.9_remind_SSP5-Base_2030': 1.09E-02,
               'ecoinvent_cutoff_3.9_remind_SSP5-Base_2040': 7.60E-03,
               'ecoinvent_cutoff_3.9_remind_SSP5-Base_2050': 4.99E-03,}

        self.pb_pop_grandf_energy_system_on_co2 = {
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2010': 4.54E-06,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg1150_2030': 4.14E-06,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg1150_2040': 4.18E-06,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg1150_2050': 4.26E-06,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg500_2030': 4.67E-06,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg500_2040': 5.13E-06,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg500_2050': 5.22E-06,
               'ecoinvent_cutoff_3.9_remind_SSP2-Base_2030': 3.79E-06, 
               'ecoinvent_cutoff_3.9_remind_SSP2-Base_2040': 3.54E-06,
               'ecoinvent_cutoff_3.9_remind_SSP2-Base_2050': 3.37E-06,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg1150_2030': 4.91E-06,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg1150_2040': 5.16E-06,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg1150_2050': 5.53E-06,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg500_2030': 6.89E-06,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg500_2040': 7.21E-06,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg500_2050': 7.04E-06,
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2030': 4.57E-06,
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2040': 4.48E-06,
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2050': 4.55E-06,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg1150_2030': 3.53E-06,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg1150_2040': 3.41E-06,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg1150_2050': 3.35E-06,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg500_2030': 3.77E-06,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg500_2040': 4.00E-06,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg500_2050': 4.07E-06,
               'ecoinvent_cutoff_3.9_remind_SSP5-Base_2030': 3.37E-06,
               'ecoinvent_cutoff_3.9_remind_SSP5-Base_2040': 2.90E-06,
               'ecoinvent_cutoff_3.9_remind_SSP5-Base_2050': 2.58E-06,}
  
        self.pb_pop_eu = {
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg1150_2030': 0.0657,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg1150_2040': 0.0638,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg1150_2050': 0.0621,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg500_2030': 0.0657,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg500_2040': 0.0638,
               'ecoinvent_cutoff_3.9_remind_SSP2-PkBudg500_2050': 0.0621,
               'ecoinvent_cutoff_3.9_remind_SSP2-Base_2030': 0.0657,
               'ecoinvent_cutoff_3.9_remind_SSP2-Base_2040': 0.0638,
               'ecoinvent_cutoff_3.9_remind_SSP2-Base_2050': 0.0621,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg1150_2030': 0.0675,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg1150_2040': 0.065,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg1150_2050': 0.0658,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg500_2030': 0.0675,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg500_2040': 0.065,
               'ecoinvent_cutoff_3.9_remind_SSP1-PkBudg500_2050': 0.0658,
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2030': 0.0675,
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2040': 0.065,
               'ecoinvent_cutoff_3.9_remind_SSP1-Base_2050': 0.0658,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg1150_2030': 0.0692,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg1150_2040': 0.069,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg1150_2050': 0.0692,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg500_2030': 0.0692,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg500_2040': 0.069,
               'ecoinvent_cutoff_3.9_remind_SSP5-PkBudg500_2050': 0.0692,
               'ecoinvent_cutoff_3.9_remind_SSP5-Base_2030': 0.0692,
               'ecoinvent_cutoff_3.9_remind_SSP5-Base_2040': 0.069,
               'ecoinvent_cutoff_3.9_remind_SSP5-Base_2050': 0.0692,}               

        self.scenario = scenario
        for db in self.dbs:
            for index,method in enumerate(self.pb_names):
                self.pb_limits_pros[method + db.split('remind')[1]] = self.pb_limits[index]
                
                if self.scenario == 'europe_10mil':
     
                    self.pb_limits_pros[method + db.split('remind')[1]] *= self.pb_pop_eu[db]
                    self.pb_limits_pros[method + db.split('remind')[1]] *= 1e9
                
                elif self.scenario == 'gdp':

                    self.pb_limits_pros[method + db.split('remind')[1]] *= self.pb_gdp[db]
                
                elif self.scenario == 'energy':
                
                    self.pb_limits_pros[method + db.split('remind')[1]] *= self.pb_pop_grandf_energy_system_on_co2[db]
                                
                self.pb_limits_pros[method + db.split('remind')[1]] *= 1e-9

        if self.scenario == 'europe_10mil':
                    
            # scenario 0    
            self.pb_limits *= self.pb_pop_eu['ecoinvent_cutoff_3.9_remind_SSP2-Base_2030'] # roughly 6.7% of population will live in Europe in 2030
            self.pb_limits *=  1e9 # to counter the general 1e-9 a bit further in the code
        
        elif self.scenario == 'gdp':

            # scenario 1: GDP 

            self.pb_limits *= self.pb_gdp['ecoinvent_cutoff_3.9_remind_SSP2-Base_2030'] # roughly 6.7% of population will live in Europe in 2030

        elif self.scenario == 'energy':

            # scenario 2: 
            self.pb_limits *= self.pb_pop_grandf_energy_system_on_co2['ecoinvent_cutoff_3.9_remind_SSP2-Base_2030'] # roughly 6.7% of population will live in Europe in 2030

        self.pb_limits *= 1e-9 # budget for production of 1 kWh

                

        # the solar irradiance and ambient temperature are shifted based on starting point

        sol_irr = np.roll(sol_irr,-int(self.par['start']))
        t_amb = np.roll(t_amb,-int(self.par['start']))
        wind = np.roll(wind,-int(self.par['start']))

        self.sol_irr = np.concatenate((sol_irr, sol_irr[:int(self.par['delay_h'])])) * self.par['u_sol_irr']
        self.t_amb = np.concatenate((t_amb, t_amb[:int(self.par['delay_h'])])) + self.par['u_t_amb']
        self.wind = np.concatenate((wind, wind[:int(self.par['delay_h'])])) * self.par['u_wind']

        # the result dictionary
        self.res = {}
        
        self.par['n_pv'] *= 1e9
        self.par['n_wind'] *= 1e9
        self.par['n_pemel'] *= 1e9
        self.par['n_bat'] *= 1e9
        self.par['n_tank_h2'] *= 1e6

        self.m_h2_max = self.par['n_tank_h2']

        # initialize the operating hours of the electrolyzer array
        self.res['running_hours_pemel'] = 0.
        
            
    #############################
    # photovoltaic array module #
    #############################

    def photovoltaic(self):
        """

        The hourly photovoltaic power is quantified via CFs.

        """

        p_pv = np.zeros(len(self.sol_irr))

        # power determination for each hour in the timeframe
        for i, irr in enumerate(self.sol_irr):
            if irr > 0.:
                p_pv[i] = irr * self.par['n_pv'] # Wh

        # store the hourly pv power in the result dictionary
        self.res['p_pv'] = p_pv

    def wind_turbine(self):
        """

        The hourly wind profile is calculated by CFs.

        """

        p_wind = np.zeros(len(self.sol_irr))

        for i, wind in enumerate(self.wind):
            if wind > 0.:
                p_wind[i] = wind * self.par['n_wind']  # Wh

        # store the hourly pv power in the result dictionary
        self.res['p_wind'] = p_wind
        
    def p_demand_nominal(self):
    
        lhv_ch4 = 50. * 0.27777 * 1e-9 # TWh/kg
        
        # ammonia per hour required
        h_to_supply_demand = 8760.
        m_ch4 = self.par['demand_TWh']/h_to_supply_demand/lhv_ch4 # kg/h
        
        self.par['n_meth'] = m_ch4/3600. * 50.0e6 #W
        n_ch4 = m_ch4 / 16e-3 # mol
        
        # h2 per hour required
        n_h2 = n_ch4 * 4.
        m_h2 = n_h2 * 2e-3
        
        # co2 per hour required
        n_co2 = n_ch4
        m_co2 = n_co2 * 44e-3
        
        self.par['n_dac'] = m_co2 * 8760. / 1e3
        
        p_pem_nom = self.par['e_pem'] * m_h2 * 1e3 #Wh

        p_dac_nom = self.par['e_dac'] * m_co2 * 1e3 #Wh electric
        
        
        p_dac_nom += self.par['q_dac'] / self.par['hp_cop'] * m_co2 * 1e3 #Wh heat

        self.par['n_hp'] = self.par['q_dac'] * m_co2 * 1e3
        
        e_comp_co2 = self.compressor_co2(1.) #Wh/kg_co2
        e_comp_co2 *= m_co2/m_ch4 #Wh/kg_ch4
        self.par['e_meth'] = e_comp_co2/1e3 

        e_comp_ch4_rec = self.compressor_ch4_recycle(1.) * 7. #Wh/kg_nh3 
        self.par['e_meth'] += e_comp_ch4_rec/1e3 
        
        self.par['n_comp_co2'] = e_comp_co2 * m_ch4  #W

        self.par['n_comp_ch4_rec'] = e_comp_ch4_rec * m_ch4 #W
        
        self.par['n_liq'] = m_ch4 #W

        p_meth_nom = self.par['e_meth'] * m_ch4 * 1e3 # m_nh3 * 1e3 #Wh # to be checked that it is in function of NH3
        
        #add liquefaction energy
        p_meth_nom += self.par['e_liq'] * m_ch4 * 1e3

        
        p_desal_nom = self.par['e_desal'] * m_h2 * self.par['h2o_to_h2_ratio_pem'] * 1e3 #Wh

        p_store_h2_nom = 0. #self.par['e_h2_compr'] * m_h2 * 1e3
        
        return p_pem_nom, p_dac_nom, p_desal_nom, p_meth_nom, p_store_h2_nom
    

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
            'p_in': 30.,
            'p_out': 200.,
            'eta_c': 0.85,
            'R': 4.124,
            'n': 1.609,
        }
        
        stages = 3.

        r = (par_c['p_out'] / par_c['p_in'])**(1./stages)
        
        power = (m_h2 *
                 par_c['n'] *
                 par_c['R'] *
                 par_c['T_in'] *
                 ((r)**((par_c['n'] -
                                     1.) /
                                    par_c['n']) -
                     1.) *
                 1000. /
                 (par_c['eta_c'] *
                     (par_c['n'] -
                      1.))) * stages

        return power

    def compressor_co2(self, m_in):
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
        m_in *= 1. / 3600.

        par_c = {
            'T_in': 373.,
            'p_in': 1.,
            'p_out': 10.,
            'eta_c': 0.85,
            #'R': 4.124,
            'n': 1.31,
        }
        
        R_h2o = 0.4615 #kJ/kgK
        R_co2 = 0.1889 #kJ/kgK

        M_co2 = 44e-3
        M_h2o = 18e-3
        n_co2 = 1.
        n_h2o = 4. #1.5
        m_h2o = M_h2o * n_h2o
        m_co2 = M_co2 * n_co2
        
        m_tot = m_co2 + m_h2o
        par_c['R'] = m_co2/m_tot * R_co2 + m_h2o/m_tot * R_h2o
        
        stages = 3.

        r = (par_c['p_out'] / par_c['p_in'])**(1./stages)
        
        power = (m_in *
                 par_c['n'] *
                 par_c['R'] *
                 par_c['T_in'] *
                 ((r)**((par_c['n'] -
                                     1.) /
                                    par_c['n']) -
                     1.) *
                 1000. /
                 (par_c['eta_c'] *
                     (par_c['n'] -
                      1.))) * stages

        return power

    def compressor_ch4_recycle(self, m_in):
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
        m_in *= 1. / 3600.

        par_c = {
            'T_in': 773.,
            'p_in': 9.5,
            'p_out': 10.,
            'eta_c': 0.85,
            #'R': 0.4882,
            'n': 1.4,
        }

        R_h2 = 4.124 #kJ/kgK
        R_h2o = 0.4615 #kJ/kgK
        R_co2 = 0.1889 #kJ/kgK
        R_ch4 = 0.5182 #kJ/kgK

        M_co2 = 44e-3
        M_h2o = 18e-3
        M_h2 = 2e-3
        M_ch4 = 16e-3
        m_co2 = 0.226
        m_h2o = 0.524 #1.5
        m_h2 = 0.041 #1.5
        m_ch4 = 0.209 #1.5
        
        m_tot = m_co2 + m_h2o + m_h2 + m_ch4
        par_c['R'] = m_co2/m_tot * R_co2 + m_h2o/m_tot * R_h2o + m_h2/m_tot * R_h2 + m_ch4/m_tot * R_ch4
        
        stages = 1.

        r = (par_c['p_out'] / par_c['p_in'])**(1./stages)
        
        power = (m_in *
                 par_c['n'] *
                 par_c['R'] *
                 par_c['T_in'] *
                 ((r)**((par_c['n'] -
                                     1.) /
                                    par_c['n']) -
                     1.) *
                 1000. /
                 (par_c['eta_c'] *
                     (par_c['n'] -
                      1.))) * stages

        return power


    def charge_pemel(self, p_pemel, p_pem_nom):

        # the operating bounds
        op_lower_lim = self.par['n_pemel'] * 0.
        op_upper_lim = self.par['n_pemel'] * 1.
        
        p_pemel *= self.par['e_pem']/(self.par['e_pem'] + self.par['e_h2_compr'])

        # check if power is higher than the lowest operating point
        # note that the nominal power is required to keep it running
        # and we are adding an extra excess power
        if p_pemel + p_pem_nom > op_lower_lim and self.m_h2[self.t-1] < self.m_h2_max:

            # check if the power exceeds the upper bound
            if p_pemel + p_pem_nom > op_upper_lim:
            
                # then we can only run the extra power allowed,
                # considering we are already running at nominal power
                # to provide the required h2 to keep things running
                p_pemel_applied = op_upper_lim - p_pem_nom

            else:
                # add some extra power to make hydrogen
                p_pemel_applied = p_pemel

            m_h2 = p_pemel_applied / ( self.par['e_pem'] * 1e3)

            # produced hydrogen is added to the tank
            self.m_h2[self.t] = self.m_h2[self.t-1] + m_h2

            if self.m_h2[self.t] > self.m_h2_max:
            
                self.m_h2[self.t] = self.m_h2_max

                # define power that results in a full storage tank
                m_h2_real =  self.m_h2_max - self.m_h2[self.t-1]
                p_pemel_applied = m_h2_real * self.par['e_pem'] * 1e3

            # increase the operating hours by 1
            self.res['running_hours_pemel'] += 1.
            
        # no hydrogen production when the power falls outside the operating
        # bounds
        else:
            p_pemel_applied = 0.
            self.m_h2[self.t] = self.m_h2[self.t-1]

            # increase the operating hours by 1 because we still run the p_pem_nom
            self.res['running_hours_pemel'] += 1.

        p_pemel_applied *= (self.par['e_pem'] + self.par['e_h2_compr']) / self.par['e_pem']
        
        return p_pemel_applied

    def charge_asu(self, p_asu, p_asu_nom):

        # the operating bounds
        op_lower_lim = self.par['n_asu'] * 0.
        op_upper_lim = self.par['n_asu'] * 1.

        if p_asu + p_asu_nom > op_lower_lim:

            # check if the power exceeds the upper bound
            if p_asu + p_asu_nom > op_upper_lim:
            
                # then we can only run the extra power allowed,
                # considering we are already running at nominal power
                # to provide the required N2 to keep things running
                p_asu_applied = op_upper_lim - p_asu_nom

            else:
                # add some extra power to make nitrogen
                p_asu_applied = p_asu

            m_n2 = p_asu_applied / (self.par['e_asu'] * 1e3)

            # produced nitrogen is added to the tank
            self.m_n2[self.t] = self.m_n2[self.t-1] + m_n2

            if self.m_n2[self.t] > self.m_n2_max:
            
                self.m_n2[self.t] = self.m_n2_max

                # define power that results in a full storage tank
                m_n2_real =  self.m_n2_max - self.m_n2[self.t-1]
                p_asu_applied = m_n2_real * self.par['e_asu'] * 1e3

            # increase the operating hours by 1
            self.res['running_hours_asu'] += 1.
            
        # no hydrogen production when the power falls outside the operating
        # bounds
        else:
            p_asu_applied = 0.
            self.m_n2[self.t] = self.m_n2[self.t-1]

            # increase the operating hours by 1 because we still run the p_asu_nom
            self.res['running_hours_asu'] += 1.

        return p_asu_applied

    def charge_bat(self, p_bat):
        
        p_supply_bat = self.par['n_bat'] * self.par['bat_frac_power_to_capacity_applied_one_hour'] # maximum energy that can be added in one hour
        
        p_bat_applied = min(p_supply_bat, p_bat)
        
        # add power to battery, but do not add more than possible
        self.e_bat[self.t] = min(self.e_bat[self.t-1] + p_bat_applied*self.par['eff_bat'],
                                 self.par['n_bat'])

    def discharge_bat(self, p_bat):
    
        p_extract_bat = self.par['n_bat'] * self.par['bat_frac_power_to_capacity_extracted_one_hour'] # maximum energy that can be discharged in one hour
    
        #p_bat_applied = min(p_extract_bat, p_bat)

        p_bat_applied = p_bat
        
        # add power to battery, but do not add more than possible
        self.e_bat[self.t] = self.e_bat[self.t-1] - p_bat_applied/self.par['eff_bat']
                                 
        if p_extract_bat < p_bat or self.e_bat[self.t] < self.par['n_bat'] * self.par['bat_DOD']: 
            missing = self.par['n_bat'] * self.par['bat_DOD'] - self.e_bat[self.t]
            self.e_bat[self.t] = self.par['n_bat'] * self.par['bat_DOD']
            
            return True, missing
        else:
            return False, 0.

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
        
        self.final_check = True 


        # get the nominal power domand
        p_pem_nom, p_dac_nom, p_desal_nom, p_meth_nom, p_store_h2_nom = self.p_demand_nominal()
        #print(p_pem_nom, p_dac_nom, p_desal_nom, p_meth_nom, p_store_h2_nom)
                
        self.m_h2 = np.zeros(len(self.sol_irr))
        self.e_bat = np.ones(len(self.sol_irr)) * self.par['n_bat']
        self.p_desal = np.zeros(len(self.sol_irr))
        
        checker = []
        
        self.res['p_dem_nom'] = p_pem_nom + p_dac_nom + p_desal_nom + p_meth_nom
        
        
        # when we want to store the excess power in different forms (power, h2 and n2),
        # we divide it based on the relative demand of the different components.
        # h2 and n2 are produced instantaneously and stored in a tank, we assume it is 
        # for free to extract the h2 and n2 from the tank when we need it in a next hour.
        # storing power in a battery is subject to losses when putting it in the battery
        # and when extracting it out of the battery. Thus the power stored to run equipment
        # should overcome the efficiency losses at the time of producing NH3 when there is 
        # no renewable power 
        p_division = p_pem_nom + p_store_h2_nom + p_desal_nom + p_meth_nom / self.par['eff_bat']**2.  + p_dac_nom  / self.par['eff_bat']**2. 
        p_pem_frac = (p_pem_nom + p_store_h2_nom)/p_division
        p_desal_frac = p_desal_nom/p_division
        p_meth_frac = p_meth_nom / (self.par['eff_bat']**2.) /p_division
        p_dac_frac = p_dac_nom / (self.par['eff_bat']**2.) /p_division

        # when there is not enough renewable power, this division checks how the power is
        # divided over the different components, to calculate the missing h2, n2 and power
        # here, we don't go through the battery, so this efficiency scaling is gone
        p_division_2 = p_pem_nom + p_dac_nom + p_desal_nom + p_meth_nom 
        p_pem_frac_2 = p_pem_nom/p_division_2
        p_desal_frac_2 = p_desal_nom/p_division_2
        p_meth_frac_2 = p_meth_nom /p_division_2
        p_dac_frac_2 = p_dac_nom /p_division_2
                
        # get the hourly photovoltaic array power
        self.photovoltaic()
        
        self.wind_turbine()
        
        self.res['p_renew'] = self.res['p_pv'] + self.res['p_wind']

        p_net = self.res['p_renew'] - self.res['p_dem_nom']

        # evaluate hourly hydrogen production
        
        pause_counter = 0
        n_pauses = 0.

        self.par['n_ch4_turb'] = 0.
        self.res['ch4_needed'] = 0.

        for t, p_ren in enumerate(self.res['p_renew']):
            self.t = t
        
            if t < self.par['delay_h'] or pause_counter > 0.1:
            
                pause_counter -= 1
                pause_counter = max(pause_counter,0)

                p_dem = p_ren - 0. #self.res['p_dem_nom']
            
            else:
        
                p_dem = p_ren - self.res['p_dem_nom']
        
            if p_dem > 0.:
                        
                # make some extra hydrogen
                p_pemel_applied = self.charge_pemel(p_dem * p_pem_frac, p_pem_nom)
                
                # store the amount of power supplied to the desalination
                # for sizing purposes
                if p_pemel_applied > 1e-3:
                    self.p_desal[self.t] = p_desal_nom + p_pemel_applied / (self.par['e_pem']*1e3) * self.par['h2o_to_h2_ratio_pem'] * self.par['e_desal'] * 1e3
                else:
                    self.p_desal[self.t] = p_desal_nom
                    
                # store remaining power in battery                
                self.charge_bat(p_dem * p_meth_frac + 
                                p_dem * p_dac_frac + 
                                p_dem * p_pem_frac - p_pemel_applied 
                                )
                                            
            else:
                
                # supply remaining H2
                p_pem_still_needed = p_pem_nom - self.res['p_renew'][self.t] * p_pem_frac_2
                m_h2_still_needed = p_pem_still_needed / ( self.par['e_pem'] * 1e3)

                if self.m_h2[self.t-1] - m_h2_still_needed > 0.:

                    self.m_h2[t] = self.m_h2[self.t-1] - m_h2_still_needed
                    p_pem_needed = 0.
                    p_desal_needed = 0.
                    self.p_desal[self.t] = 0.
                    
                    # so we covered the remaining H2 with what was left in the tank,
                    # to check if the PEM was used to cover part of the H2 demand,
                    # we check if there was some renewable power to give to the PEM
                    if self.res['p_renew'][self.t] > 1e-3:
                        self.res['running_hours_pemel'] += 1.
                    
                else:
                    
                    p_pem_needed = (m_h2_still_needed - self.m_h2[self.t-1]) * self.par['e_pem'] * 1e3
                    p_desal_needed = (m_h2_still_needed - self.m_h2[self.t-1]) * self.par['h2o_to_h2_ratio_pem'] * self.par['e_desal'] * 1e3
                    self.m_h2[self.t] = 0.
                    self.res['running_hours_pemel'] += 1.
                    self.p_desal[self.t] = p_desal_needed
                    
                    
                # supply remaining power for Haber Bosch
                p_meth_still_needed = p_meth_nom - self.res['p_renew'][self.t] * p_meth_frac_2
                p_dac_still_needed = p_dac_nom - self.res['p_renew'][self.t] * p_dac_frac_2
                
                fail, missing = self.discharge_bat(p_pem_needed + p_desal_needed + p_meth_still_needed + p_dac_still_needed)
                checker.append(fail)

                if fail:

                    ch4_needed = missing / self.par['eff_ch4_turb'] #Wh of h2 needed
                    n_turb = missing
                    if n_turb > self.par['n_ch4_turb']:
                        self.par['n_ch4_turb'] = n_turb
                        
                    self.res['ch4_needed'] += ch4_needed

            self.e_bat[self.t] -= (self.par['e_loss_day']/24. * self.par['n_bat'])
            self.e_bat[self.t] = max(self.e_bat[self.t], self.par['n_bat'] * self.par['bat_DOD'])

        self.par['demand_TWh'] *= 1. - n_pauses * self.par['delay_restart'] / 8760.
        
        if n_pauses * self.par['delay_restart'] > 0 or self.res['ch4_needed']/1e12*1.05 > self.par['demand_TWh']:
            self.final_check = False

        self.n_pauses = n_pauses
                        
        self.res['n_desal'] = max(self.p_desal)
                             
        # determine the electrolyzer lifetime
        self.lifetime()

        # determine the system cost and levelized cost of hydrogen
        self.cost()
        
        self.calc_lca()

    def battery_lifetime(self,SOCarray):
        import rainflow as rf

        cycles = rf.count_cycles(SOCarray)

        init = self.par['life_bat']
        y_ref = [1e6,2e5,1e5,7e4,3.5e4,1.5e4,1e4,8e3,6e3,5e3,4.5e3] #http://imistorage.blob.core.windows.net/imidocs/7420p009%20intensium%20flex.pdf
        x_ref = [0.035,0.1,0.135,0.2,0.3,0.4,0.42,0.5,0.6,0.7,0.8]
        
        z_ref = np.polyfit(x_ref,y_ref,8)
        Nmax = np.poly1d(z_ref)
        
        bat_aging = 0.
        for i in cycles:
            bat_aging += i[1]/Nmax(i[0])
        
        if bat_aging == 0.:
            bat_life = 1e8
        else:
            bat_life = 1./bat_aging
        
        return bat_life

    def lifetime(self):
        """

        The lifetime method determines the lifetime of
        the electrolyzer array, based on the number of
        operating hours during the evaluated year.

        """

        self.lifeBat = self.par['life_bat']


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
        components_cost = pv_cost

        wind_cost = self.par['n_wind'] * (crf * self.par['capex_wind'] +
                                      self.par['opex_wind'])
        components_cost += wind_cost

        bat_cost = self.par['n_bat'] * (crf * self.par['capex_bat'] +
                                      self.par['opex_bat'])
        components_cost += bat_cost
        
        # annual cost of electrolyzer array
        pemel_cost = self.par['n_pemel'] * (self.par['capex_pemel'] *
                                            (crf + self.par['opex_pemel']))
        components_cost += pemel_cost

        h2tank_cost = self.par['n_tank_h2'] * (self.par['capex_tank_h2'] *
                                            (crf + self.par['opex_tank_h2']))
        components_cost += h2tank_cost

        meth_cost = self.par['n_meth'] * (crf * self.par['capex_meth'] + self.par['opex_meth'] * self.par['capex_meth'])
        components_cost += meth_cost

        hp_cost = self.par['n_hp'] * (crf * self.par['capex_hp'] + self.par['opex_hp'])
        components_cost += hp_cost

        dac_cost = self.par['n_dac'] * (self.par['capex_dac'] * (crf + self.par['opex_dac']))
        components_cost += dac_cost

        desal_cost = self.res['n_desal'] / (self.par['e_desal'] * 1e3) * 8760. * 1e-3 * (crf * self.par['capex_desal'] +
                                      self.par['opex_desal']* self.par['capex_desal'])
        components_cost += desal_cost

        turb_cost = self.par['n_ch4_turb'] * (crf * self.par['capex_turb'] +
                                      self.par['opex_turb'] * self.par['capex_turb'])
        components_cost += turb_cost

        comp_cost_co2 = self.par['n_comp_co2'] * (self.par['capex_comp'] *
                                            (crf + self.par['opex_comp']))
        components_cost += comp_cost_co2

        liq_cost = self.par['n_liq'] * (self.par['capex_liq'] *
                                            (crf + self.par['opex_liq']))
        components_cost += liq_cost

        comp_cost_ch4 = self.par['n_comp_ch4_rec'] * (self.par['capex_comp'] *
                                            (crf + self.par['opex_comp']))
        components_cost += comp_cost_ch4

        speed = self.par['boat_knots'] * 1.852 #km/h
        time_trip = self.par['trip_length'] / speed #h

        times_op_en_af = 8760. / (time_trip * 2. + 48.) #including 24h loading and unloading time
        
        time_at_sea = time_trip * 2. / (time_trip * 2. + 48.)
        time_at_poarts = 48. / (time_trip * 2. + 48.)
                        
        boat_cap_twh = self.par['demand_TWh'] / times_op_en_af  #TWh
        
        boat_cap_m3 = boat_cap_twh * 1e6 * 3600. / 50.0 / 430. #m3

        boat_cap_tonne = boat_cap_twh * 1e6 * 3600. / 50.0 /1e3 #tonne

        self.n_ship = self.par['trip_length'] * boat_cap_tonne * times_op_en_af #vol
                      # self.par['trip_length'] * boat_cap_tonne * 0.05 * times_op_en_af #leeg
                      
        self.n_ship *= self.par['life_sys']              
        
        boat_cap_mwh = boat_cap_twh * 1e6
        cost_boat = boat_cap_mwh * (crf * self.par['capex_boat'] + self.par['opex_boat'])        
        
        ch4_produced = self.par['demand_TWh']*1e6
        ch4_produced -= self.res['ch4_needed']/1e6  
        
        self.m_ch4_for_backup = self.res['ch4_needed']/1e3 / (55./3.6) * self.par['life_sys']

        boat_cap_t_ref = 45000. #tonne
        e_cons_boat = self.par['e_boat'] #MJ/km
        
        fuel_needed = boat_cap_tonne / boat_cap_t_ref * e_cons_boat/ 3600. * self.par['trip_length'] * times_op_en_af #MWh fuel needed per year
        
        self.m_ch4_boiled_off_for_fuel = fuel_needed / (55./3600.) * self.par['life_sys'] #kg ch4 needed for boat in lifetime
        bog_fuel = (1. - (1.-self.par['bog'])**(time_trip/24.)) * ch4_produced #MWh boil-off per year. only 1 part of the trip, het naartoe gaan, is considered as fuel

        self.mwh_produced = ch4_produced - fuel_needed
        self.fuel_needed = fuel_needed
        
        if fuel_needed < bog_fuel:
            frac_bog_fuel = (bog_fuel - fuel_needed)/bog_fuel        

            refr_kg_h = boat_cap_tonne*1e3 * self.par['bog']/24. * frac_bog_fuel
            self.refr_kg_h = refr_kg_h

            refr_cost = refr_kg_h * (self.par['capex_liq'] * (crf + self.par['opex_liq']))
                
            refr_kg_year = boat_cap_tonne*1e3 * (1. - (1.-self.par['bog'])**(time_trip/24.)) * frac_bog_fuel * times_op_en_af 
                
            h2_needed_for_refr = refr_kg_year * self.par['e_liq'] # kWh power needed
            h2_needed_for_refr /= (55./3.6 * self.par['eff_ch4_turb']) #kWh power/ (kWh power per kg in H2) = kg H2 needed
            h2_needed_for_refr *= 55./3600. #MWh H2 needed
            
            self.mwh_produced -= h2_needed_for_refr
            self.h2_needed_for_refr = h2_needed_for_refr
            
        else:
            refr_cost = 0.
            self.refr_kg_h = 0.
            self.h2_needed_for_refr = 0.
        
        components_cost += refr_cost
        
               
        if self.smr:
            
            self.mwh_produced /= self.par['e_ratio_ng_h2']
            

        if self.res['ch4_needed']/1e12 > 0.10:
            self.final_check = False

        # the levelized cost of hydrogen
        cost = components_cost + cost_boat

        extra_backup = cost / (1. - self.res['ch4_needed']/1e12 - self.fuel_needed/1e6 - self.h2_needed_for_refr/1e6) - cost
        cost += extra_backup
        extra_smr = cost * self.par['e_ratio_ng_h2'] - cost
        cost += extra_smr
        
        self.res['lcoe'] = cost / 1e6
        self.res['lcoe_tonne'] = cost / (1e6 * 3600. / 120. /1e3)
        
        if not self.final_check:
            self.res['lcoe'] = 1e8
            self.par['n_pemel'] = 1e20
            self.par['n_bat'] = 1e20
            self.par['n_tank_h2'] = 1e12
            self.par['n_tank_n2'] = 1e12
        
        self.res['cost_frac_pv'] = pv_cost / cost
        self.res['cost_frac_wind'] = wind_cost / cost
        self.res['cost_frac_bat'] = bat_cost / cost
        self.res['cost_frac_hp'] = hp_cost / cost
        self.res['cost_frac_h2_tank'] = h2tank_cost / cost
        self.res['cost_frac_pem'] = (pemel_cost) / cost
        self.res['cost_frac_dac'] = dac_cost / cost
        self.res['cost_frac_meth'] = meth_cost / cost
        self.res['cost_frac_desal'] = desal_cost / cost
        self.res['cost_frac_refr'] = refr_cost / cost
        self.res['cost_frac_boat'] = cost_boat / cost
        self.res['cost_frac_comps'] = (comp_cost_co2 + comp_cost_ch4) / cost
        self.res['cost_frac_liq'] = liq_cost / cost
        self.res['cost_frac_smr'] = extra_smr / cost

    def calc_lca(self):



        n_bat = self.par['n_bat']/1e6 / (self.lifeBat/self.par['life_sys'])
        n_dac = self.par['n_dac']/8760. / (self.par['life_dac']/self.par['life_sys'])
        n_meth = self.par['n_meth'] /1e6 / (self.par['life_meth']/self.par['life_sys'])
        n_tank_h2 = self.par['n_tank_h2']/1e3 / (self.par['life_tank_h2']/self.par['life_sys'])
        n_pem = self.par['n_pemel']/1e6 / (self.par['life_pemel']/self.par['life_sys'])
        n_pv = self.par['n_pv']/1e6 / (self.par['life_pv']/self.par['life_sys'])
        n_wind = self.par['n_wind']/1e6 / (self.par['life_wind']/self.par['life_sys'])
        n_ship = self.n_ship / (self.par['life_boat']/self.par['life_sys'])
        n_desal = self.res['n_desal']/1e6 * 24. / self.par['e_desal'] / (self.par['life_desal']/self.par['life_sys'])
        n_hp = self.par['n_hp']/1e6 / (self.par['life_hp']/self.par['life_sys'])
        n_liq = self.par['n_liq'] / (self.par['life_liq']/self.par['life_sys'])
        n_refr = self.refr_kg_h / (self.par['life_liq']/self.par['life_sys'])
        
        m_ch4_produced = self.mwh_produced/(55./3600.) * self.par['life_sys'] 
                
        n_turb = self.par['n_ch4_turb']/1e6


        for index,method in enumerate(self.pb_names):

            self.res[method] = self.dict_lcia_pb['pblcia_battery_SSP1-Base_reference'][index] * n_bat 
            self.res[method] += self.dict_lcia_pb['pblcia_pem_SSP1-Base_reference'][index] * n_pem 
            self.res[method] += self.dict_lcia_pb['pblcia_pv_SSP1-Base_reference'][index] * n_pv
            self.res[method] += self.dict_lcia_pb['pblcia_wind_SSP1-Base_reference'][index] * n_wind 
            self.res[method] += self.dict_lcia_pb['pblcia_h2storage_SSP1-Base_reference'][index] * n_tank_h2 
            self.res[method] += self.dict_lcia_pb['pblcia_dac_SSP1-Base_reference'][index] * n_dac
            self.res[method] += self.dict_lcia_pb['pblcia_meth_SSP1-Base_reference'][index] * n_meth
            self.res[method] += self.dict_lcia_pb['pblcia_lng_ship_SSP1-Base_reference'][index] * n_ship 
            self.res[method] += self.dict_lcia_pb['pblcia_desal_SSP1-Base_reference'][index] * n_desal 
            self.res[method] += self.dict_lcia_pb['pblcia_hp_SSP1-Base_reference'][index] * n_hp 
            self.res[method] += self.dict_lcia_pb['pblcia_ch4liq_SSP1-Base_reference'][index] * n_liq 
            self.res[method] += self.dict_lcia_pb['pblcia_ch4liq_SSP1-Base_reference'][index] * n_refr 
            self.res[method] += self.dict_lcia_pb['pblcia_turb_SSP1-Base_reference'][index] * n_turb 
            self.res[method] += self.dict_lcia_pb['pblcia_ng_comb_real_SSP1-Base_reference'][index] * (self.m_ch4_boiled_off_for_fuel) # + m_ch4_produced * (self.par['e_ratio_ng_h2'] - 1.)/self.par['e_ratio_ng_h2']
            extra_backup = self.res[method] / (1. - self.res['ch4_needed']/1e12 - self.fuel_needed/1e6 - self.h2_needed_for_refr/1e6) - self.res[method]
            self.res[method] += extra_backup
            if self.smr:
                extra_smr = self.res[method] * self.par['e_ratio_ng_h2'] - self.res[method]
            else:
                extra_smr = 0.
            self.res[method] += extra_smr
            
            self.res['%s_frac_bat' %method] = self.dict_lcia_pb['pblcia_battery_SSP1-Base_reference'][index] * n_bat / self.res[method]
            self.res['%s_frac_pem' %method] = self.dict_lcia_pb['pblcia_pem_SSP1-Base_reference'][index] * n_pem / self.res[method]
            self.res['%s_frac_pv' %method] = self.dict_lcia_pb['pblcia_pv_SSP1-Base_reference'][index] * n_pv / self.res[method]
            self.res['%s_frac_wind' %method] = self.dict_lcia_pb['pblcia_wind_SSP1-Base_reference'][index] * n_wind / self.res[method]
            self.res['%s_frac_h2storage' %method] = self.dict_lcia_pb['pblcia_h2storage_SSP1-Base_reference'][index] * n_tank_h2 / self.res[method]
            self.res['%s_frac_dac' %method] = self.dict_lcia_pb['pblcia_dac_SSP1-Base_reference'][index] * n_dac / self.res[method]
            self.res['%s_frac_meth' %method] = self.dict_lcia_pb['pblcia_meth_SSP1-Base_reference'][index] * n_meth / self.res[method]
            self.res['%s_frac_lng_ship' %method] = self.dict_lcia_pb['pblcia_lng_ship_SSP1-Base_reference'][index] * n_ship / self.res[method]
            self.res['%s_frac_desal' %method] = self.dict_lcia_pb['pblcia_desal_SSP1-Base_reference'][index] * n_desal / self.res[method]
            self.res['%s_frac_hp' %method] = self.dict_lcia_pb['pblcia_hp_SSP1-Base_reference'][index] * n_hp / self.res[method]
            self.res['%s_frac_ch4liq' %method] = self.dict_lcia_pb['pblcia_ch4liq_SSP1-Base_reference'][index] * n_liq / self.res[method]
            self.res['%s_frac_refr' %method] = self.dict_lcia_pb['pblcia_ch4liq_SSP1-Base_reference'][index] * n_refr / self.res[method]
            self.res['%s_frac_turb' %method] = self.dict_lcia_pb['pblcia_turb_SSP1-Base_reference'][index] * n_turb / self.res[method]
            self.res['%s_frac_ng_comb_real' %method] = self.dict_lcia_pb['pblcia_ng_comb_real_SSP1-Base_reference'][index] * (self.m_ch4_boiled_off_for_fuel + self.m_ch4_for_backup) / self.res[method]
            self.res['%s_frac_co2_capture' %method] = self.dict_lcia_pb['pblcia_co2_capture_SSP1-Base_reference'][index] * m_ch4_produced * 2.75 / self.res[method]
            self.res['%s_frac_smr' %method] = extra_smr / self.res[method]  #* n_dac / self.res[method]
            
            self.res['%s_frac_co2' %method] = self.dict_lcia_pb['pblcia_co2_SSP1-Base_reference'][index] * m_ch4_produced * 2.75 / 175198850.69553944 / self.res[method]  #* n_dac / self.res[method]

            self.res[method] /=  1e9 #(self.mwh_produced * 1e3) # scaled to a 1 kWh/y plant
            
            self.res[method] /= self.par['life_sys'] #scaled to have annual contribution to PBs

            if self.scenario == 'europe_10mil':
                self.res[method] *= 330e9 # scaled to a 330 TWh = 10 million tonne / year  

            self.res['%s_%s' %(method, 'SOS_occupation')] = self.res[method] / self.pb_limits[index]
        
            if not self.final_check:
                self.res['lcoe'] = 1e8
                self.par['n_pemel'] = 1e20
                self.par['n_bat'] = 1e20
                self.par['n_tank_h2'] = 1e12
                self.res[method] = 1e20
                self.res['%s_%s' %(method, 'SOS_occupation')] = 1e20

        if self.all_dbs:

            for db in self.dbs:
                
                for index,method in enumerate(self.pb_names):

                    method += db.split('remind')[1]
                    
                    self.res[method] = self.dict_lcia_pb['pblcia_battery' + db.split('remind')[1]][index] * n_bat
                    self.res[method] += self.dict_lcia_pb['pblcia_pem' + db.split('remind')[1]][index] * n_pem
                    self.res[method] += self.dict_lcia_pb['pblcia_pv' + db.split('remind')[1]][index] * n_pv
                    self.res[method] += self.dict_lcia_pb['pblcia_wind' + db.split('remind')[1]][index] * n_wind
                    self.res[method] += self.dict_lcia_pb['pblcia_h2storage' + db.split('remind')[1]][index] * n_tank_h2
                    self.res[method] += self.dict_lcia_pb['pblcia_dac' + db.split('remind')[1]][index] * n_dac
                    self.res[method] += self.dict_lcia_pb['pblcia_hp' + db.split('remind')[1]][index] * n_hp 
                    self.res[method] += self.dict_lcia_pb['pblcia_meth' + db.split('remind')[1]][index] * n_meth
                    self.res[method] += self.dict_lcia_pb['pblcia_lng_ship' + db.split('remind')[1]][index] * n_ship 
                    self.res[method] += self.dict_lcia_pb['pblcia_desal' + db.split('remind')[1]][index] * n_desal 
                    self.res[method] += self.dict_lcia_pb['pblcia_ch4liq' + db.split('remind')[1]][index] * n_liq 
                    self.res[method] += self.dict_lcia_pb['pblcia_ch4liq' + db.split('remind')[1]][index] * n_refr 
                    self.res[method] += self.dict_lcia_pb['pblcia_turb' + db.split('remind')[1]][index] * n_turb 
                    self.res[method] += self.dict_lcia_pb['pblcia_ng_comb_real' + db.split('remind')[1]][index] * (self.m_ch4_boiled_off_for_fuel + m_ch4_produced * (self.par['e_ratio_ng_h2'] - 1.)/self.par['e_ratio_ng_h2'] + self.m_ch4_for_backup)
                    
                    self.res['%s_frac_bat' %method] = self.dict_lcia_pb['pblcia_battery' + db.split('remind')[1]][index] * n_bat / self.res[method]
                    self.res['%s_frac_pem' %method] = self.dict_lcia_pb['pblcia_pem' + db.split('remind')[1]][index] * n_pem / self.res[method]
                    self.res['%s_frac_pv' %method] = self.dict_lcia_pb['pblcia_pv' + db.split('remind')[1]][index] * n_pv / self.res[method]
                    self.res['%s_frac_wind' %method] = self.dict_lcia_pb['pblcia_wind' + db.split('remind')[1]][index] * n_wind / self.res[method]
                    self.res['%s_frac_h2storage' %method] = self.dict_lcia_pb['pblcia_h2storage' + db.split('remind')[1]][index] * n_tank_h2 / self.res[method]
                    self.res['%s_frac_dac' %method] = self.dict_lcia_pb['pblcia_dac' + db.split('remind')[1]][index] * n_dac / self.res[method]
                    self.res['%s_frac_meth' %method] = self.dict_lcia_pb['pblcia_meth' + db.split('remind')[1]][index] * n_meth / self.res[method]
                    self.res['%s_frac_lng_ship' %method] = self.dict_lcia_pb['pblcia_lng_ship' + db.split('remind')[1]][index] * n_ship / self.res[method]
                    self.res['%s_frac_desal' %method] = self.dict_lcia_pb['pblcia_desal' + db.split('remind')[1]][index] * n_desal / self.res[method]
                    self.res['%s_frac_hp' %method] = self.dict_lcia_pb['pblcia_hp' + db.split('remind')[1]][index] * n_hp / self.res[method]
                    self.res['%s_frac_ch4liq' %method] = self.dict_lcia_pb['pblcia_ch4liq' + db.split('remind')[1]][index] * n_liq / self.res[method]
                    self.res['%s_frac_refr' %method] = self.dict_lcia_pb['pblcia_ch4liq' + db.split('remind')[1]][index] * n_refr / self.res[method]
                    self.res['%s_frac_turb' %method] = self.dict_lcia_pb['pblcia_turb' + db.split('remind')[1]][index] * n_turb / self.res[method]


                    self.res[method] /= (self.mwh_produced * 1e3) # scaled to a 1 kWh/y plant
                    
                    self.res[method] /= self.par['life_sys'] #scaled to have annual contribution to PBs

                    if self.scenario == 'europe_10mil':
                        self.res[method] *= 330e9 # scaled to a 330 GWh = 10 million kg / year  

                    self.res['%s_%s' %(method, 'SOS_occupation')] = self.res[method] / self.pb_limits_pros[method]
                    if not self.final_check:
                        self.res['%s_%s' %(method, 'SOS_occupation')] = 1e20


            
    def print_results(self):
        """

        This method prints the levelized cost of hydrogen,
        the hydrogen production, the annual energy produced
        by the photovoltaic array and the energy consumed by
        the electrolyzer array.

        """

        print('outputs:')
        print('energy produced:'.ljust(50) + '%.2f TWh' % (self.mwh_produced/1e6) )
        print('backup CH4 consumed:'.ljust(50) + '%.3f MWh' % (self.res['ch4_needed']/1e6) )
        print('number of shutdowns:'.ljust(50) + '%i ' % (self.n_pauses) )
        print('TWh delivered:'.ljust(30) + '%.2f TWh' % (self.mwh_produced * self.par['life_sys']/1e6))
        print('LCOE:'.ljust(30) + '%.2f euro/MWh' % self.res['lcoe'])
        print('LCOE:'.ljust(30) + '%.2f euro/tonne' % self.res['lcoe_tonne'])
        print('cost due to PV:'.ljust(30) + '%.1f %%' % (self.res['cost_frac_pv']*1e2))
        print('cost due to WT:'.ljust(30) + '%.1f %%' % (self.res['cost_frac_wind']*1e2))
        print('cost due to PEM:'.ljust(30) + '%.1f %%' % (self.res['cost_frac_pem']*1e2))
        print('cost due to DAC:'.ljust(30) + '%.1f %%' % (self.res['cost_frac_dac']*1e2))
        print('cost due to HP:'.ljust(30) + '%.1f %%' % (self.res['cost_frac_hp']*1e2))
        print('cost due to METH:'.ljust(30) + '%.1f %%' % (self.res['cost_frac_meth']*1e2))
        print('cost due to BAT:'.ljust(30) + '%.1f %%' % (self.res['cost_frac_bat']*1e2))
        print('cost due to H2 tank:'.ljust(30) + '%.1f %%' % (self.res['cost_frac_h2_tank']*1e2))
        print('cost due to desalination:'.ljust(30) + '%.1f %%' % (self.res['cost_frac_desal']*1e2))
        print('cost due to compressors:'.ljust(30) + '%.1f %%' % (self.res['cost_frac_comps']*1e2))
        print('cost due to liquefaction:'.ljust(30) + '%.1f %%' % (self.res['cost_frac_liq']*1e2))
        print('cost due to boat:'.ljust(30) + '%.1f %%' % (self.res['cost_frac_boat']*1e2))
        print('cost due to refrigeration on boat:'.ljust(30) + '%.1f %%' % (self.res['cost_frac_refr']*1e2))
        print('cost due to SMR:'.ljust(30) + '%.1f %%' % (self.res['cost_frac_smr']*1e2))

        print('PBLCIA:')
        for method in self.pb_names[:1]:
            print(method.ljust(50) + '%.8f' % (self.res[method]))
        print('')
        print('PBLCIA safe operating space taken:')
        for method in self.pb_names[:1]:
            print(method.ljust(50) + '%.1f %%' % (self.res['%s_%s' %(method, 'SOS_occupation')]*1e2))
        print('')
        for method in self.pb_names[:1]:
            print(method + ' due to PV:'.ljust(50) + '%.1f %%' % (self.res['%s_frac_pv' %method]*1e2))
            print(method + ' due to WT:'.ljust(50) + '%.1f %%' % (self.res['%s_frac_wind' %method]*1e2))
            print(method + ' due to PEM:'.ljust(50) + '%.1f %%' % (self.res['%s_frac_pem' %method]*1e2))
            print(method + ' due to LIQ:'.ljust(50) + '%.1f %%' % (self.res['%s_frac_ch4liq' %method]*1e2))
            print(method + ' due to BAT:'.ljust(50) + '%.1f %%' % (self.res['%s_frac_bat' %method]*1e2))
            print(method + ' due to DAC:'.ljust(50) + '%.1f %%' % (self.res['%s_frac_dac' %method]*1e2))
            print(method + ' due to HP:'.ljust(50) + '%.1f %%' % (self.res['%s_frac_hp' %method]*1e2))
            print(method + ' due to METH:'.ljust(50) + '%.1f %%' % (self.res['%s_frac_meth' %method]*1e2))
            print(method + ' due to H2 Tank:'.ljust(50) + '%.1f %%' % (self.res['%s_frac_h2storage' %method]*1e2))
            print(method + ' due to desalination:'.ljust(50) + '%.1f %%' % (self.res['%s_frac_desal' %method]*1e2))
            print(method + ' due to boat:'.ljust(50) + '%.1f %%' % (self.res['%s_frac_lng_ship' %method]*1e2))
            print(method + ' due to refr on boat:'.ljust(50) + '%.1f %%' % (self.res['%s_frac_refr' %method]*1e2))
            print(method + ' due to smr:'.ljust(50) + '%.1f %%' % (self.res['%s_frac_smr' %method]*1e2))
            print(method + ' due to NG emissions:'.ljust(50) + '%.1f %%' % (self.res['%s_frac_ng_comb_real' %method]*1e2))
            print('')

        print('size wind:'.ljust(30) + '%.2f MW' % (self.par['n_wind']/1e6))
        print('size PV:'.ljust(30) + '%.2f MW' % (self.par['n_pv']/1e6))
        print('size PEM:'.ljust(30) + '%.2f MW' % (self.par['n_pemel']/1e6))
        print('size BAT:'.ljust(30) + '%.2f MWh' % (self.par['n_bat']/1e6))
        print('life BAT:'.ljust(30) + '%.2f year' % (self.lifeBat))
        print('size HP:'.ljust(30) + '%.2f MW' % (self.par['n_hp']/1e6))
        print('size DAC:'.ljust(30) + '%.2f t/h' % (self.par['n_dac']/8760.))
        print('size methanation:'.ljust(30) + '%.2f MW' % (self.par['n_meth']/1e6))
        print('size desalination:'.ljust(30) + '%.2f MW' % (self.res['n_desal']/1e6))
        print('size desalination:'.ljust(30) + '%.2f m3/d' % (self.res['n_desal']/1e6 * 24. / self.par['e_desal']))
        print('size compressors:'.ljust(30) + '%.2f MW' % ((self.par['n_comp_ch4_rec'] + self.par['n_comp_co2'])/1e6 ))
        print('size liquefaction:'.ljust(30) + '%.2f MW' % ((self.par['n_liq'] * self.par['e_liq'])/1e3 ))
        print('size liquefaction:'.ljust(30) + '%.2f kg/h' % ((self.par['n_liq']) ))
        print('size H2 tank:'.ljust(30) + '%.2f tonne (%.2f GWh)' % (self.par['n_tank_h2']/1e3, self.par['n_tank_h2']*33.3/1e6))
        print('number of ships:'.ljust(30) + '%.2f' % (self.n_ship))
        
        

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
