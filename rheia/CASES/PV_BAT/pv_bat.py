'''
@author: Diederik Coppitters

@date: 14 Sep 2020

Short description ::
-----------------

Library for the techno-economic evaluation of 
a photovoltaic-battery system coupled to an electric demand

required packages:
- python2.7
- numpy
- scipy
- matplotlib
- rainflow

'''

import os
import numpy as np
import scipy as sp
from scipy import optimize
import pandas as pd
import matplotlib.pyplot as plt
import rainflow as rf


class ReadData:
    '''

    This class includes functions which enable to read data from the different .csv files
    
    '''

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
    '''

    This class includes:
    - component models
    - hourly evaluation of components behavior
    - quantification of LCOE and SSR

    '''

    def __init__(self,sol_irr,temp_amb,load_elec,sell_grid,parameters):
        self.sol_irr = sol_irr
        self.temp_amb = temp_amb
        self.sell_grid = sell_grid
        self.parameters = parameters
        self.load_elec = load_elec*self.parameters['uq_load_e']

        self.length = len(sol_irr)
        self.grid_electricity_array = np.zeros(self.length)
        self.sold_electricity_array = np.zeros(self.length)
        self.dcac_capacity_array = np.zeros(self.length)
        self.soc_array = np.zeros(self.length)
        
        self.dcdc_bat_capacity = [0.]
        self.dcdc_pv_capacity = 0.
        self.grid_cost = 0.
        self.grid_sold = 0.
        self.LIFETIME_SYSTEM = 20. # [year]
        #self.n_pv = self.parameters['n_pv'] / 0.240 # [# photovoltaic panels in array]
        self.n_pv = self.parameters['n_pv'] *1e3 # [# photovoltaic panels in array]

        # battery parameters
        self.C_NOM = 110. #nominal capacity [Ah]
        self.BAT_VOLTAGE = 2.*self.parameters['n_bat'] / (self.C_NOM*2.)*1e3
        self.SOC_MIN = 0.2
        self.SOC_MAX = 0.95
        self.soc = 0.0*(self.SOC_MAX - self.SOC_MIN) + self.SOC_MIN

    def elec_profiles(self):
        """
        Set the grid electricity price for buying and selling electricity. 
            
        Attributes:
            elec_profile (array): grid electricity buying price array 
            elec_profile_sale (array): grid electricity selling price array 
        """

        self.elec_profile = np.ones(self.length) * ((self.parameters['elec_cost']+20.)
                            * (1./0.3)) / 1e6
        self.elec_profile_sale = np.ones(self.length) * self.parameters['elec_cost'] / 1e6 

    def photovoltaic(self):
        """
        Determine the hourly power produced by the photovoltaic array  
            
        Attributes:
            dcdc_pv_capacity (float): the capacity of the DC-DC converter connected to the photovoltaic array

        Returns:
            pv_power (array): hourly photovoltaic array power [W]

        """
        
        '''
        pv_power = np.zeros(self.length)

        # get the specific photovoltaic panel characteristics
        pv_database = pvlib.pvsystem.retrieve_sam('CECmod')
        pv_system = pv_database.SunPower_SPR_X19_240_BLK
        #pv_system = pv_database.Sunpower_SPR_X19_240_BLK
        
        
        for i in range(self.length):
            if self.sol_irr[i] > 0.:
                pv_inputs = pvlib.pvsystem.calcparams_desoto(self.sol_irr[i]*self.parameters['uq_g'],
                                                            self.temp_amb[i]+self.parameters['uq_t'],
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
                pv_outputs = pvlib.pvsystem.max_power_point(pv_inputs[0],
                                                            pv_inputs[1],
                                                            pv_inputs[2],
                                                            pv_inputs[3],
                                                            pv_inputs[4],
                                                            method='newton')
                pv_power[i] = min( (1.+self.parameters['power_tol_pv']/100.) * pv_outputs['p_mp'] * self.n_pv, self.parameters['pv_dcdc']*1e3 )
            else:
                pv_power[i] = 0.
        
        '''
        pv_power = np.zeros(self.length)
        for i in range(self.length):
            if self.sol_irr[i] > 0.:
                pv_power[i] = min((1.+self.parameters['power_tol_pv']/100.) * self.sol_irr[i]/1e3*0.18 * self.n_pv/0.25, self.parameters['pv_dcdc']*1e3 )

        #self.dcdc_pv_capacity = max(pv_power)
        self.dcdc_pv_capacity = self.parameters['pv_dcdc']*1e3

        return pv_power

    def net_power(self):
        """
        Determine the hourly net power.
        This corresponds to the net power available (or still required) after extracting
        the hourly load from the photovoltaic power, considering DC-DC converter and
        DC-AC inverter efficiency.
        
        Attributes:
            pv_power (array): hourly photovoltaic power [W]

        Returns:
            p_net (array): the hourly net power [W]

        """

        p_net = np.zeros(self.length)
        self.pv_power = self.photovoltaic()

        for i in range(self.length):
            p_req = self.load_elec[i] / (self.parameters['eta_dcdc']*self.parameters['eta_dcac'])
            if self.pv_power[i] >= p_req:
                p_net[i] = (self.pv_power[i]-p_req) * self.parameters['eta_dcdc']
            else:
                p_net[i] = (self.pv_power[i]*self.parameters['eta_dcdc']*self.parameters['eta_dcac']
                            - self.load_elec[i]) / self.parameters['eta_dcac']
        #print(self.pv_power,p_net)
        return p_net

    ##########################################
    ## battery module
    ##########################################

    def battery(self,current,soc_init):
        """
        determine the energy stored/extracted from the battery stack
        and the new state of charge (SOC)
        
        Arguments:
            current (float): current provided/extracted [A]
            soc_init (float): initial state of charge, before applying the current 

        Returns:
            e_bat (float): energy stored/extracted in the battery [Wh]
            soc (float): updated state of charge

        """

        soc = soc_init
        EFF_BAT = 0.8
        
        if current > 0:
            soc += current*EFF_BAT/self.C_NOM
        else:
            soc += current/EFF_BAT/self.C_NOM

        e_bat = abs(soc_init-soc)*self.C_NOM*self.BAT_VOLTAGE

        return e_bat,soc

    def battery_lifetime(self,soc_array):
        """
        determine the battery lifetime
        
        Arguments:
            soc_array (array): the hourly state of charge over the year

        Returns:
            bat_life (float): battery lifetime [year]
            
        """

        # use rainflow method to determine hourly changes in state of charge
        cycles = rf.count_cycles(soc_array)
        # create polynomial which approximates the lifetime (in cycles) 
        # in function of the hourly depth of discharge
        init = self.parameters['life_bat']
        x_ref = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
        y_ref = np.array([init+6750.,init+2750.,init+1750.,init+1250.,
                        init+750.,init+450.,init+225.,init])
        poly_coeff = np.polyfit(x_ref,y_ref,5)
        f_lifetime = np.poly1d(poly_coeff)

        bat_aging = 0.
        for i in cycles:
            bat_aging += i[1]/f_lifetime(i[0])
            #print(bat_aging)

        #bat_aging *= 365.
        if bat_aging == 0.:
            bat_life = 1e8
        else:
            bat_life = 1./bat_aging

        return bat_life

    def f_battery_full(self,inp,*arguments):
        """
        Iterable function used to find the current which results in 
        a state of charge equal to the maximum state of charge.
        
        Arguments:
            inp (float): current [A]

        Returns:
            soc_res - self.SOC_MAX: difference between actual and maximum state of charge
            
        """

        soc_init = arguments
        soc_res = self.battery(inp,soc_init)[1]

        return soc_res - self.SOC_MAX

    def f_battery_empty(self,inp,*arguments):
        """
        Iterable function used to find the current which results in 
        a state of charge equal to the minimum state of charge.
        
        Arguments:
            inp (float): current [A]

        Returns:
            soc_res - self.SOC_MIN: difference between actual and minimum state of charge
            
        """

        soc_init = arguments
        soc_res = self.battery(inp,soc_init)[1]

        return soc_res - self.SOC_MIN

    def f_battery_current(self,inp,*arguments):
        """
        Iterable function used to find the current which results in 
        a battery power equal to the applied power.
        
        Arguments:
            inp (float): current [A]

        Returns:
            p_res - p_net: difference between actual and applied power
            
        """

        soc_init,p_net = arguments
        
        if abs(inp) > self.C_NOM/10.*3.:
            return 1e8
        
        p_res = self.battery(inp,soc_init)[0]
        
        return p_res - p_net

    def discharge_battery(self,load_to_cover):
        """
        Discharging the battery stack.
        
        Arguments:
            load_to_cover (float): the energy required to cover the load [Wh]

        Attributes:
            dcdc_bat_capacity (list): the hourly DC-DC converter capacity required to operate the battery
            soc (float): the actual battery state of charge
            
        Returns:
            p_bat: power discharged from the battery [W]
            
        """

        soc_init = self.soc
        soc = round(self.soc,6)

        if soc > self.SOC_MIN and self.parameters['n_bat'] > 1e-3:
            #print(load_to_cover)
            load_to_cover = min( load_to_cover, self.parameters['bat_dcdc']*1e3 ) * 1./self.parameters['eta_dcdc']
            #print(load_to_cover,'corr')
            current_init = -load_to_cover/self.BAT_VOLTAGE
            
            # Determine the required current extracted from the battery
            if abs(current_init) < self.C_NOM/10.*3.:
                arguments = soc_init,load_to_cover
                current_res = sp.optimize.fsolve(self.f_battery_current,
                                                                    current_init,
                                                                    args = arguments,
                                                                    xtol=1e-4,
                                                                    full_output = True)[0]
                current = current_res[0]
            else:
                current = current_init

            # operate at nominal current, if required current exceeds maximum discharge current
            if abs(current) > 3.*self.C_NOM/10.:
                current = -self.C_NOM/10.

            power,soc = self.battery(current,soc_init)
            p_bat = power

            # determine actual discharge power, 
            # when required discharge power exceeds minimum state of charge
            if soc < self.SOC_MIN:
                arguments = soc_init
                current_real = sp.optimize.fsolve(self.f_battery_empty,
                                                                    -0.1,
                                                                    args=arguments,
                                                                    xtol=1e-6,
                                                                    full_output = True)[0]
                p_real,soc = self.battery(current_real[0],soc_init)
                p_bat = abs(p_real)

            # determine required DC-DC converter capacity
            self.dcdc_bat_capacity.append(p_bat)
        else:
            p_bat = 0.
            current = 0.

        self.soc = soc

        return p_bat

    def charge_battery(self,e_left):
        """
        Discharging the battery stack.
        
        Arguments:
            e_left (float): the energy available to charge the battery [Wh]

        Attributes:
            dcdc_bat_capacity (list): the hourly DC-DC converter capacity required to operate the battery
            soc (float): the actual battery state of charge
            
        Returns:
            p_bat: power charged to the battery [W]
            
        """

        soc_init = self.soc
        soc = self.soc

        if soc < self.SOC_MAX and self.parameters['n_bat'] > 1e-3:
            e_left = min(e_left, self.parameters['bat_dcdc']*1e3) * self.parameters['eta_dcdc']
            current = min(e_left/self.BAT_VOLTAGE,self.C_NOM/10.)
            soc = self.battery(current,soc_init)[1]
            p_bat = current*self.BAT_VOLTAGE

            # determine actual charge power, 
            # when required charge power exceeds maximum state of charge
            if soc > self.SOC_MAX:
                current_real = sp.optimize.fsolve(self.f_battery_full,
                                                                    1.,
                                                                    args=soc_init,
                                                                    xtol=1e-4,
                                                                    full_output = True)[0]
                soc = self.battery(current_real[0],soc_init)[1]
                p_bat = current_real[0]*self.BAT_VOLTAGE

            # determine required DC-DC converter capacity
            self.dcdc_bat_capacity.append(p_bat*self.parameters['eta_dcdc'])
        else:
            p_bat = 0.

        self.soc = soc

        return p_bat

    ##########################################
    ## model evaluation
    ##########################################

    def evaluation(self):
        """
        Hourly evaluation of the system.
        
        Attributes:
            e_grid (float): hourly energy bought from the grid
            e_grid_sold (float): hourly energy sold to the grid
            dcac_capacity_array (array): array of hourly capacity required of the DC-AC inverter
            grid_cost (float): annual cost of buying electricity from the grid [euro]
            grid_electricity_array (array): array of hourly energy bought from the grid
            sold_electricity_array (array): array of hourly energy sold to the grid
            soc (float): state of charge 
            soc_array (array): array of hourly state of charge over the year 
            life_bat (float): battery stack lifetime [year]
            ssr (float): self-sufficiency ratio
            
        """

        # generate the grid electricity price profiles
        self.elec_profiles()
        
        # acquire the hourly net power available/required
        p_net = self.net_power()
        
        # hourly evaluation of the system
        for i in range(self.length):
            self.e_grid = 0.
            e_grid_sold = 0.

            if p_net[i] > 0.:
                p_excess = p_net[i]
                #print(p_excess)
                # determine power charged to the battery
                p_bat = self.charge_battery(p_excess)

                # sell remaining energy to the grid
                p_excess = p_excess - p_bat/self.parameters['eta_dcdc']
                e_grid_sold = p_excess*self.parameters['eta_dcac']
                self.grid_sold += e_grid_sold * self.elec_profile_sale[i]
                self.dcac_capacity_array[i] += self.load_elec[i] + e_grid_sold
            else:
                p_rem = abs(p_net[i])

                # determine power discharged from the battery
                p_bat = self.discharge_battery(p_rem)
                p_bat *= self.parameters['eta_dcdc']

                # buy remaining energy from the grid
                p_rem = p_rem - p_bat
                self.e_grid += p_rem*self.parameters['eta_dcac']
                self.dcac_capacity_array[i] += (self.load_elec[i] - self.e_grid)
                self.grid_cost += self.e_grid * self.elec_profile[i]

            self.grid_electricity_array[i] = self.e_grid
            self.sold_electricity_array[i] = e_grid_sold
            self.soc -= self.parameters['self_discharge']/24.
            self.soc = max(self.soc,self.SOC_MIN)
            self.soc_array[i] = self.soc
        self.soc_array = np.tile(self.soc_array,365)
        # determine battery lifetime
        self.life_bat = self.battery_lifetime(self.soc_array)

        # determine self-sufficiency ratio
        self.ssr = round((1. - sum(self.grid_electricity_array)/(sum(self.load_elec))) *100.,2)

        # determine the system cost
        self.cost()

    def cost(self):
        """
        Determination of the system cost and the levelized cost of electricity.
        
        Attributes:
            lcoe (float): levelized cost of electricity
            
        """

        inv_rate = (self.parameters['disc_rate']-0.02)/(1. + 0.02)
        crf = (((1.+inv_rate)**self.LIFETIME_SYSTEM-1.)
              / (inv_rate*(1.+inv_rate)**self.LIFETIME_SYSTEM))**(-1)
        
        # cost related to photovoltaic array
        p_pv = self.parameters['n_pv']
        #p_dcdcmppt = self.dcdc_pv_capacity
        p_dcdcmppt = self.parameters['pv_dcdc']*1e3
        pv_cost = p_pv * (crf*self.parameters['capex_pv'] + self.parameters['opex_pv'])
        pv_dcdc_cost = p_dcdcmppt * (crf*self.parameters['capex_dcdc']
                        + self.parameters['opex_dcdc']*self.parameters['capex_dcdc'])
        components_cost = pv_cost + pv_dcdc_cost
        
        # cost related to battery stack
        e_bat = self.parameters['n_bat']
        p_dcdc_bat = self.parameters['bat_dcdc']*1e3
        bat_cost = e_bat * (crf*self.parameters['capex_bat'] + self.parameters['opex_bat'])
        bat_cost_repl = crf * sum([(1.+inv_rate)**(-(i+1.)*self.life_bat)*e_bat*self.parameters['repl_bat']
                        for i in range(int(self.LIFETIME_SYSTEM/self.life_bat))])
        bat_dcdc_cost = p_dcdc_bat * (crf*self.parameters['capex_dcdc_bat']
                        + self.parameters['opex_dcdc']*self.parameters['capex_dcdc_bat'])
        components_cost += bat_cost + bat_dcdc_cost
        arc = bat_cost_repl # annual replacement cost
        
        # cost related to DC-AC inverter
        p_dcac = max(self.dcac_capacity_array)
        dcac_cost = p_dcac * (crf*self.parameters['capex_dcac']
                    + self.parameters['opex_dcac']*self.parameters['capex_dcac'])
        components_cost += dcac_cost
        
        # if selling to the grid is not considered,
        # set the gain from selling electricity equal to 0.
        if not self.sell_grid:
            self.grid_sold = 0.
        
        cost = arc + components_cost + self.grid_cost*365. - self.grid_sold*365.
        self.lcoe = cost/(sum(self.load_elec)*365.)*1e6
        #cost = arc + components_cost + self.grid_cost - self.grid_sold
        #self.lcoe = cost/(sum(self.load_elec))*1e6

    ##########################################
    ## generate results
    ##########################################

    def print_results(self):
        """
        Print the results
           
        """

        print('outputs:')
        print('LCOE:'.ljust(30) + str(round(self.lcoe,2)) + ' euro/MWh')
        print('SSR:'.ljust(30) + str(round(self.ssr,2)) + " %")
        print('PV electricity generated:'.ljust(30) + str(round(sum(self.pv_power)*365./1e6,2)) + ' MWh')
        print('grid electricity:'.ljust(30) + str(round(sum(self.grid_electricity_array)*365./1e6,2)) + ' MWh')
        print('electricity sold:'.ljust(30) + str(round(sum(self.sold_electricity_array)*365./1e6,2)) + ' MWh')
        print("life battery:".ljust(30) + str(round(self.life_bat,2)) + ' year')
        plt.plot((self.soc_array-self.SOC_MIN)/(self.SOC_MAX-self.SOC_MIN),'-g')
        plt.xlabel("time [s]")
        plt.ylabel("SOC")
        plt.show()

    def get_objectives(self):
        """
        Return the main objectives of the system
        
        Returns:
            lcoe (float): levelized cost of electricity [euro/MWh]
            ssr (float): self-sufficiency ratio            
           
        """

        return self.lcoe,self.ssr
