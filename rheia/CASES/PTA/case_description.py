"""
The :py:mod:`case_description` module contains a function
to read and store the fixed parameters for model evaluation and a function
to evaluate the system model.
"""

import os
import rheia.CASES.PTA.pta_backup as lb
import numpy as np
import pandas as pd

def set_params():
    """
    Set the fixed parameters for each model evaluation.

    """

    path = os.path.dirname(os.path.abspath(__file__))

    filename_climate = os.path.join(os.path.abspath(
                                    os.path.join(path,
                                                 os.pardir)),
                                    'DATA',
                                    'climate',
                                    'climate_oman_cf_pv.csv')

    # get the solar irradiance and ambient temperature data
    my_data = lb.ReadData(filename_climate)
    sol_irr, t_amb, wind_pow = my_data.load_climate()
    
    start = 'wind' #'solar' #'wind'

    if start is 'solar':
        print('start at max value solar')
        max_value = max(sol_irr)
        index = list(sol_irr).index(max_value) - 24

    else:
        print('start at max value wind')
        max_value = max(wind_pow)
        index = list(wind_pow).index(max_value) - 24
    
    bog = False

    csv_names = ['pblcia_battery_SSP1-Base_reference.csv',
                 'pblcia_hb_SSP1-Base_reference.csv',
                 'pblcia_h2storage_SSP1-Base_reference.csv',
                 'pblcia_pem_SSP1-Base_reference.csv',
                 'pblcia_pv_SSP1-Base_reference.csv',
                 'pblcia_wind_SSP1-Base_reference.csv',
                 'pblcia_ship_SSP1-Base_reference.csv',
                 'pblcia_desal_SSP1-Base_reference.csv',
                 'pblcia_nh3_bog_SSP1-Base_reference.csv',
                 'pblcia_hfo_refr_SSP1-Base_reference.csv',
                 'pblcia_asu_SSP1-Base_reference.csv',
                 'pblcia_compr_SSP1-Base_reference.csv',
                 'pblcia_turb_SSP1-Base_reference.csv',
                 'pblcia_nh3_comb_SSP1-Base_reference.csv',
                ]

    dict_lcia_pb = {}
    for index, name in enumerate(csv_names):

        filename_csv = os.path.join(os.path.abspath(
                                        os.path.join(path,
                                                     os.pardir)),
                                                    'PTA',
                                                    csv_names[index])

        with open(filename_csv, 'r') as file:
            lines = file.readlines()

            for i, line in enumerate(lines):
                line_split = line.split(",")
                x_int = [float(el) for el in line_split[:-1]]
    
        dict_lcia_pb[name[:-4]] = x_int

    scenario = 'europe_10mil'
    bog_flared = False
        
    # store data in params list
    params = [sol_irr, t_amb, wind_pow, bog_flared, False, {}, dict_lcia_pb, False, scenario]
    
    return params

def evaluate(x_in, params=[]):
    '''
    Evaluation of the system objectives for one given design.

    Parameters
    ----------
    x_in : tuple
        An enumerate object for the input sample.
        The first element of x
        - the index of the input sample in the list of samples -
        can be used for multiprocessing purposes of executable files
        with input and output text files.
        The second element of x - the input sample -
        is a dictionary with the names and values for the model parameters
        and design variables.
    params : list, optional
        List with fixed data, used during model evaluation. The default is [].

    Returns
    -------
    lcoh : float
        the levelized cost of hydrogen
    mh2: float
        hydrogen production
    '''

    arguments = params + [x_in[1]]

    # Evaluation object
    my_evaluation = lb.Evaluation(*arguments)

    # evaluate system model
    my_evaluation.evaluation()

    # get values for lcoh and mh2
    lcoe = my_evaluation.res['lcoe']
    
    pb_lcia_names = [
        'climate change CO2 concentration',
        'climate change energy imbalance',
        'stratospheric ozone depletion',
        'ocean acidification',
        'biogeochemical flows P',
        'biogeochemical flows N',
        'land system change global',
        'freshwater use global',
        'biosphere integrity',]

    res_pb_lcia = []
    
    for ij in pb_lcia_names:
        res_pb_lcia.append(my_evaluation.res['%s' %ij])
        
    pb_lcia_ccco2 = res_pb_lcia[0]
    pb_lcia_ccei = res_pb_lcia[1]
    pb_lcia_sod = res_pb_lcia[2]
    pb_lcia_oa = res_pb_lcia[3]
    pb_lcia_bfp = res_pb_lcia[4]
    pb_lcia_bfn = res_pb_lcia[5]
    pb_lcia_lscg = res_pb_lcia[6]
    pb_lcia_fug = res_pb_lcia[7]
    pb_lcia_bi = res_pb_lcia[8]
    
    return lcoe, pb_lcia_ccco2, pb_lcia_ccei, pb_lcia_oa, pb_lcia_bi, pb_lcia_bfn