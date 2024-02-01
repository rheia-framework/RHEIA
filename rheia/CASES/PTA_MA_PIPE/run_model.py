import os
import rheia

import rheia.CASES.PTH2.pth2_h2_backup as lb
import rheia.CASES.PTH2_MA.pth2_h2_backup as lb_ma
import rheia.CASES.PTH2_MA_PIPE.pth2_h2_backup_pipeline as lb_ma_pipe

import rheia.CASES.PTA.pta_backup as lb_pta
import rheia.CASES.PTA_MA.pta_backup as lb_pta_ma
import rheia.CASES.PTA_MA_PIPE.pta_backup_pipeline as lb_pta_ma_pipe

import rheia.CASES.PTCH4.ptch4_backup as lb_ptch4
import rheia.CASES.PTCH4_MA.ptch4_backup as lb_ptch4_ma
import rheia.CASES.PTCH4_MA_PIPE.ptch4_backup_pipeline as lb_ptch4_ma_pipe

import numpy as np
import pandas as pd

case = 'PTH2'

path = os.path.dirname(rheia.__file__)

# the climate file considered
filename_climate = os.path.join(path,
                                'CASES',
                                'DATA',
                                'climate',
                                'climate_oman_cf_pv.csv')

# the object to read in the data
my_data = lb.ReadData(filename_climate)

# get the solar irradiance and ambient temperature
sol_irr, t_amb, wind_pow = my_data.load_climate()

max_value = max(wind_pow)
index = list(wind_pow).index(max_value)  -24

# retrieve the deterministic values for the model parameters
parameters = my_data.load_parameters()

lengthh = 8760

path = os.path.dirname(os.path.abspath(__file__))

csv_names = ['pblcia_battery_SSP1-Base_reference.csv',
             'pblcia_h2liq_SSP1-Base_reference.csv',
             'pblcia_h2storage_SSP1-Base_reference.csv',
             'pblcia_pem_SSP1-Base_reference.csv',
             'pblcia_pv_SSP1-Base_reference.csv',
             'pblcia_wind_SSP1-Base_reference.csv',
             'pblcia_ship_SSP1-Base_reference.csv',
             'pblcia_desal_SSP1-Base_reference.csv',
             'pblcia_h2_bog_SSP1-Base_reference.csv',
             'pblcia_hfo_refr_SSP1-Base_reference.csv',
             'pblcia_ch4_pipeline_SSP1-Base_reference.csv',
             'pblcia_hb_SSP1-Base_reference.csv',
             'pblcia_nh3_bog_SSP1-Base_reference.csv',
             'pblcia_asu_SSP1-Base_reference.csv',
             'pblcia_compr_SSP1-Base_reference.csv',
             'pblcia_turb_SSP1-Base_reference.csv',
             'pblcia_nh3_comb_SSP1-Base_reference.csv',
             'pblcia_dac_SSP1-Base_reference.csv',
             'pblcia_meth_SSP1-Base_reference.csv',
             'pblcia_lng_ship_SSP1-Base_reference.csv',
             'pblcia_lng_fuel_SSP1-Base_reference.csv',
             'pblcia_hp_SSP1-Base_reference.csv',
             'pblcia_ch4liq_SSP1-Base_reference.csv',
             'pblcia_pcch4_SSP1-Base_reference.csv',
             'pblcia_co2_SSP1-Base_reference.csv',
             'pblcia_co2_capture_SSP1-Base_reference.csv',
             'pblcia_ng_comb_real_SSP1-Base_reference.csv',
            ]

dict_lcia = {}
dict_lcia_pb = {}
for index, name in enumerate(csv_names):

    filename_csv = os.path.join(os.path.abspath(
                                    os.path.join(path,
                                                 os.pardir)),
                                                case,
                                                csv_names[index])

    with open(filename_csv, 'r') as file:
        lines = file.readlines()

        for i, line in enumerate(lines):
            line_split = line.split(",")
            x_int = [float(el) for el in line_split[:-1]]

    dict_lcia_pb[name[:-4]] = x_int

    
scenario = 'europe_10mil'


# define the design to be tested
inputs = {'n_pv': 0.876, #GWp
          'n_wind': 0.3098, #GWp
          'n_pemel': 0.4085, #GW
          'n_bat': 0.606, #GWh
          'n_tank_h2': 0.174, #kt
          'start': 4021 #starting hour in the profile
          }

if 'PTH2' in case:

    if '_MA' in case:

        if '_PIPE' in case:
            # instantiate from the Evaluation class
            my_evaluation = lb_ma_pipe.Evaluation(sol_irr[:lengthh], t_amb[:lengthh], wind_pow[:lengthh], False, False, {}, dict_lcia_pb, False, scenario, {**parameters, **inputs})
        else:
            # instantiate from the Evaluation class
            my_evaluation = lb_ma.Evaluation(sol_irr[:lengthh], t_amb[:lengthh], wind_pow[:lengthh], False, False, {}, dict_lcia_pb, False, scenario, {**parameters, **inputs})
    else:
        # instantiate from the Evaluation class
        my_evaluation = lb.Evaluation(sol_irr[:lengthh], t_amb[:lengthh], wind_pow[:lengthh], False, False, {}, dict_lcia_pb, False, scenario, {**parameters, **inputs})

if 'PTA' in case:

    if '_MA' in case:

        if '_PIPE' in case:
            # instantiate from the Evaluation class
            my_evaluation = lb_pta_ma_pipe.Evaluation(sol_irr[:lengthh], t_amb[:lengthh], wind_pow[:lengthh], False, False, {}, dict_lcia_pb, False, scenario, {**parameters, **inputs})
        else:
            # instantiate from the Evaluation class
            my_evaluation = lb_pta_ma.Evaluation(sol_irr[:lengthh], t_amb[:lengthh], wind_pow[:lengthh], False, False, {}, dict_lcia_pb, False, scenario, {**parameters, **inputs})
    else:
        # instantiate from the Evaluation class
        my_evaluation = lb_pta.Evaluation(sol_irr[:lengthh], t_amb[:lengthh], wind_pow[:lengthh], False, False, {}, dict_lcia_pb, False, scenario, {**parameters, **inputs})

if 'PTCH4' in case:

    if '_MA' in case:

        if '_PIPE' in case:
            # instantiate from the Evaluation class
            my_evaluation = lb_ptch4_ma_pipe.Evaluation(sol_irr[:lengthh], t_amb[:lengthh], wind_pow[:lengthh], False, False, {}, dict_lcia_pb, False, scenario, {**parameters, **inputs})
        else:
            # instantiate from the Evaluation class
            my_evaluation = lb_ptch4_ma.Evaluation(sol_irr[:lengthh], t_amb[:lengthh], wind_pow[:lengthh], False, False, {}, dict_lcia_pb, False, scenario, {**parameters, **inputs})
    else:
        # instantiate from the Evaluation class
        my_evaluation = lb_ptch4.Evaluation(sol_irr[:lengthh], t_amb[:lengthh], wind_pow[:lengthh], False, False, {}, dict_lcia_pb, False, scenario, {**parameters, **inputs})
            

# evaluate the system
my_evaluation.evaluation()

# print the results
my_evaluation.print_results()
    