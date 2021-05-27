import pytest
import numpy as np
import rheia.UQ.pce as pce
from rheia.CASES.determine_stoch_des_space import load_case

@pytest.fixture
def run_dict():
    dict_uq = {'case':                  'H2_FUEL',
               'pol order':             2,
               'objective names':       ['lcoh','mh2'],
               'objective of interest': 'lcoh',
               'sampling method': 'SOBOL',
               'n jobs': 1,
               'draw pdf cdf':          [True, 1e5],
               'results dir':           'opt_design_tutorial'
               }                
        
    return dict_uq

@pytest.fixture
def case(run_dict):

    design_space = 'design_space_tutorial_uq'
    space_obj, eval_func, params = load_case(
        run_dict, design_space, uq_bool=True)
    run_dict['evaluate'] = eval_func
    return space_obj, params

@pytest.fixture
def input_case_data(run_dict, case):

    data_obj = pce.Data(run_dict, case[0])
    return data_obj

@pytest.fixture
def input_RandomExperiment(input_case_data):
    input_case_data.read_stoch_parameters()
    input_case_data.create_samples_file()
    random_exp_obj = pce.RandomExperiment(input_case_data, 0)
    return random_exp_obj

def input_pce(input_RandomExperiment):
    pce_obj = pce.PCE(input_RandomExperiment)
    
    return pce_obj




def test_Data_inputs(run_dict, input_case_data):
    assert input_case_data.inputs == run_dict 

def test_Data_space_obj(input_case_data, case):
    assert input_case_data.space_obj == case[0] 

def test_read_stoch_parameters(input_case_data):
    input_case_data.read_stoch_parameters()
    input_case_data.stoch_data['l_b'] = [round(num, 3) for num in input_case_data.stoch_data['l_b']]
    input_case_data.stoch_data['u_b'] = [round(num, 3) for num in input_case_data.stoch_data['u_b']]
    input_case_data.stoch_data['mean'] = [round(num, 3) for num in input_case_data.stoch_data['mean']]
    input_case_data.stoch_data['deviation'] = [round(num, 3) for num in input_case_data.stoch_data['deviation']]
    assert input_case_data.stoch_data['names'] == ['u_sol_irr', 'u_t_amb', 'capex_pv', 'opex_pv', 'capex_pemel', 'opex_pemel', 'life_pemel', 'repl_pemel', 'capex_dcdc', 'opex_dcdc', 'int_rate', 'infl_rate']
    assert input_case_data.stoch_data['types'] == ['Uniform', 'Uniform', 'Uniform', 'Uniform', 'Uniform', 'Uniform', 'Uniform', 'Uniform', 'Uniform', 'Uniform', 'Uniform', 'Uniform']
    assert input_case_data.stoch_data['l_b'] == [0.9, -0.4, 350.0, 16.0, 1400.0, 0.03, 60000.0, 0.15, 100.0, 0.01, 0.04, 0.01]
    assert input_case_data.stoch_data['u_b'] == [1.1, 0.4, 600.0, 19.0, 2100.0, 0.05, 100000.0, 0.2, 200.0, 0.05, 0.08, 0.03]
    assert input_case_data.stoch_data['mean'] == [1.0, 0.0, 475.0, 17.5, 1750.0, 0.04, 80000.0, 0.175, 150.0, 0.03, 0.06, 0.02]
    assert input_case_data.stoch_data['deviation'] == [0.1, 0.4, 125.0, 1.5, 350.0, 0.01, 20000.0, 0.025, 50.0, 0.02, 0.02, 0.01]

def test_RandomExperiment_my_data(input_RandomExperiment, input_case_data):
    assert input_RandomExperiment.my_data == input_case_data

def test_RandomExperiment_obj_position(input_RandomExperiment):
    assert input_RandomExperiment.objective_position == 0

def test_RandomExperiment_dimension(input_RandomExperiment):
    assert input_RandomExperiment.dimension == 12

def test_RandomExperiment_n_terms(input_RandomExperiment):
    input_RandomExperiment.n_terms()
    assert input_RandomExperiment.n_samples == 182

def test_read_prev_samples(input_RandomExperiment):
    input_RandomExperiment.read_previous_samples(False)
    assert len(input_RandomExperiment.x_prev) == 182
    assert len(input_RandomExperiment.x_prev[0]) == 12
    assert len(input_RandomExperiment.y_prev) == 182
    assert input_RandomExperiment.y_prev[0][0] == 7.598744

def test_create_distr(input_RandomExperiment):
    input_RandomExperiment.create_distributions()
    assert input_RandomExperiment.polytypes == ['Legendre']*12

def test_create_samples(input_RandomExperiment):
    input_RandomExperiment.read_previous_samples(False)
    input_RandomExperiment.create_distributions()
    input_RandomExperiment.create_samples()
    assert input_RandomExperiment.size == 182
    assert len(input_RandomExperiment.x_u) == 182
    assert len(input_RandomExperiment.x_u_scaled) == 182

