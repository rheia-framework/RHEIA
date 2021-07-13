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

@pytest.fixture
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

def test_multindices(input_pce):
    assert input_pce.multindices(range(0,2)) == [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                                 (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)]

def test_ols(input_pce):
    a = np.array([[ 1.,-0.96263],
                  [ 1.,-0.94701],
                  [ 1.,0.0529],
                  [ 1.,0.5529],
                  ])
    b = np.array([[7.598744],
                  [7.586986],
                  [7.298572],
                  [6.390994]])

    test_res = input_pce.ols(a,b)[0][0]
    assert round(test_res,3) == 6.996
    
def test_calc_a(input_RandomExperiment, input_pce):
    input_RandomExperiment.read_previous_samples(False)
    input_RandomExperiment.create_distributions()
    input_RandomExperiment.create_samples()
    multindices = input_pce.multindices(range(0,91))
    test_res = input_pce.calc_a(multindices)
    assert round(test_res[0][1],3) == -0.963
    assert round(test_res[3][12],3) == -0.533
    
def test_run(input_RandomExperiment, input_pce):
    input_RandomExperiment.n_terms()
    input_RandomExperiment.read_previous_samples(False)
    input_RandomExperiment.create_distributions()
    input_RandomExperiment.create_samples()
    input_RandomExperiment.y = input_RandomExperiment.y_prev
    input_pce.run()
    coeff = input_pce.coefficients
    assert round(coeff[0][0],3) == 7.730
    assert round(coeff[3][0],3) == 0.559
    assert round(input_pce.moments['mean'],3) == 7.730 
    assert round(input_pce.moments['variance'],3) == 0.742

def test_calc_sobol(input_RandomExperiment, input_pce):
    input_RandomExperiment.n_terms()
    input_RandomExperiment.read_previous_samples(False)
    input_RandomExperiment.create_distributions()
    input_RandomExperiment.create_samples()
    input_RandomExperiment.y = input_RandomExperiment.y_prev
    input_pce.run()
    input_pce.calc_sobol()
    assert round(input_pce.sensitivity['s_i'][0],3) == 0.059
    assert round(input_pce.sensitivity['s_i'][6],3) == 0.012    
    assert round(input_pce.sensitivity['s_tot_i'][2],3) == 0.142
    assert round(input_pce.sensitivity['s_tot_i'][9],3) == 0.002    

def test_calc_loo(input_RandomExperiment, input_pce):
    input_RandomExperiment.n_terms()
    input_RandomExperiment.read_previous_samples(False)
    input_RandomExperiment.create_distributions()
    input_RandomExperiment.create_samples()
    input_RandomExperiment.y = input_RandomExperiment.y_prev
    input_pce.run()
    input_pce.calc_loo()
    assert round(input_pce.loo,5) == 0.00595
