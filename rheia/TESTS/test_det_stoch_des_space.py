import os, sys
import shutil
import pytest

from rheia.CASES.determine_stoch_des_space import StochasticDesignSpace, load_case
    
@pytest.fixture
def run_dict():
    dict_opt = {'case':                'H2_FUEL',
                'objectives':          {'DET': (-1, 1)},
                'stop':                ('BUDGET', 2000),
                'n jobs':              1,
                'population size':     20,
                'results dir':         'run_tutorial',
                }
                
    
    
    return dict_opt

@pytest.fixture
def input_case(run_dict):

    design_space = 'design_space'
    opt_type = list(run_dict['objectives'].keys())[0]
    case = run_dict['case']
    test_class = StochasticDesignSpace(opt_type, case, design_space)
    
    return test_class

def test_case_name(input_case):
    assert input_case.case == 'H2_FUEL'

def test_design_space_name(input_case):       
    assert input_case.design_space == 'design_space'

def test_opt_type(input_case):       
    assert input_case.opt_type == 'DET'

def test_n_dim(input_case):       
    assert input_case.n_dim == 2

def test_n_par(input_case):       
    assert input_case.n_par == 13

def test_par_dict(input_case):       
    assert input_case.par_dict == {'n_pv': 5.,
                                   'u_sol_irr': 1,
                                   'u_t_amb': 0,
                                   'capex_pv': 475,
                                   'opex_pv': 17.5,
                                   'capex_pemel': 1750,
                                   'opex_pemel': 0.04,
                                   'life_pemel': 80000,
                                   'repl_pemel': 0.175,
                                   'capex_dcdc': 150,
                                   'opex_dcdc': 0.03,
                                   'int_rate': 0.06,
                                   'infl_rate': 0.02}
                                   
def test_var_dict(input_case):       
    assert input_case.var_dict == {'n_dcdc_pv': [0., 7.], 
                                   'n_pemel': [0., 7.]}

def test_lb(input_case):       
    assert input_case.l_b == [0., 0.]

def test_ub(input_case):       
    assert input_case.u_b == [7., 7.]

def test_upar_dict(input_case):       
    assert input_case.upar_dict == {}
    
def test_attach_objectives(run_dict, input_case):
    input_case.attach_objectives(run_dict['objectives']['DET'])
    assert input_case.obj == (-1, 1)
   
                
    

    
