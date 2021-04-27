import os, sys
import shutil
import pytest

path = os.path.split(
    os.path.dirname(
        os.path.abspath(__file__)))[0]
sys.path.insert(0, os.path.join(path,'CASES'))
from determine_stoch_des_space import stochastic_design_space, load_case

@pytest.fixture
def define_case():

    case = 'CASE1'
    path = os.path.split(
        os.path.dirname(
            os.path.abspath(__file__)))[0]

    path_case = os.path.join(path,'CASES',case)
    if not os.path.exists(path_case):
        os.makedirs(path_case)
        with open(os.path.join(path_case, 'design_space'), "w") as f:
            f.write("var_1 var 0 10 \n")
            f.write("par_1 par 1")

        with open(os.path.join(path_case, "stochastic_space"), "w") as f:
            f.write("var_1 absolute Uniform 1 \n")
            f.write("par_1 relative Gaussian 0.2")
            
        shutil.copy(os.path.join(path,'CASES','REF','case_description.py'),path_case)
    
@pytest.fixture
def input_case(define_case):

    design_space = 'design_space'
    opt_type = 'ROB'
    case = 'CASE1'
    test_class = stochastic_design_space(opt_type, case, design_space)

    path = os.path.split(
        os.path.dirname(
            os.path.abspath(__file__)))[0]

    path_case = os.path.join(path,'CASES',case)
    shutil.rmtree(path_case)
    
    return test_class

def test_load_case_uq(define_case):
    run_dict = {'case': 'CASE1'}
    test_space_obj, test_eval_func, test_params = load_case(run_dict, 'design_space', uq = True)
    assert test_space_obj.case == 'CASE1'
    #assert test_obj.lb = []
    #assert test_obj.ub = []
    

def test_load_case_det(define_case):
    run_dict = {'case': 'CASE1',
                'objectives': {'DET': (-1, 1)}
               }
    test_space_obj, test_eval_func, test_params = load_case(run_dict, 'design_space')
    assert test_space_obj.case == 'CASE1'

def test_load_case_rob(define_case):
    run_dict = {'case': 'CASE1',
                'objectives': {'ROB': (-1, 1)}
               }
    test_space_obj, test_eval_func, test_params = load_case(run_dict, 'design_space')
    assert test_space_obj.case == 'CASE1'


def test_case_name(input_case):
    assert input_case.case == 'CASE1'

def test_design_space_name(input_case):       
    assert input_case.design_space == 'design_space'

def test_opt_type(input_case):       
    assert input_case.opt_type == 'ROB'

def test_n_dim(input_case):       
    assert input_case.n_dim == 1

def test_n_par(input_case):       
    assert input_case.n_par == 1

def test_par_dict(input_case):       
    assert input_case.par_dict == {'par_1': 1}

def test_var_dict(input_case):       
    assert input_case.var_dict == {'var_1': [0., 10.]}

def test_lb(input_case):       
    assert input_case.lb == [0.]

def test_ub(input_case):       
    assert input_case.ub == [10.]

def test_upar_dict(input_case):       
    assert input_case.upar_dict == {'par_1': ['relative', 'Gaussian', 0.2], 'var_1': ['absolute', 'Uniform', 1.]}
    
def test_attach_objectives(input_case):
    input_case.attach_objectives((1,1,-1))
    assert input_case.obj == (1,1,-1)
   
                
    

    
