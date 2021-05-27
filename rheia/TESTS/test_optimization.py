import pytest
import numpy as np
import rheia.OPT.optimization as opt

@pytest.fixture
def run_dict():
    dict_opt = {'case':                'H2_FUEL',
                'objectives':          {'DET': (-1, 1)},
                'stop':                ('BUDGET', 2000),
                #'n jobs':              int(mp.cpu_count() / 2),
                'population size':     20,
                'results dir':         'run_tutorial',
                }
    
    return dict_opt

@pytest.fixture
def input_case(run_dict):

    design_space = 'design_space'
    space_obj, eval_func, params = opt.load_case(run_dict, design_space)

    return space_obj, eval_func, params
    
    #test_class = stochastic_design_space(input_dict['objectives'].keys[0], input_dict['case'], design_space)


def test_parse_available_opt():
    methods = opt.parse_available_opt()
    assert methods[0][0] == 'NSGA2'

def test_load_optimizer():
    opt_obj = opt.load_optimizer('NSGA2')
    assert opt_obj.__name__ == 'NSGA2'

def test_check_existing_results(run_dict, input_case):
    bool_test = opt.check_existing_results(run_dict, input_case[0])
    assert bool_test

def test_check_existing_results_false(run_dict, input_case):
    run_dict['results dir'] = 'run_test_no_results'
    bool_test = opt.check_existing_results(run_dict, input_case[0])
    assert bool_test == False

def test_scale_samples_to_design_space(input_case):
    nondim_doe = np.array([[0.5, 0.5],[0.1,0.9]])
    dim_doe = opt.scale_samples_to_design_space(nondim_doe, input_case[0])
    assert round(dim_doe[0][0],1) == 3.5
    assert round(dim_doe[0][1],1) == 3.5
    assert round(dim_doe[1][0],1) == 0.7
    assert round(dim_doe[1][1],1) == 6.3