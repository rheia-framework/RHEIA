'''
Module to test the det_stoch_des_space module.
'''

import pytest
import numpy as np
import rheia.OPT.optimization as opt


@pytest.fixture
def run_dict():
    """
    Generate a generic input dictionary for optimization.

    Returns
    -------
    dict_opt : dict
        Input dictionary for optimization.

    """

    dict_opt = {'case': 'H2_FUEL',
                'objectives': {'DET': (-1, 1)},
                'stop': ('BUDGET', 2000),
                'population size': 20,
                'results dir': 'run_tutorial',
                }

    return dict_opt


@pytest.fixture
def input_case(run_dict):
    """
    Instantiate a StochasticDesignSpace object, get the function
    to evaluate the model and the list of fixed parameters.

    Parameters
    ----------
    run_dict : dict
        Input dictionary for optimization.

    Returns
    -------
    space_obj : object
        StochasticDesignSpace object.
    params : list
        List of fixed parameters for evaluation of the model.
    eval_func : function
        Function for evaluation of the model.

    """

    design_space = 'design_space'
    space_obj, eval_func, params = opt.load_case(run_dict, design_space)

    return space_obj, eval_func, params


def test_parse_available_opt():
    """
    Assert if the NSGA2 method is available in the package.
    """

    methods = opt.parse_available_opt()
    assert methods[0][0] == 'NSGA2'


def test_load_optimizer():
    """
    Assert the name for a test object of the NSGA2 class.
    """

    opt_obj = opt.load_optimizer('NSGA2')
    assert opt_obj.__name__ == 'NSGA2'


def test_check_existing_results(run_dict, input_case):
    """
    Assert function that checks if existing results are present.

    Parameters
    ----------
    run_dict : dict
        Input dictionary for optimization.
    input_case : list
        List that includes the StochasticDesignSpace object,
        list with fixed inputs for model evaluation and the function to
        evaluate the model.

    """

    bool_test = opt.check_existing_results(run_dict, input_case[0])
    assert bool_test


def test_check_existing_results_false(run_dict, input_case):
    """
    Assert function that checks if existing results are present.

    Parameters
    ----------
    run_dict : dict
        Input dictionary for optimization.
    input_case : list
        List that includes the StochasticDesignSpace object,
        list with fixed inputs for model evaluation and the function to
        evaluate the model.

    """

    run_dict['results dir'] = 'run_test_no_results'
    bool_test = opt.check_existing_results(run_dict, input_case[0])
    assert bool_test is False


def test_scale_samples_to_design_space(input_case):
    """
    Assert function that scales samples to the design space.

    Parameters
    ----------
    run_dict : dict
        Input dictionary for optimization.
    input_case : list
        List that includes the StochasticDesignSpace object,
        list with fixed inputs for model evaluation and the function to
        evaluate the model.

    """

    nondim_doe = np.array([[0.5, 0.5], [0.1, 0.9]])
    dim_doe = opt.scale_samples_to_design_space(nondim_doe, input_case[0])
    assert round(dim_doe[0][0], 1) == 3.5
    assert round(dim_doe[0][1], 1) == 3.5
    assert round(dim_doe[1][0], 1) == 0.7
    assert round(dim_doe[1][1], 1) == 6.3
