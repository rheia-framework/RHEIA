'''
Module to test the det_stoch_des_space module.
'''

import pytest
import rheia.OPT.optimization as opt
import rheia.OPT.genetic_algorithms as ga
from deap import creator, base


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
                'n jobs': 1,
                'population size': 20,
                'results dir': 'run_tutorial',
                }

    return dict_opt


@pytest.fixture
def input_case(run_dict):
    """
    Instantiate a NSGA2 object.

    Parameters
    ----------
    run_dict : dict
        Input dictionary for optimization.

    Returns
    -------
    nsga_obj : object
        NSGA2 object.
    space_obj : object
        StochasticDesignSpace object.
    params : list
        List of fixed inputs for evaluation of the model.
    eval_func : function
        Function for evaluation of the model.

    """

    design_space = 'design_space'
    space_obj, eval_func, params = opt.load_case(run_dict, design_space)
    run_dict['evaluate'] = eval_func
    nsga_obj = ga.NSGA2(run_dict, space_obj, True, params)
    return nsga_obj, space_obj, params, eval_func


def test_nsga_obj_run_dict(run_dict, input_case):
    """
    Assert the input dictionary for optimization in NSGA2 object.

    Parameters
    ----------
    run_dict : dict
        Input dictionary for optimization.
    input_case : list
        List that includes the NSGA2 object, StochasticDesignSpace object,
        list with fixed inputs for model evaluation and the function to
        evaluate the model.

    """

    assert input_case[0].run_dict == run_dict


def test_nsga_obj_space_obj(input_case):
    """
    Assert the allocation of the StochasticDesignSpace object in NSGA2 object.

    Parameters
    ----------
    input_case : list
        List that includes the NSGA2 object, StochasticDesignSpace object,
        list with fixed inputs for model evaluation and the function to
        evaluate the model.

    """

    assert input_case[0].space_obj == input_case[1]


def test_nsga_obj_start_from_last_gen(input_case):
    """
    Assert the boolean for starting from the last generation in NSGA2 object.

    Parameters
    ----------
    input_case : list
        List that includes the NSGA2 object, StochasticDesignSpace object,
        list with fixed inputs for model evaluation and the function to
        evaluate the model.

    """

    assert input_case[0].start_from_last_gen


def test_nsga_obj_params(run_dict, input_case):
    """
    Assert the list with fixed parameters in NSGA2 object.

    Parameters
    ----------
    input_case : list
        List that includes the NSGA2 object, StochasticDesignSpace object,
        list with fixed inputs for model evaluation and the function to
        evaluate the model.

    """

    assert input_case[0].params == input_case[2]


def test_parse_status(input_case):
    """
    Assert the status in the STATUS file in NSGA2 object.

    Parameters
    ----------
    input_case : list
        List that includes the NSGA2 object, StochasticDesignSpace object,
        list with fixed inputs for model evaluation and the function to
        evaluate the model.

    """

    assert input_case[0].parse_status() == (110, 2000)


def test_define_samples_to_eval(input_case):
    """
    Assert the creation of a population in NSGA2 object.

    Parameters
    ----------
    input_case : list
        List that includes the NSGA2 object, StochasticDesignSpace object,
        list with fixed inputs for model evaluation and the function to
        evaluate the model.

    """

    pop = [[2., 5.], [3., 6.], [1., 2.], [7., 2.]]
    assert input_case[0].define_samples_to_eval(pop) == ([[5.0,
                                                           1.0,
                                                           0.0,
                                                           475.0,
                                                           17.5,
                                                           1750.0,
                                                           0.04,
                                                           80000.0,
                                                           0.175,
                                                           150.0,
                                                           0.03,
                                                           0.06,
                                                           0.02,
                                                           2.0,
                                                           5.0],
                                                          [5.0,
                                                           1.0,
                                                           0.0,
                                                           475.0,
                                                           17.5,
                                                           1750.0,
                                                           0.04,
                                                           80000.0,
                                                           0.175,
                                                           150.0,
                                                           0.03,
                                                           0.06,
                                                           0.02,
                                                           3.0,
                                                           6.0],
                                                          [5.0,
                                                           1.0,
                                                           0.0,
                                                           475.0,
                                                           17.5,
                                                           1750.0,
                                                           0.04,
                                                           80000.0,
                                                           0.175,
                                                           150.0,
                                                           0.03,
                                                           0.06,
                                                           0.02,
                                                           1.0,
                                                           2.0],
                                                          [5.0,
                                                           1.0,
                                                           0.0,
                                                           475.0,
                                                           17.5,
                                                           1750.0,
                                                           0.04,
                                                           80000.0,
                                                           0.175,
                                                           150.0,
                                                           0.03,
                                                           0.06,
                                                           0.02,
                                                           7.0,
                                                           2.0]],
                                                         [])


def test_evaluate_samples(input_case):
    """
    Assert the result from evaluation a generic sample in NSGA2 object.

    Parameters
    ----------
    input_case : list
        List that includes the NSGA2 object, StochasticDesignSpace object,
        list with fixed inputs for model evaluation and the function to
        evaluate the model.

    """

    samples = [[5.0, 1.0, 0.0, 475.0, 17.5, 1750.0, 0.04,
                80000.0, 0.175, 150.0, 0.03, 0.06, 0.02, 2.0, 5.0]]
    fitness = input_case[0].evaluate_samples(samples)
    assert round(fitness[0][0], 0) == 14
    assert round(fitness[0][1], 0) == 95


def test_assign_fitness_to_population(input_case):
    """
    Assert the assignment of the model output to the fitness in NSGA2 object.

    Parameters
    ----------
    input_case : list
        List that includes the NSGA2 object, StochasticDesignSpace object,
        list with fixed inputs for model evaluation and the function to
        evaluate the model.

    """

    weigth = list(input_case[0].space_obj.obj.values())[0]
    creator.create('Fitness', base.Fitness, weights=weigth)
    creator.create('Individual', list, fitness=creator.Fitness)
    samples = [creator.Individual([5.0,
                                   1.0,
                                   0.0,
                                   475.0,
                                   17.5,
                                   1750.0,
                                   0.04,
                                   80000.0,
                                   0.175,
                                   150.0,
                                   0.03,
                                   0.06,
                                   0.02,
                                   2.0,
                                   5.0])]
    fitness = [[14., 95.]]
    test_pop = input_case[0].assign_fitness_to_population(samples, fitness, [])
    assert test_pop[0].fitness.values == (14., 95.)


def test_return_opt_methods():
    """
    Assert if the NSGA2 method is available in the package.
    """

    methods = ga.return_opt_methods()
    assert methods == ['NSGA2']


def test_return_opt_obj():
    """
    Assert the name for a test object of the NSGA2 class.
    """

    test_object = ga.return_opt_obj('NSGA2')
    assert test_object.__name__ == 'NSGA2'
