import pytest
import numpy as np
import rheia.OPT.optimization as opt
import rheia.OPT.genetic_algorithms as ga
from deap import creator, base

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
    space_obj, eval_func, params = opt.load_case(run_dict, design_space)
    run_dict['evaluate'] = eval_func
    nsga_obj = ga.NSGA2(run_dict, space_obj, True, params)
    return nsga_obj, space_obj, params, eval_func


def test_nsga_obj_run_dict(run_dict, input_case):
    assert input_case[0].run_dict == run_dict

def test_nsga_obj_space_obj(run_dict, input_case):
    assert input_case[0].space_obj == input_case[1]

def test_nsga_obj_start_from_last_gen(run_dict, input_case):
    assert input_case[0].start_from_last_gen == True

def test_nsga_obj_params(run_dict, input_case):
    assert input_case[0].params == input_case[2]

def test_parse_status(input_case):
    assert input_case[0].parse_status() == (110,2000)

def test_define_samples_to_eval(input_case):
    pop = [[2.,5.],[3.,6.],[1.,2.],[7.,2.]]
    samples = input_case[0].define_samples_to_eval(pop)
    assert input_case[0].define_samples_to_eval(pop) == ([[5.0, 1.0, 0.0, 475.0, 17.5, 1750.0, 0.04, 80000.0, 0.175, 150.0, 0.03, 0.06, 0.02, 2.0, 5.0], 
                                                     [5.0, 1.0, 0.0, 475.0, 17.5, 1750.0, 0.04, 80000.0, 0.175, 150.0, 0.03, 0.06, 0.02, 3.0, 6.0],
                                                     [5.0, 1.0, 0.0, 475.0, 17.5, 1750.0, 0.04, 80000.0, 0.175, 150.0, 0.03, 0.06, 0.02, 1.0, 2.0],
                                                     [5.0, 1.0, 0.0, 475.0, 17.5, 1750.0, 0.04, 80000.0, 0.175, 150.0, 0.03, 0.06, 0.02, 7.0, 2.0]], [])

def test_evaluate_samples(input_case):
    samples = [[5.0, 1.0, 0.0, 475.0, 17.5, 1750.0, 0.04, 80000.0, 0.175, 150.0, 0.03, 0.06, 0.02, 2.0, 5.0]]
    fitness = input_case[0].evaluate_samples(samples) 
    assert round(fitness[0][0],0) == 14
    assert round(fitness[0][1],0) == 95

def test_assign_fitness_to_population(input_case):
    weigth = list(input_case[0].space_obj.obj.values())[0]
    creator.create('Fitness', base.Fitness, weights=weigth)
    creator.create('Individual', list, fitness=creator.Fitness)
    samples = [creator.Individual([5.0, 1.0, 0.0, 475.0, 17.5, 1750.0, 0.04, 80000.0, 0.175, 150.0, 0.03, 0.06, 0.02, 2.0, 5.0])]
    fitness = [[14.,95.]] 
    test_pop = input_case[0].assign_fitness_to_population(samples, fitness, [])
    assert test_pop[0].fitness.values == (14.,95.)

def test_return_opt_methods():
    methods = ga.return_opt_methods()
    assert methods == ['NSGA2']

def test_return_opt_obj():
    object = ga.return_opt_obj('NSGA2')
    assert object.__name__ == 'NSGA2'

