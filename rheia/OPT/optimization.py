"""
The :py:mod:`optimization` module contains functions to load the optimizer,
create the starting samples and run the optimization.
"""

import os
import imp
from shutil import copyfile
from pyDOE import lhs
from rheia.CASES.determine_stoch_des_space import load_case, check_dictionary
import numpy as np
import pandas as pd
#######################
# optimizer functions #
#######################


def parse_available_opt():
    """

    Parse all available optimizers. In this version,
    only NSGA-II is available.

    Returns
    -------
    methods : list of tuples
        The available optimizers with the corresponding classes.

    """

    opt_dir = os.path.join(os.path.split(
        os.path.dirname(
            os.path.abspath(__file__)))[0],
        'OPT')

    # get the optimizer modules
    files = [f for f in os.listdir(opt_dir) if f.endswith('algorithms.py')]
    methods = []
    for file in files:
        obj = imp.load_source(file.split('.')[0], os.path.join(opt_dir, file))
        tmp = obj.return_opt_methods()
        for method in tmp:
            methods.append((method, obj))

    return methods


def load_optimizer(optimizer):
    """

    Load the selected optimizer.

    Parameters
    ----------
    optimizer : string
        The name of the optimizer (currently 'NSGA2').

    Returns
    -------
    opt_obj : object
        The optimization object.

    """

    opt_list = parse_available_opt()
    optimizers = [elem[0] for elem in opt_list]

    # check if optimizer exists and load optimization module
    if optimizer not in optimizers:
        raise KeyError('Optimizer is not available!')

    # get object from optimizer class
    for opt in opt_list:
        if optimizer == opt[0]:
            opt_obj = opt[1].return_opt_obj(optimizer)

    return opt_obj

####################
# starting samples #
####################


def check_existing_results(run_dict, space_obj):
    """
    Check if previously generated results exists in the provided results
    directory.

    Parameters
    ----------
    run_dict : dict
        The dictionary with user-defined values for
        the characterization of the optimization.
    space_obj : object
        The object with information on the design space and stochastic space

    Returns
    -------
    start_from_last_gen : bool
        Boolean that indicates if results exist in the considered result
        directory.

    """

    # path of the population file
    pop_file = os.path.join(os.path.split(os.path.dirname(
        os.path.abspath(__file__)))[0],
        'RESULTS',
        space_obj.case,
        list(run_dict['objectives'].keys())[0], run_dict['results dir'],
        'population.csv')

    # path of the fitness file
    fitness_file = os.path.join(os.path.split(os.path.dirname(
        os.path.abspath(__file__)))[0],
        'RESULTS',
        space_obj.case,
        list(run_dict['objectives'].keys())[0], run_dict['results dir'],
        'fitness.csv')

    # check if the population and fitness file exist
    start_from_last_gen = False
    if os.path.isfile(fitness_file) and os.path.isfile(pop_file):

        start_from_last_gen = True

    return start_from_last_gen


def scale_samples_to_design_space(nondim_doe, space_obj):
    """
    Scales the starting sample to the given design space.

    Parameters
    ----------
    nondim_doe : ndarray
        The non-dimensionalized design of experiment.
    space_obj : object
        The object with information on the design space and stochastic space

    Returns
    -------
    dim_doe : ndarray
        The design of experiment scaled up to the design space.

    """

    dim_doe = np.zeros(nondim_doe.shape)
    for j in range(space_obj.n_dim):
        dim_doe[:, j] = ((space_obj.u_b[j] - space_obj.l_b[j]) *
                         nondim_doe[:, j] + space_obj.l_b[j])

    return dim_doe


def write_starting_samples(doe, filename):
    """

    Writes the starting samples to file.

    Parameters
    ----------
    doe : list
        The set of starting samples.
    filename : string
        The path to the filename for the starting samples.

    """
    with open(filename, 'w') as file:
        for x_in in doe:
            if not isinstance(x_in, list):
                x_in = x_in.tolist()
            for item in x_in:
                file.write('%.8f, ' % item)
            file.write('\n')


def create_starting_samples(run_dict, space_obj, start_from_last_gen):
    """

    Load the starting samples for the optimization run.

    Parameters
    ----------
    run_dict : dict
        The dictionary with user-defined values for
        the characterization of the optimization.
    space_obj : object
        The object with information on the design space and stochastic space
    start_from_last_gen : bool
        Boolean that determines if the starting samples
        start from a previously generated population.

    """

    # define folder to put starting samples file
    doe_path = os.path.join(os.path.split(
        os.path.dirname(
            os.path.abspath(__file__)))[0],
        'OPT',
        'INPUTS',
        space_obj.case,
        '%iD' %
        space_obj.n_dim)

    # check if folder exists in INPUTS
    # if not, create one
    if not os.path.exists(doe_path):
        os.makedirs(doe_path)

    # create the doe set of samples
    doe_filename = os.path.join(doe_path,
                                'DOE_n%i.csv' % run_dict['population size'])

    if not start_from_last_gen:
        # if the starting population needs to be created
        if 'AUTO' in run_dict['x0'][0]:

            # generate the starting population randomly
            if 'RANDOM' in run_dict['x0'][1]:
                ddoe = np.random.random(
                    (run_dict['population size'], space_obj.n_dim))

            # generate the starting population
            # based on Latin Hypercube Sampling
            if 'LHS' in run_dict['x0'][1]:
                ddoe = lhs(space_obj.n_dim,
                           samples=run_dict['population size'])

            # scale the starting samples
            doe = scale_samples_to_design_space(ddoe, space_obj)

            # Write starting samples to file
            write_starting_samples(doe, doe_filename)

        else:
            # the starting population is provided in a custom file
            doe_custom = os.path.join(os.path.split(
                os.path.dirname(
                    os.path.abspath(__file__)))[0],
                'CASES',
                space_obj.case,
                run_dict['x0'][1])

            # check if the custom file exists
            if not os.path.isfile(doe_custom):
                raise NameError(
                    """The initial population file %s
                    is not found in the case folder.""" %
                    os.path.basename(doe_custom))

            copyfile(doe_custom, doe_filename)

    else:
        # the starting population is the final generation from the existing
        # population in the result directory
        doe_custom = os.path.join(os.path.split(
            os.path.dirname(
                os.path.abspath(__file__)))[0],
            'RESULTS',
            space_obj.case,
            list(run_dict['objectives'].keys())[0],
            run_dict['results dir'],
            'population.csv')
            
        df = pd.read_csv(doe_custom, header=None)
        df.drop(df.tail(1).index,inplace=True)

        rows = df.tail(run_dict['population size']).to_numpy()
        

        doe = []
        for line in rows:
            dummy = []
            for l in line:
                dummy.append(float(l))
            doe.append(dummy)
            
        # check if the size of the existing population matches with the
        # provided one in the optimization dictionary
        if len(doe) == run_dict['population size']:
            write_starting_samples(doe, doe_filename)
        else:
            raise ValueError(""" The defined population size does not match
                                 the population size from a previous run in
                                 this result directory.""")

########################
# run the optimization #
########################


def run_opt(run_dict, design_space='design_space.csv'):
    """
    This function runs the optimization pipeline.
    First, the case, optimizer and configuration
    are loaded. Thereafter, the starting samples
    are created, the specific optimization class
    instantiated and the :py:meth:`run_optimizer`
    method is called.

    Parameters
    ----------
    run_dict : dict
        The dictionary with user-defined values for
        the characterization of the optimization.
    design_space : string, optional
        The design_space filename. The default is 'design_space.csv'.

    """

    # evaluate if the optimization dictionary is properly characterized
    check_dictionary(run_dict)

    # load the object on the design space, the evaluation function
    # and the params provided for each model evaluation
    space_obj, eval_func, params = load_case(run_dict, design_space)

    # add the evaluation function to the optimization dictionary
    run_dict['evaluate'] = eval_func

    # load optimizer class
    opt_class = load_optimizer('NSGA2')

    # check if previous results in this result directory exist
    start_from_last_gen = check_existing_results(run_dict, space_obj)

    # define starting samples
    create_starting_samples(run_dict, space_obj, start_from_last_gen)

    # optimizer object
    res = opt_class(run_dict, space_obj, start_from_last_gen, params)

    # run the optimizer
    res.run_optimizer()
