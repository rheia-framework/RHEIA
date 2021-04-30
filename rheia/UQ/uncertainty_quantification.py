"""
The :mod:`uncertainty_quantification` module provides functions to
execute uncertainty quantification.

"""

import os
from pyDOE import lhs
from rheia.CASES.determine_stoch_des_space import load_case, check_dictionary
import rheia.UQ.pce as uq


def get_design_variables(case):
    """
    This function loads the design variable names and bounds
    out of the :file:`design_space` file.

    Parameters
    ----------
    case : string
        The name of the case.


    Returns
    -------
    var_dict : dict
        A dictionary which includes the design variables and their bounds.

    """
    var_dict = {}

    path = os.path.dirname(os.path.abspath(__file__))
    path_to_read = os.path.join(
        os.path.abspath(
            os.path.join(
                path,
                os.pardir)),
        'CASES',
        case,
        'design_space')

    # read in the design variable bounds
    with open(path_to_read, 'r') as file:
        for line in file:
            tmp = line.split()
            if tmp[1] == 'var':
                var_dict[tmp[0]] = [float(tmp[2]), float(tmp[3])]

    return var_dict


def set_design_samples(var_dict, n_samples):
    """
    Based on the design variable characteristics,
    a set of design samples is created through
    Latin Hypercube Sampling.

    Parameters
    ----------
    var_dict : dict
        A dictionary which includes the design variables and their bounds.
    n_samples : int
        The number of design samples to be created.

    Returns
    -------
    samples : array
        The generated design samples.

    """

    # generate the samples through Latin Hypercube Sampling
    samples = lhs(len(var_dict), samples=n_samples)

    bounds = list(var_dict.values())

    # scale the samples based on the design variable bounds
    for i, bound in enumerate(bounds):
        samples[:, i] *= (bound[1] - bound[0])
        samples[:, i] += bound[0]

    return samples


def write_design_space(case, iteration, var_dict, sample, ds = 'design_space'):
    """
    A new design space file is created. In this file,
    the model parameters are copied from the original file,
    i.e. file:`design_space`. The design variable names are copied,
    but the bounds are loaded out of the array `sample`.
    This function is of interest when evaluating the LOO error
    or Sobol' indices for several design samples.

    Parameters
    ----------
    case : string
        The name of the case.
    iteration : int
        The index of the design sample
        out of the collection of generated design samples.
    var_dict : dict
        A dictionary which includes the design variables and their bounds.
    sample : array
        The design sample out of the collection of generated design samples.
    ds : string, optional
        The design_space filename. The default is 'design_space'.

    """

    path = os.path.dirname(os.path.abspath(__file__))

    des_var_file = os.path.join(os.path.abspath(os.path.join(path, os.pardir)),
                                'CASES',
                                case,
                                'design_space',
                                )

    new_des_var_file = os.path.join(
        os.path.abspath(
            os.path.join(
                path,
                os.pardir)),
        'CASES',
        case,
        '%s_%i' % (ds, iteration)
    )

    # write the new design_space file if it does not exist already
    if not os.path.isfile(new_des_var_file):
        with open(des_var_file, 'r') as file:
            text = []
            for line in file.readlines():
                found = False
                tmp = line.split()
                for index, name in enumerate(list(var_dict.keys())):

                    if name == tmp[0]:
                        text.append('%s par %f \n' % (name, sample[index]))
                        found = True
                if not found:
                    text.append(line)

        with open(new_des_var_file, 'w') as file:
            for item in text:
                file.write("%s" % item)


def run_uq(run_dict, design_space='design_space'):
    """
    This function is the main to run uncertainty quantification.
    First, the input distributions are created,
    followed by the reading of previously evaluated samples.
    Thereafter, the new samples are created and evaluated
    in the system model when desired. Finally, the PCE is
    constructed, the statistical moments printed and
    the distributions generated (when desired) for the
    quantity of interest.

    Parameters
    ----------
    run_dict : dict
        The dictionary with information on the uncertainty quantification.
    design_space : string, optional
        The design_space filename. The default is 'design_space'.

    """

    # check if the UQ dictionary is properly characterized
    check_dictionary(run_dict, uq_bool=True)

    objective_position = run_dict['objective names'].index(
        run_dict['objective of interest'])

    # load the object on the design space, the evaluation function
    # and the params provided for each model evaluation
    space_obj, eval_func, params = load_case(
        run_dict, design_space, uq_bool=True,
        create_only_samples=run_dict['create only samples'])

    my_data = uq.Data(run_dict, space_obj)

    # acquire information on stochastic parameters
    my_data.read_stoch_parameters()

    # create result csv file to capture all input-output of the samples
    my_data.create_samples_file()

    # create experiment object
    my_experiment = uq.RandomExperiment(my_data, objective_position)

    # create uniform/gaussian distributions and corresponding orthogonal
    # polynomials
    my_experiment.create_distributions()

    # calculate the number of terms in the PCE according to the
    # truncation scheme
    my_experiment.n_terms()

    # read in the previously generated samples
    my_experiment.read_previous_samples(run_dict['create only samples'])

    # create a design of experiment for the remaining samples
    # to be evaluated
    my_experiment.create_samples(
        size=my_experiment.n_samples - len(my_experiment.x_prev))

    # check if the samples need to be evaluated or not
    my_experiment.create_only_samples(run_dict['create only samples'])

    # when the PCE needs to be constructed
    if not run_dict['create only samples']:

        # evaluate the samples remaining to reach the required
        # number of samples for the PCE
        if my_experiment.n_samples > len(my_experiment.x_prev):
            my_experiment.evaluate(eval_func, params)
        elif my_experiment.n_samples == len(my_experiment.x_prev):
            my_experiment.y = my_experiment.y_prev
        else:
            my_experiment.y = my_experiment.y_prev[:my_experiment.n_samples]

        # create PCE object
        my_pce = uq.PCE(my_experiment)

        # evaluate the PCE
        my_pce.run()

        # calculate the LOO error
        my_pce.calc_loo()

        # calculate the Sobol' indices
        my_pce.calc_sobol()

        # extract and print results
        my_pce.print_res()

        # generate the pdf and cdf when desired
        if run_dict['draw pdf cdf'][0]:
            my_pce.draw(int(run_dict['draw pdf cdf'][1]))
