"""
The :mod:`uncertainty_quantification` module provides functions to
execute uncertainty quantification.

"""

import os
import pandas as pd
from scipy.stats import qmc
from rheia.CASES.determine_stoch_des_space import (
    DESIGN_SPACE_COLUMNS, load_case, check_dictionary)
import rheia.UQ.pce as uq


def get_design_variables(case):
    """
    This function loads the design variable names and bounds
    out of the :file:`design_space.csv` file.

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
        'design_space.csv')

    data = pd.read_csv(path_to_read, dtype=str, keep_default_na=False,
                       na_filter=False)
    missing_columns = [column for column in DESIGN_SPACE_COLUMNS
                       if column not in data.columns]
    if missing_columns:
        raise ValueError(
            """The design_space file should contain the columns: %s."""
            % ', '.join(DESIGN_SPACE_COLUMNS))

    data = data[DESIGN_SPACE_COLUMNS]
    data = data.apply(lambda col: col.str.strip())

    for row in data.itertuples(index=False):
        if row.type == 'var':
            var_dict[row.name] = [float(row.value), float(row.upper_bound)]

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

    # generate unit-hypercube samples through SciPy's maintained LHS sampler
    sampler = qmc.LatinHypercube(d=len(var_dict))
    samples = sampler.random(n=n_samples)

    bounds = list(var_dict.values())

    # scale the samples based on the design variable bounds
    for i, bound in enumerate(bounds):
        samples[:, i] *= (bound[1] - bound[0])
        samples[:, i] += bound[0]

    return samples


def write_design_space(case, iteration, var_dict, sample, ds='design_space.csv'):
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
                                'design_space.csv',
                                )

    new_des_var_file = os.path.join(
        os.path.abspath(
            os.path.join(
                path,
                os.pardir)),
        'CASES',
        case,
        '%s_%i%s' % (ds[:-4], iteration, ds[-4:])
    )

    # write the new design_space file if it does not exist already
    if not os.path.isfile(new_des_var_file):
        data = pd.read_csv(des_var_file, dtype=str, keep_default_na=False,
                           na_filter=False)
        missing_columns = [column for column in DESIGN_SPACE_COLUMNS
                           if column not in data.columns]
        if missing_columns:
            raise ValueError(
                """The design_space file should contain the columns: %s."""
                % ', '.join(DESIGN_SPACE_COLUMNS))

        data = data[DESIGN_SPACE_COLUMNS]
        for index, name in enumerate(var_dict):
            mask = data['name'] == name
            data.loc[mask, 'type'] = 'par'
            data.loc[mask, 'value'] = '%f' % sample[index]
            data.loc[mask, 'upper_bound'] = ''

        data.to_csv(new_des_var_file, index=False, lineterminator='\n')


def run_uq(run_dict, design_space='design_space.csv'):
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

    # Validate and complete optional UQ settings before any objects are built.
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

    # calculate the full basis size and the required number of samples
    my_experiment.n_terms()

    # read in the previously generated samples
    my_experiment.read_previous_samples(run_dict['create only samples'])

    n_existing = len(my_experiment.x_prev)
    n_to_create = my_experiment.n_samples - n_existing

    # In create-only mode, only append new input samples to samples.csv.
    if run_dict['create only samples']:
        my_experiment.create_samples(size=n_to_create)
        my_experiment.create_only_samples(True)
        return

    # Evaluate only the new samples needed to reach the requested DOE size.
    # If enough previous samples exist, select the required subset instead.
    if n_to_create > 0:
        my_experiment.evaluate(eval_func, params)
    else:
        my_experiment.create_samples(size=n_to_create)
        my_experiment.y = my_experiment.y_prev[:my_experiment.n_samples]

    # create and fit PCE object
    my_pce = uq.PCE(my_experiment)
    my_pce.run(run_dict['uq method'])

    # calculate final accuracy and sensitivity metrics
    my_pce.calc_loo()
    my_pce.calc_sobol()

    # print a summary and write result files
    my_pce.print_res()

    # generate the pdf and cdf when desired
    if run_dict['draw pdf cdf'][0]:
        my_pce.draw(int(run_dict['draw pdf cdf'][1]))

