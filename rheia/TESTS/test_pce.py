'''
Module to test the det_stoch_des_space module.
'''

import pytest
import numpy as np
import rheia.UQ.pce as pce
from rheia.CASES.determine_stoch_des_space import load_case


@pytest.fixture
def run_dict():
    """
    Generate a generic input dictionary for uncertainty quantification.

    Returns
    -------
    dict_uq : dict
        Input dictionary for uncertainty quantification.

    """

    dict_uq = {'case': 'H2_FUEL',
               'pol order': 2,
               'objective names': ['lcoh', 'mh2'],
               'objective of interest': 'lcoh',
               'sampling method': 'SOBOL',
               'n jobs': 1,
               'draw pdf cdf': [True, 1e5],
               'results dir': 'opt_design_tutorial'
               }

    return dict_uq


@pytest.fixture
def case(run_dict):
    """
    Instantiate a StochasticDesignSpace object and
    the list of fixed parameters.

    Parameters
    ----------
    run_dict : dict
        Input dictionary for uncertainty quantification.

    Returns
    -------
    space_obj : object
        StochasticDesignSpace object.
    params : list
        List of fixed parameters for evaluation of the model.

    """

    design_space = 'design_space_tutorial_uq.csv'
    space_obj, eval_func, params = load_case(
        run_dict, design_space, uq_bool=True)
    run_dict['evaluate'] = eval_func
    return space_obj, params


@pytest.fixture
def input_case_data(run_dict, case):
    """
    Instantiate a Data object.

    Parameters
    ----------
    run_dict : dict
        Input dictionary for uncertainty quantification.
    case : list
        List that includes the StochasticDesignSpace object and
        list with fixed inputs for model evaluation.

    Returns
    -------
    data_obj : object
        Data object.
    """

    data_obj = pce.Data(run_dict, case[0])
    return data_obj


@pytest.fixture
def input_random_experiment(input_case_data):
    """
    Instantiate a RandomExperiment object.

    Parameters
    ----------
    input_case_data : object
        Data object.

    Returns
    -------
    random_exp_obj : object
        RandomExperiment object.
    """

    input_case_data.read_stoch_parameters()
    input_case_data.create_samples_file()
    random_exp_obj = pce.RandomExperiment(input_case_data, 0)
    return random_exp_obj


@pytest.fixture
def input_pce(input_random_experiment):
    """
    Instantiate a PCE object.

    Parameters
    ----------
    input_random_experiment : object
        RandomExperiment object.

    Returns
    -------
    pce_obj : object
        PCE object.
    """

    pce_obj = pce.PCE(input_random_experiment)

    return pce_obj


def test_data_inputs(run_dict, input_case_data):
    """
    Assert the input dictionary in the Data object.

    Parameters
    ----------
    run_dict : dict
        Input dictionary for uncertainty quantification.
    input_case_data : object
        Data object.
    """

    assert input_case_data.inputs == run_dict


def test_data_space_obj(input_case_data, case):
    """
    Assert the input dictionary in the Data object.

    Parameters
    ----------
    input_case_data : object
        Data object.
    case : list
        List that includes the StochasticDesignSpace object,
        list with fixed inputs for model evaluation.
    """

    assert input_case_data.space_obj == case[0]


def test_read_stoch_parameters(input_case_data):
    """
    Assert the stochastic data dictionary in the Data object.

    Parameters
    ----------
    input_case_data : object
        Data object.
    """

    input_case_data.read_stoch_parameters()
    input_case_data.stoch_data['l_b'] = [
        round(num, 3) for num in input_case_data.stoch_data['l_b']]
    input_case_data.stoch_data['u_b'] = [
        round(num, 3) for num in input_case_data.stoch_data['u_b']]
    input_case_data.stoch_data['mean'] = [
        round(num, 3) for num in input_case_data.stoch_data['mean']]
    input_case_data.stoch_data['deviation'] = [
        round(num, 3) for num in input_case_data.stoch_data['deviation']]
    assert input_case_data.stoch_data['names'] == [
        'u_sol_irr',
        'u_t_amb',
        'capex_pv',
        'opex_pv',
        'capex_pemel',
        'opex_pemel',
        'life_pemel',
        'repl_pemel',
        'capex_dcdc',
        'opex_dcdc',
        'int_rate',
        'infl_rate']
    assert input_case_data.stoch_data['types'] == [
        'Uniform',
        'Uniform',
        'Uniform',
        'Uniform',
        'Uniform',
        'Uniform',
        'Uniform',
        'Uniform',
        'Uniform',
        'Uniform',
        'Uniform',
        'Uniform']
    assert input_case_data.stoch_data['l_b'] == [
        0.9, -0.4, 350.0, 16.0, 1400.0, 0.03, 60000.0,
        0.15, 100.0, 0.01, 0.04, 0.01]
    assert input_case_data.stoch_data['u_b'] == [
        1.1, 0.4, 600.0, 19.0, 2100.0, 0.05, 100000.0,
        0.2, 200.0, 0.05, 0.08, 0.03]
    assert input_case_data.stoch_data['mean'] == [
        1.0, 0.0, 475.0, 17.5, 1750.0, 0.04, 80000.0,
        0.175, 150.0, 0.03, 0.06, 0.02]
    assert input_case_data.stoch_data['deviation'] == [
        0.1, 0.4, 125.0, 1.5, 350.0, 0.01, 20000.0,
        0.025, 50.0, 0.02, 0.02, 0.01]


def test_random_experiment_my_data(input_random_experiment, input_case_data):
    """
    Assert the Data object in the RandomExperiment object.

    Parameters
    ----------
    input_random_experiment : object
        RandomExperiment object.
    input_case_data : object
        Data object.
    """

    assert input_random_experiment.my_data == input_case_data


def test_random_experiment_obj_position(input_random_experiment):
    """
    Assert the objective position in the RandomExperiment object.

    Parameters
    ----------
    input_random_experiment : object
        RandomExperiment object.
    """

    assert input_random_experiment.objective_position == 0


def test_random_experiment_dimension(input_random_experiment):
    """
    Assert the stochastic dimension in the RandomExperiment object.

    Parameters
    ----------
    input_random_experiment : object
        RandomExperiment object.
    """

    assert input_random_experiment.dimension == 12


def test_random_experiment_n_terms(input_random_experiment):
    """
    Assert the number of terms in the PCE in the RandomExperiment object.

    Parameters
    ----------
    input_random_experiment : object
        RandomExperiment object.
    """

    input_random_experiment.n_terms()
    assert input_random_experiment.n_samples == 182


def test_read_prev_samples(input_random_experiment):
    """
    Assert reading the existing samples in the RandomExperiment object.

    Parameters
    ----------
    input_random_experiment : object
        RandomExperiment object.
    """

    input_random_experiment.read_previous_samples(False)
    assert len(input_random_experiment.x_prev) == 182
    assert len(input_random_experiment.x_prev[0]) == 12
    assert len(input_random_experiment.y_prev) == 182
    assert input_random_experiment.y_prev[0][0] == 7.598744


def test_create_distr(input_random_experiment):
    """
    Assert the list of polynomial families in the RandomExperiment object.

    Parameters
    ----------
    input_random_experiment : object
        RandomExperiment object.
    """

    input_random_experiment.create_distributions()
    assert input_random_experiment.polytypes == ['Legendre'] * 12


def test_create_samples(input_random_experiment):
    """
    Assert the samples size in the RandomExperiment object.

    Parameters
    ----------
    input_random_experiment : object
        RandomExperiment object.
    """

    input_random_experiment.read_previous_samples(False)
    input_random_experiment.create_distributions()
    input_random_experiment.create_samples()
    assert input_random_experiment.size == 182
    assert len(input_random_experiment.x_u) == 182
    assert len(input_random_experiment.x_u_scaled) == 182


def test_multindices(input_pce):
    """
    Assert the generic multindices in the PCE object.

    Parameters
    ----------
    input_pce : object
        PCE object.
    """

    assert input_pce.multindices(range(0, 2)) == [
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)]


def test_ols(input_pce):
    """
    Assert the least-square regression method in the PCE object.

    Parameters
    ----------
    input_pce : object
        PCE object.
    """
    a = np.array([[1., -0.96263],
                  [1., -0.94701],
                  [1., 0.0529],
                  [1., 0.5529],
                  ])
    b = np.array([[7.598744],
                  [7.586986],
                  [7.298572],
                  [6.390994]])

    test_res = input_pce.ols(a, b)[0][0]
    assert round(test_res, 3) == 6.996


def test_calc_a(input_random_experiment, input_pce):
    """
    Assert the calculation of the A matrix in the PCE object.

    Parameters
    ----------
    input_random_experiment : object
        RandomExperiment object.
    input_pce : object
        PCE object.
    """

    input_random_experiment.read_previous_samples(False)
    input_random_experiment.create_distributions()
    input_random_experiment.create_samples()
    multindices = input_pce.multindices(range(0, 91))
    test_res = input_pce.calc_a(multindices)
    assert round(test_res[0][1], 3) == -0.963
    assert round(test_res[3][12], 3) == -0.533


def test_run(input_random_experiment, input_pce):
    """
    Assert the construction of the PCE in the PCE object.

    Parameters
    ----------
    input_random_experiment : object
        RandomExperiment object.
    input_pce : object
        PCE object.
    """

    input_random_experiment.n_terms()
    input_random_experiment.read_previous_samples(False)
    input_random_experiment.create_distributions()
    input_random_experiment.create_samples()
    input_random_experiment.y = input_random_experiment.y_prev
    input_pce.run()
    coeff = input_pce.coefficients
    assert round(coeff[0][0], 3) == 7.730
    assert round(coeff[3][0], 3) == 0.559
    assert round(input_pce.moments['mean'], 3) == 7.730
    assert round(input_pce.moments['variance'], 3) == 0.742


def test_calc_sobol(input_random_experiment, input_pce):
    """
    Assert the calculation of the Sobol' indices in the PCE object.

    Parameters
    ----------
    input_random_experiment : object
        RandomExperiment object.
    input_pce : object
        PCE object.
    """

    input_random_experiment.n_terms()
    input_random_experiment.read_previous_samples(False)
    input_random_experiment.create_distributions()
    input_random_experiment.create_samples()
    input_random_experiment.y = input_random_experiment.y_prev
    input_pce.run()
    input_pce.calc_sobol()
    assert round(input_pce.sensitivity['s_i'][0], 3) == 0.059
    assert round(input_pce.sensitivity['s_i'][6], 3) == 0.012
    assert round(input_pce.sensitivity['s_tot_i'][2], 3) == 0.142
    assert round(input_pce.sensitivity['s_tot_i'][9], 3) == 0.002


def test_calc_loo(input_random_experiment, input_pce):
    """
    Assert the calculation of the Leave-One-Out error in the PCE object.

    Parameters
    ----------
    input_random_experiment : object
        RandomExperiment object.
    input_pce : object
        PCE object.
    """

    input_random_experiment.n_terms()
    input_random_experiment.read_previous_samples(False)
    input_random_experiment.create_distributions()
    input_random_experiment.create_samples()
    input_random_experiment.y = input_random_experiment.y_prev
    input_pce.run()
    input_pce.calc_loo()
    assert round(input_pce.loo, 5) == 0.00595
