'''
Module to test the post_process module.
'''

import pytest

from rheia.POST_PROCESS.post_process import PostProcessOpt, PostProcessUQ


@pytest.fixture
def input_case_uq():
    """
    Instantiate a PostProcessUQ object.

    Returns
    -------
    test_obj : object
        PostProcessUQ object.

    """
    case = 'H2_FUEL'
    test_obj = PostProcessUQ(case, 1)

    return test_obj


def test_get_fitness_population():
    """
    Assert the population and fitness samples imported from fitness and
    population files.
    """

    case = 'H2_FUEL'
    eval_type = 'DET'
    test_obj = PostProcessOpt(case, eval_type)

    x, y = test_obj.get_fitness_population('run_tutorial')

    assert round(x[0][0], 3) == 7.796
    assert round(x[0][12], 3) == 11.117
    assert round(x[1][7], 3) == 106.956
    assert round(y[0][0], 3) == 1.933
    assert round(y[0][12], 3) == 4.382
    assert round(y[1][7], 3) == 2.959


def test_get_pdf(input_case_uq):
    """
    Assert the data points for the Probability Density Function, stored in a
    UQ results folder.

    Parameters
    ----------
    input_case_uq : object
        PostProcessUQ object.

    """

    result_dir = 'opt_design_tutorial'
    objective = 'lcoh'
    x, y = input_case_uq.get_pdf(result_dir, objective)

    assert round(x[0], 3) == 4.912
    assert round(x[5], 3) == 5.195
    assert round(x[7], 3) == 5.308
    assert round(y[23], 3) == 0.102
    assert round(y[29], 3) == 0.201
    assert round(y[35], 3) == 0.280


def test_get_cdf(input_case_uq):
    """
    Assert the data points for the Cumulative Distribution Function, stored in
    a UQ results folder.

    Parameters
    ----------
    input_case_uq : object
        PostProcessUQ object.

    """

    result_dir = 'opt_design_tutorial'
    objective = 'lcoh'
    x, y = input_case_uq.get_cdf(result_dir, objective)

    assert round(x[0], 3) == 4.912
    assert round(x[5], 3) == 5.195
    assert round(x[7], 3) == 5.308
    assert round(y[23], 3) == 0.038
    assert round(y[29], 3) == 0.092
    assert round(y[35], 3) == 0.177


def test_get_loo(input_case_uq):
    """
    Assert the Leave-One-Out error, stored in a UQ results folder.

    Parameters
    ----------
    input_case_uq : object
        PostProcessUQ object.

    """

    result_dir = 'opt_design_tutorial'
    objective = 'lcoh'
    loo = input_case_uq.get_loo(result_dir, objective)

    assert round(loo, 3) == 0.063


def test_get_sobol(input_case_uq):
    """
    Assert the Sobol' indices, stored in a UQ results folder.

    Parameters
    ----------
    input_case_uq : object
        PostProcessUQ object.

    """

    result_dir = 'opt_design_tutorial'
    objective = 'lcoh'
    names, sobol = input_case_uq.get_sobol(result_dir, objective)

    assert names[1] == 'capex_pemel'
    assert names[4] == 'opex_pemel'
    assert round(sobol[1], 3) == 0.323
    assert round(sobol[4], 3) == 0.059
