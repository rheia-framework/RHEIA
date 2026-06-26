'''
Module to test the post_process module.
'''

import numpy as np
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
    test_obj = PostProcessUQ(case, 2)

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


def test_get_hypervolume(tmp_path):
    """
    Assert the hypervolume is calculated for every generation.
    """

    result_path = tmp_path / 'run_hv'
    result_path.mkdir()
    with open(result_path / 'fitness.csv', 'w') as file:
        file.write('1,4\n')
        file.write('3,2\n')
        file.write('4,1\n')
        file.write('-\n')
        file.write('1,3\n')
        file.write('2,2\n')
        file.write('3,1\n')
    with open(result_path / 'population.csv', 'w') as file:
        file.write('0,0\n')
        file.write('0,0\n')
        file.write('0,0\n')
        file.write('-\n')
        file.write('0,0\n')
        file.write('0,0\n')
        file.write('0,0\n')

    test_obj = PostProcessOpt('dummy', 'DET')
    test_obj.result_path = str(tmp_path)

    generations, hypervolume = test_obj.get_hypervolume(
        'run_hv', reference_point=[5., 5.])

    assert np.array_equal(generations, np.array([1, 2]))
    assert np.allclose(hypervolume, np.array([9., 13.]))


def test_get_hypervolume_with_maximization_objective(tmp_path):
    """
    Assert maximization objectives are transformed before hypervolume.
    """

    result_path = tmp_path / 'run_hv_max'
    result_path.mkdir()
    with open(result_path / 'fitness.csv', 'w') as file:
        file.write('1,1\n')
        file.write('3,3\n')
    with open(result_path / 'population.csv', 'w') as file:
        file.write('0,0\n')
        file.write('0,0\n')

    test_obj = PostProcessOpt('dummy', 'DET')
    test_obj.result_path = str(tmp_path)

    generations, hypervolume = test_obj.get_hypervolume(
        'run_hv_max',
        reference_point=[5., 0.],
        objective_weights=[-1., 1.])

    assert np.array_equal(generations, np.array([1]))
    assert np.allclose(hypervolume, np.array([8.]))


def test_get_hypervolume_rejects_bad_reference_point(tmp_path):
    """
    Assert the reference point must be worse than the Pareto points.
    """

    result_path = tmp_path / 'run_hv_bad_reference'
    result_path.mkdir()
    with open(result_path / 'fitness.csv', 'w') as file:
        file.write('1,4\n')
        file.write('3,2\n')
    with open(result_path / 'population.csv', 'w') as file:
        file.write('0,0\n')
        file.write('0,0\n')

    test_obj = PostProcessOpt('dummy', 'DET')
    test_obj.result_path = str(tmp_path)

    with pytest.raises(ValueError, match="reference point"):
        test_obj.get_hypervolume('run_hv_bad_reference',
                                 reference_point=[2., 5.])


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

    assert round(x[0], 3) == 5.187
    assert round(x[5], 3) == 5.481
    assert round(x[7], 3) == 5.598
    assert round(y[23], 3) == 0.204
    assert round(y[29], 3) == 0.324
    assert round(y[35], 3) == 0.415


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

    assert round(x[0], 3) == 5.187
    assert round(x[5], 3) == 5.481
    assert round(x[7], 3) == 5.598
    assert round(y[23], 3) == 0.084
    assert round(y[29], 3) == 0.180
    assert round(y[35], 3) == 0.315


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

    assert round(loo, 5) == 0.00595
