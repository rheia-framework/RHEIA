"""
Workflow tests for public optimization and UQ entry points.
"""

import csv
import os
import shutil

import numpy as np
import pandas as pd

import rheia.OPT.optimization as rheia_opt
from rheia.CASES.determine_stoch_des_space import load_case
from rheia.POST_PROCESS.post_process import PostProcessUQ
import rheia.UQ.uncertainty_quantification as rheia_uq


PACKAGE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))


def _results_path(*parts):
    return os.path.join(PACKAGE_ROOT, 'RESULTS', *parts)


def _read_csv_rows(path):
    with open(path, newline='') as file:
        return list(csv.reader(file))


def _cleanup(path):
    if os.path.isdir(path):
        shutil.rmtree(path)


def _lognormal_to_normal(mean, deviation):
    sigma = np.sqrt(np.log(1. + (deviation / mean) ** 2))
    mu = np.log(mean) - 0.5 * sigma ** 2
    return mu, sigma


def _four_bar_uq_dict(results_dir, n_samples, create_only_samples=True):
    return {
        'case': 'FOUR_BAR_TRUSS',
        'pol order': 2,
        'uq method': 'sparse',
        'n samples': n_samples,
        'objective names': ['volume', 'displacement'],
        'objective of interest': 'displacement',
        'sampling method': 'SOBOL',
        'n jobs': 1,
        'results dir': results_dir,
        'draw pdf cdf': [False],
        'create only samples': create_only_samples,
    }


def test_uq_no_model_creates_pce_training_samples():
    """
    Assert NO_MODEL can generate input-only samples for a PCE.
    """

    results_dir = 'pytest_no_model_create_only_samples'
    result_path = _results_path('NO_MODEL', 'UQ', results_dir)
    _cleanup(result_path)

    run_dict = {
        'case': 'NO_MODEL',
        'pol order': 2,
        'uq method': 'full',
        'objective names': ['output_1'],
        'objective of interest': 'output_1',
        'sampling method': 'SOBOL',
        'n jobs': 1,
        'results dir': results_dir,
        'draw pdf cdf': [False],
        'create only samples': True,
    }

    try:
        rheia_uq.run_uq(run_dict)

        samples_path = os.path.join(result_path, 'samples.csv')
        sample_rows = _read_csv_rows(samples_path)

        assert len(sample_rows) == 13
        assert sample_rows[0] == ['var_1', 'par_1', 'output_1']
        assert all(len(row) == 2 for row in sample_rows[1:])
        assert os.listdir(result_path) == ['samples.csv']

    finally:
        _cleanup(result_path)


def test_uq_create_only_samples_accepts_lognormal_distribution():
    """
    Assert stochastic_space.csv can define a lognormal input distribution.
    """

    case = 'PYTEST_LOGNORMAL'
    cases_path = os.path.join(PACKAGE_ROOT, 'CASES')
    case_path = os.path.join(cases_path, case)
    results_dir = 'pytest_lognormal_create_only_samples'
    result_path = _results_path(case, 'UQ', results_dir)

    _cleanup(case_path)
    _cleanup(result_path)

    try:
        os.makedirs(case_path)
        with open(os.path.join(case_path, '__init__.py'), 'w'):
            pass
        with open(os.path.join(case_path, 'case_description.py'), 'w'):
            pass
        with open(os.path.join(case_path, 'design_space.csv'), 'w') as file:
            file.write('name,type,value,upper_bound\n')
            file.write('x,par,10,\n')
        with open(os.path.join(case_path, 'stochastic_space.csv'), 'w') as file:
            file.write('name,relation,distribution,deviation\n')
            file.write('x,absolute,Lognormal,2\n')

        run_dict = {
            'case': case,
            'pol order': 1,
            'uq method': 'full',
            'objective names': ['output_1'],
            'objective of interest': 'output_1',
            'sampling method': 'SOBOL',
            'n jobs': 1,
            'results dir': results_dir,
            'draw pdf cdf': [False],
            'create only samples': True,
        }

        rheia_uq.run_uq(run_dict)

        sample_rows = _read_csv_rows(os.path.join(result_path, 'samples.csv'))

        assert len(sample_rows) == 5
        assert sample_rows[0] == ['x', 'output_1']
        assert all(len(row) == 1 for row in sample_rows[1:])
        assert all(float(row[0]) > 0. for row in sample_rows[1:])

    finally:
        _cleanup(case_path)
        _cleanup(result_path)


def test_uq_lognormal_matches_equivalent_gaussian_latent_model():
    """
    Assert Lognormal PCE matches an equivalent latent-Gaussian model.
    """

    lognormal_case = 'PYTEST_LOGNORMAL_EQUIV'
    gaussian_case = 'PYTEST_GAUSSIAN_EQUIV'
    results_dir = 'pytest_lognormal_gaussian_equivalence'
    cases_path = os.path.join(PACKAGE_ROOT, 'CASES')
    lognormal_case_path = os.path.join(cases_path, lognormal_case)
    gaussian_case_path = os.path.join(cases_path, gaussian_case)
    lognormal_result_case_path = _results_path(lognormal_case)
    gaussian_result_case_path = _results_path(gaussian_case)
    lognormal_result_path = _results_path(lognormal_case, 'UQ', results_dir)
    gaussian_result_path = _results_path(gaussian_case, 'UQ', results_dir)

    for path in [lognormal_case_path, gaussian_case_path,
                 lognormal_result_case_path, gaussian_result_case_path]:
        _cleanup(path)

    mean_x = 10.
    std_x = 2.
    mu_x, sigma_x = _lognormal_to_normal(mean_x, std_x)

    lognormal_model = """
def evaluate(x_in, params=[]):
    _, sample = x_in
    u_1 = sample['u_1']
    u_2 = sample['u_2']
    x = sample['x']
    output = 2. + 0.4 * u_1 - 0.2 * u_2 + 0.3 * x + 0.01 * x * u_1
    return (output,)
"""

    gaussian_model = """
import math


def evaluate(x_in, params=[]):
    _, sample = x_in
    u_1 = sample['u_1']
    u_2 = sample['u_2']
    x = math.exp(sample['x'])
    output = 2. + 0.4 * u_1 - 0.2 * u_2 + 0.3 * x + 0.01 * x * u_1
    return (output,)
"""

    try:
        for case_path, model in [(lognormal_case_path, lognormal_model),
                                 (gaussian_case_path, gaussian_model)]:
            os.makedirs(case_path)
            with open(os.path.join(case_path, '__init__.py'), 'w'):
                pass
            with open(os.path.join(case_path, 'case_description.py'),
                      'w') as file:
                file.write(model)

        with open(os.path.join(lognormal_case_path, 'design_space.csv'),
                  'w') as file:
            file.write('name,type,value,upper_bound\n')
            file.write('u_1,par,1,\n')
            file.write('u_2,par,2,\n')
            file.write('x,par,%s,\n' % mean_x)
        with open(os.path.join(lognormal_case_path, 'stochastic_space.csv'),
                  'w') as file:
            file.write('name,relation,distribution,deviation\n')
            file.write('u_1,absolute,Uniform,0.3\n')
            file.write('u_2,absolute,Uniform,0.4\n')
            file.write('x,absolute,Lognormal,%s\n' % std_x)

        with open(os.path.join(gaussian_case_path, 'design_space.csv'),
                  'w') as file:
            file.write('name,type,value,upper_bound\n')
            file.write('u_1,par,1,\n')
            file.write('u_2,par,2,\n')
            file.write('x,par,%s,\n' % mu_x)
        with open(os.path.join(gaussian_case_path, 'stochastic_space.csv'),
                  'w') as file:
            file.write('name,relation,distribution,deviation\n')
            file.write('u_1,absolute,Uniform,0.3\n')
            file.write('u_2,absolute,Uniform,0.4\n')
            file.write('x,absolute,Gaussian,%s\n' % sigma_x)

        base_run_dict = {
            'pol order': 3,
            'uq method': 'full',
            'objective names': ['output'],
            'objective of interest': 'output',
            'sampling method': 'SOBOL',
            'n jobs': 1,
            'results dir': results_dir,
            'draw pdf cdf': [False],
        }

        for case in [lognormal_case, gaussian_case]:
            run_dict = dict(base_run_dict, case=case)
            rheia_uq.run_uq(run_dict)

        lognormal_samples = np.asarray(
            _read_csv_rows(os.path.join(lognormal_result_path, 'samples.csv'))[
                1:], dtype=float)
        gaussian_samples = np.asarray(
            _read_csv_rows(os.path.join(gaussian_result_path, 'samples.csv'))[
                1:], dtype=float)

        assert np.allclose(lognormal_samples[:, :2], gaussian_samples[:, :2])
        assert np.allclose(lognormal_samples[:, 2],
                           np.exp(gaussian_samples[:, 2]))
        assert np.allclose(lognormal_samples[:, 3], gaussian_samples[:, 3])

        lognormal_post = PostProcessUQ(lognormal_case, 3)
        gaussian_post = PostProcessUQ(gaussian_case, 3)

        lognormal_mean_std = lognormal_post.get_mean_std(
            results_dir, 'output')
        gaussian_mean_std = gaussian_post.get_mean_std(
            results_dir, 'output')
        assert np.allclose(lognormal_mean_std, gaussian_mean_std)

        lognormal_names, lognormal_sobol = lognormal_post.get_sobol(
            results_dir, 'output')
        gaussian_names, gaussian_sobol = gaussian_post.get_sobol(
            results_dir, 'output')

        assert lognormal_names == gaussian_names
        assert np.allclose(lognormal_sobol, gaussian_sobol)

    finally:
        pass


def test_uq_create_only_samples_can_continue_with_sobol():
    """
    Assert create-only UQ appends only new Sobol samples on continuation.
    """

    continued_dir = 'pytest_four_bar_uq_create_only_continue'
    one_shot_dir = 'pytest_four_bar_uq_create_only_one_shot'
    continued_path = _results_path('FOUR_BAR_TRUSS', 'UQ', continued_dir)
    one_shot_path = _results_path('FOUR_BAR_TRUSS', 'UQ', one_shot_dir)

    _cleanup(continued_path)
    _cleanup(one_shot_path)

    try:
        rheia_uq.run_uq(
            _four_bar_uq_dict(continued_dir, 4),
            design_space='design_space_uq.csv')
        samples_path = os.path.join(continued_path, 'samples.csv')
        first_rows = _read_csv_rows(samples_path)

        assert len(first_rows) == 5
        assert len(first_rows[0]) == 7
        assert all(len(row) == 5 for row in first_rows[1:])
        assert os.listdir(continued_path) == ['samples.csv']

        rheia_uq.run_uq(
            _four_bar_uq_dict(continued_dir, 6),
            design_space='design_space_uq.csv')
        continued_rows = _read_csv_rows(samples_path)

        rheia_uq.run_uq(
            _four_bar_uq_dict(one_shot_dir, 6),
            design_space='design_space_uq.csv')
        one_shot_rows = _read_csv_rows(
            os.path.join(one_shot_path, 'samples.csv'))

        assert len(continued_rows) == 7
        assert continued_rows[:5] == first_rows
        assert continued_rows == one_shot_rows

    finally:
        _cleanup(continued_path)
        _cleanup(one_shot_path)


def test_uq_reuses_existing_evaluated_samples_for_pce():
    """
    Assert UQ builds result files from existing evaluated samples without appending.
    """

    results_dir = 'pytest_four_bar_uq_reuse_evaluated'
    result_path = _results_path('FOUR_BAR_TRUSS', 'UQ', results_dir)
    _cleanup(result_path)

    try:
        create_dict = _four_bar_uq_dict(results_dir, 12)
        rheia_uq.run_uq(create_dict, design_space='design_space_uq.csv')

        samples_path = os.path.join(result_path, 'samples.csv')
        sample_rows = _read_csv_rows(samples_path)
        header = sample_rows[0]
        input_rows = sample_rows[1:]

        space_obj, eval_func, params = load_case(
            create_dict, 'design_space_uq.csv', uq_bool=True)

        evaluated_rows = []
        for index, row in enumerate(input_rows):
            values = [float(value) for value in row]
            sample_values = list(space_obj.par_dict.values())
            for name, value in zip(header[:5], values):
                sample_values[list(space_obj.par_dict.keys()).index(name)] = value
            sample = space_obj.convert_into_dictionary(sample_values)
            volume, displacement = eval_func((index, sample), params=params)
            evaluated_rows.append(values + [volume, displacement])

        with open(samples_path, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(header)
            writer.writerows(evaluated_rows)

        before = _read_csv_rows(samples_path)
        run_dict = _four_bar_uq_dict(results_dir, 12, create_only_samples=False)
        rheia_uq.run_uq(run_dict, design_space='design_space_uq.csv')
        after = _read_csv_rows(samples_path)

        assert after == before
        assert os.path.isfile(os.path.join(
            result_path,
            'sparse_pce_order_2_displacement_n_samples_12.txt'))
        assert os.path.isfile(os.path.join(
            result_path,
            'sparse_pce_order_2_displacement_n_samples_12_Sobol_indices.csv'))

    finally:
        _cleanup(result_path)


def test_optimization_four_bar_workflow_creates_result_files():
    """
    Assert a small FOUR_BAR_TRUSS optimization run creates valid output files.
    """

    results_dir = 'pytest_four_bar_optimization'
    result_path = _results_path('FOUR_BAR_TRUSS', 'DET', results_dir)
    _cleanup(result_path)

    run_dict = {
        'case': 'FOUR_BAR_TRUSS',
        'objectives': {'DET': (-1, -1)},
        'stop': ('BUDGET', 10),
        'n jobs': 1,
        'population size': 4,
        'results dir': results_dir,
    }

    try:
        rheia_opt.run_opt(run_dict)

        status_path = os.path.join(result_path, 'STATUS.txt')
        population_path = os.path.join(result_path, 'population.csv')
        fitness_path = os.path.join(result_path, 'fitness.csv')

        assert os.path.isfile(status_path)
        assert os.path.isfile(population_path)
        assert os.path.isfile(fitness_path)

        status = pd.read_csv(status_path, sep='\\s+')
        population_rows = _read_csv_rows(population_path)
        fitness_rows = _read_csv_rows(fitness_path)

        assert status['EVALS'].iloc[-1] >= run_dict['stop'][1]
        assert any(row == ['-'] for row in population_rows)
        assert any(row == ['-'] for row in fitness_rows)

    finally:
        _cleanup(result_path)
