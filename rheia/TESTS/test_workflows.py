"""
Workflow tests for public optimization and UQ entry points.
"""

import csv
import os
import shutil

import pandas as pd

import rheia.OPT.optimization as rheia_opt
from rheia.CASES.determine_stoch_des_space import load_case
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
