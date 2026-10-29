"""
Module to test the optional 'seed' key, which makes a run reproducible.

Omitting the key has to leave every sampler behaving as it did before the
key existed, so the tests here only assert what seeding adds: the same seed
reproduces a run, and a different seed changes it.
"""

import os
import shutil

import numpy as np
import pytest

import rheia.OPT.optimization as rheia_opt
import rheia.UQ.uncertainty_quantification as rheia_uq
from rheia.CASES.determine_stoch_des_space import check_dictionary

PACKAGE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir))
RESULTS = os.path.join(PACKAGE_ROOT, 'RESULTS')
INPUTS = os.path.join(PACKAGE_ROOT, 'OPT', 'INPUTS')


def _wipe(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
    elif os.path.isfile(path):
        os.remove(path)


def _capture(directory):
    """Return the exact contents of every file in a result directory."""
    captured = {}
    for name in sorted(os.listdir(directory)):
        full = os.path.join(directory, name)
        if os.path.isfile(full):
            with open(full) as file:
                captured[name] = file.read()
    return captured


# --------------------------------------------------------------------------
# validation
# --------------------------------------------------------------------------

def _uq_dict(overrides=None):
    run_dict = {'case': 'FOUR_BAR_TRUSS',
                'pol order': 1,
                'objective names': ['volume', 'displacement'],
                'objective of interest': 'displacement',
                'sampling method': 'RANDOM',
                'n jobs': 1,
                'results dir': 'pytest_seed',
                'draw pdf cdf': [False]}
    if overrides:
        run_dict.update(overrides)
    return run_dict


def test_check_dictionary_accepts_a_seed():
    check_dictionary(_uq_dict({'seed': 0}), uq_bool=True)
    check_dictionary(_uq_dict({'seed': 12345}), uq_bool=True)
    # None means "leave the run unseeded" and is explicitly allowed
    check_dictionary(_uq_dict({'seed': None}), uq_bool=True)


def test_check_dictionary_does_not_invent_a_seed():
    """Absence of the key must stay absent, so the samplers keep their
    previous behaviour rather than silently becoming seeded."""
    run_dict = _uq_dict()
    check_dictionary(run_dict, uq_bool=True)

    assert 'seed' not in run_dict


@pytest.mark.parametrize('bad', [-1, 1.5, 'abc', True])
def test_check_dictionary_rejects_an_invalid_seed(bad):
    with pytest.raises((TypeError, ValueError)):
        check_dictionary(_uq_dict({'seed': bad}), uq_bool=True)


def test_seed_is_accepted_for_optimization_too():
    check_dictionary({'case': 'FOUR_BAR_TRUSS',
                      'objectives': {'DET': (-1, -1)},
                      'stop': ('BUDGET', 20),
                      'population size': 10,
                      'results dir': 'pytest_seed',
                      'n jobs': 1,
                      'seed': 3})


# --------------------------------------------------------------------------
# set_design_samples
# --------------------------------------------------------------------------

def test_set_design_samples_is_reproducible_with_a_seed():
    var_dict = {'a': [0., 1.], 'b': [10., 20.]}

    first = rheia_uq.set_design_samples(var_dict, 8, seed=4)
    second = rheia_uq.set_design_samples(var_dict, 8, seed=4)
    other = rheia_uq.set_design_samples(var_dict, 8, seed=5)

    assert np.array_equal(first, second)
    assert not np.array_equal(first, other)

    # the samples still respect the bounds
    assert np.all(first[:, 0] >= 0.) and np.all(first[:, 0] <= 1.)
    assert np.all(first[:, 1] >= 10.) and np.all(first[:, 1] <= 20.)


# --------------------------------------------------------------------------
# uncertainty quantification
# --------------------------------------------------------------------------

def _run_uq_with_seed(seed, label, sampling):
    results_dir = 'pytest_seed_%s' % label
    path = os.path.join(RESULTS, 'FOUR_BAR_TRUSS', 'UQ', results_dir)
    _wipe(path)

    run_dict = _uq_dict({'results dir': results_dir,
                         'sampling method': sampling,
                         'seed': seed})
    rheia_uq.run_uq(run_dict, design_space='design_space_uq.csv')

    captured = _capture(path)
    _wipe(path)
    return captured


@pytest.mark.parametrize('sampling', ['RANDOM', 'SOBOL'])
def test_uq_run_is_reproducible_with_a_seed(sampling):
    first = _run_uq_with_seed(7, 'a_%s' % sampling, sampling)
    second = _run_uq_with_seed(7, 'b_%s' % sampling, sampling)
    other = _run_uq_with_seed(99, 'c_%s' % sampling, sampling)

    assert first == second
    assert first != other


# --------------------------------------------------------------------------
# optimization
# --------------------------------------------------------------------------

def _run_opt_with_seed(seed, label, x0):
    results_dir = 'pytest_seed_%s' % label
    path = os.path.join(RESULTS, 'FOUR_BAR_TRUSS', 'DET', results_dir)
    population_size = 8
    doe_file = os.path.join(INPUTS, 'FOUR_BAR_TRUSS', '4D',
                            'DOE_n%i.csv' % population_size)
    doe_existed = os.path.isfile(doe_file)

    _wipe(path)
    if not doe_existed:
        _wipe(doe_file)

    rheia_opt.run_opt({'case': 'FOUR_BAR_TRUSS',
                       'objectives': {'DET': (-1, -1)},
                       'stop': ('BUDGET', 24),
                       'population size': population_size,
                       'n jobs': 1,
                       'results dir': results_dir,
                       'x0': ('AUTO', x0),
                       'seed': seed})

    captured = _capture(path)
    _wipe(path)
    if not doe_existed:
        _wipe(doe_file)
    return captured


@pytest.mark.parametrize('x0', ['LHS', 'RANDOM'])
def test_optimization_is_reproducible_with_a_seed(x0):
    """Covers the starting population and the crossover and mutation
    operators, which draw from different generators."""
    first = _run_opt_with_seed(7, 'a_%s' % x0, x0)
    second = _run_opt_with_seed(7, 'b_%s' % x0, x0)
    other = _run_opt_with_seed(99, 'c_%s' % x0, x0)

    assert first == second
    assert first != other
