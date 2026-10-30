"""
Regression tests for previously reported defects.

Each test here reproduces a concrete failure that users could hit, so that
the fix cannot silently regress.
"""

import csv
import os
import shutil

import pytest

import rheia.OPT.optimization as rheia_opt
import rheia.UQ.uncertainty_quantification as rheia_uq
from rheia.CASES.determine_stoch_des_space import (
    StochasticDesignSpace, check_dictionary)
from rheia.OPT.genetic_algorithms import NSGA2

PACKAGE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir))


def _base_opt_dict():
    return {'case': 'FOUR_BAR_TRUSS',
            'objectives': {'DET': (-1, -1)},
            'stop': ('BUDGET', 20),
            'population size': 10,
            'results dir': 'pytest_regression',
            'n jobs': 1}


# --------------------------------------------------------------------------
# A run dictionary must stay reusable after a run.
#
# run_opt() attaches the evaluation function for internal use. It used to
# attach it to the caller's own dictionary, which check_dictionary() then
# rejected on the next call, so a script that ran two optimizations with the
# same dictionary crashed on the second one.
# --------------------------------------------------------------------------

def test_run_opt_twice_with_the_same_dictionary():
    """The real reproduction: two runs from one dictionary, as a user
    comparing two result directories in a single script would write it."""
    run_dict = _base_opt_dict()
    results_root = os.path.join(PACKAGE_ROOT, 'RESULTS', 'FOUR_BAR_TRUSS',
                                'DET')
    # only the starting-sample file this run creates; the INPUTS folder also
    # holds files that ship with the package
    doe_file = os.path.join(PACKAGE_ROOT, 'OPT', 'INPUTS', 'FOUR_BAR_TRUSS',
                            '4D', 'DOE_n%i.csv' % run_dict['population size'])
    doe_existed = os.path.isfile(doe_file)
    created = []

    try:
        for name in ('pytest_regression_run1', 'pytest_regression_run2'):
            run_dict['results dir'] = name
            created.append(os.path.join(results_root, name))
            rheia_opt.run_opt(run_dict)

        assert 'evaluate' not in run_dict
        for path in created:
            assert os.path.isfile(os.path.join(path, 'fitness.csv'))
    finally:
        for path in created:
            if os.path.isdir(path):
                shutil.rmtree(path)
        if not doe_existed and os.path.isfile(doe_file):
            os.remove(doe_file)


def test_check_dictionary_still_rejects_unknown_keys():
    """The fix above must not weaken validation of genuine typos."""
    run_dict = _base_opt_dict()
    run_dict['populaton size'] = 10  # typo

    with pytest.raises(ValueError):
        check_dictionary(run_dict)


# --------------------------------------------------------------------------
# A starting population sitting exactly on a design variable bound must load.
#
# read_doe() nudges every sample by +2e-8. The bounded crossover and mutation
# operators return samples that sit exactly on a bound, so on a restart that
# nudge pushed them out of the design space and raised.
# --------------------------------------------------------------------------

@pytest.fixture
def nsga2_four_bar(tmp_path):
    space_obj = StochasticDesignSpace(
        'DET', 'FOUR_BAR_TRUSS', 'design_space.csv')
    run_dict = _base_opt_dict()
    run_dict['population size'] = 1
    return NSGA2(run_dict, space_obj, False, []), space_obj


def _write_doe(path, rows):
    with open(path, 'w') as file:
        for row in rows:
            file.write(''.join('%.8f, ' % value for value in row) + '\n')


def test_read_doe_accepts_sample_on_the_upper_bound(nsga2_four_bar, tmp_path):
    optimizer, space_obj = nsga2_four_bar
    path = os.path.join(str(tmp_path), 'doe_upper.csv')
    _write_doe(path, [space_obj.u_b])

    doe = optimizer.read_doe(path)

    assert len(doe) == 1
    for value, upper in zip(doe[0], space_obj.u_b):
        assert value <= upper          # clipped back into the design space
        assert value == pytest.approx(upper)


def test_read_doe_accepts_sample_on_the_lower_bound(nsga2_four_bar, tmp_path):
    optimizer, space_obj = nsga2_four_bar
    path = os.path.join(str(tmp_path), 'doe_lower.csv')
    _write_doe(path, [space_obj.l_b])

    doe = optimizer.read_doe(path)

    for value, lower in zip(doe[0], space_obj.l_b):
        assert value >= lower
        assert value == pytest.approx(lower)


def test_read_doe_keeps_the_offset_for_interior_samples(nsga2_four_bar,
                                                        tmp_path):
    """Interior samples must come back exactly as before: value + 2e-8."""
    optimizer, space_obj = nsga2_four_bar
    interior = [(low + up) / 2. for low, up in zip(space_obj.l_b,
                                                   space_obj.u_b)]
    path = os.path.join(str(tmp_path), 'doe_interior.csv')
    _write_doe(path, [interior])

    doe = optimizer.read_doe(path)

    for value, original in zip(doe[0], interior):
        assert value == pytest.approx(original + 2e-8, abs=1e-15)


def test_read_doe_still_rejects_out_of_bounds_samples(nsga2_four_bar,
                                                      tmp_path):
    """Clipping must not swallow a genuinely invalid starting population."""
    optimizer, space_obj = nsga2_four_bar
    too_high = [value + 1. for value in space_obj.u_b]
    path = os.path.join(str(tmp_path), 'doe_high.csv')
    _write_doe(path, [too_high])

    with pytest.raises(TypeError):
        optimizer.read_doe(path)

    too_low = [value - 1. for value in space_obj.l_b]
    path = os.path.join(str(tmp_path), 'doe_low.csv')
    _write_doe(path, [too_low])

    with pytest.raises(TypeError):
        optimizer.read_doe(path)


# --------------------------------------------------------------------------
# The shipped examples must describe the case they actually run.
# --------------------------------------------------------------------------

def test_h2_fuel_example_objective_names_match_the_model():
    """
    The H2_FUEL example used to carry the four-bar-truss objective names, so
    it reported hydrogen production under the label 'displacement'. The
    counts matched, so nothing raised and the mislabelling was silent.
    """
    examples = os.path.join(PACKAGE_ROOT, os.pardir, 'examples')
    path = os.path.abspath(os.path.join(examples, 'run_uq_h2_fuel.py'))
    if not os.path.isfile(path):
        pytest.skip('examples directory is not part of the installation')

    namespace = {}
    with open(path) as file:
        source = file.read()
    exec(compile(source, path, 'exec'), namespace)  # noqa: S102

    dict_uq = namespace['dict_uq']
    assert dict_uq['case'] == 'H2_FUEL'
    # H2_FUEL.case_description.evaluate returns (lcoh, mh2)
    assert dict_uq['objective names'] == ['lcoh', 'mh2']
    assert dict_uq['objective of interest'] in dict_uq['objective names']


# --------------------------------------------------------------------------
# The samples file is read by column name, not by column position.
#
# The header was written but never read back, so the quantity of interest
# was taken from a fixed column index. A samples file filled in externally
# with its objective columns in another order was then read silently from
# the wrong column.
# --------------------------------------------------------------------------

UQ_RESULTS = os.path.join(PACKAGE_ROOT, 'RESULTS', 'FOUR_BAR_TRUSS', 'UQ')


def _uq_run_dict(results_dir):
    return {'case': 'FOUR_BAR_TRUSS',
            'pol order': 1,
            'objective names': ['volume', 'displacement'],
            'objective of interest': 'displacement',
            'sampling method': 'SOBOL',
            'n jobs': 1,
            'results dir': results_dir,
            'draw pdf cdf': [False]}


def _mean_from(results_dir):
    path = os.path.join(UQ_RESULTS, results_dir,
                        'full_pce_order_1_displacement.txt')
    with open(path) as file:
        for line in file:
            if line.strip().startswith('mean'):
                return float(line.split()[-1])
    raise AssertionError('no mean in %s' % path)


def _wipe_dir(path):
    if os.path.isdir(path):
        shutil.rmtree(path)


def test_samples_file_objectives_are_read_by_name_not_position():
    """A samples file whose objective columns are swapped, header
    included, must give the same answer as the original file."""
    reference, swapped = 'pytest_cols_ref', 'pytest_cols_swapped'
    _wipe_dir(os.path.join(UQ_RESULTS, reference))
    _wipe_dir(os.path.join(UQ_RESULTS, swapped))

    try:
        rheia_uq.run_uq(_uq_run_dict(reference),
                        design_space='design_space_uq.csv')

        with open(os.path.join(UQ_RESULTS, reference, 'samples.csv')) as file:
            rows = list(csv.reader(file))

        os.makedirs(os.path.join(UQ_RESULTS, swapped))
        with open(os.path.join(UQ_RESULTS, swapped, 'samples.csv'), 'w',
                  newline='') as file:
            writer = csv.writer(file, lineterminator='\n')
            for row in rows:
                # swap the two objective columns, header included
                writer.writerow(row[:-2] + [row[-1], row[-2]])

        rheia_uq.run_uq(_uq_run_dict(swapped),
                        design_space='design_space_uq.csv')

        assert _mean_from(swapped) == pytest.approx(_mean_from(reference))
    finally:
        _wipe_dir(os.path.join(UQ_RESULTS, reference))
        _wipe_dir(os.path.join(UQ_RESULTS, swapped))


def test_samples_file_with_unmatched_header_falls_back_and_warns():
    """An older file whose header does not name the objectives keeps
    working, read in column order, but the user is told."""
    reference, renamed = 'pytest_cols_ref2', 'pytest_cols_renamed'
    _wipe_dir(os.path.join(UQ_RESULTS, reference))
    _wipe_dir(os.path.join(UQ_RESULTS, renamed))

    try:
        rheia_uq.run_uq(_uq_run_dict(reference),
                        design_space='design_space_uq.csv')

        with open(os.path.join(UQ_RESULTS, reference, 'samples.csv')) as file:
            rows = list(csv.reader(file))

        os.makedirs(os.path.join(UQ_RESULTS, renamed))
        with open(os.path.join(UQ_RESULTS, renamed, 'samples.csv'), 'w',
                  newline='') as file:
            writer = csv.writer(file, lineterminator='\n')
            # header uses different spellings, values keep their order
            writer.writerow(rows[0][:-2] + ['VOLUME', 'DISPLACEMENT'])
            writer.writerows(rows[1:])

        with pytest.warns(UserWarning, match='order they appear'):
            rheia_uq.run_uq(_uq_run_dict(renamed),
                            design_space='design_space_uq.csv')

        # positional fallback reproduces the previous behaviour exactly
        assert _mean_from(renamed) == pytest.approx(_mean_from(reference))
    finally:
        _wipe_dir(os.path.join(UQ_RESULTS, reference))
        _wipe_dir(os.path.join(UQ_RESULTS, renamed))


def test_appending_to_a_reordered_samples_file_is_refused():
    """New rows are written in run-dictionary order, so appending to a
    reordered file would leave two column orders in one file."""
    reference, swapped = 'pytest_cols_ref3', 'pytest_cols_append'
    _wipe_dir(os.path.join(UQ_RESULTS, reference))
    _wipe_dir(os.path.join(UQ_RESULTS, swapped))

    try:
        rheia_uq.run_uq(_uq_run_dict(reference),
                        design_space='design_space_uq.csv')

        with open(os.path.join(UQ_RESULTS, reference, 'samples.csv')) as file:
            rows = list(csv.reader(file))

        os.makedirs(os.path.join(UQ_RESULTS, swapped))
        with open(os.path.join(UQ_RESULTS, swapped, 'samples.csv'), 'w',
                  newline='') as file:
            writer = csv.writer(file, lineterminator='\n')
            # swapped columns, and one row short so a new one is needed
            writer.writerow(rows[0][:-2] + [rows[0][-1], rows[0][-2]])
            for row in rows[1:-1]:
                writer.writerow(row[:-2] + [row[-1], row[-2]])

        run_dict = _uq_run_dict(swapped)
        with pytest.raises(ValueError, match='different order'):
            rheia_uq.run_uq(run_dict, design_space='design_space_uq.csv')
    finally:
        _wipe_dir(os.path.join(UQ_RESULTS, reference))
        _wipe_dir(os.path.join(UQ_RESULTS, swapped))
