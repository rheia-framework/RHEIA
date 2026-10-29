"""
Module to test the characterize_uncertainty module.
"""

import csv
import math
import os
import shutil
from statistics import NormalDist

import matplotlib
import pytest
from scipy import stats

import rheia.UQ.characterize_uncertainty as char_uq
import rheia.UQ.uncertainty_quantification as rheia_uq
from rheia.UQ.pce import _lognormal_to_normal

# headless: plot_fit()'s show=True must not block or pop up a window here
matplotlib.use('Agg')

PACKAGE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir))
CASES_PATH = os.path.join(PACKAGE_ROOT, 'CASES')

HEADER = 'name,distribution,x_low,p_low,x_high,p_high\n'


def _cleanup(path):
    if os.path.isdir(path):
        shutil.rmtree(path)


def _design(case_path):
    return os.path.join(case_path, 'design_space.csv')


def _stoch(case_path):
    return os.path.join(case_path, 'stochastic_space.csv')


def _char(case_path):
    return os.path.join(case_path, char_uq.CHARACTERIZATION_FILE)


def _read_csv_rows(path):
    with open(path, newline='') as file:
        return list(csv.DictReader(file))


def _make_case(case, design_rows=None, stoch_rows=None,
               characterization=None):
    """Create a bare case folder, optionally with design/stochastic
    spaces and a characterization file."""
    case_path = os.path.join(CASES_PATH, case)
    os.makedirs(case_path)
    with open(os.path.join(case_path, '__init__.py'), 'w'):
        pass
    with open(os.path.join(case_path, 'case_description.py'), 'w'):
        pass

    if design_rows is not None:
        with open(_design(case_path), 'w') as file:
            file.write('name,type,value,upper_bound\n')
            file.writelines(design_rows)
    if stoch_rows is not None:
        with open(_stoch(case_path), 'w') as file:
            file.write('name,relation,distribution,deviation\n')
            file.writelines(stoch_rows)
    if characterization is not None:
        with open(_char(case_path), 'w') as file:
            file.write(HEADER)
            file.writelines(characterization)

    return case_path


# --------------------------------------------------------------------------
# fit_distribution: each fitted distribution must reproduce the exact two
# percentiles it was fitted from.
# --------------------------------------------------------------------------

@pytest.mark.parametrize('x_low,p_low,x_high,p_high', [
    (10., 0.1, 90., 0.9),
    (-5., 0.05, 5., 0.95),
])
def test_fit_gaussian_reproduces_quantiles(x_low, p_low, x_high, p_high):
    mean, std = char_uq.fit_distribution(
        'Gaussian', x_low, p_low, x_high, p_high)
    dist = stats.norm(loc=mean, scale=std)

    assert dist.cdf(x_low) == pytest.approx(p_low)
    assert dist.cdf(x_high) == pytest.approx(p_high)


def test_fit_uniform_reproduces_quantiles():
    x_low, p_low, x_high, p_high = 10., 0.2, 40., 0.8
    mean, deviation = char_uq.fit_distribution(
        'Uniform', x_low, p_low, x_high, p_high)
    dist = stats.uniform(loc=mean - deviation, scale=2. * deviation)

    assert dist.cdf(x_low) == pytest.approx(p_low)
    assert dist.cdf(x_high) == pytest.approx(p_high)


def test_fit_lognormal_reproduces_quantiles():
    x_low, p_low, x_high, p_high = 500., 0.50, 3000., 0.99
    mean, deviation = char_uq.fit_distribution(
        'Lognormal', x_low, p_low, x_high, p_high)

    mu, sigma = _lognormal_to_normal(mean, deviation)
    dist = stats.lognorm(s=sigma, scale=math.exp(mu))

    assert dist.cdf(x_low) == pytest.approx(p_low)
    assert dist.cdf(x_high) == pytest.approx(p_high)


def test_fit_lognormal_matches_cooks_method_reference_values():
    """
    Cross-check against Cook's method applied directly (mu/sigma fit on
    ln(x)), for a declared value at p50 and a design value at p99.
    """
    norm = NormalDist()

    x_low, p_low, x_high, p_high = 500., 0.50, 3000., 0.99
    z_1, z_2 = norm.inv_cdf(p_low), norm.inv_cdf(p_high)
    x_1, x_2 = math.log(x_low), math.log(x_high)
    expected_sigma = (x_2 - x_1) / (z_2 - z_1)
    expected_mu = (x_1 * z_2 - x_2 * z_1) / (z_2 - z_1)

    mean, deviation = char_uq.fit_distribution(
        'Lognormal', x_low, p_low, x_high, p_high)
    mu, sigma = _lognormal_to_normal(mean, deviation)

    assert mu == pytest.approx(expected_mu)
    assert sigma == pytest.approx(expected_sigma)


def test_fit_distribution_rejects_bad_quantile_order():
    with pytest.raises(ValueError):
        char_uq.fit_distribution('Gaussian', 10., 0.9, 90., 0.1)
    with pytest.raises(ValueError):
        char_uq.fit_distribution('Gaussian', 90., 0.1, 10., 0.9)


def test_fit_distribution_rejects_unknown_distribution():
    with pytest.raises(ValueError):
        char_uq.fit_distribution('Triangular', 10., 0.1, 90., 0.9)


def test_fit_lognormal_rejects_nonpositive_quantiles():
    with pytest.raises(ValueError):
        char_uq.fit_distribution('Lognormal', -1., 0.1, 90., 0.9)


# --------------------------------------------------------------------------
# read_characterization_csv: parsing and validating the percentile file.
# --------------------------------------------------------------------------

def test_read_characterization_csv_parses_rows():
    case = 'PYTEST_CHAR_READ'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        _make_case(case, characterization=[
            'th_cond,Lognormal,500,0.5,3000,0.99\n',
            ' capex , Uniform , 400 , 0.05 , 550 , 0.95 \n',
            '\n'])  # the trailing blank line is ignored

        params = char_uq.read_characterization_csv(_char(case_path))

        assert len(params) == 2
        assert params[0] == {'name': 'th_cond',
                             'distribution': 'Lognormal',
                             'x_low': 500., 'p_low': 0.5,
                             'x_high': 3000., 'p_high': 0.99}
        # surrounding whitespace is stripped from every field
        assert params[1]['name'] == 'capex'
        assert params[1]['distribution'] == 'Uniform'
        assert params[1]['x_high'] == 550.
    finally:
        _cleanup(case_path)


def test_read_characterization_csv_rejects_missing_file():
    with pytest.raises(FileNotFoundError):
        char_uq.read_characterization_csv(
            os.path.join(CASES_PATH, 'PYTEST_DOES_NOT_EXIST', 'nope.csv'))


def test_read_characterization_csv_rejects_missing_columns():
    case = 'PYTEST_CHAR_COLUMNS'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        _make_case(case)
        with open(_char(case_path), 'w') as file:
            file.write('name,distribution,x_low,p_low\n')
            file.write('x,Gaussian,10,0.1\n')

        with pytest.raises(ValueError, match='should contain the columns'):
            char_uq.read_characterization_csv(_char(case_path))
    finally:
        _cleanup(case_path)


def test_read_characterization_csv_rejects_incomplete_row():
    case = 'PYTEST_CHAR_INCOMPLETE'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        _make_case(case, characterization=['x,Gaussian,10,0.1,90,\n'])
        with pytest.raises(ValueError, match='incomplete'):
            char_uq.read_characterization_csv(_char(case_path))
    finally:
        _cleanup(case_path)


def test_read_characterization_csv_rejects_non_numeric_value():
    case = 'PYTEST_CHAR_NAN'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        _make_case(case, characterization=['x,Gaussian,low,0.1,90,0.9\n'])
        with pytest.raises(ValueError, match='not a number'):
            char_uq.read_characterization_csv(_char(case_path))
    finally:
        _cleanup(case_path)


def test_read_characterization_csv_rejects_unknown_distribution():
    case = 'PYTEST_CHAR_DIST'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        _make_case(case, characterization=['x,Triangular,10,0.1,90,0.9\n'])
        with pytest.raises(ValueError, match='unknown distribution'):
            char_uq.read_characterization_csv(_char(case_path))
    finally:
        _cleanup(case_path)


def test_read_characterization_csv_rejects_empty_template():
    """The shipped template holds only a header -- running it unfilled
    must say so rather than silently doing nothing."""
    case = 'PYTEST_CHAR_EMPTY'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        _make_case(case, characterization=[])
        with pytest.raises(ValueError, match='No parameters found'):
            char_uq.read_characterization_csv(_char(case_path))
    finally:
        _cleanup(case_path)


@pytest.mark.parametrize('case', ['NO_MODEL', 'REF'])
def test_shipped_templates_are_valid_headers(case):
    """The templates shipped in NO_MODEL and REF must carry exactly the
    expected columns, since users copy them into their own case."""
    path = os.path.join(CASES_PATH, case, char_uq.CHARACTERIZATION_FILE)
    assert os.path.isfile(path)

    with open(path, newline='') as file:
        reader = csv.DictReader(file)
        assert reader.fieldnames == char_uq.CHARACTERIZATION_COLUMNS


# --------------------------------------------------------------------------
# characterize_from_csv: the CSV-driven workflow.
# --------------------------------------------------------------------------

def test_characterize_from_csv_creates_spaces_when_absent():
    """A case folder holding only a characterization file gets both space
    files created from it."""
    case = 'PYTEST_CHAR_CSV_CREATE'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        _make_case(case, characterization=['x,Gaussian,10,0.1,90,0.9\n'])

        fitted = char_uq.characterize_from_csv(case)

        assert set(fitted) == {'x'}

        design_rows = _read_csv_rows(_design(case_path))
        assert len(design_rows) == 1
        assert design_rows[0]['name'] == 'x'
        assert design_rows[0]['type'] == 'par'
        assert design_rows[0]['upper_bound'] == ''
        assert float(design_rows[0]['value']) == pytest.approx(fitted['x'][0])

        stoch_rows = _read_csv_rows(_stoch(case_path))
        assert len(stoch_rows) == 1
        assert stoch_rows[0]['name'] == 'x'
        assert stoch_rows[0]['relation'] == 'absolute'
        assert stoch_rows[0]['distribution'] == 'Gaussian'
        assert float(
            stoch_rows[0]['deviation']) == pytest.approx(fitted['x'][1])
    finally:
        _cleanup(case_path)


def test_characterize_from_csv_adds_and_replaces():
    """New parameters are appended, already-characterized ones replaced,
    and untouched rows left alone."""
    case = 'PYTEST_CHAR_CSV_MERGE'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        _make_case(case,
                   design_rows=['keep_me,par,42,\n', 'x,par,1,\n'],
                   stoch_rows=['keep_me,absolute,Uniform,1\n',
                               'x,absolute,Uniform,0.5\n'],
                   characterization=['x,Gaussian,10,0.1,90,0.9\n',
                                     'brand_new,Uniform,1,0.25,3,0.75\n'])

        fitted = char_uq.characterize_from_csv(case)

        design_rows = _read_csv_rows(_design(case_path))
        assert [row['name'] for row in design_rows] == [
            'keep_me', 'x', 'brand_new']
        # the untouched row survives verbatim
        assert design_rows[0]['value'] == '42'
        # the existing one is replaced in place, keeping its position
        assert float(design_rows[1]['value']) == pytest.approx(fitted['x'][0])

        stoch_rows = _read_csv_rows(_stoch(case_path))
        assert [row['name'] for row in stoch_rows] == [
            'keep_me', 'x', 'brand_new']
        assert stoch_rows[0]['distribution'] == 'Uniform'
        assert stoch_rows[1]['distribution'] == 'Gaussian'
        assert stoch_rows[2]['distribution'] == 'Uniform'
    finally:
        _cleanup(case_path)


def test_characterize_from_csv_reports_added_and_updated(capsys):
    case = 'PYTEST_CHAR_CSV_REPORT'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        _make_case(case,
                   design_rows=['x,par,1,\n'],
                   stoch_rows=['x,absolute,Uniform,0.5\n'],
                   characterization=['x,Gaussian,10,0.1,90,0.9\n',
                                     'brand_new,Uniform,1,0.25,3,0.75\n'])

        char_uq.characterize_from_csv(case)

        out = capsys.readouterr().out
        assert 'added:   brand_new' in out
        assert 'updated: x' in out
    finally:
        _cleanup(case_path)


def test_characterize_from_csv_respects_overwrite_false():
    case = 'PYTEST_CHAR_CSV_NO_OVERWRITE'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        _make_case(case,
                   design_rows=['x,par,1,\n'],
                   stoch_rows=['x,absolute,Uniform,0.5\n'],
                   characterization=['x,Gaussian,10,0.1,90,0.9\n'])

        with pytest.raises(ValueError, match='already exist'):
            char_uq.characterize_from_csv(case, overwrite=False)

        # nothing was written
        stoch_rows = _read_csv_rows(_stoch(case_path))
        assert stoch_rows[0]['distribution'] == 'Uniform'
    finally:
        _cleanup(case_path)


def test_characterize_from_csv_write_false_touches_nothing():
    case = 'PYTEST_CHAR_CSV_PREVIEW'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        _make_case(case, characterization=['x,Gaussian,10,0.1,90,0.9\n'])

        fitted = char_uq.characterize_from_csv(case, write=False)

        assert set(fitted) == {'x'}
        assert not os.path.exists(_design(case_path))
        assert not os.path.exists(_stoch(case_path))
    finally:
        _cleanup(case_path)


def test_characterize_from_csv_requires_existing_case():
    case = 'PYTEST_CHAR_CSV_MISSING'
    _cleanup(os.path.join(CASES_PATH, case))

    with pytest.raises(FileNotFoundError):
        char_uq.characterize_from_csv(case)


def test_characterize_from_csv_rejects_duplicate_names():
    case = 'PYTEST_CHAR_CSV_DUP'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        _make_case(case, characterization=['x,Gaussian,10,0.1,90,0.9\n',
                                           'x,Uniform,1,0.25,3,0.75\n'])
        with pytest.raises(ValueError, match='Duplicate parameter name'):
            char_uq.characterize_from_csv(case)
    finally:
        _cleanup(case_path)


def test_characterize_from_csv_produces_a_runnable_uq_case():
    """The case built from a characterization file must be usable by
    RHEIA's own UQ pipeline, not just match the CSV schema."""
    case = 'PYTEST_CHAR_CSV_RUNNABLE'
    case_path = os.path.join(CASES_PATH, case)
    case_results_path = os.path.join(PACKAGE_ROOT, 'RESULTS', case)
    _cleanup(case_path)
    _cleanup(case_results_path)

    try:
        _make_case(case,
                   characterization=['x,Gaussian,10,0.1,90,0.9\n',
                                     'y,Lognormal,500,0.5,3000,0.99\n'])
        char_uq.characterize_from_csv(case)

        results_dir = 'pytest_char_csv_runnable'
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

        result_path = os.path.join(case_results_path, 'UQ', results_dir)
        assert os.path.isfile(os.path.join(result_path, 'samples.csv'))
    finally:
        _cleanup(case_path)
        _cleanup(case_results_path)


def test_characterize_from_csv_creates_case_description_when_absent():
    """A case folder holding only a characterization file still becomes
    loadable, via the NO_MODEL case_description template."""
    case = 'PYTEST_CHAR_CSV_BARE'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        os.makedirs(case_path)
        with open(_char(case_path), 'w') as file:
            file.write(HEADER)
            file.write('x,Gaussian,10,0.1,90,0.9\n')

        char_uq.characterize_from_csv(case)

        assert os.path.isfile(os.path.join(case_path, '__init__.py'))
        with open(os.path.join(case_path, 'case_description.py')) as file:
            case_description = file.read()
        with open(char_uq.NO_MODEL_CASE_DESCRIPTION) as file:
            assert case_description == file.read()
    finally:
        _cleanup(case_path)


# --------------------------------------------------------------------------
# characterize_parameters: the programmatic counterpart.
# --------------------------------------------------------------------------

def test_characterize_parameters_new_case_creates_loadable_case():
    case = 'PYTEST_CHARACTERIZE_NEW'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        fitted = char_uq.characterize_parameters(
            case,
            [{'name': 'x', 'distribution': 'Gaussian',
              'x_low': 10., 'p_low': 0.1, 'x_high': 90., 'p_high': 0.9}],
            new_case=True)

        assert set(fitted) == {'x'}

        assert os.path.isfile(os.path.join(case_path, '__init__.py'))
        with open(os.path.join(case_path, 'case_description.py')) as file:
            case_description = file.read()
        with open(char_uq.NO_MODEL_CASE_DESCRIPTION) as file:
            assert case_description == file.read()

        design_rows = _read_csv_rows(_design(case_path))
        assert len(design_rows) == 1
        assert design_rows[0]['name'] == 'x'
        assert design_rows[0]['type'] == 'par'
        assert float(design_rows[0]['value']) == pytest.approx(fitted['x'][0])

        stoch_rows = _read_csv_rows(_stoch(case_path))
        assert len(stoch_rows) == 1
        assert stoch_rows[0]['name'] == 'x'
        assert stoch_rows[0]['relation'] == 'absolute'
        assert stoch_rows[0]['distribution'] == 'Gaussian'
        assert float(
            stoch_rows[0]['deviation']) == pytest.approx(fitted['x'][1])

        # The produced case must be usable by RHEIA's own UQ pipeline.
        results_dir = 'pytest_characterize_new_case'
        case_results_path = os.path.join(PACKAGE_ROOT, 'RESULTS', case)
        result_path = os.path.join(case_results_path, 'UQ', results_dir)
        _cleanup(case_results_path)
        try:
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
            assert os.path.isfile(os.path.join(result_path, 'samples.csv'))
        finally:
            _cleanup(case_results_path)

    finally:
        _cleanup(case_path)


def test_characterize_parameters_new_case_rejects_existing_folder():
    case = 'PYTEST_CHARACTERIZE_EXISTS'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)
    os.makedirs(case_path)

    try:
        with pytest.raises(FileExistsError):
            char_uq.characterize_parameters(
                case,
                [{'name': 'x', 'distribution': 'Uniform',
                  'x_low': 1., 'p_low': 0.1, 'x_high': 2., 'p_high': 0.9}],
                new_case=True)
    finally:
        _cleanup(case_path)


def test_characterize_parameters_merges_without_disturbing_rows():
    case = 'PYTEST_CHARACTERIZE_MERGE'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        _make_case(case,
                   design_rows=['keep_me,par,42,\n'],
                   stoch_rows=['keep_me,absolute,Uniform,1\n'])

        char_uq.characterize_parameters(
            case,
            [{'name': 'new_param', 'distribution': 'Gaussian',
              'x_low': 10., 'p_low': 0.1, 'x_high': 90., 'p_high': 0.9}],
            new_case=False)

        design_rows = _read_csv_rows(_design(case_path))
        assert [row['name'] for row in design_rows] == [
            'keep_me', 'new_param']
        assert design_rows[0]['value'] == '42'

        stoch_rows = _read_csv_rows(_stoch(case_path))
        assert [row['name'] for row in stoch_rows] == ['keep_me', 'new_param']
        assert stoch_rows[0]['distribution'] == 'Uniform'
        assert stoch_rows[1]['distribution'] == 'Gaussian'

    finally:
        _cleanup(case_path)


def test_characterize_parameters_rejects_collision_without_overwrite():
    case = 'PYTEST_CHARACTERIZE_COLLISION'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    params = [{'name': 'x', 'distribution': 'Gaussian',
               'x_low': 10., 'p_low': 0.1, 'x_high': 90., 'p_high': 0.9}]

    try:
        _make_case(case,
                   design_rows=['x,par,1,\n'],
                   stoch_rows=['x,absolute,Uniform,0.5\n'])

        with pytest.raises(ValueError):
            char_uq.characterize_parameters(
                case, params, new_case=False, overwrite=False)

        char_uq.characterize_parameters(
            case, params, new_case=False, overwrite=True)

        assert len(_read_csv_rows(_design(case_path))) == 1
        stoch_rows = _read_csv_rows(_stoch(case_path))
        assert len(stoch_rows) == 1
        assert stoch_rows[0]['distribution'] == 'Gaussian'

    finally:
        _cleanup(case_path)


def test_characterize_parameters_merge_requires_existing_case():
    case = 'PYTEST_CHARACTERIZE_MISSING'
    _cleanup(os.path.join(CASES_PATH, case))

    with pytest.raises(FileNotFoundError):
        char_uq.characterize_parameters(
            case,
            [{'name': 'x', 'distribution': 'Uniform',
              'x_low': 1., 'p_low': 0.1, 'x_high': 2., 'p_high': 0.9}],
            new_case=False)


def test_characterize_parameters_rejects_duplicate_names():
    with pytest.raises(ValueError, match='Duplicate parameter name'):
        char_uq.characterize_parameters(
            'PYTEST_CHARACTERIZE_DUPLICATE',
            [{'name': 'x', 'distribution': 'Uniform',
              'x_low': 1., 'p_low': 0.1, 'x_high': 2., 'p_high': 0.9},
             {'name': 'x', 'distribution': 'Gaussian',
              'x_low': 1., 'p_low': 0.1, 'x_high': 2., 'p_high': 0.9}],
            write=False)


# --------------------------------------------------------------------------
# write=False: preview while iterating on the percentiles, without
# touching disk.
# --------------------------------------------------------------------------

def test_characterize_parameters_write_false_creates_no_new_case():
    case = 'PYTEST_CHARACTERIZE_PREVIEW_NEW'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        fitted = char_uq.characterize_parameters(
            case,
            [{'name': 'x', 'distribution': 'Gaussian',
              'x_low': 10., 'p_low': 0.1, 'x_high': 90., 'p_high': 0.9}],
            new_case=True, write=False)

        assert set(fitted) == {'x'}
        assert not os.path.exists(case_path)  # nothing created at all
    finally:
        _cleanup(case_path)


def test_characterize_parameters_write_false_leaves_case_untouched():
    case = 'PYTEST_CHARACTERIZE_PREVIEW_MERGE'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        _make_case(case,
                   design_rows=['keep_me,par,42,\n'],
                   stoch_rows=['keep_me,absolute,Uniform,1\n'])

        design_before = _read_csv_rows(_design(case_path))
        stoch_before = _read_csv_rows(_stoch(case_path))

        # Even a name collision with an existing parameter is fine in
        # preview mode -- nothing is written, so there's nothing to
        # silently clobber.
        fitted = char_uq.characterize_parameters(
            case,
            [{'name': 'keep_me', 'distribution': 'Gaussian',
              'x_low': 10., 'p_low': 0.1, 'x_high': 90., 'p_high': 0.9}],
            new_case=False, write=False)

        assert set(fitted) == {'keep_me'}
        assert _read_csv_rows(_design(case_path)) == design_before
        assert _read_csv_rows(_stoch(case_path)) == stoch_before
    finally:
        _cleanup(case_path)


def test_characterize_parameters_write_false_plot_show_does_not_raise():
    """plot=True/show=True must render and close without blocking (Agg
    backend, see module top) and without writing anything."""
    case = 'PYTEST_CHARACTERIZE_PREVIEW_PLOT'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        char_uq.characterize_parameters(
            case,
            [{'name': 'x', 'distribution': 'Lognormal',
              'x_low': 500., 'p_low': 0.50,
              'x_high': 3000., 'p_high': 0.99}],
            new_case=True, write=False, plot=True, show=True)

        assert not os.path.exists(case_path)
    finally:
        _cleanup(case_path)


def test_characterize_parameters_write_false_then_true_agree():
    """Previewing must not change what actually gets written once
    write=True -- same inputs, same fit."""
    case = 'PYTEST_CHARACTERIZE_PREVIEW_THEN_WRITE'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    params = [{'name': 'x', 'distribution': 'Lognormal',
               'x_low': 500., 'p_low': 0.50,
               'x_high': 3000., 'p_high': 0.99}]

    try:
        previewed = char_uq.characterize_parameters(
            case, params, new_case=True, write=False)
        written = char_uq.characterize_parameters(
            case, params, new_case=True, write=True)

        assert previewed == written
    finally:
        _cleanup(case_path)


def test_characterize_parameters_saves_plot_when_writing():
    """plot=True saves a '<name>_fit.png' next to the written spaces."""
    case = 'PYTEST_CHARACTERIZE_PLOT_WRITE'
    case_path = os.path.join(CASES_PATH, case)
    _cleanup(case_path)

    try:
        char_uq.characterize_parameters(
            case,
            [{'name': 'x', 'distribution': 'Gaussian',
              'x_low': 10., 'p_low': 0.1, 'x_high': 90., 'p_high': 0.9}],
            new_case=True, plot=True)

        assert os.path.isfile(os.path.join(case_path, 'x_fit.png'))
    finally:
        _cleanup(case_path)


def test_plot_fit_saves_file_for_each_distribution(tmp_path):
    for distribution in ['Uniform', 'Gaussian', 'Lognormal']:
        value, deviation = char_uq.fit_distribution(
            distribution, 500., 0.5, 3000., 0.99)
        path = os.path.join(str(tmp_path), '%s.png' % distribution)
        char_uq.plot_fit(distribution, value, deviation,
                         500., 0.5, 3000., 0.99, path=path)
        assert os.path.isfile(path)


def test_plot_fit_rejects_unknown_distribution():
    with pytest.raises(ValueError):
        char_uq.plot_fit('Triangular', 1., 1., 0.5, 0.1, 1.5, 0.9)
