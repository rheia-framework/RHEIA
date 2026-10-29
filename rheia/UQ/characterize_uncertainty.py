"""
The :py:mod:`characterize_uncertainty` module characterizes an uncertain
parameter's distribution from two percentiles (Cook's method), and writes
the result straight into a case's :file:`design_space.csv` and
:file:`stochastic_space.csv`.

Two percentiles pin down a two-parameter distribution exactly, with no other
data needed -- e.g. a declared value at p50 and a design value at p99 from
a standard or specification. This covers the three distributions RHEIA's PCE
supports (Uniform, Gaussian, Lognormal); for a Lognormal input the fit
follows Cook, J.D. (2010), "Determining distribution parameters from
quantiles".

The percentile estimates are collected in a
:file:`uncertainty_characterization.csv` file in the case folder, with one
row per uncertain parameter::

    name,distribution,x_low,p_low,x_high,p_high
    th_cond,Lognormal,500,0.5,3000,0.99

:py:func:`characterize_from_csv` reads that file, fits every row, and merges
the result into the case's :file:`design_space.csv` and
:file:`stochastic_space.csv`. See
:file:`examples/characterize_uncertainty.py` for a worked example, and
:py:func:`characterize_parameters` to pass the estimates programmatically
instead of through a file.
"""

import csv
import math
import os
from statistics import NormalDist

import numpy as np

from rheia.CASES.determine_stoch_des_space import (
    DESIGN_SPACE_COLUMNS, STOCHASTIC_SPACE_COLUMNS)
from rheia.UQ.pce import _lognormal_to_normal

_NORMAL = NormalDist()

CASES_DIR = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, 'CASES'))
NO_MODEL_CASE_DESCRIPTION = os.path.join(CASES_DIR, 'NO_MODEL',
                                         'case_description.py')

CHARACTERIZATION_FILE = 'uncertainty_characterization.csv'
CHARACTERIZATION_COLUMNS = ['name', 'distribution',
                            'x_low', 'p_low', 'x_high', 'p_high']


def _validate_quantiles(x_low, p_low, x_high, p_high):
    if not 0. < p_low < p_high < 1.:
        raise ValueError(
            "require 0 < p_low < p_high < 1 (got p_low=%r, p_high=%r)"
            % (p_low, p_high))
    if not x_low < x_high:
        raise ValueError(
            "require x_low < x_high (got x_low=%r, x_high=%r)"
            % (x_low, x_high))


def _fit_gaussian(x_low, p_low, x_high, p_high):
    """Two percentiles of a Gaussian pin down its mean/std exactly:
    x = mean + std*z, where z = Phi^-1(p)."""
    z_low, z_high = _NORMAL.inv_cdf(p_low), _NORMAL.inv_cdf(p_high)

    std = (x_high - x_low) / (z_high - z_low)
    mean = (x_low * z_high - x_high * z_low) / (z_high - z_low)

    return mean, std


def _fit_lognormal(x_low, p_low, x_high, p_high):
    """Cook's method: fit a Gaussian to ln(x), the same way
    _fit_gaussian() does, then convert the resulting latent (mu, sigma) to
    the physical mean/std RHEIA's stochastic_space.csv expects for a
    Lognormal row -- the mathematical inverse of
    rheia.UQ.pce._lognormal_to_normal(), so constructing the PCE from
    this fit's (mean, deviation) reproduces the same (mu, sigma) here."""
    if x_low <= 0. or x_high <= 0.:
        raise ValueError("a lognormal requires strictly positive percentiles")

    mu, sigma = _fit_gaussian(math.log(x_low), p_low, math.log(x_high), p_high)

    mean = math.exp(mu + sigma ** 2 / 2.)
    deviation = mean * math.sqrt(math.exp(sigma ** 2) - 1.)

    return mean, deviation


def _fit_uniform(x_low, p_low, x_high, p_high):
    """A Uniform's CDF is linear, so the two percentiles interpolate the
    bounds directly: x = a + p*(b-a)."""
    width = (x_high - x_low) / (p_high - p_low)
    a = x_low - p_low * width
    b = a + width

    mean = (a + b) / 2.
    deviation = (b - a) / 2.

    return mean, deviation


_FIT_FUNCTIONS = {
    'Uniform': _fit_uniform,
    'Gaussian': _fit_gaussian,
    'Lognormal': _fit_lognormal,
}


def fit_distribution(distribution, x_low, p_low, x_high, p_high):
    """
    Fit `distribution` to two percentiles (x_low at p_low, x_high at p_high).

    Parameters
    ----------
    distribution : str
        One of 'Uniform', 'Gaussian', 'Lognormal' -- the same
        case-sensitive names stochastic_space.csv itself uses.
    x_low, x_high : float
        The two elicited values, x_low < x_high.
    p_low, p_high : float
        The percentiles (in (0, 1)) those values represent, p_low < p_high.

    Returns
    -------
    value : float
        The design_space.csv value for this parameter (the distribution's
        mean).
    deviation : float
        The stochastic_space.csv deviation for this parameter, under
        relation='absolute' (see fit function docstrings for what
        'deviation' means per distribution).
    """
    if distribution not in _FIT_FUNCTIONS:
        raise ValueError(
            "Unknown distribution '%s'. Choose between %s."
            % (distribution, ', '.join(_FIT_FUNCTIONS)))

    _validate_quantiles(x_low, p_low, x_high, p_high)

    return _FIT_FUNCTIONS[distribution](x_low, p_low, x_high, p_high)


def _read_csv_rows(path, columns):
    with open(path, newline='') as file:
        return [{col: row[col] for col in columns}
                for row in csv.DictReader(file)]


def _write_csv_rows(path, columns, rows):
    with open(path, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=columns, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def read_characterization_csv(path):
    """
    Read a :file:`uncertainty_characterization.csv` file into the list of
    parameter dictionaries that :py:func:`characterize_parameters` expects.

    Parameters
    ----------
    path : str
        Path to the characterization file. It should contain the columns
        'name', 'distribution', 'x_low', 'p_low', 'x_high' and 'p_high',
        with one uncertain parameter per row.

    Returns
    -------
    parameters : list of dict
        One dictionary per row, with 'name' and 'distribution' as strings
        and 'x_low', 'p_low', 'x_high', 'p_high' as floats.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(
            'Missing file: "%s" not found. Create it with the columns: %s.'
            % (path, ', '.join(CHARACTERIZATION_COLUMNS)))

    with open(path, newline='') as file:
        reader = csv.DictReader(file)
        fieldnames = reader.fieldnames or []
        missing_columns = [column for column in CHARACTERIZATION_COLUMNS
                           if column not in fieldnames]
        if missing_columns:
            raise ValueError(
                'The characterization file %s should contain the columns: %s.'
                % (path, ', '.join(CHARACTERIZATION_COLUMNS)))
        raw_rows = list(reader)

    parameters = []
    for index, raw in enumerate(raw_rows):
        values = {column: (raw[column] or '').strip()
                  for column in CHARACTERIZATION_COLUMNS}

        # skip fully empty lines, so a trailing newline is not an error
        if not any(values.values()):
            continue

        empty = [column for column, value in values.items() if value == '']
        if empty:
            raise ValueError(
                'Row %d of %s is incomplete: %s has no value.'
                % (index + 2, path, ', '.join(empty)))

        param = {'name': values['name'],
                 'distribution': values['distribution']}
        for column in ['x_low', 'p_low', 'x_high', 'p_high']:
            try:
                param[column] = float(values[column])
            except ValueError:
                raise ValueError(
                    'Row %d of %s: "%s" is not a number for %s.'
                    % (index + 2, path, values[column], column)) from None

        if param['distribution'] not in _FIT_FUNCTIONS:
            raise ValueError(
                "Row %d of %s: unknown distribution '%s' for %s. "
                "Choose between %s."
                % (index + 2, path, param['distribution'], param['name'],
                   ', '.join(_FIT_FUNCTIONS)))

        parameters.append(param)

    if not parameters:
        raise ValueError(
            'No parameters found in %s. Add one row per uncertain parameter, '
            'below the header line.' % path)

    return parameters


def plot_fit(distribution, value, deviation, x_low, p_low, x_high, p_high,
             path=None, show=False):
    """Plot the fitted PDF with the two elicited percentiles marked, as a
    sanity check on the fit.

    Parameters
    ----------
    distribution : str
        One of 'Uniform', 'Gaussian', 'Lognormal'.
    value, deviation : float
        The fitted mean and deviation, as returned by
        :py:func:`fit_distribution`.
    x_low, p_low, x_high, p_high : float
        The two elicited percentiles, marked on the plot.
    path : str, optional
        Where to save the plot (e.g. ending in .png or .pdf). Not saved
        when None (the default) -- useful together with show=True while
        still iterating on x_low/p_low/x_high/p_high, before a case
        directory necessarily exists to save into.
    show : bool, optional
        If True, display the plot interactively
        (matplotlib.pyplot.show(), which blocks until the window is
        closed). The default is False.
    """
    import matplotlib.pyplot as plt
    from scipy import stats

    if distribution == 'Uniform':
        dist = stats.uniform(loc=value - deviation, scale=2. * deviation)
    elif distribution == 'Gaussian':
        dist = stats.norm(loc=value, scale=deviation)
    elif distribution == 'Lognormal':
        mu, sigma = _lognormal_to_normal(value, deviation)
        dist = stats.lognorm(s=sigma, scale=math.exp(mu))
    else:
        raise ValueError(
            "Unknown distribution '%s'. Choose between %s."
            % (distribution, ', '.join(_FIT_FUNCTIONS)))

    x = np.linspace(dist.ppf(0.001), dist.ppf(0.999), 500)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(x, dist.pdf(x), color='black')
    for x_q, p_q in ((x_low, p_low), (x_high, p_high)):
        ax.axvline(x_q, color='#2a78d6', linestyle='--', linewidth=1)
        ax.annotate('p%g' % (p_q * 100), (x_q, 0), xytext=(3, 6),
                    textcoords='offset points', color='#2a78d6', fontsize=9)
    ax.set_xlabel('value')
    ax.set_ylabel('probability density')
    ax.set_title('%s fit: mean=%.4g, deviation=%.4g'
                 % (distribution, value, deviation))
    fig.tight_layout()
    if path is not None:
        fig.savefig(path)
    if show:
        plt.show()
    plt.close(fig)


def _fit_all(parameters):
    """Fit every parameter, after checking the names are unique."""
    names = [param['name'] for param in parameters]
    if len(names) != len(set(names)):
        raise ValueError(
            "Duplicate parameter name(s): %s"
            % ', '.join(sorted(name for name in set(names)
                               if names.count(name) > 1)))

    return {param['name']: fit_distribution(
        param['distribution'], param['x_low'], param['p_low'],
        param['x_high'], param['p_high']) for param in parameters}


def _create_case_folder(case_dir):
    """Create a case folder with an empty __init__.py and a
    case_description.py copied from the NO_MODEL template, so the case is
    immediately loadable with 'create only samples': True."""
    os.makedirs(case_dir, exist_ok=True)

    init_path = os.path.join(case_dir, '__init__.py')
    if not os.path.isfile(init_path):
        open(init_path, 'w').close()

    description_path = os.path.join(case_dir, 'case_description.py')
    if not os.path.isfile(description_path):
        with open(NO_MODEL_CASE_DESCRIPTION) as src:
            template = src.read()
        with open(description_path, 'w') as dst:
            dst.write(template)


def _merge_into_case(case_dir, parameters, fitted, overwrite, plot, show):
    """Merge the fitted parameters into the case's design_space.csv and
    stochastic_space.csv, creating either file when it does not exist yet.

    Returns (added, updated): the parameter names that were appended and
    the ones that replaced an existing characterization.
    """
    design_space_path = os.path.join(case_dir, 'design_space.csv')
    stochastic_space_path = os.path.join(case_dir, 'stochastic_space.csv')

    if os.path.isfile(design_space_path):
        design_rows = _read_csv_rows(design_space_path, DESIGN_SPACE_COLUMNS)
    else:
        design_rows = []
    if os.path.isfile(stochastic_space_path):
        stoch_rows = _read_csv_rows(stochastic_space_path,
                                    STOCHASTIC_SPACE_COLUMNS)
    else:
        stoch_rows = []

    design_by_name = {row['name']: row for row in design_rows}
    stoch_by_name = {row['name']: row for row in stoch_rows}

    names = [param['name'] for param in parameters]
    colliding = set(design_by_name) & set(names)
    if colliding and not overwrite:
        raise ValueError(
            "Parameter(s) %s already exist in %s -- pass overwrite=True to "
            "update them, or rename."
            % (', '.join(sorted(colliding)), design_space_path))

    added, updated = [], []

    for param in parameters:
        name = param['name']
        distribution = param['distribution']
        value, deviation = fitted[name]

        design_row = {'name': name, 'type': 'par',
                      'value': repr(value), 'upper_bound': ''}
        stoch_row = {'name': name, 'relation': 'absolute',
                     'distribution': distribution,
                     'deviation': repr(deviation)}

        if name in design_by_name:
            design_by_name[name].update(design_row)
            updated.append(name)
        else:
            design_rows.append(design_row)
            design_by_name[name] = design_row
            added.append(name)

        if name in stoch_by_name:
            stoch_by_name[name].update(stoch_row)
        else:
            stoch_rows.append(stoch_row)
            stoch_by_name[name] = stoch_row

        if plot:
            plot_fit(distribution, value, deviation,
                     param['x_low'], param['p_low'],
                     param['x_high'], param['p_high'],
                     path=os.path.join(case_dir, '%s_fit.png' % name),
                     show=show)

    _write_csv_rows(design_space_path, DESIGN_SPACE_COLUMNS, design_rows)
    _write_csv_rows(stochastic_space_path, STOCHASTIC_SPACE_COLUMNS,
                    stoch_rows)

    return added, updated


def _report(case, case_dir, added, updated):
    """Print which parameters were appended and which replaced an existing
    characterization, so a replacement is never silent."""
    print("Characterized %d parameter(s) in case '%s' (%s):"
          % (len(added) + len(updated), case, case_dir))
    if added:
        print('  added:   %s' % ', '.join(added))
    if updated:
        print('  updated: %s' % ', '.join(updated))
    print('  written to design_space.csv and stochastic_space.csv')


def characterize_from_csv(case, filename=CHARACTERIZATION_FILE,
                          overwrite=True, write=True, plot=False, show=False,
                          verbose=True):
    """
    Fit every uncertain parameter listed in a case's
    :file:`uncertainty_characterization.csv` from its two elicited
    percentiles, and write the result into that case's
    :file:`design_space.csv` and :file:`stochastic_space.csv`.

    The case folder has to exist and contain the characterization file;
    :file:`design_space.csv` and :file:`stochastic_space.csv` do not --
    they are created when absent, and merged into when present. A
    parameter that is not in those files yet is appended; one that is
    already there is replaced, so re-running after updating an estimate
    simply refreshes it. Rows in those files that the characterization
    file does not mention are left untouched.

    Parameters
    ----------
    case : str
        The case name (the folder name under CASES).
    filename : str, optional
        The characterization file to read from the case folder. The
        default is 'uncertainty_characterization.csv'. A path to a file
        outside the case folder is also accepted.
    overwrite : bool, optional
        Whether a parameter that already exists in the case may be
        replaced. The default is True, since the characterization file is
        the source of truth for the parameters it lists; the replacement
        is reported rather than silent. Pass False to raise on a name that
        already exists instead.
    write : bool, optional
        If True (the default), write the fit into the case's
        :file:`design_space.csv` and :file:`stochastic_space.csv`. If
        False, only fit -- and, if plot or show is set, render -- the
        result; nothing is written to disk, so the same call can be
        repeated freely while tuning the percentiles.
    plot : bool, optional
        If True, render the fitted probability density function with the
        two elicited percentiles marked (see :py:func:`plot_fit`) for
        every parameter. Saved as '<name>_fit.png' in the case folder when
        write=True; not saved when write=False (pass show=True to see it
        instead).
    show : bool, optional
        If True, display each parameter's plot interactively (blocks on
        matplotlib.pyplot.show() until the window is closed). Independent
        of whether the plot is also saved to disk.
    verbose : bool, optional
        If True (the default), print which parameters were added and which
        were updated.

    Returns
    -------
    dict
        parameter name -> (value, deviation), the fitted values.
    """
    case_dir = os.path.join(CASES_DIR, case)
    path = filename if os.path.isabs(filename) else os.path.join(case_dir,
                                                                 filename)

    if not os.path.isdir(case_dir):
        raise FileNotFoundError(
            "Case '%s' does not exist at %s. Create the case folder and put "
            "a %s file in it (see CASES/NO_MODEL/%s for a template)."
            % (case, case_dir, CHARACTERIZATION_FILE, CHARACTERIZATION_FILE))

    parameters = read_characterization_csv(path)
    fitted = _fit_all(parameters)

    if not write:
        if plot or show:
            for param in parameters:
                value, deviation = fitted[param['name']]
                plot_fit(param['distribution'], value, deviation,
                         param['x_low'], param['p_low'],
                         param['x_high'], param['p_high'],
                         path=None, show=show)
        if verbose:
            print('write=False -- nothing written to disk. '
                  'Set write=True once the fit looks right.')
        return fitted

    _create_case_folder(case_dir)
    added, updated = _merge_into_case(case_dir, parameters, fitted,
                                      overwrite, plot, show)

    if verbose:
        _report(case, case_dir, added, updated)

    return fitted


def characterize_parameters(case, parameters, new_case=False, overwrite=False,
                            write=True, plot=False, show=False, verbose=True):
    """
    Fit each parameter in `parameters` from its two elicited percentiles,
    and, once satisfied with the fit, write the result into `case`'s
    :file:`design_space.csv` and :file:`stochastic_space.csv`.

    This is the programmatic counterpart of
    :py:func:`characterize_from_csv`, taking the percentile estimates
    directly instead of reading them from a file.

    Picking good percentiles is usually iterative -- call this a few times
    with write=False (and plot=True, show=True) to see how the fit reacts
    to x_low/p_low/x_high/p_high before committing anything to disk; the
    fitted values are always returned regardless of write, so they can be
    inspected either way.

    Parameters
    ----------
    case : str
        The case name (the folder name under CASES).
    parameters : list of dict
        One dict per uncertain parameter, each with keys 'name',
        'distribution' ('Uniform', 'Gaussian' or 'Lognormal'), 'x_low',
        'p_low', 'x_high', 'p_high'. Every parameter is written as a model
        parameter ('par' in design_space.csv), not a design variable.
    new_case : bool, optional
        If True, create a new case folder (CASES/<case>) containing only
        these parameters, plus a case_description.py copied from the
        NO_MODEL template (no model connected -- works immediately with
        'create only samples': True; fill in evaluate()/set_params() later
        to connect a real model). Raises if the case folder already
        exists. If False (the default), the case folder and both CSVs
        must already exist, and the parameters are merged into them.
        Ignored when write=False (nothing is created either way).
    overwrite : bool, optional
        When merging into an existing case (new_case=False, write=True),
        whether a parameter name that's already present should be updated
        in place. Raises on a name collision when False (the default), to
        avoid silently discarding an existing characterization.
    write : bool, optional
        If True (the default), write the fit into the case's
        :file:`design_space.csv` and :file:`stochastic_space.csv`
        (creating the case first when new_case=True). If False, only fit
        -- and, if plot or show is set, render -- the result; nothing is
        written to disk and the case does not need to exist yet.
    plot : bool, optional
        If True, render the fitted probability density function with the
        two elicited percentiles marked (see :py:func:`plot_fit`) for
        every parameter. Saved as '<name>_fit.png' in the case folder when
        write=True; not saved when write=False (pass show=True to see it
        instead, since there may be no case directory yet to save into).
    show : bool, optional
        If True, display each parameter's plot interactively (blocks on
        matplotlib.pyplot.show() until the window is closed). Independent
        of whether the plot is also saved to disk.
    verbose : bool, optional
        If True (the default), print which parameters were added and which
        were updated.

    Returns
    -------
    dict
        parameter name -> (value, deviation), the fitted values.
    """
    fitted = _fit_all(parameters)

    if not write:
        if plot or show:
            for param in parameters:
                value, deviation = fitted[param['name']]
                plot_fit(param['distribution'], value, deviation,
                         param['x_low'], param['p_low'],
                         param['x_high'], param['p_high'],
                         path=None, show=show)
        return fitted

    case_dir = os.path.join(CASES_DIR, case)
    design_space_path = os.path.join(case_dir, 'design_space.csv')
    stochastic_space_path = os.path.join(case_dir, 'stochastic_space.csv')

    if new_case:
        if os.path.exists(case_dir):
            raise FileExistsError(
                "Case '%s' already exists at %s -- pass new_case=False to "
                "merge into it instead." % (case, case_dir))
        _create_case_folder(case_dir)
    else:
        if not (os.path.isfile(design_space_path) and
                os.path.isfile(stochastic_space_path)):
            raise FileNotFoundError(
                "Case '%s' (or its design_space.csv/stochastic_space.csv) "
                "does not exist at %s -- pass new_case=True to create it."
                % (case, case_dir))

    added, updated = _merge_into_case(case_dir, parameters, fitted,
                                      overwrite, plot, show)

    if verbose:
        _report(case, case_dir, added, updated)

    return fitted
