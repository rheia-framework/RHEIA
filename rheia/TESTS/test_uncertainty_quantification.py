'''
Module to test the uncertainty_quantification module.
'''

import os
import pytest
import pandas as pd
from rheia.CASES.determine_stoch_des_space import (
    DESIGN_SPACE_COLUMNS, check_dictionary)
import rheia.UQ.uncertainty_quantification as uq


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


def test_get_design_variables(run_dict):
    """
    Assert the design variable dictionary.

    Parameters
    ----------
    run_dict : dict
        Input dictionary for uncertainty quantification.

    """

    var_dict = uq.get_design_variables(run_dict['case'])
    assert var_dict == {'n_pemel': [0., 7.],
                        'n_dcdc_pv': [0., 7.]}


def test_set_design_samples(run_dict):
    """
    Assert the design samples.

    Parameters
    ----------
    run_dict : dict
        Input dictionary for uncertainty quantification.

    """

    var_dict = uq.get_design_variables(run_dict['case'])
    samples = uq.set_design_samples(var_dict, 4)
    assert len(samples) == 4
    assert len(samples[0]) == 2


def test_draw_pdf_cdf_accepts_integral_float(run_dict):
    """
    Assert scientific notation floats are accepted for draw sample counts.
    """

    run_dict.update({
        'uq method': 'sparse',
        'n samples': 20,
        'draw pdf cdf': [True, 1e5],
    })

    check_dictionary(run_dict, uq_bool=True)

    assert run_dict['draw pdf cdf'][1] == 100000


def test_draw_pdf_cdf_rejects_non_integral_float(run_dict):
    """
    Assert fractional sample counts are rejected.
    """

    run_dict.update({
        'uq method': 'sparse',
        'n samples': 20,
        'draw pdf cdf': [True, 1.5],
    })

    with pytest.raises(TypeError):
        check_dictionary(run_dict, uq_bool=True)


def test_write_design_space_keeps_header(run_dict):
    """
    Assert that generated design_space files use the named schema.
    """

    var_dict = uq.get_design_variables(run_dict['case'])
    ds = 'design_space_test_generated.csv'
    path = os.path.join(os.path.dirname(__file__), os.pardir, 'CASES',
                        run_dict['case'], 'design_space_test_generated_0.csv')

    if os.path.isfile(path):
        os.remove(path)

    try:
        uq.write_design_space(run_dict['case'], 0, var_dict, [1., 2.], ds=ds)
        data = pd.read_csv(path, dtype=str, keep_default_na=False,
                           na_filter=False)

        assert list(data.columns) == DESIGN_SPACE_COLUMNS
        assert data.loc[data['name'] == 'n_dcdc_pv', 'type'].iloc[0] == 'par'
        assert data.loc[data['name'] == 'n_dcdc_pv', 'upper_bound'].iloc[0] == ''
    finally:
        if os.path.isfile(path):
            os.remove(path)
