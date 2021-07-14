'''
Module to test the uncertainty_quantification module.
'''

import pytest
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
