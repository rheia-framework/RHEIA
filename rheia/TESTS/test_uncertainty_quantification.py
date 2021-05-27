import pytest
import numpy as np
import rheia.UQ.uncertainty_quantification as uq

@pytest.fixture
def run_dict():
    dict_opt = {'case':                'H2_FUEL',
                'objectives':          {'DET': (-1, 1)},
                'stop':                ('BUDGET', 2000),
                'n jobs':              1,
                'population size':     20,
                'results dir':         'run_tutorial',
                }
                
    
    
    return dict_opt

def test_get_design_variables(run_dict):
    var_dict = uq.get_design_variables(run_dict['case'])
    assert var_dict == {'n_pemel': [0.,7.],
                        'n_dcdc_pv': [0.,7.]}

def test_set_design_samples(run_dict):
    var_dict = uq.get_design_variables(run_dict['case'])
    samples = uq.set_design_samples(var_dict, 4)
    assert len(samples) == 4
    assert len(samples[0]) == 2
