"""
The :py:mod:`case_description` module contains a function
to evaluate the four bar truss system model.
"""

import rheia.CASES.ENERGYSCOPE.uq_estd as es
import energyscope as es2
import os
import yaml
import rheia

def evaluate(x_in, params=[]):
    '''
    Evaluation of the system objectives for one given design.

    Parameters
    ----------
    x_in : tuple
        An enumerate object for the input sample.
        The first element of x
        - the index of the input sample in the list of samples -
        can be used for multiprocessing purposes of executable files
        with input and output text files.
        The second element of x - the input sample -
        is a dictionary with the names and values for the model parameters
        and design variables.
    params : list, optional
        List with fixed data, used during model evaluation. The default is [].

    Returns
    -------
    total_cost : float
        cost for the sample [euro]
    '''

    path_rheia = os.path.dirname(rheia.__file__)    
    
    path_es = os.path.dirname(os.path.dirname(es2.__file__))

    # the config file considered
    config_file = os.path.join(path_rheia,
                               'CASES',
                               'ENERGYSCOPE',
                               'config_ref_UQ.yaml')

    with open(config_file, 'r') as f:
        lines = yaml.load(f, Loader=yaml.Loader)
        
    # the config file considered
    config_file2 = os.path.join(path_rheia,
                               'CASES',
                               'ENERGYSCOPE',
                               'config_ref_UQ.yaml')

    lines['case_study'] = '%s/Run_%s' %(params[0]['results dir'],x_in[0])

    with open(config_file2, 'w') as file:
        documents = yaml.dump(lines, file)
    
    print(x_in)
    
    # evaluate the energyscope model
    total_cost = es.run_ESTD_UQ([x_in, params[0]['results dir']])
    

    return total_cost