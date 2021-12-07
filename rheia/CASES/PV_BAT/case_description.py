"""
The :py:mod:`case_description` module contains a function
to read and store the fixed parameters for model evaluation and a function
to evaluate the system model.
"""

import os
import rheia.CASES.PV_BAT.pv_bat as pv_bat


def set_params():
    """
    Set the fixed parameters for each model evaluation.

    """

    path = os.path.dirname(os.path.abspath(__file__))

    filename_climate = os.path.join(os.path.abspath(
                                    os.path.join(path,
                                                 os.pardir)),
                                    'DATA',
                                    'climate',
                                    'climate_Brussels_1day.csv')

    filename_demand = os.path.join(os.path.abspath(
                                   os.path.join(path,
                                                os.pardir)),
                                   'DATA',
                                   'demand',
                                   'load_Brussels_dwelling_1day.csv')

    # get the solar irradiance, ambient temperature data and electric demand
    my_data = pv_bat.ReadData(filename_climate, filename_demand)
    sol_irr, t_amb = my_data.load_climate()
    load_elec = my_data.load_demand()

    # store data in params list
    params = [sol_irr, t_amb, load_elec]

    return params


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
    lcoe : float
        the levelized cost of electricity
    ssr: float
        the self-sufficiency ratio
    '''

    arguments = params + [False] + [x_in[1]]

    # Evaluation object
    my_evaluation = pv_bat.Evaluation(*arguments)

    # evaluate system model
    my_evaluation.evaluation()

    # get values for objectives
    lcoe = my_evaluation.lcoe
    ssr = my_evaluation.ssr

    return lcoe, ssr
