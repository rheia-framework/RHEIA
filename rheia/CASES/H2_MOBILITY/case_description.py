"""
The :py:mod:`case_description` module contains a function
to read and store the fixed parameters for model evaluation and a function
to evaluate the system model.
"""

import os
import rheia.CASES.H2_MOBILITY.h2_mobility as lb


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
                                    'climate_Brussels.csv')

    # get the solar irradiance and ambient temperature data
    my_data = lb.ReadData(filename_climate)
    sol_irr, t_amb = my_data.load_climate()

    # store data in params list
    params = [sol_irr, t_amb]

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
    lcom : float
        the levelized cost of driving
    lco2 : float
        the carbon intensity
    grid_e_bought : float
        the annual grid electricity bought
    '''

    arguments = params + [x_in[1]]

    # Evaluation object
    my_evaluation = lb.Evaluation(*arguments)

    # evaluate system model
    succes = my_evaluation.evaluation()
    if succes:

        # when constrained satisfied, get values for objectives
        lcom = my_evaluation.res['lcom']
        lco2 = my_evaluation.res['lco2']
        grid_e_bought = sum(my_evaluation.res['grid_e_buy'])/1e6
    else:

        # when constrained violated, return large value for objectives
        lcom = 1e8
        lco2 = 1e8
        grid_e_bought = 1e8

    return lcom, lco2, grid_e_bought
