"""
The :py:mod:`case_description` module contains a function
to read and store the fixed parameters for model evaluation and a function
to evaluate the system model.
"""


def set_params():
    """
    Set the fixed parameters for each model evaluation.

    """

    params = []

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
    y_1 : float
        model output 1
    y_2 : float
        model output 2
    '''

    # provide your system model and pass x_in as an argument
    y_1 = 1.
    y_2 = 1.

    # return the objectives
    return y_1, y_2
