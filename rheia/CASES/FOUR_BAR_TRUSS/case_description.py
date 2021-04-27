"""
The :py:mod:`case_description` module contains a function
to evaluate the four bar truss system model.
"""

from rheia.CASES.FOUR_BAR_TRUSS.four_bar_truss import four_bar_truss

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
    vol : float
        truss volume
    disp : float
        displacement of outer node
    '''

    # evaluate four bar truss model
    vol, disp = four_bar_truss(x_in[1])

    return vol, disp
