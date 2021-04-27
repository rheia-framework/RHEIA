"""
The :py:mod:`case_description` module contains a function
to evaluate the EnergyPLAN system model.
"""

from rheia.CASES.ENERGYPLAN.run_energyplan import energyplan


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
    co2 : float
        total CO2-emission
    fuel : float
        primary energy consumption
    '''

    # evaluate EnergyPLAN model
    co2, fuel = energyplan(x_in)

    return co2, fuel
