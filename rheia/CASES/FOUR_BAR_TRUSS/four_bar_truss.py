"""
The :py:mod:`four_bar_truss` module contains a function
to evaluate the four bar truss system model.
"""

def four_bar_truss(x_i):
    """

    This function evaluates the volume and displacement
    of a four-bar truss system.

    Parameters
    ----------
    x_i : array
        The input sample.

    Returns
    -------
    vol : float
        The volume of the truss.
    disp : float
        The displacement of the node.

    """

    # the volume
    vol = x_i['L'] * (2. * x_i['A_1'] + 2.**(0.5) *
                       x_i['A_2'] + x_i['A_3']**(0.5) + x_i['A_4'])

    # the displacement of the outer node
    disp = x_i['F'] * x_i['L'] * (2. / (x_i['A_1'] * x_i['E_1']) +
                                  2. * 2**(0.5) / (x_i['A_2'] * x_i['E_2']) -
                                  2. * 2**(0.5) / (x_i['A_3'] * x_i['E_3']) +
                                  2. / (x_i['A_4'] * x_i['E_4']))

    return vol, disp
