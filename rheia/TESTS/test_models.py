'''
Module to test the hydrogen-based energy system models.
'''

import os
import pytest
import rheia
import rheia.CASES.H2_FUEL.h2_fuel as lb
import rheia.CASES.H2_POWER.h2_power as lb2
import rheia.CASES.H2_MOBILITY.h2_mobility as lb3


def test_h2_fuel():
    """
    Test the LCOH and m_h2 for a generic design in the photovoltaic-hydrogen
    system model.
    """

    path = os.path.dirname(rheia.__file__)

    # the climate file considered
    filename_climate = os.path.join(path,
                                    'CASES',
                                    'DATA',
                                    'climate',
                                    'climate_Brussels.csv')

    # the object to read in the data
    my_data = lb.ReadData(filename_climate)

    # get the solar irradiance and ambient temperature
    sol_irr, t_amb = my_data.load_climate()

    # retrieve the deterministic values for the model parameters
    parameters = my_data.load_parameters()

    # define the design to be tested
    inputs = {'n_dcdc_pv': 5.,
              'n_pemel': 4.}

    # instantiate from the Evaluation class
    my_evaluation = lb.Evaluation(sol_irr, t_amb, {**parameters, **inputs})

    # evaluate the system
    my_evaluation.evaluation()

    assert round(my_evaluation.res['lcoh'], 3) == 10.172
    assert round(my_evaluation.res['m_h2'], 3) == 115.441


def test_h2_power():
    """
    Test the LCOE and SSR for a generic design in the power-to-power
    system model.
    """

    path = os.path.dirname(rheia.__file__)

    # the climate file considered
    filename_climate = os.path.join(path,
                                    'CASES',
                                    'DATA',
                                    'climate',
                                    'climate_Brussels.csv')

    filename_demand = os.path.join(path,
                                   'CASES',
                                   'DATA',
                                   'demand',
                                   'load_Brussels_dwelling.csv')

    # the object to read in the data
    my_data = lb2.ReadData(filename_climate, filename_demand)

    # get the solar irradiance and ambient temperature
    sol_irr, t_amb = my_data.load_climate()

    # get the electric load
    load_elec = my_data.load_demand()

    # retrieve the deterministic values for the model parameters
    parameters = my_data.load_parameters()

    # define the design to be tested
    inputs = {'n_pv': 10.,
              'n_pemel': 2.,
              'n_pemfc': 1.,
              'n_tank': 100.,
              }

    # instantiate from the Evaluation class
    my_evaluation = lb2.Evaluation(
        sol_irr, t_amb, load_elec, {
            **parameters, **inputs})

    # evaluate the system
    my_evaluation.evaluation()

    assert round(my_evaluation.res['lcoe'], 3) == 484.926
    assert round(my_evaluation.res['ssr'], 3) == 0.773


def test_h2_mobility():
    """
    Test the LCOD and CI for a generic design in the power-to-mobility
    system model.
    """

    path = os.path.dirname(rheia.__file__)

    # the climate file considered
    filename_climate = os.path.join(path,
                                    'CASES',
                                    'DATA',
                                    'climate',
                                    'climate_Brussels.csv')

    # the object to read in the data
    my_data = lb3.ReadData(filename_climate)

    # get the solar irradiance and ambient temperature
    sol_irr, t_amb = my_data.load_climate()

    # retrieve the deterministic values for the model parameters
    parameters = my_data.load_parameters()

    # define the design to be tested
    inputs = {'n_pv': 1000.,
              'n_pemel': 2000.,
              'n_tank': 10000.,
              'n_h2_bus': 1.,
              }

    # instantiate from the Evaluation class
    my_evaluation = lb3.Evaluation(sol_irr, t_amb, {**parameters, **inputs})

    # evaluate the system
    my_evaluation.evaluation()

    assert round(my_evaluation.res['lcom'], 3) == 1.357
    assert round(my_evaluation.res['ci'], 3) == 1.318
