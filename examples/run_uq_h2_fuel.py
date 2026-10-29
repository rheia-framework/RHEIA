"""
Run a full-PCE uncertainty quantification for the H2_FUEL case.

The H2_FUEL model returns the levelized cost of hydrogen and the hydrogen
production, in that order, so 'objective names' lists them in that order.
Here the hydrogen production is the quantity of interest; set 'objective of
interest' to 'lcoh' to analyse the levelized cost instead.
"""

import rheia.UQ.uncertainty_quantification as rheia_uq


dict_uq = {
    'case': 'H2_FUEL',
    'pol order': 3,
    'uq method': 'full',
    'objective names': ['lcoh', 'mh2'],
    'objective of interest': 'mh2',
    'n jobs': 4,
    'results dir': 'example_uq_h2_fuel',
    'draw pdf cdf': [False],
}


if __name__ == '__main__':
    rheia_uq.run_uq(dict_uq, design_space='design_space_tutorial_0.csv')
