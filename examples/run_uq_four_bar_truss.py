"""
Run a small sparse-PCE uncertainty quantification for FOUR_BAR_TRUSS.
"""

import rheia.UQ.uncertainty_quantification as rheia_uq


dict_uq = {
    'case': 'FOUR_BAR_TRUSS',
    'pol order': 2,
    'uq method': 'sparse',
    'n samples': 12,
    'objective names': ['volume', 'displacement'],
    'objective of interest': 'displacement',
    'sampling method': 'SOBOL',
    'n jobs': 1,
    'results dir': 'example_uq',
    'draw pdf cdf': [False],
}


if __name__ == '__main__':
    rheia_uq.run_uq(dict_uq, design_space='design_space_uq.csv')
