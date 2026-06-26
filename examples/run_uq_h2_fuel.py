"""
Run a small sparse-PCE uncertainty quantification for FOUR_BAR_TRUSS.
"""

import rheia.UQ.uncertainty_quantification as rheia_uq


dict_uq = {
    'case': 'H2_FUEL',
    'pol order': 3,
    'uq method': 'full',
    # 'n samples': 12,
    'objective names': ['volume', 'displacement'],
    'objective of interest': 'displacement',
    'n jobs': 4,
    'results dir': 'run_123',
    'draw pdf cdf': [False],
}


if __name__ == '__main__':
    rheia_uq.run_uq(dict_uq, design_space='design_space_tutorial_0.csv')
