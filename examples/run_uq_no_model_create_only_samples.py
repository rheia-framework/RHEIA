"""
Generate PCE training samples without a connected model.

The NO_MODEL case defines the uncertain inputs in design_space.csv and
stochastic_space.csv, while create only samples=True prevents model
evaluation and writes input-only samples.csv rows.
"""

import rheia.UQ.uncertainty_quantification as rheia_uq


dict_uq = {
    'case': 'NO_MODEL',
    'pol order': 2,
    'uq method': 'full',
    'objective names': ['output_1'],
    'objective of interest': 'output_1',
    'sampling method': 'SOBOL',
    'n jobs': 1,
    'results dir': 'example_create_only_samples',
    'create only samples': True,
    'draw pdf cdf': [False],
}


if __name__ == '__main__':
    rheia_uq.run_uq(dict_uq)
