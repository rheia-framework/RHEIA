"""
Run a small deterministic optimization for the FOUR_BAR_TRUSS case.
"""

import multiprocessing as mp

import rheia.OPT.optimization as rheia_opt


dict_opt = {
    'case': 'FOUR_BAR_TRUSS',
    'objectives': {'DET': (-1, -1)},
    'stop': ('BUDGET', 40),
    'n jobs': 1,
    'population size': 10,
    'results dir': 'example_optimization',
}


if __name__ == '__main__':
    # For multiprocessing, keep this call inside the main guard on Windows.
    # Example: dict_opt['n jobs'] = max(1, int(mp.cpu_count() / 2))
    rheia_opt.run_opt(dict_opt)
