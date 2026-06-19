import rheia.UQ.uncertainty_quantification as rheia_uq
import multiprocessing as mp


objectives = ['lcoh', 'mh2']

for i in objectives:
    dict_uq = {'case':                  'H2_FUEL',
               'pol order':             2,
               # 'uq method':             'full',
               'uq method':             'sparse',
               'n samples':             20,
               'objective names':       objectives,
               'objective of interest': i,
               'results dir':           'run_2_sparse',
               'draw pdf cdf':          [True, 1e5],
               # 'create only samples': True,
               }

    if __name__ == '__main__':
        rheia_uq.run_uq(dict_uq, design_space = 'design_space_tutorial_uq.csv')