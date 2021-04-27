import os, sys
import shutil
import pytest
import multiprocessing as mp 
import gc

path = os.path.split(
    os.path.dirname(
        os.path.abspath(__file__)))[0]
sys.path.insert(0, os.path.join(path,'OPT'))
from optimization import run_opt

def del_results(case, opt_type):
    path = os.path.split(
        os.path.dirname(
            os.path.abspath(__file__)))[0]

    dir_results = os.path.join(path,
                                'RESULTS',
                                case,
                                opt_type,
                                'testdir'
                                )

    shutil.rmtree(dir_results)

def test_pvh2():
    case = 'PV_H2'
    opt_type = 'DET'
    dict_opt = {'case': case,
                'objectives':          {opt_type: (-1, 1)}, 
                'stop':                ('BUDGET', 12),
                'population size':     3,
                'results dir':      'testdir',
               }

    run_opt(dict_opt)
    del_results(case, opt_type)

def test_pvh2_rob():
    case = 'PV_H2'
    opt_type = 'ROB'
    dict_opt = {'case': case,
                'objectives':          {opt_type: (-1, -1)}, 
                'stop':                ('BUDGET', 30),
                'population size':     2,
                'results dir':      'testdir',
                'pol order': 1,
                'objective names': ['lcoe','ssr'],
                'objective of interest': ['lcoe'],
                'n jobs':              int(mp.cpu_count()/2), 
               }
    
    run_opt(dict_opt)
    del_results(case, opt_type)
