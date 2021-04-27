import os, sys
import shutil
import pytest
import multiprocessing as mp 
import numpy as np
import shutil

path = os.path.split(
    os.path.dirname(
        os.path.abspath(__file__)))[0]
sys.path.insert(0, os.path.join(path,'UQ'))
from uncertainty_quantification import run_uq

def create_design_space_test(case):
    path = os.path.split(
        os.path.dirname(
            os.path.abspath(__file__)))[0]

    des_var_file = os.path.join(path,
                                'CASES',
                                case,
                                'design_space',
                                )

    new_des_var_file = os.path.join(path,
                                'CASES',
                                case,
                                'design_space_test',
                                )

    if not os.path.isfile(new_des_var_file):
        with open(des_var_file, 'r') as f:
            text = []
            for line in f.readlines():
                found = False
                tmp = line.split()
                if tmp[1] == 'var':
                        text.append('%s par %f \n' %(tmp[0], np.mean([float(tmp[2]), float(tmp[3])])))
                        found = True
                if not found:
                    text.append(line)

        with open(new_des_var_file, 'w') as f:
            for item in text:
                f.write("%s" % item)

def del_results(case):
    path = os.path.split(
        os.path.dirname(
            os.path.abspath(__file__)))[0]

    dir_results = os.path.join(path,
                                'RESULTS',
                                case,
                                'UQ',
                                'testdir'
                                )

    new_des_var_file = os.path.join(path,
                                'CASES',
                                case,
                                'design_space_test',
                                )

    shutil.rmtree(dir_results)
    os.remove(new_des_var_file)

def test_pv_h2():
    case = 'PV_H2'
    dict_uq = {'case': case,
               'pol order': 1,
               'objective names': ['obj1', 'obj2'],
               'objective of interest': 'obj1',
               'results dir': 'testdir', 
               'n jobs': int(mp.cpu_count()/2), 
              }  
    create_design_space_test(case)
    run_uq(dict_uq, design_space = 'design_space_test')
    del_results(case)