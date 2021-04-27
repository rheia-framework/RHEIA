import subprocess as sp
import os

path = os.path.dirname(os.path.abspath(__file__))
path_cmd = '"' + path + '"'

#os.system('cd %s' %path)


tests = [f for f in os.listdir(path) if f.startswith('test_')]

for test in tests:
    cmd = 'pytest %s' %os.path.join(path_cmd, test)
    os.system(cmd)
    
