import os
import sys
path = os.path.dirname(os.path.abspath(__file__))


sys.path.insert(0, os.path.join(path,'POST_PROCESS'))

sys.path.insert(1, os.path.join(path,'OPT'))

sys.path.insert(2, os.path.join(path,'UQ'))