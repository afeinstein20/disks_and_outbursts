import os
import numpy as np

files = os.listdir('.')
paths = np.sort([i for i in files if i.endswith('.out')] )
#paths = [i for i in paths if 'noflare' in i]

for fn in paths:
    print(fn)
    os.system('python dynchem.py {0}'.format(fn))
