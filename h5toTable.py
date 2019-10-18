import os
import sys
import h5py
import numpy as np
from tqdm import tqdm
from astropy.table import Table


directory=str(sys.argv[1])
path_file=str(sys.argv[2])
output_fn=str(sys.argv[3])

files = os.listdir(directory)
files = np.sort([os.path.join(directory, i) for i in 
                 files if i.endswith('.h5')])

times = np.loadtxt(path_file)[0:,0] + 1e-5

f = h5py.File(files[0], 'r')

names = []
for i in range(len(f['Species'])):
    l = f['Species'][i].decode('utf-8')
    names.append(l)
names = np.append(['time'], names)

t = Table(names=names)

for i, fn in tqdm(enumerate(files)):
    f = h5py.File(fn, 'r')
    data = f['Abundances'][0][0]
    data = np.append([times[i]], data)
    t.add_row(data)


t.write(os.path.join(directory, output_fn), format='ascii')

