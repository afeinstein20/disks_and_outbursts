import numpy as np
import h5py as h5
import glob
import sys
import os
import time
from astropy.table import Table

def split_file(fn):
    """
    Opens a give filename and splits by line.
    """
    return open(fn).read().split('\n')   

def run(i, snooze=False, sleep_time=10):
    """
    Sends commands to run astrochem and move output file.
    sleep_time : delay run by x [s]
    """
    global h5file
    os.system('astrochem input.ini')
    if snooze is True:
        time.sleep(sleep_time)
    os.system('mv astrochem_output.h5 {}'.format(h5file.format(i)))
    return

def append_table(fn, index, t=None):
    """
    Adds row to track file. If index is 0, creates table.
    Returns: astropy.table.Table
    """
    global times
    
    f=h5.File(h5file.format(i),"r")

    if index == 0:
        names = []
        for n in f['Species']:
            names.append(n.decode('utf-8'))
        names = np.append(['time'], names)

        t = Table(names=names)

    data = f['Abundances'][0][1]
    data = np.append([times[i]], data)
    t.add_row(data)
    return t

def create_source_file(i):
    """
    Reads the modsrc.txt file and creates source.mdl file for 
    index passed in.
    """
    global taus, rhos, tgs

    fsrc = open('source.mdl','w')
    contents = split_file('modsrc.txt')
    contents.append("%i    %8.2f     %8.2e    %8.2f   %8.2f" %(0,taus[i],
                                                               rhos[i],
                                                               tgs[i],
                                                               tgs[i]))
    for i,c in enumerate(contents):
        if c != '':
            if i != len(contents)-1:
                c += '\n'
            fsrc.write(c)
    fsrc.close()
    return


def create_input_file(i):
    """
    Creates the new input.ini file.
    """
    global times, h5file, UVflux

    finp=open('input.ini','w')
    contents = split_file('modinput.txt')

    input_ini = []

    new_additions = False
    
    # Takes modinput.txt --> input.ini
    for j in range(len(contents)):
        if contents[j] != '':
            row = contents[j].split(' ')
            row[-1] += '\n'

            non_dupes = ['[output]\n', 'time_steps',
                         'abundances']

            # If extra chemicals were added from a previous run,
            # avoids duplicates in the input file
            if i>0 and new_additions is True:
                for n in non_dupes:
                    if n in row:
                        finp.write(' '.join(str(e) for e in row))

            else:
                if 'chi' in row:
                    finp.write('chi = {0}\n'.format(UVflux[i]))
                else:
                    finp.write(' '.join(str(e) for e in row))
    
                    # Sets final time and error on time
                    if 'ti' in row:
                        if i == 0:
                            finp.write('tf = 1e6\n')
                            finp.write('rel_err = 1e-10\n')
                        else:
                            diff = times[i]-times[i-1]
                            tf = str('%12.6e' % diff)
                            finp.write('tf = {}\n'.format(tf))
                            finp.write('rel_err = 1e-8\n')

                    # Loads in abundances from previous run
                    if '[abundances]\n' in row and i != 0:
                        f=h5.File(h5file.format(i-1),"r")
                        for k in range(len(f['Species'])):
                            line = [f['Species'][k].decode('utf-8'), 
                                    "=", str('%10.4e' % f['Abundances'][0,1,k])+'\n']
                            finp.write(' '.join(str(e) for e in line))
                        new_additions = True

    finp.close()
    return


# Path file
path = str(sys.argv[1])
strength = int(sys.argv[2])
path_name = path+'.out'
xray_name = path+'_xray{0}.txt'.format(strength)

# Track file
track_name = str(sys.argv[3])

try:
    snooze = bool(sys.argv[4])
    sleep_time = int(sys.argv[5])
except:
    snooze = False
    sleep_time=0


# Sets output name
h5file='ast_out{0:05d}.h5'

# Loading in the changing parameters
data=np.loadtxt(path_name, skiprows=1)
UVflux=np.loadtxt(xray_name, skiprows=1)

times=data[:,0] + 1.e-5 # Time
tgs=data[:,4] # Gas Temperature
rhos=2*data[:,3]/2.3/1.67e-24 # Density
taus=data[:,6]*1.086 # Optical Depth

for i in range(5):
    create_input_file(i)
    create_source_file(i)
    run(i, snooze=snooze, sleep_time=sleep_time)
    
    if i == 0:
        t = append_table(h5file.format(i), i)
    else:
        t = append_table(h5file.format(i), i, t=t)

# Saves the tracks to a file
t.write(track_name, format='ascii')
