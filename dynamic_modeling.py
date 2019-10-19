import os
import time
import h5py as h5
import numpy as np
from datetime import date
from astropy.table import Table


__all__ = ['dynamical_astrochem']

class dynamical_astrochem(object):
    
    def __init__(self, model=None, path_fn=None, track_fn=None,
                 out_fn=None, snooze=False, sleep_time=10):
        """
        Users can either pass in the disk_model class
        or the files they wish to use to update the 
        input.ini files.
        
        Parameters
        ----------
        model : disk_model.disk_model
        path_fn : .out filename
        track_fn : track filename
        out_fn : what to rename the .h5 files to
        snooze : pauses between each astrochem run
        sleep_time : number of seconds to pause between runs
        """
        if model is not None:
            self.path_fn = model.path_fn
            self.time = model.table['time'].data
            self.tgas = model.table['tgas'].data
            self.rhos = model.table['n'].data
            self.taus = model.table['opt_depth'].data
            self.chi  = model.table['UV_lum'].data
        else:
            self.path_fn = path_fn
            try:
                data = np.loadtxt(self.path_fn)
                self.time = data[:,0] + 1.e-5
                self.tgas = data[:,4]
                self.rhos = 2*data[:,3]/2.3/1.67e-24
                self.taus = data[:,6]*1.086
                self.chi  = np.full(self.time.shape, 100)
            except:
                data = Table.read(self.path_fn, format='ascii')
                self.time = data['time'].data
                self.tgas = data['tgas'].data
                self.rhos = data['n'].data
                self.taus = data['opt_depth'].data
                self.chi  = data['UV_lum'].data


        if track_fn is not None:
            self.track_fn = track_fn
        else:
            self.track_fn = 'track.txt'

        if out_fn is not None:
            self.out_fn = out_fn
        else:
            today = date.today().strftime("%d%m%Y")
            self.out_fn = today+ '_ast_out{0:06d}.h5'

        self.snooze = snooze
        self.sleep_time = sleep_time


    def run(self):
        """
        Runs the script.
        """
        for i in range(len(self.time)):
            self.create_input_file(i)
            self.create_source_file(i)
            self.astrochem(i)

            if i == 0:
                t = self.append_table(i)
            else:
                t = self.append_table(i, t=t)

        t.write(self.track_fn, format='ascii')
                              

    def split_file(self, fn):
        """
        Opens a given filename and splits by line.
        """
        return open(fn).read().split('\n')

    def astrochem(self, i):
        """
        Sends commands to run astrochem and move output file.
        """
        os.system('astrochem input.ini')

        if self.snooze:
            time.sleep(self.sleep_time)

        os.system('mv astrochem_output.h5 {}'.format(self.out_fn.format(i)))
        return

    def append_table(self, index, t=None):
        """
        Adds row to track file. If index is 0, creates table.
        Returns: astropy.table.Table
        """
        f = h5.File(self.out_fn.format(index), "r")
                     
        if index == 0:
            names = []
            for n in f['Species']:
                names.append(n.decode('utf-8'))
            names = np.append(['time'], names)
            t = Table(names=names)

        data = f['Abundances'][0][1]
        data = np.append( [self.time[index]], data)
        t.add_row(data)
        return t

    def create_source_file(self, i):
        """   
        Reads the modsrc.txt file and creates source.mdl file for
        index passed in.
        """
        fsrc = open('source.mdl','w')
        contents = self.split_file('modsrc.txt')
        contents.append("%i    %8.2f     %8.2e    %8.2f   %8.2f" %(0,self.taus[i],
                                                                   self.rhos[i],
                                                                   self.tgas[i],
                                                                   self.tgas[i]))
        for i,c in enumerate(contents):
            if c != '':
                if i != len(contents)-1:
                    c += '\n'
            fsrc.write(c)
        fsrc.close()
        return

    def create_input_file(self, i):
        """
        Creates the new input.ini file.
        """
        finp=open('input.ini','w')
        contents = self.split_file('modinput.txt')
        
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
                    finp.write('chi = {0}\n'.format(self.chi[i]))
                else:
                    finp.write(' '.join(str(e) for e in row))

                    # Sets final time and error on time
                    if 'ti' in row:
                        if i == 0:
                            finp.write('tf = 1e6\n')
                            finp.write('rel_err = 1e-10\n')
                        else:
                            diff = self.time[i]-self.time[i-1]
                            tf = str('%12.6e' % diff)
                            finp.write('tf = {}\n'.format(tf))
                            finp.write('rel_err = 1e-8\n')

                    # Loads in abundances from previous run
                    if '[abundances]\n' in row and i != 0:
                        f=h5.File(self.out_fn.format(i-1),"r")
                        for k in range(len(f['Species'])):
                            line = [f['Species'][k].decode('utf-8'),
                                    "=", str('%10.4e' % f['Abundances'][0,1,k])+'\n']
                            finp.write(' '.join(str(e) for e in line))
                        new_additions = True

        finp.close()
        return
