import numpy as np
import h5py as h5
import glob
import sys
import os
import time

path_file = sys.argv[1]

data=np.loadtxt(path_file, skiprows=1)

ftime=open('stat_times.txt','w')
fphys=open('stat_phys.txt','w')
ftrack=open(path_file[:-4]+'.track','w')


times=data[:,0]
times=times+1.e-5

uv = data[:,2]

tgs=data[:,6]
#tgs[0]=tgs[1]

rhos=data[:,5]
#rhos[0]=rhos[1]

taus=data[0:,8]

ns=len(rhos)


h5file = 'ast_out{0:06d}.h5'

for i in range(0,ns):
    print i
    fsrc=open('source.mdl','w')
    finp=open('input.ini','w')
    fim=open('modinput.txt','r')
    fsm=open('modsrc.txt','r')

    if i == 0:
#create input file
        for j in range(0,12):
            line=fim.readline()
            if 'chi' in line:
                finp.write('chi = {}\n'.format(uv[i]))
            else:
                finp.write(line)
        finp.write('tf = 1e6\n')
        finp.write('rel_err = 1e-10\n')
        line=fim.readline()
        for j in range(13,26):
            line=fim.readline()
            finp.write(line)
        finp.close()
#create src file
        for j in range(0,3):
            line=fsm.readline()
            fsrc.write(line)
        fsrc.write("%i    %8.2f     %8.2e    %8.2f   %8.2f\n" %(0,taus[i]*1.086,2.*rhos[i]/2.3/1.67e-24,tgs[i],tgs[i]))
        fsrc.close()

        os.system('astrochem input.ini')
        os.system('mv astrochem_output.h5 {}'.format(h5file.format(i)))

        f=h5.File(h5file.format(i),"r")
        abuns=f['Abundances']
        spec=f['Species']
        tsteps=f['TimeSteps']
        nspec=spec.size

        for j in range(0,nspec):
            ftrack.write(spec[j]+'\t')
        ftrack.write('\n')
        for j in range(0,nspec):
            ab=str('%10.4e' % abuns[0,1,j])
            ftrack.write(str(ab)+'\t')
        ftrack.write('\n')
        f.close()


    if i > 0:
    #create input file
        f=h5.File(h5file.format(i-1),"r")
        abuns=f['Abundances']
        spec=f['Species']
        tsteps=f['TimeSteps']
        nspec=spec.size
        for j in range(0,12):
            line=fim.readline()
            if 'chi' in line:
                print(uv[i])
                finp.write('chi = {}\n'.format(uv[i]))
            else:
                finp.write(line)
            diffs=times[i]-times[i-1]
        ts=str('%12.6e' % diffs)
        finp.write('tf = '+ts+'\n')
        finp.write('rel_err = 1e-8\n')
        line=fim.readline()
        line=fim.readline()
        finp.write(line)
        line=fim.readline()
        finp.write(line)
        line=fim.readline()
        line=fim.readline()
        for j in range(0,5):
            line=fim.readline()
        for j in range(0,nspec):
            ab=str('%10.4e' % abuns[0,1,j])
            finp.write(spec[j]+' = '+ab+'\n')
        for j in range(0,4):
            line=fim.readline()
            finp.write(line)
        finp.close()
        f.close()

       #create src file
        for j in range(0,3):
            line=fsm.readline()
            fsrc.write(line)
        fsrc.write("%i    %8.2f     %8.2e    %8.2f   %8.2f\n" %(0,taus[i]*1.086,2.*rhos[i]/2.3/1.67e-24,tgs[i],tgs[i]))
        fsrc.close()
        if i < 100000:
            ct=str(i)
        if i < 10000:
            ct='0'+str(i)
        if i < 1000:
            ct='00'+str(i)
        if i < 100:
            ct='000'+str(i)
        if i < 10:
            ct='0000'+str(i)
        os.system('astrochem input.ini')
        os.system('astrochem input.ini')
        os.system('astrochem input.ini')
        os.system('mv astrochem_output.h5 {}'.format(h5file.format(i)))

        f=h5.File(h5file.format(i),"r")
        abuns=f['Abundances']
        spec=f['Species']
        tsteps=f['TimeSteps']
        nspec=spec.size
        
        for j in range(0,nspec):
            ab=str('%10.4e' % abuns[0,1,j])
            ftrack.write(str(ab)+'\t')
        ftrack.write('\n')
        f.close()


ftime.close()
fphys.close()
fsrc.close()
finp.close()
fim.close()
fsm.close()
ftrack.close()
