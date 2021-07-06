import numpy as np
import os
import sys
from constants import *
import model_dust as mod_dust
import model_gas as mod_gas
import yaml
import plot_summary as ps
import save_spec as ssp

PATH = "/Users/jenniferbergner/Research/adina_flares/disks_and_outbursts/radmc/"
sys.path.append(PATH + 'scripts/')

models = ['fiducial']
modpath = PATH + "models/"

for mod in models:
    if not os.path.exists('%s%s/' %(modpath, mod)):
        os.system('mkdir %s%s/' %(modpath, mod))

os.system('python save_model_params.py')

for mod in models:
    print(mod)
    os.chdir(modpath + mod)

    with open('model_inputs.yaml') as infile:
        mi = yaml.load(infile, Loader=yaml.CLoader)

    ssp.save_spectrum(mi, PATH, modpath + mod)
    os.system('cp %s/dustkappa_atmosphere.inp .' %(PATH))
    os.system('cp %s/dustkappa_midplane.inp .' %(PATH))

    mod_dust.setup_radmc(mi)
    os.system('rm -rf nohup.out')
    os.system('nohup radmc3d mctherm')
    print('completed thermal mc')
    mod_gas.save_gasdisk(mi, mod)

    ps.plot_model(mi, mod)
