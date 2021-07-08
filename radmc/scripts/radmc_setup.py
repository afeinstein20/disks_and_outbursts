import numpy as np
import os
import sys
from .constants import *
from .model_dust import setup_radmc
from .model_gas import save_gasdisk
import yaml
from .save_model_params import save_mod_params
from .plot_summary import plot_model
from .save_spec import save_spectrum
from .model_dust_inputstar import setup_uv

__all__ = ['setup']

def setup(PATH, models=['fiducial'], star_params={}, uv=True):

#    sys.path.append(PATH + 'scripts/')
    modpath = os.path.join(PATH, "models/")
    
    for mod in models:
        if not os.path.exists('%s%s/' %(modpath, mod)):
            os.system('mkdir %s%s/' %(modpath, mod))

    save_mod_params(mod_name=models[0], 
                    star_dict=star_params, 
                    mod_path=os.path.join(modpath, models[0])
                   )

    for mod in models:
        print(mod)
        os.chdir(modpath + mod)
        
        with open('model_inputs.yaml') as infile:
            mi = yaml.load(infile, Loader=yaml.CLoader)

        ## Adds in UV component
        save_spectrum(mi, PATH, os.path.join(modpath, models[0]))

        os.system('cp %s/dustkappa_atmosphere.inp .' %(PATH))
        os.system('cp %s/dustkappa_midplane.inp .' %(PATH))


        if uv == True:
            setup_uv(mi)
        else:
            setup_radmc(mi)


        os.system('rm -rf nohup.out')
        os.system('nohup ./radmc3d mctherm')
        print('completed thermal mc')
        save_gasdisk(mi, mod)
        
        plot_model(mi, mod)
