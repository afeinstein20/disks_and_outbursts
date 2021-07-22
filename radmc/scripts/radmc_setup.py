import numpy as np
import os
import sys
from .constants import *
from .model_dust import setup_radmc
from .model_gas import save_gasdisk
import yaml
from .save_model_params import save_mod_params
from .save_model_params_herbig import save_mod_params_herbig
from .plot_summary import plot_model, plot_ind_pop
from .save_spec import save_spectrum
from .model_dust_inputstar import setup_uv

__all__ = ['setup']

def setup(PATH, models=['fiducial'], star_params={}, disk_params={},
          grid_params={}, uv=True, run=True, disk_type='TW Hydra'):

    modpath = os.path.join(PATH, "models/")
    
    for mod in models:
        if not os.path.exists('%s%s/' %(modpath, mod)):
            os.system('mkdir %s%s/' %(modpath, mod))

    if disk_type.lower() == 'tw hydra':
        save_mod_params(mod_name=models[0], 
                        star_dict=star_params, 
                        disk_dict=disk_params,
                        grid_dict=grid_params,
                        mod_path=os.path.join(modpath, models[0]),
                        )
    elif disk_type.lower() == 'herbig':
        save_mod_params_herbig(mod_name=models[0],
                               star_dict=star_params,
                               disk_dict=disk_params,
                               grid_dict=grid_params,
                               mod_path=os.path.join(modpath, models[0]),
                               )
    else:
        print("Disk type not implemented yet.")
        sys.exit()

    for mod in models:
        print(mod)
        os.chdir(modpath + mod)
        
        with open('model_inputs.yaml') as infile:
            mi = yaml.load(infile, Loader=yaml.CLoader)

        ## Adds in UV component
        if run == True:
            save_spectrum(mi, PATH, os.path.join(modpath, models[0]))

            os.system('cp %s/dustkappa_atmosphere.inp .' %(PATH))
            os.system('cp %s/dustkappa_midplane.inp .' %(PATH))


            if uv == True:
                setup_uv(mi)
            else:
                setup_radmc(mi)

#        plot_model(mi, mod)
#        plot_ind_pop(mi, mod)

            os.system('rm -rf nohup.out')
            os.system('nohup ./radmc3d mctherm')
            print('completed thermal mc')
            save_gasdisk(mi, mod)

        plot_model(mi, mod)
        plot_ind_pop(mi, mod)

