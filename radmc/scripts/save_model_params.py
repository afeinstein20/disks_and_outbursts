import numpy as np
import os, sys
from .constants import *

__all__ = ['save_mod_params']

def save_mod_params(mod_name, mod_path, star_dict={}, disk_dict={}, grid_dict={}):
####### Fiducial disk parameters (TW Hya-like) #######
    disk_keys = np.array(list(disk_dict.keys()))
    star_keys = np.array(list(star_dict.keys()))
    grid_keys = np.array(list(grid_dict.keys()))

    # disk parameters
    # Huang+ 2018
    if 'r_in' in disk_keys: # Inner radius (AU)
        r_in = disk_dict['r_in'] + 0.0
    else:
        r_in   = 0.05  

    if 'r_out' in disk_keys: # Outer (gas disk) radius (AU)
        r_out = disk_dict['r_out'] + 0.0
    else:
        r_out  = 200

    if 'r_peb' in disk_keys: # Outer (pebble disk) radius (AU)
        r_peb = disk_dict['r_peb'] + 0.0
    else:
        r_peb  = 200

    if 'r_snow' in disk_keys: # Snow line radius (AU)
        r_snow = disk_dict['r_snow'] + 0.0
    else:
        r_snow = 0
    
    # Cleeves+ 2015
    H_c = 15        # Scale height at R_H AU (AU)
    R_H = 150       # Characteristic radius for scale height (AU)
    hh = 1.3        # Scale height gradient
    sigma_c = 0.04   # Characteristic *dust* surface density at R_c AU (g cm^-2)
    R_c = 150       # Characteristic radius for surface density (AU)
    gam = 1.       # Surface density gradient
    XH_mid = 0.2   # scale height midplane to atmosphere
    XR_mid = 0.9  # density fraction midplane to atmosphere'
    
    # star parameters
    # Cleeves+ 2015
    if 'T_star' in star_keys:
        T_star = star_dict['T_star'] + 0.0
    else:
        T_star = 4110     # Stellar effective temperature (K)

    if 'L_star' in star_keys:
        L_star = star_dict['L_star'] + 0.0
    else:
        L_star = 0.54      # Stellar bolometric luminosity (Lsol)
    
    if 'L_uv_star' in star_keys:
        L_uv_star = star_dict['L_uv_star']/lsol + 0.0
    else:
        L_uv_star = 2.9e31/lsol  # Stellar UV luminosity (Lsol)

    if 'R_star' in star_keys:
        R_star = star_dict['R_star'] + 0.0
    else:
        R_star = 1.04      # Stellar radius (Rsol)

    if 'M_star' in star_keys:
        M_star = star_dict['M_star'] + 0.0
    else:
        M_star = 0.8        # Stellar mass (Msol)
        
    if 'dpc' in star_keys:
        dpc = star_dict['dpc'] + 0.0
    else:
        dpc = 60           # distance to star (pc)

    # gas parameters
    # from Huang+2018
    Tc_atm = 125    # Characteristic atmosphere temperature (K) at R_T AU
    R_T = 10        # Characteristic radius for temperature (AU)
    q_atm = 0.47     # Atmosphere temperature gradient
    delta = 2.0     # Shape of vertical temperature gradient
    zq = 4          # temperature coupling scale height
    
    # grid parameters
    if 'nr' in grid_keys:
        nr = grid_dict['nr']
    else:
        nr = 90

    if 'ntheta' in grid_keys:
        ntheta = grid_dict['ntheta']
    else:
        ntheta = 80


    nphot = 1000000

    ff = os.path.join(mod_path, 'model_inputs.yaml')
    

    inps = ['nr: %1.0f' %(nr), 'ntheta: %1.0f' %(ntheta), 'nphi: 1', 'ped: 0.1',
            'rin: %1.7f  # AU; inner radius' %(r_in),
            'rout: %1.3f  # AU; small dust radius' %(r_out),
            'rsnow: %1.1f # AU; snow line radius' %(r_snow),
            'r_peb: %1.1f  # AU; mm dust radius' %(r_peb),
            'XH_mid: %1.1f  # scale height midplane to atmosphere' %(XH_mid),
            'XR_mid: %1.2f  # density fraction midplane to atmosphere' %(XR_mid),
            'eps_gas: 0.01  # dust to gas ratio',
            'gam: %1.2f  # surface density gradient' %(gam),
            'H_c: %1.1f  # AU; scale height' %(H_c),
            'R_H: %1.1f  # AU; characteristic radius for scale height' %(R_H),
            'hh: %1.2f  # scale height gradient' %(hh),
            'sigma_c: %1.4f  # g cm^-2; characteristic surface density' %(sigma_c),
            'r_c: %1.1f  # AU; characteristic radius for surface density' %(R_c),
            'dpc: %1.1f  # pc; distance to source' %(dpc),
            'T_star: %4.1f  # K; star temperature' %(T_star),
            'L_star: %2.2f  # Lsol; star luminosity' %(L_star),
            'Luv_star: %1.3e  # Lsol; star UV luminosity' %(L_uv_star),
            'R_star: %2.2f  # Rsol; star radius' %(R_star),
            'M_star: %2.2f  # Msol; star mass' %(M_star),
            'Tc_atm: %2.1f  # K; atmosphere temp at R_T' %(Tc_atm),
            'R_T: %1.1f  # AU; characteristic radius for atm temp' %(R_T),
            'q_atm: %1.2f  # atmospehre temperature gradient' %(q_atm),
            'delta: %1.1f  # shape of vertical temp gradient' %(delta),
            'zq: %1.1f   # temperature coupling scale height' %(zq),
            'nphot: %1.0f' %(nphot)
            ]
    
    print(ff)
    with open(ff, 'w+') as wfile:
        wfile.write('# %s model inputs\n' %(mod_name))
        for ii in inps:
            wfile.write('%s \n' %(ii))
