import numpy as np
import csv
from scipy.interpolate import interp1d
from constants import *

def BB(lam,T):
    '''Blackbody at wavelength (microns)'''
    v = 3e10/(lam*1e-4)
    return 2*hh*v**3/cc**2/(np.exp(hh*v/(kB*T)) - 1) # erg/s/cm^2/Hz/sr

def save_spectrum(mi, specpath, savepath):

    rstar = mi['R_star']*rsol # cm
    teff = mi['T_star']       # K
    luv_star = mi['Luv_star']*lsol # erg/s

    ## Read in TW Hya UV spectrum
    uwave, ufl = [],[]
    with open(specpath + 'twhya_uv.csv', 'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            uwave.append(float(row[0]))
            ufl.append(float(row[1]))

    uwave_A, jj = np.unique(uwave, return_index = True)  # Angstroms
    ufl_A = np.array(ufl)[jj]                            # photons/cm^2/s/A

    ## Convert from A^-1 to Hz^-1
    ufr_hz = cc/(uwave_A*1e-8) # frequencies in Hz
    ufl_hz = ufl_A*uwave_A/ufr_hz # fluxes in ph/cm^s/s/Hz
    ff = np.argsort(ufr_hz) # resort by frequency
    ufr_hz = ufr_hz[ff]
    ufl_hz = ufl_hz[ff]

    ## Convert from ph to erg
    ufl_erg = ufl_hz*(hh*ufr_hz) # erg/cm^2/s/Hz

    uwave_um = uwave_A[ff]/1e4 # wavelengths in um

    ## scale to star's measured luminosity
    luv = np.trapz(ufl_erg, ufr_hz)*4*np.pi*rstar**2 # erg/s
    facu = luv_star / luv
    ufl_scale = ufl_erg * facu

    ## downsample uv spectrum
    ulam_out = np.concatenate([uwave_um[:340:15], uwave_um[340:400:4], uwave_um[400:466:15], uwave_um[466:470], uwave_um[470::4]])
    ufl_out = np.concatenate([ufl_scale[:340:15], ufl_scale[340:400:4], ufl_scale[400:466:15], ufl_scale[466:470], ufl_scale[470::4]])

    ## sample longer wavelengths (um)
    lam1, lam2, lam3, lam4 = 3e-1, 7.0e0, 25.0e1, 1.0e4
    nlam = 15
    lam12    = np.logspace(np.log10(lam1),np.log10(lam2),nlam,endpoint=False)
    lam23    = np.logspace(np.log10(lam2),np.log10(lam3),nlam,endpoint=False)
    lam34    = np.logspace(np.log10(lam3),np.log10(lam4),nlam,endpoint=True)
    vlam      = np.concatenate([lam12,lam23,lam34])

    ## Add blackbody, uv spectra
    bb = BB(vlam, teff)*4*np.pi  #erg/s/cm^2/Hz
    all_lam = np.concatenate([ulam_out, vlam])
    all_fl = np.concatenate([ufl_out, bb])
    fl_pc = all_fl*(rstar/pc)**2
    ii = np.argsort(all_lam)
    all_lam = all_lam[ii]
    fl_pc = fl_pc[ii]

    np.savetxt(savepath + '/starspec.txt',
           [all_lam, fl_pc])

    print('saved star spectrum')
