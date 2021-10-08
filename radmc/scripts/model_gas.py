import radmc3dPy as rmc
from scipy.interpolate import griddata
import numpy as np
import sys
import os
import csv
from .constants import *
import dill as pickle
from .disk_eqns import T_rho_gas

__all__ = ['save_gasdisk']

def save_gasdisk(mi, mod):

    # Load dust densities & temperatures
    data = rmc.analyze.readData(ddens = True, dtemp = True)

    # Grid coordinate system
    grid = data.grid
    rx = grid.x
    theta = grid.y
    phi = grid.z
    ncells = grid.nx*grid.ny*grid.nz

    #Transform spherical to cartesian coordinates
    zz = np.zeros(ncells)
    rr = np.zeros_like(zz)
    tdustsm = np.zeros_like(zz)
    tdustlg = np.zeros_like(zz)
    ddustsm = np.zeros_like(zz)
    ddustlg = np.zeros_like(zz)
    tgas = np.zeros_like(zz)
    dgas = np.zeros_like(zz)

    ## Calculate gas densities and intensity field
    kk = 0
    for jj in range(grid.ny):
        for ii in range(grid.nx):
            zz[kk] = rx[ii]*np.cos(theta[jj])
            rr[kk] = rx[ii]*np.sin(theta[jj])

            tdustsm[kk] = data.dusttemp[ii,jj,0,0]
            tdustlg[kk] = data.dusttemp[ii,jj,0,1]
            ddustsm[kk] = data.rhodust[ii,jj,0,0]
            ddustlg[kk] =  data.rhodust[ii,jj,0,1]
            tgas[kk], dgas[kk] = T_rho_gas(rr[kk], zz[kk], tdustsm[kk], mi)

            kk += 1

    # R, Z in AU units and sorted low to high
    rr = rr/au
    zz = zz/au
    rsort = rr[np.argsort(rr)]
    zsort = zz[np.argsort(zz)]
    nr = mi['nr']
    nz = mi['ntheta']

    # higher-res sampling across r and z boundaries
    re, ze = np.mgrid[np.min(rr):np.max(rr):8*(nr+1)*1j, np.min(zz):np.max(zz):8*(nz+1)*1j]
    # add 1 row mirrored below the midplane to enable linear interpolation down to z=0
    r_re = rr.reshape(nz,nr)
    z_re = zz.reshape(nz,nr)
    r_add = np.vstack([r_re[-1], r_re])
    z_add = np.vstack([-z_re[-1], z_re])
    r_rav = r_add.ravel()
    z_rav = z_add.ravel()

    dlist = ['ddustsm', 'ddustlg', 'tdustsm', 'tdustlg', 'tgas', 'dgas']
    datadict = {}
    for dd,data in enumerate([ddustsm, ddustlg, tdustsm, tdustlg, tgas, dgas]):
        dname = dlist[dd]
        data_re = data.reshape(nz, nr)
        data_add = np.vstack([data_re[-1], data_re])
        data_rav = data_add.ravel()

        # interpolate onto linear grid
        data_i = griddata((r_rav,z_rav), data_rav, (re,ze), method='linear')
        if dname=='dgas':
            data_i = data_i/(2.37*mH) # gas density to cm^-3
        datadict[dname] = data_i

    with open('diskdata_interpolated.pkl', 'wb') as savefile:
        pickle.dump([datadict['ddustsm'], datadict['ddustlg'], datadict['tdustsm'], datadict['tdustlg'],
            datadict['dgas'], datadict['tgas'], re, ze], savefile)

    with open('diskdata_raw.pkl', 'wb') as savefile:
        pickle.dump([ddustsm, ddustlg, tdustsm, tdustlg, tgas, dgas, rr, zz], savefile)
    print('Saved gas disk data')
