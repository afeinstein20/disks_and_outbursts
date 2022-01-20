import yaml
from scipy.ndimage import map_coordinates
import pickle
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
import matplotlib as mpl
from astropy import constants, units
from mpl_toolkits.axes_grid1 import make_axes_locatable

import warnings
warnings.filterwarnings("ignore")

__all__ = ['open_raw', 'open_interp', 'pressure_profile', 'scale_height',
           'plot_contour', 'make_line_profiles']

def open_raw(path):
    """ 
    Opens the non-interpolated radmc output file. Takes the path
    to where the files are stored.
    """
    with open(os.path.join(path, 'diskdata_raw.pkl'), 'rb') as infile:
        diskinp = pickle.load(infile)
        ddustsm, ddustlg, tdustsm, tdustlg, dgas, tgas, re, ze = diskinp
    return ddustsm, ddustlg, tdustsm, tdustlg, dgas, tgas, re, ze

def open_interp(path):
    """
    Opens the interpolated radmc output file. Takes the path to
    where the files are stored.
    """
    with open(os.path.join(path,'diskdata_interpolated.pkl'), 'rb') as infile:
                  diskinp = pickle.load(infile)
                  ddustsm, ddustlg, tdustsm, tdustlg, dgas, tgas, re, ze = diskinp
    return ddustsm, ddustlg, tdustsm, tdustlg, dgas, tgas, re, ze

def pressure_profile(rho, T, mu=2.3):
    """
    Calculates the pressure profile. Takes the density and temperature 2D arrays.
    """
    rho = rho * units.g / units.cm**3
    T = T * units.K
    
    kB = 1.3806*10**-23 * units.Joule / units.K
    mH = 1.67*10**-24 * units.g
    
    cs = (kB * T) / (mu * mH)
    P = rho * cs
    return P.to(units.Joule/units.cm**3)

def scale_height(r, midplane, mstar, mu=2.3):
    """
    Calculates the scale height. Takes r (radius, in AU), the midplane temperature (K),
    and the mass of the star (solar masses).
    """
    Tm = np.nanmean(midplane)*units.K
    r = r * units.AU
    mstar = mstar * units.Msun
    
    kB = 1.3806*10**-23 * units.Joule / units.K
    mH = 1.67*10**-24 * units.g
    
    num = kB * Tm * r**3
    denom = constants.G * mstar * mu * mH
    return np.sqrt(num/denom).to(units.AU)


def plot_contour(data, cmap, ticks, ax, label, index=5, contourcolor=['k'], 
                 skipticks=0, logNorm=False, xlim=None, ylim=None, 
                 pressure=False, mu=2.3, add_cbar=True, add_xlabel=True):
    """
    Plots the contour maps.
    Takes:
       data : radmc interpolated output files (all)
       cmap : the colormap you wish to use
       ticks : the colorbar ticks to scale the contour by
       label : what to label the colorbar with
       index : the index for the data array that you want to plot
    """
    cmap = mpl.cm.get_cmap(cmap)
    cmap.set_under('w')
    cmap.set_over(cmap(1.0))
    
    if pressure == False:
        dat = data[index] + 0.0
    else:
        dat = pressure_profile(data[4], data[5], mu=mu).value
    
    if logNorm == False:
        im = ax.pcolormesh(data[-2], data[-1], dat, cmap=cmap,
                           vmin=np.nanmin(ticks), vmax=np.nanmax(ticks))#,
                       #norm=mpl.colors.LogNorm(vmin=np.nanmin(ticks), vmax=np.nanmax(ticks)))
    else:
        im = ax.pcolormesh(data[-2], data[-1], dat, cmap=cmap,
                           norm=mpl.colors.LogNorm(vmin=np.nanmin(ticks), 
                                                   vmax=np.nanmax(ticks)))
        
    ax.contour(data[-2], data[-1], dat, ticks, colors=contourcolor)
    div = make_axes_locatable(ax)
    if add_cbar:
        cax = div.append_axes('top', size='6%', pad=0.1)
        cbar = plt.colorbar(im, cax=cax, orientation='horizontal', ticks=ticks[skipticks:])
        cbar.set_label(label, fontsize=16)
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.xaxis.set_label_position('top')
        cbar.ax.tick_params(labelsize = 14)
    if add_xlabel:
        ax.set_xlabel('R [AU]')
    
    if xlim is not None:
        ax.set_xlim(xlim[0], xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])
    
    ax.set_rasterized(True)
    
    return


def make_line_profiles(data, index, zrs=[0,0.1,0.2,0.3,0.4], p=None,
                       ntheta=90, nr=110):
    """
    Plots the line countours at a given z/r location
    Takes:
       data : the raw radmc output data
       index : which index to access within the data structure
    """
    
    r2d = data[-2].reshape(ntheta,nr)
    z2d = data[-1].reshape(ntheta,nr)

    if p is None:
        var2d = data[index].reshape(ntheta,nr)
    else:
        var2d = p.reshape(ntheta,nr)

    zvals = []
    rvals = []
    var_slices = []
    for nn in zrs:
        vv = np.argmin(np.abs(z2d/r2d-nn),axis=0)[0]
        var_slices.append(var2d[vv])
        rvals.append(r2d[vv])
        zvals.append(z2d[vv])
        
    return rvals, zvals, var_slices
