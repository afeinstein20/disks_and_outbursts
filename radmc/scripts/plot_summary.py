import numpy as np
import os
import dill as pickle
import sys
import csv
from astropy import units, constants
from scipy.interpolate import griddata
from .constants import *
import radmc3dPy as rmc
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

mpl.rcParams['font.size'] = 13
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['ytick.minor.size'] = 3

__all__ = ['plot_model', 'pressure_profile']

def plot_model(mi, mod):

    # load disk data
    with open('diskdata.pkl', 'rb') as infile:
        diskinp = pickle.load(infile)
        ddustsm, ddustlg, tdustsm, tdustlg, dgas, tgas, re, ze = diskinp
    pgas = pressure_profile(dgas*units.g/units.cm**3,
                                tgas*units.K).value

    # plot setup
    f,ax = plt.subplots(1,5, figsize = (24,4), sharex = True, sharey = True)
    f.subplots_adjust(wspace = 0.13)
    for xx in range(1,4):
        plt.setp(ax[xx].get_xticklabels(), visible = False)
        plt.setp(ax[xx].get_yticklabels(), visible = False)
    ax[0].set_xlabel("R (AU)")
    ax[0].set_ylabel("Z (AU)", labelpad = 1)
    ax[0].set_xlim((0,mi['rout']))
    ax[0].set_ylim((-0.1,1.0))#20))

    dcmap = mpl.cm.get_cmap('BuPu')
    dcmap.set_under('w')
    dcmap.set_over(dcmap(1.0))

    tcmap = mpl.cm.get_cmap('RdYlBu_r')
    tcmap.set_under(tcmap(0.0))
    tcmap.set_over(tcmap(1.0))
    
    pmap = mpl.cm.get_cmap('viridis')
    pmap.set_under(pmap(0.0))
    pmap.set_over(pmap(1.0))

    # plot formatting
    ddust_ticks = np.logspace(-22, -14, 5)
    dgas_ticks = np.logspace(5,14,6) #np.logspace(5, 10, 6)
    tdust_ticks = np.linspace(30,200,9)#np.linspace(10, 90, 9)
    tgas_ticks = np.linspace(30,200,9)#np.linspace(10, 90, 9)
    pgas_ticks = np.linspace(-15,20,9)

    ticklist = [ddust_ticks, dgas_ticks, tdust_ticks, tgas_ticks, pgas_ticks]
    labels = ["Total dust density (g cm$^{-3}$)", "Gas density (cm$^{-3}$)", "Dust temperature (K)",
              "Gas temperature (K)", r"log$_{10}$(P [J/cm$^3$])"]
    cmaps = [dcmap, dcmap, tcmap, tcmap, pmap]

    for ii, data in enumerate([ddustsm+ddustlg, dgas, tdustsm, tgas, np.log10(pgas)]):
        ticks = ticklist[ii]
        tickspan = np.max(ticks) - np.min(ticks)
        if ii in [0,1]:
            aa = ax[ii].pcolormesh(re, ze, data, cmap = cmaps[ii], norm =
                mpl.colors.LogNorm(vmin = np.min(ticks),
                vmax = np.max(ticks)))
        else:
            aa = ax[ii].pcolormesh(re, ze, data, cmap = cmaps[ii], norm =
                mpl.colors.Normalize(vmin = np.min(ticks) - 0.05*tickspan,
                vmax = np.max(ticks) + 0.05*tickspan))
        ax[ii].contour(re, ze, data, ticks, colors = ['k'])

        div = make_axes_locatable(ax[ii])
        cax = div.append_axes('top', size ='6%', pad = 0.1)
        cbar = plt.colorbar(aa, cax = cax, orientation = 'horizontal', ticks = ticks)
        cbar.set_label(labels[ii])
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.xaxis.set_label_position('top')
        cbar.ax.tick_params(labelsize = 12)

    f.set_rasterized(True)
    f.savefig('%s_summary.png' %(mod), bbox_inches = 'tight')
    print('plotted %s summary figure' %(mod))


# Value for mu comes from https://arxiv.org/pdf/1509.06382.pdf
def pressure_profile(rho, T, mu=2.3):
    kB = 1.3806*10**-23 * units.Joule / units.K
    mH = 1.67*10**-24 * units.g
    cs = (kB * T) / (mu * mH)
    P = rho * cs
    return P.to(units.Joule/units.cm**3)
