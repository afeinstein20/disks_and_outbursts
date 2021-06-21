import numpy as np
from constants import *
import radmc3dPy as rmc
from scipy.integrate import odeint, simps, quad
from scipy.interpolate import interp1d

def Tdust(rr, zz):
    ''''Dust temperature (output from radmc3d dust radiative transfer) at rr, zz based on nearest
       r, theta value'''

    if not given:
        cwd = os.getcwd()
        os.chdir(radmc_dir)
        data = rmc.analyze.readData(dtemp = True, binary = False)
        os.chdir(cwd)

        rpoints = data.grid.x
        tpoints = data.grid.y
        tdust_radmc = data.dusttemp[:,:,0,0]

        # spherical r, theta coords from cartesian r,z coords
        r_sp = np.sqrt(rr**2 + zz**2)
        theta = np.acos(zz/r_sp)

        ii_r = np.argmin(np.abs(rpoints - r_sp))
        ii_t = np.argmin(np.abs(tpoints - theta))

        return tdust_radmc[ii_r, ii_t]

def sigma_dust(rr, sigma_c, r_c, gam):
    '''Dust surface density at radius rr'''
    sig_dust = sigma_c*(rr/r_c)**(-gam)*np.exp(-(rr/r_c))**(2-gam)
    return sig_dust

def rho_dust_2(rr, zz, mi):
    '''Dust density for a 2-component dust structure.  Input cartesian r,z coordinates, characteristic
       surface density (g/cm^2), characteristic radius, surface density power law, 10AU scale height,
       scale height power law'''
    delt = 1.0
    # cut off density beyond inner & outer radii
    if rr < mi["rin"]*au:
        delt = 0.0
    if rr > mi['rout']*au:
        delt = 0.0

    sig_dust = sigma_dust(rr, mi["sigma_c"], mi["r_c"]*au, mi["gam"])*delt
    H_atm = mi["H_c"]*au*(rr/(mi["R_H"]*au))**mi["hh"]
    H_mid = mi["XH_mid"]*H_atm

    rho_atm = (1-mi["XR_mid"])*sig_dust/(np.sqrt(2*np.pi)*H_atm)*np.exp(-0.5*(zz/H_atm)**2)
    rho_mid = mi["XR_mid"]*sig_dust/(np.sqrt(2*np.pi)*H_mid)*np.exp(-0.5*(zz/H_mid)**2)

    if rr > mi["r_peb"]*au:
        rho_mid = 0.0

    return rho_atm, rho_mid

def sigma_gas(rr, mi):
    '''Gas surface density at radius rr, assuming vertically integrated gas:dust ratio of 100:1
       at each radius'''
    return sigma_dust(rr, mi["sigma_c"], mi["r_c"]*au, mi["gam"])/mi["eps_gas"]

def T_rho_gas(rr, zz, Tdust, mi):
    '''Gas temperature & density, parameterized as in Dartois+2003'''
    # Gas temp calculation
    delta = mi["delta"]
    Tgas_mid = Tdust
    Tgas_atm = mi["Tc_atm"]*(rr/(mi["R_T"]*au))**(-mi["q_atm"])
    Omega = np.sqrt(G0*mi["M_star"]*msol/rr**3)
    c_s = np.sqrt(kB*Tgas_mid/(mu_gas*mH))
    H_gas = c_s/Omega
    z_q = mi['zq']*H_gas

    if zz < z_q:
        Tg = Tgas_atm + (Tgas_mid - Tgas_atm)*(np.cos(np.pi*zz/(2*z_q)))**(2*delta)
    else:
        Tg = Tgas_atm

    if Tg < Tdust:
        Tg = Tdust

    if Tg > 1e5:
        Tg = 1e5

    # Gas density calculation
    if zz > z_q:
        logtgrad = 0
    else:
        logtgrad = -1./Tg*np.sin(2*np.pi*zz/(2*z_q))*np.pi/(2*z_q)*(Tgas_mid - Tgas_atm)

    y0 = 10**80
    maxz = 10*z_q  # integration outer bound

    def intfunc(rho, zz):
        '''diff eq describing vertical hydrostatic eq for gas temp distribution'''
        ai = G0 * mi["M_star"]*msol * zz / (rr**2 + zz**2)**1.5
        bi = mu_gas * mH / (kB * Tg)
        ci = logtgrad
        return -rho * (ai * bi + ci)

    if zz > maxz:
        rhogas = 0  # assume densities ~ 0 past max z value
        return Tg, rhogas

    zvals = np.logspace(np.log10(0.1), np.log10(maxz + 0.1), 1000) - 0.1
    sols = np.squeeze(odeint(intfunc, y0, zvals))

    bounds = sigma_gas(rr, mi)
    scale = bounds/(2*simps(sols, zvals))

    if zz == 0:
        rhogas = scale*y0
        return Tg, rhogas

    else:
        ff = interp1d(zvals, sols)
        rhogas = float(ff(zz))*scale
        return Tg, rhogas
