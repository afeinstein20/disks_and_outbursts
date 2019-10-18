import numpy as np
import scipy as sc
from tqdm import tqdm
from datetime import date
from astropy import units as u
import scipy.optimize as optimize
from astropy import constants as const
from scipy.interpolate import interp1d
from astropy.table import Table, Column

__all__ = ['disk_model']

class disk_model(object):
    
    def __init__(self, time=np.linspace(0, 1e6, 1000),
                 r=50.0, z=4.5, mid_T=21, 
                 M_star=1, kappa=1, sigma=40,
                 zr=None):
        """
        Creates a .out file to run through astrochem.

        Parameters
        ----------
        time : time array [years]
        r : radial distance [AU]
        z : vertical distance [AU]
        mid_T : midplane temperature [K]
        M_star : mass of star [M_sun]
        kappa : opacity [cm^2 g]
        sigma : surface density [g cm^-2]
        zr : z/r, optional
        """
        self.time = time
        if type(time) != u.quantity.Quantity:
            self.time *= u.year

        self.t_gas_mid = mid_T * u.K
        self.t_gas_atm = 3.0 * mid_T * u.K

        self.sigma = sigma * u.g / u.cm**2
        self.kappa = kappa * u.g * u.cm**2

        self.M_star = (M_star * const.M_sun).to(u.g)

        self.L_star = 0.23 * M_star**2.3 * const.L_sun
        self.L_star_time = np.full(self.time.shape, self.L_star.value)
        self.L_star_time_flare = None

        self.fn = None

        self.r = (r * u.AU).to(u.cm)

        if zr is not None:
            self.z = zr * self.r
        else:
            self.z = (z * u.AU).to(u.cm)

        self.scale_height()
        self.density_init()
        self.density_rad()

        tvals = self.z / (self.H * np.sqrt(2))
        self.optical_depth(sc.special.erfc(tvals))

        self.cosmic_rays()


    def scale_height(self):
        """
        Attributes
        ----------
        H : scale height [cm]
        """
        num = np.sqrt( const.k_B * self.t_gas_mid / (2.3*const.m_p))
        num = num.to(u.m/u.s)
        denom = np.sqrt(const.G * self.M_star / self.r**3)
        H = num/denom
        self.H = H.to(u.cm)

    def density_init(self):
        """
        Attributes
        ----------
        rho_o : initial density [g cm^-3]
        """
        rho_o = self.sigma/ (np.sqrt(2 * np.pi) * self.H)
        self.rho_o = rho_o.to(u.g/u.cm**3)
    
    def density_rad(self):
        """
        Attributes
        ----------
        rho : density [g cm^-3]
        """
        exp = (-self.z**2)/(2.0 * self.H**2)
        rho = self.rho_o * np.exp(exp.value)
        self.rho = rho
    
    def optical_depth(self, taus):
        """
        Attributes
        ----------
        tau : optical depth [g^2]
        """
        od = np.sqrt(np.pi * self.H**2 * 0.5) * self.kappa * self.rho_o * taus
        self.tau = od
    
    def cosmic_rays(self):
        """
        Attributes
        ----------
        cr : cosmic ray??
        """
        cosmic_constant = 5.0e-17
        exp_denom = 96.0 * u.g / u.cm**2
        self.cr = cosmic_constant * np.exp(-self.sigma/exp_denom)
        

    def inject_flare(self, amp=100, t0_ind=200,
                     rise=0.000003, fall=0.02, 
                     xray=False, base_xray=10):
        """
        Injects a Gaussian decay in the luminosity of the star.

        Parameters
        ---------- 
        amp : float, amplitude of the increase in brightness.
        t0_ind : int, where the flare occurs in the time array.
        rise : float, the steepness of the rise of the flare.
        fall : float, the steepness of the decay of the flare.
        xray : bool, returns a flare to feed into the chi value.
        base_xray : float, base x-ray flux value.
        """
        def gauss_rise(time, amp, t0, rise):
            return amp * np.exp( -(time - t0)**2.0 / (2.0*rise**2.0) )

        def exp_decay(time, amp, t0, fall):
            return amp * np.exp( -(time - t0) / fall )

        time = self.time.value

        if xray:
            baseline = base_xray
        else:
            baseline = self.L_star.value

        growth = np.where(time <= time[t0_ind])[0]
        decay  = np.where(time >  time[t0_ind])[0]

        rise = gauss_rise(time[growth], amp*baseline,
                          time[t0_ind], rise)
        fall = exp_decay( time[decay] , amp*baseline,
                          time[t0_ind], fall)

        if xray is False:
            self.L_star_time_flare = self.L_star_time + np.append(rise, fall)
            self.L_star_time_flare *= u.W
        else:
            self.xray_change = np.append(rise, fall)+baseline
            self.xray_factor = baseline*amp

    def calc_temp(self, z_q=2, d=3, flare=True):
        """
        Calculates the temperature in the disk.

        Parameters
        ---------- 
        z_q : 
        d : 
        flare : bool, tells which luminosity to use.

        Attributes
        ---------- 
        temp : temperature [K]
        """
        temp = np.zeros(len(self.time))

        if flare is True:
            if self.L_star_time_flare is None:
                raise ValueError("Please call disk_model.inject_flare() first.")
            else:
                L_star = self.L_star_time_flare
        else:
            L_star = self.L_star_time

        for i in tqdm(range(len(self.time))):
            del_lum = (L_star[i]/const.L_sun)**0.25

            if self.z < z_q*self.H:
                t1 = self.t_gas_atm - self.t_gas_mid
                t2 = np.sin( ( (np.pi*self.z.to(u.AU))/(2.0*z_q*self.H.to(u.AU)) ).value)**(2.0*d)
                t = self.t_gas_mid + t1*t2
            else:
                t = self.t_gas_atm
                
            temp[i] = t.to(u.K).value * del_lum

        self.temp = temp


    def default_file(self, xray):
        """
        Creates a default file name if one is not passed in.
        """
        today = date.today().strftime("%d%m%Y")
        fn = '{0}_Mstar{3}_r{1}_z{2}'.format(today,
                                             self.r.to(u.AU).value,
                                             self.z.to(u.AU).value,
                                             np.round(self.M_star.to(const.M_sun).value,2))
        self.fn = fn
        return
        

    def write_path_file(self, fn=None, alpha=0.001):
        """
        Creates the .out file with parameters.
        Columns are: time, r, z, n, tgas, alpha, opt_depth.
        """
        if fn is None:
            self.default_file(False)
            fn = self.fn + '.out'
        else:
            self.fn = fn.split('.')[0]

        t = Table()
        t.add_column(Column(data=self.time, name='time'))
        t.add_column(Column(data=np.full(self.time.shape, self.r.to(u.AU).value), name='r'))
        t.add_column(Column(data=np.full(self.time.shape, self.z.to(u.AU).value), name='z'))
        t.add_column(Column(data=np.full(self.time.shape, self.rho.to(u.g/u.cm**3).value), name='n'))
        t.add_column(Column(data=self.temp, name='tgas'))
        t.add_column(Column(data=np.full(self.time.shape, alpha), name='alpha'))
        t.add_column(Column(data=np.full(self.time.shape, self.tau.value), name='opt_depth'))

        t.write(fn, format='ascii')

        return


    def write_chi_file(self):
        """
        Creates a file with the change in UV flux.
        """
        if self.fn is None:
            raise ValueError("Please create .out file to go with this x-ray change.")
        
        t = Table()
        t.add_column(Column(data=self.xray_change, name='chi'))
        t.write(self.fn+'_xray{0}.txt'.format(self.xray_factor), format='ascii')
