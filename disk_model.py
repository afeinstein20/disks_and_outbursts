import numpy as np
from datetime import date
from scipy.special import erf
from astropy import units as u
from astropy import constants as c
from astropy.table import Table, Column

__all__ = ['disk_model']

class disk_model(object):

    def __init__(self, r=1, z=0.03,
                 Mstar=1, gas_dist=2000,
                 opacity=0.01,
                 time=np.linspace(0,5,1000)):
        """
        Parameters
        ---------- 
        r : radial distance [AU]
        z : vertical distance from midplane [AU]
        Mstar : mass of the host star [M_sun]
        gas_dist : distribution of gas for surface density [g cm^-2]
        opacity : opacity of the disk [cm^2 g]
        time : array of times [days]

        Atributes
        ----------
        r : radial distance [AU]
        z : vertical distance [AU]
        Mstar : mass of the host star [M_sun]
        Lstar : luminosity of the host star [L_sun]
        time : array of times [years, unless units otherwise specified]
        """
        self.r = r*u.AU
        self.z = z*u.AU
        self.opacity = opacity * u.cm**2 * u.g

        self.midplane_temp()
        self.Tatm = 3*self.Tmid

        self.Mstar = Mstar*c.M_sun
        self.Lstar = 0.23 * (self.Mstar/c.M_sun).value**2.3 * c.L_sun

        self.time = time
        if type(time) != u.quantity.Quantity:
            self.time *= u.day

        self.scale_height()
        self.init_density(gas_dist * u.g / u.cm**2)
        self.density()
        self.optical_depth()
        self.temperature()

        self.lum_flare = None
        self.flare = None
        self.lum_UV = None

        return


    def midplane_temp(self):
        """
        Defines the initial midplane temperature.

        Attributes
        ---------- 
        Tmid : [K]
        """
        Tmid = 150*u.K / np.sqrt(self.r.to(u.AU).value)
        self.Tmid = Tmid


    def scale_height(self):
        """
        Calculates the scale height at disk radius r.

        Attributes
        ---------- 
        H : scale height [cm]
        """
        m_hyd   = c.m_p
        k_boltz = c.k_B

        num   = k_boltz * self.Tmid * self.r**3
        denom = 2.3 * m_hyd * c.G * self.Mstar

        self.H = np.sqrt(num/denom).to(u.cm)
        

    def init_density(self, gas_dist):
        """
        Calculates the initial density in the disk.

        Attributes
        ---------- 
        rho_o : [g cm^-3]
        """
        sigma_o = self.surface_density(r=0*u.AU, gas_dist=gas_dist)
        rho_o = sigma_o / (np.sqrt(2*np.pi) * self.H)
        self.rho_o = rho_o.to(u.g / u.cm**3)


    def density(self):
        """
        Calculates the density at a given distance from the midplane.

        Attributes
        ---------- 
        rho : [g cm^-3]
        """
        fraction = ((self.z/self.H).to(u.AU/u.AU))**2
        rho = self.rho_o * np.exp( -0.5 * fraction )
        self.rho = rho.to(u.g / u.cm**3)
        

    def surface_density(self, r, gas_dist):
        """
        Calculates the surface gas density at a given
        radial distance, r.

        Parameters
        ----------
        gas_dist : float, optional
             The distribution of surface gas density. 
             Default is 2000 [g cm^-2].

        Returns
        ---------- 
        surface density : [g cm^-2]
        """
        r = r.to(u.AU)

        if r.value == 0:
            return gas_dist
        else:
            return gas_dist / r.value


    def optical_depth(self):
        """
        Calculates the optical depth.

        Parameters
        ---------- 
        opacity : opacity of the disk [cm^2 g]

        Attributes
        ---------- 
        tau : [g^2]
        """
        term1   = np.sqrt( (np.pi * self.H**2) / 2)
        errfunc = erf(self.z / (np.sqrt(2) * self.H.to(u.AU)) )
        tau     = term1 * self.opacity * self.rho_o * errfunc
        self.tau = tau.to(u.g**2)


    def temperature(self, zq=3.0, delta=2.0):
        """
        Describes the temperature at some radial distance, r, and
        some vertical distance, z, by a piece-wise function.

        Parameters
        ---------- 
        zq : modifying parameter for the atmosphere of the disk.
             Default = 3.
        delta : factor in defining the temperature profile. 
             Default = 2.

        Attributes
        ---------- 
        T : [K]
        """
        if self.z < zq*self.H:
            fraction = (self.z/self.H).to(u.AU/u.AU)
            sin = np.sin( np.pi / (2*zq) * fraction.value )
            sin = (self.Tatm - self.Tmid) * sin**(2*delta)
            self.T = self.Tmid + sin
        else:
            self.T = self.Tatm
        

    def flare_model(self, amp, t0_ind, rise, fall):
        """
        Creates the flare model.
        """
        def gauss_rise(time, amp, t0, rise):
            return amp * np.exp( -(time - t0)**2.0 / (2.0*rise**2.0) )

        def exp_decay(time, amp, t0, fall):
            return amp * np.exp( -(time - t0) / fall )

        time = self.time.value

        growth = np.where(time <= time[t0_ind])[0]
        decay  = np.where(time >  time[t0_ind])[0]
            
        rise = gauss_rise(time[growth], amp,
                          time[t0_ind], rise)
        fall = exp_decay( time[decay] , amp,
                          time[t0_ind], fall)

        flare = np.append(rise, fall)

        return flare+1


    def luminosity_flare(self, amp=100, t0_ind=200,
                         rise=0.0002, fall=0.03):
        """
        Injects a flare-like shape into the luminosity of the star
        over time. Can pass in a list of all parameters for multiple
        flares to be injected.

        Parameters
        ---------- 
        amp : the amplitude of the flare compared to the base luminosity.
              Default is 100 times brighter during the flare.
        t0_ind : where the peak of the flare will occur. Default is 
                 at the beginning of the time array.
        rise : defines the rise timescale of the flare. Default = 0.0002.
        fall : defintes the decay timescale of the flare. Default = 0.02.
        multi : bool, allows the user to inject multiple flares.
                Will overwrite attributes set in this function.

        Attributes
        ----------
        flare : the flare model injected into the luminosity.
        lum_flare : array of luminosity with flare(s) injection.
        delta_T : changing temperature with the luminosity
        """
        if type(amp) == int or type(amp) == float:
            flare = self.flare_model(amp, t0_ind, rise, fall)
            self.flare = flare
            self.lum_flare = self.Lstar * self.flare
            
        else:
            self.flare = np.ones(self.time.shape)
            self.lum_flare = np.full(self.time.shape, self.Lstar)
            for i in range(len(amp)):
                flare = self.flare_model(amp[i]-1, t0_ind[i], rise[i], fall[i])
                self.flare *= flare

            self.lum_flare = self.Lstar * self.flare

        self.delta_T = self.T * ((self.lum_flare / c.L_sun)**0.25).value


    def UV_flare(self, base=10, factor=100):
        """
        Defines an array of chi (UV Flux) values for 
        the input.ini file.
        
        Parameters
        ---------- 
        base : the base chi value for the UV energy.
        factor : the conversion between white light flare
             energy and energy in the UV. Default = 100.

        Attributes
        ---------- 
        lum_UV : array of chi values [unitless]
        """
        if self.flare is None:
            raise ValueError("Please call disk_model.luminosity_flare() to set flares.")
        else:
            self.lum_UV = ((self.flare-1) * factor) + base

        
    def default_file(self):
        """                                                                                                     
        Creates a default file name if one is not passed in.                                                    
        """
        today = date.today().strftime("%d%m%Y")
        fn = '{0}_Mstar{3}_r{1}_z{2}'.format(today,
                                             self.r.to(u.AU).value,
                                             self.z.to(u.AU).value,
                                             np.round(self.Mstar.to(c.M_sun).value,2))
        self.fn = fn


    def write_path_file(self, fn=None):
        """
        Creates the .out file with parameters.
        Columns are: time, luminosity, UV luminosity, r, z, n, tgas, alpha, tau
        """
        if fn is None:
            self.default_file()
        else:
            self.fn = fn.split('.')[0]

        if self.lum_flare is None:
            fn = self.fn + '_noflare.out'
            luminosity = np.full(self.time.shape, self.Lstar)
            UV_lum = np.full(self.time.shape, 100)
            temp_array = np.full(self.time.shape, self.T)
        else:
            fn = self.fn + '_flare.out'
            luminosity = self.lum_flare

            if self.lum_UV is None:
                raise ValueError("You forgot to update the UV flux with your flare.")
            
            UV_lum = self.lum_UV
            temp_array = self.delta_T
                
        t = Table()
        t.add_column(Column(data=self.time, name='time'))
        t.add_column(Column(data=luminosity, name='lum'))
        t.add_column(Column(data=UV_lum, name='UV_lum'))
        t.add_column(Column(data=np.full(self.time.shape, self.r.to(u.AU).value), name='r'))
        t.add_column(Column(data=np.full(self.time.shape, self.z.to(u.AU).value), name='z'))
        t.add_column(Column(data=np.full(self.time.shape, self.rho.to(u.g/u.cm**3).value), 
                            name='n'))
        t.add_column(Column(data=temp_array, name='tgas'))
        t.add_column(Column(data=np.full(self.time.shape, self.opacity.value), name='alpha'))
        t.add_column(Column(data=np.full(self.time.shape, self.tau.value), name='opt_depth'))

        t.write(fn, format='ascii')
