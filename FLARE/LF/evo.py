# ---
from scipy.stats import linregress
import numpy as np
import scipy.integrate as cp

import astropy.cosmology
import astropy.units as u


# FOR NOW: The cosmology is forced to be the default astropy cosmology.
cosmo = astropy.cosmology.default_cosmology.get()

geo = (4. * np.pi * (100. * 10. * 3.0867 * 10 ** 16) ** 2)  # factor relating the L to M in cm^2


def L(M):
    return 10 ** (-0.4 * (M + 48.6)) * geo


def M_to_log10L(M):
    return -0.4 * (M + 48.6) + np.log10(geo)


def M(log10L):
    return -2.5 * (log10L - np.log10(geo)) - 48.6


def dVc(z):
    return cosmo.differential_comoving_volume(z).value


def _integ(x,a):
    return x**(a-1) * np.exp(-x)


def CulmPhi(log10L, sp):
    y = log10L - sp['log10L*']
    x = 10 ** y
    alpha = sp['alpha']

    gamma = cp.quad(_integ, x, np.inf, args=alpha + 1)[0]
    num = gamma * sp['phi*']

    return num


def CDF(log10L_limit, normed=True):
    log10Ls = np.arange(self.sp['log10L*'] + 5., log10L_limit - 0.01, -0.01)

    CDF = np.array([self.CulmPhi(log10L) for log10L in log10Ls])

    if normed: CDF /= CDF[-1]

    return log10Ls, CDF


def N_exact(volume, bin_edges, sp):
    # --- return the exact number of galaxies expected in each bin

    CulmN = np.array([CulmPhi(x, sp) for x in bin_edges]) * volume

    return -(CulmN[1:] - CulmN[0:-1])


def bin(log10L_sample, bins):
    # --- bins can either be the number of bins or the bin_edges

    N_sample, bin_edges = np.histogram(log10L_sample, bins=bins, normed=False)

    return N_sample



class linear:


    def __init__(self, lp, z=8):
    
        # lp is a dictionary of the parameters of the linear evolution model

        self.lp = lp


    def parameters(self,z):
    
        # use linear evolution model
        # get parameters as a function of z
        # returns a dictionary of parameters
        p = {}
        for param in self.lp:
            p[param] = self.lp[param][0]*z + self.lp[param][1]
  
        return p

        
    def N(self, area = 1., cosmo = cosmo, redshift_limits = [8., 15.], log10L_limits = [27.5, 30.], dz = 0.05, dlog10L = 0.05):

        # calculates the number of galaxies in each bin on a grid defined by redshift_limits, log10L_limits, dz, dlog10L
        # and area based on a luminosity function evolution model.
        
        area_sm = area                      # Area in square arcmin
        area_sd = area_sm / 3600.           # Area in square degrees
        area_sr = (np.pi/180.)**2 * area_sd # Area in steradian


        if not cosmo: cosmo = FLARE.default_cosmology()

        # Setting the bin edges as well as centres for later operations
        bin_edges = {'log10L': np.arange(log10L_limits[0],log10L_limits[-1]+dlog10L,dlog10L), 'z': np.arange(redshift_limits[0],redshift_limits[-1]+dz,dz)}
        bin_centres = {'log10L': np.arange(bin_edges['log10L'][0]+dlog10L/2.,bin_edges['log10L'][-1]-dlog10L/2.,dlog10L), 'z': np.arange(bin_edges['z'][0]+dz/2.,bin_edges['z'][-1]-dz/2.,dz)}

        # Using astropy.cosmology to calculate the volume in each redshift bin
        volumes = np.asarray([ cp.quad(dVc, bin_edges['z'][i-1], bin_edges['z'][i])[0] for i in range(1,len(bin_edges['z']))])

        # Initialising the output array
        N = np.zeros((len(bin_centres['log10L']), len(bin_centres['z'])))

        # Loop calculates LF for each input z (bin centres) and returns the exact numbers expected in each bin
        # (There may be a better option for generating this)
        for i in range(len(bin_centres['z'])):
            params = self.parameters(bin_centres['z'][i])

            sp = {}
            sp['alpha'] = params['alpha']
            sp['phi*'] = 10**params['log10phi*']
            sp['log10L*'] = M_to_log10L(params['M*'])

            N_ext = N_exact(volumes[i] * area_sr, bin_edges['log10L'], sp)

            for j in range(len(N_ext)):
                N[j,i] = N_ext[j]

        return bin_edges, N

'''    
    
    def sample(self, area = 1., cosmo = cosmo, redshift_limits = [8., 15.], log10L_min = 28., seed = False):
    
        # sample the LF evolution model in a given volume or area
        # if area will need to also give a cosmology 
        
        #if not cosmo: cosmo = FLARE.default_cosmology()
    
        
        return 1. #redshifts, L
    

    def bin_sample(self, redshifts, L, redshift_limits = [8., 15.], log10L_limits = [28., 32.], dz = 0.05, dlog10L = 0.05)
    
        # bin the sample
    
        return 1. #bin_edges, N
    
  
'''
  
  
class existing_model:
  
    def __init__(self, model = bluetides()):
        self.model = model
        self.lp = {}
        
    def interpolate_parameters(self, z=8.):
    
        # interpolates parameters as a function of z
        # returns a dictionary of the Schechter function parameters for given redshift(s)

        z_mod = self.model.redshifts
        alpha_mod = self.model.alpha
        log10phi_mod = self.model.phi_star
        log10M_mod = self.model.M_star
        p = {'alpha': np.interp(z, z_mod, alpha_mod), 'log10phi*': np.interp(z, z_mod, log10phi_mod), 'M*': np.interp(z, z_mod, log10M_mod)}

        return p
  
    def calculate_linear_evolution_coeffs(self):

        # Function that calculates the linear evolution coeffs
        # returns a dictionary of linear model coefficients and goodness of fit

        z_mod = self.model.redshifts
        alpha_mod = self.model.alpha
        log10phi_mod = self.model.phi_star
        M_mod = self.model.M_star

        # The output contains full linregress output (0th and 1st element contain the slope and intercept respectively)
        fit_alpha = linregress(z_mod, alpha_mod)
        fit_log10phi = linregress(z_mod, log10phi_mod)
        fit_M = linregress(z_mod, M_mod)

        self.lp =  {'alpha': fit_alpha, 'log10phi*': fit_log10phi, 'M*': fit_M}

        return self.lp
        
  
  

  
    
class bluetides: # --- based on bluetides simulation

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form
        self.redshifts = [8.0, 9.0, 10.0, 11.0, 12.0, 13.0]            # array of redshifts
        self.phi_star = [-3.92, -4.2, -4.7, -4.79, -5.09, -5.71]       # array of log10(phi_star) values
        self.M_star = [-20.93, -20.68, -20.69, -20.17, -19.92, -19.91] # array of M_star values
        self.alpha = [-2.04, -2.1, -2.27, -2.27, -2.35, -2.54]         # array of alpha values
            

class Ma2019:
    # --- LF evolution model based on Ma et al. (2019) (f_dust = 0.8)

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form
        self.redshifts = [5., 6., 7., 8.0, 9.0, 10.0]                  # array of redshifts
        self.phi_star = [-3.55, -3.44, -4.09, -3.98, -4.57, -4.74]     # array of log10(phi_star) values
        self.M_star = [-21.77, -21.34, -21.73, -20.97, -21.30, -20.90] # array of M_star values
        self.alpha = [-1.9, -1.87, -2.05, -2.08, -2.20, -2.31]         # array of alpha values

# class Mason15(existing_model): # --- based on Mason et al. (2015)
# 
#     self.redshifts = # array of redshifts
#     self.phi_star = # array of phi_star value to interpolate
#     self.M_star = #
#     self.alpha = #        
















'''
def evo_plot(bin_edges, N, cosmo = False, f_limits = False, save_file = False):

    # --- make nice plot
    
    
    
    if not cosmo: cosmo = FLARE.default_cosmology()    
    
    
    # --- draw lines of constant flux
    
    if f_limits:
    
        for f_limit in f_limits:
        
'''