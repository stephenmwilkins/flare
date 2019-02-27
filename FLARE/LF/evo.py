

# --- in this code hard-code various luminosity function evolution models. Starting with bluetides.
from scipy.stats import linregress
import numpy as np




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

        
    def N(self, area = 1., cosmo = False, redshift_limits = [8., 15.], log10L_limits = [28., 32.], dz = 0.05, dlog10L = 0.05):

        # calculate the number of galaxies in each bin on a grid defined by redshift_limits, log10L_limits, dz, dlog10L
        if not cosmo: cosmo = FLARE.default_cosmology()
        
    
        return bin_edges, N
    
    
    def sample(self, area = 1., cosmo = False, redshift_limits = [8., 15.], log10L_min = 28., seed = False):
    
        # sample the LF evolution model in a given volume or area
        # if area will need to also give a cosmology 
        
        if not cosmo: cosmo = FLARE.default_cosmology()
    
        
        return redshifts, L
    
    
    def bin_sample(self, redshifts, L, redshift_limits = [8., 15.], log10L_limits = [28., 32.], dz = 0.05, dlog10L = 0.05)
    
        # bin the sample
    
        return bin_edges, N
    
  
  
  
  
class existing_model:
  
    def __init__(self, model = bluetides()):
        self.model = model
        self.lp = {}
        
    def interpolate_parameters(self, z=8.):
    
        # get parameters as a function of z
        # returns a dictionary of the Schechter function parameters
        z_mod = self.model.redshifts
        alpha_mod = self.model.alpha
        log10phi_mod = self.model.phi_star
        log10M_mod = self.model.M_star
        p = {'alpha': np.interp(z, z_mod, alpha_mod), 'log10phi*': np.interp(z, z_mod, log10phi_mod), 'log10M*': np.interp(z, z_mod, log10M_mod)}
        return p
  
    def calculate_linear_evolution_coeffs(self):
    
        # calculate the linear evolution coeffs using z as the normalisation
        z_mod = self.model.redshifts
        alpha_mod = self.model.alpha
        log10phi_mod = self.model.phi_star
        log10M_mod = self.model.M_star

        fit_alpha = linregress(z_mod, alpha_mod)
        fit_log10phi = linregress(z_mod, log10phi_mod)
        fit_log10M = linregress(z_mod, log10M_mod)
        # add a check to quantify how linear it actually is

        self.lp =  {'alpha': fit_alpha, 'log10phi*': fit_log10phi, 'log10M*': fit_log10M}

        return self.lp
        
  
  

  
    
class bluetides: # --- based on bluetides simulation

    def __init__(self):
        self.redshifts = [8.0, 9.0, 10.0, 11.0, 12.0, 13.0]            # array of redshifts
        self.phi_star = [-3.92, -4.2, -4.7, -4.79, -5.09, -5.71]        # array of phi_star value to interpolate
        self.M_star = [-20.93, -20.68, -20.69, -20.17, -19.92, -19.91]  #
        self.alpha = [-2.04, -2.1, -2.27, -2.27, -2.35, -2.54] 
            
           
# class Mason15(existing_model): # --- based on Mason et al. (2015)
# 
#     self.redshifts = # array of redshifts
#     self.phi_star = # array of phi_star value to interpolate
#     self.M_star = #
#     self.alpha = #        

















def evo_plot(bin_edges, N, cosmo = False, f_limits = False, save_file = False):

    # --- make nice plot
    
    
    
    if not cosmo: cosmo = FLARE.default_cosmology()    
    
    
    # --- draw lines of constant flux
    
    if f_limits:
    
        for f_limit in f_limits:
        
            