

# --- in this code hard-code various luminosity function evolution models. Starting with bluetides.






class linear():


    def __init__(self, lp, z=8):
    
        # lp is a dictionary of the parameters of the linear evolution model

        self.lp = lp


    def parameters(self,z):
    
        # use linear evolution model
        # get parameters as a function of z
        # returns a dictionary of parameters
        
        return p
        
    def N(self, cosmo = False, redshift_limits = [8., 15.], log10L_limits = [28., 32.], dz = 0.05, dlog10L = 0.05):
    
        # calculate the number of galaxies in each bin on a grid defined by redshift_limits, log10L_limits, dz, dlog10L
    
        return bin_edges, N
    
    
    def sample(self, area = 1., cosmo = False, redshift_limits = [8., 15.], log10L_min = 28., seed = False):
    
        # sample the LF evolution model in a given volume or area
        # if area will need to also give a cosmology 
        
        if not cosmo: cosmo = FLARE.default_cosmology()
    
        
        return redshifts, L
    
    
    def bin_sample(self, redshifts, L, redshift_limits = [8., 15.], log10L_limits = [28., 32.], dz = 0.05, dlog10L = 0.05)
    
        # bin the sample
    
        return bin_edges, N
    
  
  
  
  
class existing_model()
  
    def __init__(self):
    
        self.calculate_linear_evolution_coeffs(self, z=8)
  
    def interpolate_parameters(self,z):
    
        # get parameters as a function of z
        # returns a dictionary of the Schechter function parameters
        
        return p
  
    def calculate_linear_evolution_coeffs(self, z=8):
    
        # calculate the linear evolution coeffs using z as the normalisation
    
        # add a check to quantify how linear it actually is

        self.lp =     
    
  
  

  
    
class bluetides(existing_model): # --- based on bluetides simulation

    self.redshifts = # array of redshifts
    self.phi_star = # array of phi_star value to interpolate
    self.M_star = #
    self.alpha = #
        
       
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
        
            