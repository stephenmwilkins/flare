


import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import cm 

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import FLARE.LF
from FLARE.LF import evo

#cosmo = FLARE.default_cosmology() # WMAP9

import astropy.cosmology
cosmo = astropy.cosmology.default_cosmology.get()

# --- simple LF examples


#p = {'alpha': -2., 'log10L*': 29., 'log10phi*': -2.0}

#LF = FLARE.LF.Schechter(p)


#LF.phi(L_bins) # --- return list of phi values (Mpc^-3 dex^-1) for a given set of L_bins


# --- LF evolution examples


evo_model = evo.linear(evo.bluetides()) # initialises the linear evolution model with the bluetides parameters


p = evo_model.parameters(z = 8.5) # return model parameters at z=8.5


# --- return grid with number of galaxies in each bin

bin_edges, N = evo_model.N(cosmo = cosmo, redshift_limits = [8., 15.], log10L_limits = [27., 30.], dz = 0.05, dlog10L = 0.05)

# <<<<< make plot


evo.evo_plot(bin_edges, N)
plt.show()

# --- return a random sample of redshifts and luminosities
'''
redshifts, L = evo_model.sample(area = 1., cosmo = cosmo, redshift_limits = [8., 15.], log10L_limits = [28., 32.])


bin_edges, N = bin_sample(redshifts, L, redshift_limits = [8., 15.], log10L_limits = [28., 32.])
'''
# <<<<< make plot