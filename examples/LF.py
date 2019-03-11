


import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import cm 

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import FLARE.LF
from FLARE.LF import evo
from FLARE.photom import m_to_flux

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

bin_edges, bin_centres, N = evo_model.N(cosmo = cosmo, redshift_limits = [8., 15.], log10L_limits = [27., 30.], dz = 0.05, dlog10L = 0.05)

# <<<<< make plot
evo.evo_plot(bin_edges, N)
plt.show()

# Plot with f_limits:
m_limits = np.arange(26.,31.,1.)
f_limits = (FLARE.photom.m_to_flux(m_limits)/1E9).tolist()
evo.evo_plot(bin_edges, N, f_limits=f_limits)
plt.show()

# Plot of binned samples area 100 sq arcmin
bin_edges, bin_centres, N_binned = evo_model.bin_sample(area=100.)

evo.evo_plot(bin_edges, N_binned)
plt.show()

# --- return a random sample of redshifts and luminosities
sample = evo_model.sample(area = 100.)

# Bin sample and plot
evo.flux_sample_bin_plot(sample)
plt.show()