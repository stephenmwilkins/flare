import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import FLARE.LF
from FLARE.LF import evo, completeness, lf_parameters, LF_plots
from FLARE.photom import m_to_flux

import astropy.cosmology

cosmo = astropy.cosmology.default_cosmology.get()


# --- LF evolution examples

# select LF parameters
m = getattr(lf_parameters, 'FLARES')()

evo_model = evo.linear(m) # initialises the linear evolution model with the bluetides parameters


p = evo_model.parameters(z = 8.5) # return model parameters at z=8.5
print(p)


# --- return grid with number of galaxies in each bin

bin_edges, bin_centres, N = evo_model.N(cosmo = cosmo, redshift_limits = [8., 15.], log10L_limits = [27., 30.], dz = 0.05, dlog10L = 0.05)

# <<<<< make plot
LF_plots.evo_plot(bin_edges, N)
plt.show()


# --- return grid with simple error function completeness model

flux_limit = 10 # in nanojansky



c = evo_model.completeness_erf(bin_centres, flux_limit, stretch=1000, cosmo=cosmo)


print(c.shape)



X, Y = np.meshgrid(bin_centres['z'], bin_centres['log10L'])
c_plot = plt.pcolormesh(X, Y, c)
c_bar = plt.colorbar(c_plot)
c_bar.set_label('completeness', rotation=90)
plt.ylabel(r"$\rm log_{10}(L_{\nu} \; / \; erg\, s^{-1}\, Hz^{-1})$")
plt.xlabel(r"$\rm z$")
plt.show()


# --- finally, convolve the two to generate the LFE model with simple completeness cut.

LFE = np.multiply(N, c)

# and plot:

LF_plots.evo_plot(bin_edges, LFE)
plt.show()


# --- Extra example for using a sampled completeness

flux_limit = 10 # in nanojansky

# bin_centres should be the same ones used for N above (this may take a while if samples at each bin is high)
c = completeness.completeness_sample(bin_centres, flux_limit, stretch=1, cosmo=cosmo, N_samples=100)

X, Y = np.meshgrid(bin_centres['z'], bin_centres['log10L'])
c_plot = plt.pcolormesh(X, Y, c)
c_bar = plt.colorbar(c_plot)
c_bar.set_label('completeness', rotation=90)
plt.ylabel(r"$\rm log_{10}(L_{\nu} \; / \; erg\, s^{-1}\, Hz^{-1})$")
plt.xlabel(r"$\rm z$")
plt.show()


# --- finally, convolve the two to generate the LFE model with simple completeness cut.

LFE = np.multiply(N, c)

# and plot:

LF_plots.evo_plot(bin_edges, LFE)
plt.show()
