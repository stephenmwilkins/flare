


import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import cm 

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import FLARE
import FLARE.LF
from FLARE.LF import evo
from FLARE.photom import m_to_flux

cosmo = FLARE.default_cosmo() # WMAP9

f_limit_deep = FLARE.photom.m_to_flux(26.)

area = 60.*60.*40.

print('area of Euclid deep: {0} arcmin2'.format(area))




evo_model = evo.linear(evo.bluetides()) # initialises the linear evolution model with the bluetides parameters
bin_edges, bin_centres, N = evo_model.N(cosmo = cosmo, redshift_limits = [8., 15.], log10L_limits = [27., 30.], dz = 0.1, dlog10L = 0.01, flux_min = f_limit_deep)
n = np.sum(N)
print('density per arcmin: {0:9.2f}'.format(n))
print('number in Eudlid deep: {0:9.1f}'.format(n*area))

evo.evo_plot(bin_edges, N, f_limits=[f_limit_deep])
plt.show()


evo_model = evo.interp(evo.bluetides()) # initialises the linear evolution model with the bluetides parameters
bin_edges, bin_centres, N = evo_model.N(cosmo = cosmo, redshift_limits = [8., 15.], log10L_limits = [27., 30.], dz = 0.1, dlog10L = 0.01, flux_min = f_limit_deep)
n = np.sum(N)
print('density per arcmin: {0:9.2f}'.format(n))
print('number in Eudlid deep: {0:9.1f}'.format(n*area))

evo.evo_plot(bin_edges, N, f_limits=[f_limit_deep])
plt.show()

