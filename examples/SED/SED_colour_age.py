
import numpy as np
import matplotlib.pyplot as plt 

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import FLARE
from FLARE.SED import models
from FLARE.SED import SFZH
import FLARE.filters

z = 8

cosmo = FLARE.default_cosmo()

SPS = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')
filters = FLARE.filters.HST + FLARE.filters.Spitzer
F = FLARE.filters.add_filters(filters, new_lam = SPS.lam * (1. + z)) # --- NOTE: need to give it the redshifted 

cxf = ['HST.WFC3.f125w', 'Spitzer.IRAC.ch1']
cyf = ['Spitzer.IRAC.ch1', 'Spitzer.IRAC.ch2']


sfzh = SFZH.simple(SPS, {'log10age': 7.05, 'log10Z': -2.})

plt.imshow(sfzh.T)
plt.show()


# 
# # ---- constant
# 
# for tau_V in [0.0, 0.1, 0.25, 0.5]:
# 
#     SPS.dust = {'model': 'simple', 'slope': -1.0, 'tau_V': tau_V}
# 
#     for a0 in [7.,8.,9.]:
# 
#         sfzh, sfr = SFZH.constant(SPS, {'log10_duration': a0, 'log10Z': -2., 'log10M*': 8.})
# 
#         SED = SPS.get_Lnu(sfzh)
# 
#         # --- create observed SED
# 
#         SED.total.get_fnu(cosmo, z=8.0) # calculate observer frame wavelength
# 
#         SED.total.get_Fnu(F) # generates Fnu (broad band fluxes)
# 
#         fluxes = SED.total.Fnu
# 
#         cx = -2.5*np.log10(fluxes[cxf[0]]/fluxes[cxf[1]])
#         cy = -2.5*np.log10(fluxes[cyf[0]]/fluxes[cyf[1]])
#         
#         print(a0,cx,cy, SED.A1500())
# 
#     
