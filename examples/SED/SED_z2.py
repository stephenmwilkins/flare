
import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import FLARE
from FLARE.SED import models
from FLARE.SED import SFZH
import FLARE.photom
import FLARE.filters
import FLARE.observatories





# SPS = models.SPS('P2/ModSalpeter_100')
SPS = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')



sfzh, sfr = SFZH.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': 9.5, 'log10Z': -2., 'log10M*': 7.5})

print('star formation rate: {0}'.format(sfr))

SED = SPS.get_Lnu(sfzh, {'fesc': 0.0}, dust = False)




# --- create observed SED

cosmo = FLARE.default_cosmo()
z = 1

SED.total.get_fnu(cosmo, z) # calculate observer frame wavelength

plt.plot(SED.total.lamz, np.log10(SED.total.fnu), zorder = 1) # plot SED


filters = ['Webb.NIRCam.F200W']

F = FLARE.filters.add_filters(filters, new_lam = SED.total.lamz) # --- NOTE: need to give it the redshifted
SED.total.get_Fnu(F) # generates Fnu (broad band fluxes)
for f in filters:
    print(f, FLARE.photom.flux_to_m(SED.total.Fnu[f]))
