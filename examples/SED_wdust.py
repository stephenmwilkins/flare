
import numpy as np
import matplotlib.pyplot as plt 

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import FLARE
from FLARE.SED import models
from FLARE.SED import SFZH
import FLARE.filters






SPS = models.SPS('P2/ModSalpeter_100')




sfzh, sfr = SFZH.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': 8., 'log10Z': -2., 'log10M*': 8.})

print('star formation rate: {0}'.format(sfr))

fesc = 0.0

SED = SPS.get_Lnu(sfzh, {'fesc': fesc, 'log10tau_V': -1.0}, dust_model = 'very_simple')


# --- create observed SED

cosmo = FLARE.default_cosmo()
z = 8.5

SED.total.get_fnu(cosmo, z) # calculate observer frame wavelength
SED.total_intrinsic.get_fnu(cosmo, z) # calculate observer frame wavelength

plt.plot(SED.total_intrinsic.lamz, np.log10(SED.total_intrinsic.fnu), zorder = 1) # plot SED
plt.plot(SED.total.lamz, np.log10(SED.total.fnu), zorder = 1) # plot SED

plt.axvline(1216*(1+z), c='k', lw=1, alpha = 0.5)

filters = FLARE.filters.NIRCam_W
F = FLARE.filters.add_filters(filters, new_lam = SED.total.lamz) # --- NOTE: need to give it the redshifted 
SED.total.get_Fnu(F) # generates Fnu (broad band fluxes)
for f in filters: plt.scatter(F[f].pivwv(), np.log10(SED.total.Fnu[f]), edgecolor = 'k', zorder = 2, label = f)



plt.xlim([5000.,50000.])

mx = np.max(np.log10(SED.total.fnu))
plt.ylim([mx-4., mx+0.3])
plt.show()



