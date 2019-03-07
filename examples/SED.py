
import numpy as np
import matplotlib.pyplot as plt 

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import FLARE
from FLARE.SED import models
from FLARE.SED import SFZH
import FLARE.filters







SPS = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')

sfzh, sfr = SFZH.constant(SPS, {'log10_duration': 8., 'log10Z': -2., 'log10M*': 8.})

print('star formation rate: {0}'.format(sfr))

# plt.imshow(sfzh.T)
# plt.show()

SED = SPS.get_Lnu(sfzh)

plt.plot(np.log10(SED.stellar.lam), np.log10(SED.stellar.lnu))
plt.plot(np.log10(SED.total.lam), np.log10(SED.total.lnu))

plt.xlim([2.7, 4.])

mx = np.max(np.log10(SED.total.lnu))
plt.ylim([mx-4., mx+0.3])
plt.show()




# --- create observed SED

cosmo = FLARE.default_cosmo()

SED.total.get_fnu(cosmo, z=8.0) # calculate observer frame wavelength

plt.plot(SED.total.lamz, np.log10(SED.total.fnu), zorder = 1) # plot SED


filters = FLARE.filters.HST + FLARE.filters.Spitzer
F = FLARE.filters.add_filters(filters, new_lam = SED.total.lamz) # --- NOTE: need to give it the redshifted 
SED.total.get_Fnu(F) # generates Fnu (broad band fluxes)
for f in filters: plt.scatter(F[f].pivwv(), np.log10(SED.total.Fnu[f]), edgecolor = 'k', zorder = 2, label = f)



plt.xlim([5000.,50000.])

mx = np.max(np.log10(SED.total.fnu))
plt.ylim([mx-4., mx+0.3])
plt.show()



