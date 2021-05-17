
import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from flare.SED import models
import flare.filters

# import flare









# --- create a simple "beta" model with a power law slope and a break at 912.

m = models.beta(np.arange(0., 2000., 1.), -2.5, 1E28, normalisation_wavelength = 1500.)


plt.plot(m.lam, m.lnu,  zorder = 1)

plt.show()

# --- now move to the observer frame

cosmo = flare.default_cosmo()

z = 8.

m.get_fnu(cosmo, z, include_IGM = True) # --- generates redshifted wavelength grid and fnu/nJy

plt.plot(m.lamz, m.fnu,  zorder = 1)
plt.savefig('beta1.png')

filters = flare.filters.HST

F = flare.filters.add_filters(filters, new_lam = m.lam * (1. + z)) # --- NOTE: need to give it the redshifted

m.get_Fnu(F) # generates Fnu (broad band fluxes)

for f in filters: plt.scatter(F[f].pivwv(), m.Fnu[f], edgecolor = 'k', zorder = 2, label = f)


plt.legend()

plt.show()
plt.savefig('beta2.png')
