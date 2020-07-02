
import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import FLARE
from FLARE.SED import dust_curves



lam = np.arange(1000.,10000.,1)

for curve in dust_curves.curves:

    plt.plot(np.log10(lam), getattr(dust_curves, curve)().tau(lam), label = curve)



plt.legend()

plt.axhline(1.0, c='k', alpha = 0.2)
plt.axvline(np.log10(5500.), c='k', alpha = 0.2)

plt.xlabel(r'$\log_{10}(\lambda/\AA)$')
plt.ylabel(r'$\tau/\tau_V$')

plt.show()
