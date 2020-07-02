
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import FLARE
from FLARE.SED import models
from FLARE.SED import SFZH


SPS = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')

N = 100

for i, log10age in enumerate(np.linspace(6.0,8.0,N)):

    print(i)

    sfzh = SFZH.simple(SPS, {'log10age': log10age, 'log10Z': -2.})

    SED = SPS.get_Lnu(sfzh)

    plt.plot(np.log10(SED.stellar.lam), np.log10(SED.stellar.lnu), c = cm.viridis(i/N), lw=1, alpha = 0.5)

    

    
mx = np.max(np.log10(SED.stellar.lnu)) # only last one
plt.xlim([2.7, 4.])
plt.ylim([mx-4., mx+3.0])
plt.show()


