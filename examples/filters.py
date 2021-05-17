

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from flare import filters




filter_set = filters.HST
filter_set = filters.NIRCam_W + filters.MIRI

filter_set = filters.Euclid + filters.HSC

F = filters.add_filters(filter_set) # --- NOTE: need to give it the redshifted


for i,f in enumerate(F['filters']):

    c = cm.jet(i/len(F['filters']))

    plt.plot(F[f].lam, F[f].T, lw = 1, zorder = 1, c=c) # --- plot filter transmission curve

    plt.scatter(F[f].pivwv(), F[f].Tpeak()*1.1, edgecolor = 'k', zorder = 3, label = f, c=[c]) # --- plot pivot wavelength

    plt.plot([F[f].pivwv()-F[f].rectw()/2., F[f].pivwv()+F[f].rectw()/2.], [F[f].Tpeak()*1.1]*2, lw =2, c='k', zorder = 2) # --- plot rectangular width of filter

    print(f,F[f].pivwv(),F[f].meanwv(), F[f].rectw())


plt.legend()

plt.savefig('filters.pdf')
