

import numpy as np


import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from flare import filters
from flare import observatories





cmap = mpl.cm.rainbow
norm = mpl.colors.Normalize(vmin=3, vmax=5.5)


Observatories = ['Hubble', 'Webb', 'Euclid', 'Roman', 'Spitzer','VISTA']


print(observatories.Webb.NIRCam.Short.filters) # a single channel (most instruments have a single channel)
print(observatories.Webb.NIRCam.filters) # includes all channels
print(observatories.Webb.filters) # includes all instruments


fig, axes = plt.subplots(len(Observatories), 1, sharex=True, figsize=(6, 6))


for Observatory, ax in zip(Observatories, axes):

    obs = observatories.observatories[Observatory]

    ax.text(-0.6, 0.7, Observatory)

    for inst in obs.instruments:

        F = filters.add_filters(obs.instrument[inst].get_filter_IDs(exclude='M'))

        for i,f in enumerate(F['filters']):

            pivwv = F[f].pivwv()
            c = cmap(norm(np.log10(pivwv)))
            ax.plot(np.log10(F[f].lam)-4, F[f].T, lw = 1, zorder = 1, c=c) # --- plot filter transmission curve
            ax.set_ylim([0,1.1])

    axes[-1].set_xlabel(r'$\rm log_{10}(\lambda/\mu m) $')


fig.savefig('all_filters.pdf')
