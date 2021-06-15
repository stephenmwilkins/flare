

import numpy as np


import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

params = {"ytick.color" : "w",
          "xtick.color" : "w",
          "axes.labelcolor" : "w",
          "axes.edgecolor" : "w"}
plt.rcParams.update(params)

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# import flare.plt
from flare import filters
from flare import observatories





cmap = mpl.cm.rainbow
norm = mpl.colors.Normalize(vmin=200, vmax=1500)


observatory = 'Hubble'




fig = plt.figure(figsize = (8, 2.5))

left  = 0.15
bottom = 0.25
height = 0.7
width = 0.8

ax = fig.add_axes((left, bottom, width, height))

obs = observatories.observatories[observatory]

for inst in obs.instruments:

    F = filters.add_filters(obs.instrument[inst].get_filter_IDs(exclude='M'))

    for i,f in enumerate(F['filters']):

        pivwv = F[f].pivwv()
        c = cmap(norm(pivwv/10))
        # ax.plot(F[f].lam/10, F[f].T, lw = 1, zorder = 1, c=c) # --- plot filter transmission curve

        ax.fill_between(F[f].lam/10, 0.0* F[f].T, F[f].T, lw = 1, zorder = 1, color='k', alpha=0.5) # --- plot filter transmission curve
        ax.plot(F[f].lam/10, F[f].T, lw = 1, zorder = 1, c='w') # --- plot filter transmission curve

ax.set_ylim([0,1.1])
ax.set_xlabel(r'$\rm wavelength/nm$')
ax.set_ylabel(r'$\rm transmission$')

fig.savefig(f'filters_{observatory}.pdf', transparent=True)
fig.savefig(f'filters_{observatory}.png', transparent=True)
