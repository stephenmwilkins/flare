

import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import flare
import flare.filters

import SynthObs.Morph.PSF as PSF


for filters, observatory in zip([flare.filters.WFC3NIR_W, flare.filters.NIRCam_W, flare.filters.Euclid_NISP],['Hubble', 'Webb', 'Euclid']):

    fig, axes = plt.subplots(1, len(filters), figsize = (len(filters)*2, 2))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.0, wspace=0.0, hspace=0.0)

    for i,f in enumerate(filters):

        psf = getattr(PSF, observatory+'PSF')(f)

        axes[i].imshow(np.log10(psf.data))
        axes[i].get_xaxis().set_ticks([])
        axes[i].get_yaxis().set_ticks([])

    plt.show()
