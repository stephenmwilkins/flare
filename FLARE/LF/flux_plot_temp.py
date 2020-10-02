import FLARE
from FLARE.photom import lum_to_flux

import matplotlib.pyplot as plt

import numpy as np


def flux_sample_bin_plot(samples, cosmo=False, bins=100, range=[[0., 15.], [0., 3.]], save_file=False):
    # --- make nice plot

    if not cosmo: cosmo = FLARE.core.default_cosmo()

    fig = plt.figure(figsize=(6, 5))

    z_array = samples['z']
    f_array = lum_to_flux(10. ** samples['log10L'], cosmo, z_array) * (10. ** 9)

    cm = plt.get_cmap('plasma')

    plt.hist2d(z_array, np.log10(f_array), bins=bins, range=range, cmap=cm)

    bar = plt.colorbar(orientation='vertical')
    bar.set_label(r'$\rm N$', rotation=90)

    plt.ylabel(r"$\rm log_{10}(f_{\nu}/ [nJy])$")
    plt.xlabel(r"$\rm z$")
    plt.xlim(range[0][0], range[0][-1])
    plt.ylim(range[1][0], range[1][-1])

    if save_file == False:
        return fig
    else:
        plt.savefig(save_file + '.png', dpi=300)