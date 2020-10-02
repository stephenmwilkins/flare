import FLARE

import matplotlib.pyplot as plt

from FLARE.photom import flux_to_L

import numpy as np


def evo_plot(bin_edges, N, cosmo=False, f_limits=False, save_file=False):
    # --- make nice plot

    if not cosmo: cosmo = FLARE.core.default_cosmo()

    fig = plt.figure(figsize=(6, 5))

    X, Y = np.meshgrid(bin_edges['z'], bin_edges['log10L'])

    cm = plt.get_cmap('plasma')
    plt.pcolormesh(X, Y, np.log10(N), cmap=cm)

    # --- draw lines of constant flux

    if type(f_limits) is list or type(f_limits) is np.ndarray or type(f_limits) is range:

        for f_limit in f_limits:
            plt.plot(bin_edges['z'], np.log10(flux_to_L(f_limit, cosmo, bin_edges['z'])), 'k--', alpha=0.8)

    if type(f_limits) is float:
        plt.plot(bin_edges['z'], np.log10(flux_to_L(f_limits, cosmo, bin_edges['z'])), 'k--', alpha=0.8)

    bar = plt.colorbar(orientation='vertical')
    bar.set_label(r'$\rm log_{10}(N \; / \; arcmin^{-2})$', rotation=90)

    plt.ylabel(r"$\rm log_{10}(L_{\nu} \; / \; erg\, s^{-1}\, Hz^{-1})$")
    plt.xlabel(r"$\rm z$")
    plt.xlim(min(bin_edges['z']), max(bin_edges['z']))
    plt.ylim(min(bin_edges['log10L']), max(bin_edges['log10L']))

    if save_file == False:
        return fig
    else:
        plt.savefig(save_file + '.png', dpi=300)


def LF(z, log10L_bins, model):
    # THIS WILL BE individual LF plot
    return


def LF():
    # THIS WILL BE LF PLOT FOR A SELECTION OF MODELS
    return
