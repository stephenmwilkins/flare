import FLARE

import matplotlib.pyplot as plt

plt.style.use('simple') # --- makes nicer plots

from FLARE.photom import flux_to_L, M_to_lum, lum_to_M
from FLARE.LF import lf_parameters

import numpy as np


def _integ(x, a):
    return x ** (a) * np.exp(-x)


def _integ2(x, a):
    return 0.4*np.log(10)*10**(-0.4*x* (a+1)) * np.exp(-10**(-0.4*x))


def _integ_dblpow(x, a, b):
    return 1 / (10 ** (x*(a+1)) + 10 ** (x*(b+1)))


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


def lf_for_z(z, log10L_bin_centres, models):
    # THIS WILL BE individual LF plot


    fig = plt.figure()

    if type(models) == str:
        models = [models]

    for model in models:

        m = getattr(lf_parameters, model)()

        if hasattr(m, 'beta'):
            s = abs(np.array(m.redshifts) - z) < 0.1

            x = 0.4 * (lum_to_M(10 ** log10L_bin_centres) - np.array(m.M_star)[s][0])
            phistar = 10 ** np.array(m.phi_star)[s][0]
            alpha = np.array(m.alpha)[s][0]
            beta = np.array(m.beta)[s][0]

            phi = phistar / (10 ** (x * (alpha + 1)) + 10 ** (x * (beta + 1)))

            plt.scatter(log10L_bin_centres, np.log10(phi))

        else:

            s = abs(np.array(m.redshifts) - z) < 0.1

            x = 10**(log10L_bin_centres - np.log10(M_to_lum(np.array(m.M_star)[s][0])))
            phistar = 10**np.array(m.phi_star)[s][0]
            alpha = np.array(m.alpha)[s][0]
            plt.scatter(log10L_bin_centres, np.log10(phistar*(x)**(alpha + 1) * np.exp(-x)))

    plt.ylabel(r'$\rm log_{10}(\phi \; / \; cMpc^{-3} \; dex^{-1})$')
    plt.xlabel(r"$\rm log_{10}(L_{UV} \; / \; erg\, s^{-1}\, Hz^{-1})$")

    return fig


def lf_params_zs(zs, log10L_bin_centres, models):
    # THIS WILL BE LF PLOT FOR A SELECTION OF MODELS

    if len(zs)%2 == 1:
        N = int(len(zs)/2) + 1
    else:
        N = int(len(zs) / 2)

    fig, axes = plt.subplots(N, 2, figsize=(6, N*2), sharex=True, sharey=True)


    plt.subplots_adjust(left=0.10, bottom=0.10, top=0.95, right=0.85, wspace=0., hspace=-0.)

    if type(models) == str:
        models = [models]

    for z, ax in zip(zs, axes.flatten()):

        for model in models:

            m = getattr(lf_parameters, model)()

            if hasattr(m, 'beta'):
                s = abs(np.array(m.redshifts) - z) < 0.5

                try:
                    x = 0.4 * (lum_to_M(10 ** log10L_bin_centres) - np.array(m.M_star)[s][0])
                    phistar = 10 ** np.array(m.phi_star)[s][0]
                    alpha = np.array(m.alpha)[s][0]
                    beta = np.array(m.beta)[s][0]

                    phi = phistar / (10 ** (x * (alpha + 1)) + 10 ** (x * (beta + 1)))

                    ax.scatter(log10L_bin_centres, np.log10(phi), s=5)
                except:
                    continue


            else:

                s = abs(np.array(m.redshifts) - z) < 0.5

                try:
                    x = 10**(log10L_bin_centres - np.log10(M_to_lum(np.array(m.M_star)[s][0])))
                    phistar = 10**np.array(m.phi_star)[s][0]
                    alpha = np.array(m.alpha)[s][0]
                    ax.scatter(log10L_bin_centres, np.log10(phistar*(x)**(alpha + 1) * np.exp(-x)), s=5)

                except:
                    continue


        ax.text(0.7, 0.9, r'$\rm z=[{0},{1})$'.format(z - 0.5, z + 0.5), fontsize=8, transform=ax.transAxes)

        ylim = np.array([-8., -4.01])

        ax.set_ylim(-8.)
        #ax.set_xlim([27.25, 30.])


    fig.text(0.5, 0.05, r'$\rm\log_{10}(L_{FUV}/erg\ s^{-1}\ Hz^{-1})$', horizontalalignment='center', verticalalignment='center')
    fig.text(0.05, 0.5, r'$\rm\log_{10}(\phi/Mpc^{-3})$', horizontalalignment='center', verticalalignment='center', rotation=90)



    return
