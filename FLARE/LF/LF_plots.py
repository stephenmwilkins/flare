import FLARE

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mpl

plt.style.use('simple') # --- makes nicer plots

from FLARE.photom import flux_to_L, M_to_lum, lum_to_M
from FLARE.LF import lf_parameters

import numpy as np


def logerr_lo(phi, sigma, ulims):
    err = np.log10(phi) - np.log10(phi - sigma)
    for i, x in enumerate(err):
        if not np.isfinite(x):
            err[i] = 100.
        if ulims[i] == True:
            err[i] = 0.5
    return err


def logerr_hi(phi, sigma):
    err = np.log10(phi + sigma) - np.log10(phi)
    for i, x in enumerate(err):
        if not np.isfinite(x):
            err[i] = 100.
    return err

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


def lf_for_z(z, log10L_bin_centres, models, legend=False):
    # THIS WILL BE individual LF plot

    colors = ['b', 'r', 'orange', 'g', 'k', 'c', 'y', 'm', 'darkslateblue', 'gray', 'tomato', 'lawngreen', 'teal', 'wheat']

    fig = plt.figure()

    if type(models) == str:
        models = [models]

    for i, model in enumerate(models):

        m = getattr(lf_parameters, model)()

        if hasattr(m, 'beta'):
            s = abs(np.array(m.redshifts) - z) < 0.1

            x = 0.4 * (lum_to_M(10 ** log10L_bin_centres) - np.array(m.M_star)[s][0])
            phistar = 10 ** np.array(m.phi_star)[s][0]
            alpha = np.array(m.alpha)[s][0]
            beta = np.array(m.beta)[s][0]

            phi = phistar / (10 ** (x * (alpha + 1)) + 10 ** (x * (beta + 1)))

            plt.scatter(log10L_bin_centres, np.log10(phi), c=colors[i], label=model)

        else:

            s = abs(np.array(m.redshifts) - z) < 0.1

            x = 10**(log10L_bin_centres - np.log10(M_to_lum(np.array(m.M_star)[s][0])))
            phistar = 10**np.array(m.phi_star)[s][0]
            alpha = np.array(m.alpha)[s][0]
            plt.scatter(log10L_bin_centres, np.log10(phistar*(x)**(alpha + 1) * np.exp(-x)), c=colors[i], label=model)

    plt.ylabel(r'$\rm log_{10}(\phi \; / \; cMpc^{-3} \; dex^{-1})$')
    plt.xlabel(r"$\rm log_{10}(L_{UV} \; / \; erg\, s^{-1}\, Hz^{-1})$")

    if legend:
        plt.legend(loc='best')

    return fig


def lf_params_zs(zs, log10L_bin_centres, models, legend=False):
    # THIS WILL BE LF PLOT FOR A SELECTION OF MODELS

    colours = ['b', 'r', 'orange', 'g', 'k', 'c', 'y', 'm', 'darkslateblue', 'gray', 'tomato', 'lawngreen', 'teal', 'wheat']

    if len(zs)%2 == 1:
        N = int(len(zs)/2) + 1
    else:
        N = int(len(zs) / 2)

    fig, axes = plt.subplots(N, 2, figsize=(6, N*2), sharex=True, sharey=True)


    plt.subplots_adjust(left=0.10, bottom=0.10, top=0.95, right=0.85, wspace=0., hspace=-0.)

    if type(models) == str:
        models = [models]

    Nxx = len(models)

    for z, ax in zip(zs, axes.flatten()):

        if legend:
            if Nxx > 7:
                if ax == axes[0, 0]:
                    try:
                        g = []
                        for i, model in enumerate(models[:7]):
                            colors_plot = colours[:7]
                            g.append(ax.scatter([-99.], [-99.], c=colors_plot[i], s=5, label=model))
                        ax.legend(handles=g, loc='lower left')

                    except:
                        continue
                if ax == axes[0, 1]:
                    try:
                        g = []
                        for i, model in enumerate(models[7:]):
                            colors_plot = colours[7:]
                            g.append(ax.scatter([-99.], [-99.], c=colors_plot[i], s=5, label=model))
                        ax.legend(handles=g, loc='lower left')

                    except:
                        continue
            else:
                if ax == axes[0, 0]:
                    try:
                        g = []
                        for i, model in enumerate(models):
                            g.append(ax.scatter([-99.], [-99.], c=colours[i], s=5, label=model))
                        ax.legend(handles=g, loc='lower left')

                    except:
                        continue


        for i, model in enumerate(models):

            m = getattr(lf_parameters, model)()

            if hasattr(m, 'beta'):
                s = abs(np.array(m.redshifts) - z) < 0.5

                try:
                    x = 0.4 * (lum_to_M(10 ** log10L_bin_centres) - np.array(m.M_star)[s][0])
                    phistar = 10 ** np.array(m.phi_star)[s][0]
                    alpha = np.array(m.alpha)[s][0]
                    beta = np.array(m.beta)[s][0]

                    phi = phistar / (10 ** (x * (alpha + 1)) + 10 ** (x * (beta + 1)))

                    ax.scatter(log10L_bin_centres, np.log10(phi), s=5, c=colours[i])
                except:
                    continue


            else:

                s = abs(np.array(m.redshifts) - z) < 0.5

                try:
                    x = 10**(log10L_bin_centres - np.log10(M_to_lum(np.array(m.M_star)[s][0])))
                    phistar = 10**np.array(m.phi_star)[s][0]
                    alpha = np.array(m.alpha)[s][0]
                    ax.scatter(log10L_bin_centres, np.log10(phistar*(x)**(alpha + 1) * np.exp(-x)), s=5, c=colours[i])

                except:
                    continue


        ax.text(0.7, 0.9, r'$\rm z=[{0:.1f}, {1:.1f})$'.format(z - 0.5, z + 0.5), fontsize=8, transform=ax.transAxes)

        ylim = np.array([-8., -4.01])

        ax.set_ylim(-8., -1.5)
        ax.set_xlim([min(log10L_bin_centres), max(log10L_bin_centres)])


    fig.text(0.5, 0.05, r'$\rm\log_{10}[L_{FUV} \; / \; (erg \; s^{-1} \; Hz^{-1})]$', horizontalalignment='center', verticalalignment='center')
    fig.text(0.05, 0.5, r'$\rm\log_{10}[\phi \; / \; Mpc^{-3}]$', horizontalalignment='center', verticalalignment='center', rotation=90)



    return fig




def lf_multi(zs, log10L_bin_centres, models, binned_lf=False, legend=False):
    # THIS WILL BE LF PLOT FOR A SELECTION OF MODELS

    markers = ['D', 's', 'h', '^', 'o', '*', 'v']

    cmap = mpl.cm.tab20
    cmap_marker = mpl.cm.Dark2_r

    if len(zs)%2 == 1:
        N = int(len(zs)/2) + 1
    else:
        N = int(len(zs) / 2)

    fig, axes = plt.subplots(N, 2, figsize=(6, N*2), sharex=True, sharey=True)


    plt.subplots_adjust(left=0.10, bottom=0.10, top=0.95, right=0.85, wspace=0., hspace=-0.)

    if type(models) == str:
        models = [models]

    Nxx = len(models)

    for z, ax in zip(zs, axes.flatten()):

        if legend:
            if Nxx > 7:
                if ax == axes[0, 0]:
                    try:
                        g = []
                        for i, model in enumerate(models[:7]):
                            label = getattr(lf_parameters, model)().label
                            g.append(Line2D([-99., -98.], [-99., -98.], color=cmap(i/19), label=label, alpha=0.6))
                        ax.legend(handles=g, loc='lower left', fontsize=6)

                    except:
                        continue
                if ax == axes[0, 1]:
                    try:
                        g = []
                        for i, model in enumerate(models[7:]):
                            i = i+7
                            label = getattr(lf_parameters, model)().label
                            g.append(Line2D([-99., -98.], [-99., -98.], color=cmap(i/19), label=label, alpha=0.6))
                        ax.legend(handles=g, loc='lower left', fontsize=6)

                    except:
                        continue
            else:
                if ax == axes[0, 0]:
                    try:
                        g = []
                        for i, model in enumerate(models):
                            label = getattr(lf_parameters, model)().label
                            g.append(Line2D([-99., -98.], [-99., -98.], color=cmap(i/19), label=label, alpha=0.6))
                        ax.legend(handles=g, loc='lower left', fontsize=6)

                    except:
                        continue

        if ax == axes[1, 0]:

            if binned_lf:

                if type(binned_lf) == list:
                    g = []

                    for l, lf in enumerate(binned_lf):
                        g.append(plt.errorbar(-99., -99., yerr=1., color=cmap_marker(l/7), linestyle='', marker=markers[l], mec='k', alpha=0.6, label=lf['label'], markersize=3, elinewidth=1, capsize=1, capthick=1))
                    ax.legend(handles=g, loc='lower left', fontsize=6)

                else:
                    g = []

                    for l in range(1):
                        g.append(
                            plt.errorbar(-99., -99., yerr=1., fmt='kD', alpha=0.6, label=binned_lf['label'], markersize=3,
                                         elinewidth=1, capsize=1, capthick=1))
                    ax.legend(handles=g, loc='lower left', fontsize=6)


            # except:
            #    continue

        for i, model in enumerate(models):

            m = getattr(lf_parameters, model)()

            if hasattr(m, 'beta'):
                s = abs(np.array(m.redshifts) - z) < 0.5

                try:
                    x = 0.4 * (lum_to_M(10 ** log10L_bin_centres) - np.array(m.M_star)[s][0])
                    phistar = 10 ** np.array(m.phi_star)[s][0]
                    alpha = np.array(m.alpha)[s][0]
                    beta = np.array(m.beta)[s][0]

                    phi = phistar / (10 ** (x * (alpha + 1)) + 10 ** (x * (beta + 1)))

                    ax.plot(log10L_bin_centres, np.log10(phi), c=cmap(i/19), alpha=0.6)
                except:
                    continue


            else:

                s = abs(np.array(m.redshifts) - z) < 0.5

                try:
                    x = 10**(log10L_bin_centres - np.log10(M_to_lum(np.array(m.M_star)[s][0])))
                    phistar = 10**np.array(m.phi_star)[s][0]
                    alpha = np.array(m.alpha)[s][0]
                    ax.plot(log10L_bin_centres, np.log10(phistar*(x)**(alpha + 1) * np.exp(-x)), c=cmap(i/19), alpha=0.6)

                except:
                    continue

        if binned_lf:
            if type(binned_lf) == list:
                for l, lf in enumerate(binned_lf):
                    try:
                        q = abs(np.array([float(k) for k in [*lf['LF'].keys()]]) - z) < 0.5
                        z_key = str(int(np.array([float(k) for k in [*lf['LF'].keys()]])[q]))
                        log10L_bin = np.log10(M_to_lum(np.array(lf['LF'][str(int(z_key))]['M'])))
                        phi_bin = np.array(lf['LF'][str(int(z_key))]['phi'])
                        err = np.array(lf['LF'][str(int(z_key))]['phi_err'])
                        uplims = lf['LF'][str(int(z_key))]['uplim']

                        if lf['both_err']:
                            if lf['log_err']:
                                phi_err = err
                            else:
                                phi_err = [logerr_lo(phi_bin, err[0], uplims), logerr_hi(phi_bin, err[1])]
                        else:
                            if lf['log_err']:
                                phi_err = err
                            else:
                                phi_err = [logerr_lo(phi_bin, err, uplims), logerr_hi(phi_bin, err)]

                        ax.errorbar(log10L_bin, np.log10(phi_bin), yerr=phi_err, color=cmap_marker(l/7), linestyle='', marker=markers[l], mec='k', alpha=0.8, zorder=50+l,
                                    markersize=3, elinewidth=1, capsize=1, capthick=1, uplims=uplims)

                    except:
                        continue


            else:
                try:
                    q = abs(np.array([float(k) for k in [*binned_lf['LF'].keys()]]) - z) < 0.5
                    z_key = str(int(np.array([float(k) for k in [*binned_lf['LF'].keys()]])[q]))
                    log10L_bin = np.log10(M_to_lum(np.array(binned_lf['LF'][str(int(z_key))]['M'])))
                    phi_bin = np.array(binned_lf['LF'][str(int(z_key))]['phi'])
                    err = np.array(binned_lf['LF'][str(int(z_key))]['phi_err'])
                    uplims = binned_lf['LF'][str(int(z_key))]['uplim']
                    if binned_lf['both_err']:
                        if binned_lf['log_err']:
                            phi_err = err
                        else:
                            phi_err = [logerr_lo(phi_bin, err[0], uplims), logerr_hi(phi_bin, err[1])]
                            print(phi_err)
                    else:
                        if binned_lf['log_err']:
                            phi_err = err
                        else:
                            phi_err = [logerr_lo(phi_bin, err, uplims) , logerr_hi(phi_bin, err)]
                    ax.errorbar(log10L_bin, np.log10(phi_bin), yerr=phi_err, fmt='kD', alpha=0.6, zorder=50, markersize=3, elinewidth=1, capsize=1, capthick=1, uplims=uplims)

                except:
                    continue


        ax.text(0.7, 0.9, r'$\rm z=[{0:.1f}, {1:.1f})$'.format(z - 0.5, z + 0.5), fontsize=8, transform=ax.transAxes)

        ylim = np.array([-8., -4.01])

        ax.set_ylim(-8., -1.5)
        ax.set_xlim([min(log10L_bin_centres), max(log10L_bin_centres)])


    fig.text(0.5, 0.05, r'$\rm\log_{10}[L_{FUV} \; / \; (erg \; s^{-1} \; Hz^{-1})]$', horizontalalignment='center', verticalalignment='center')
    fig.text(0.05, 0.5, r'$\rm\log_{10}[\phi \; / \; Mpc^{-3}]$', horizontalalignment='center', verticalalignment='center', rotation=90)



    return fig
