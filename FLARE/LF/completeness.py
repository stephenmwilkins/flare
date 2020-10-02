import FLARE

import numpy as np

from FLARE.photom import lum_to_flux, flux_to_L

import scipy.special as cps


def completeness_cut(bin_centres, flux_limit, stretch=1.0, cosmo=False):
    # calculates a grid for completeness.
    # Should be done to the same size as N, take bin_centres from N.
    # takes flux_limit and background (in nJy), and snr_target
    # stretch is variable.

    dz = bin_centres['z'][1] - bin_centres['z'][0]
    dlog10L = bin_centres['log10L'][1] - bin_centres['log10L'][0]

    if not cosmo: cosmo = FLARE.core.default_cosmo()

    c = np.zeros((len(bin_centres['z']), len(bin_centres['log10L'])))

    for i, z in enumerate(bin_centres['z']):
        for j, log10L in enumerate(bin_centres['log10L']):
            f = lum_to_flux(10 ** log10L, cosmo, z)
            c[i, j] = 0.5 * (1.0 + cps.erf(stretch * (f - flux_limit)))

    return c.T


def completeness_sample(bin_centres, flux_limit, stretch=1.0, cosmo=False, N_samples=100):
    # calculates a grid for completeness.
    # Should be done to the same size as N, take bin_centres from N.
    # takes flux_limit and background (in nJy), and snr_target
    # stretch is variable.

    dz = bin_centres['z'][1] - bin_centres['z'][0]
    dlog10L = bin_centres['log10L'][1] - bin_centres['log10L'][0]

    if not cosmo: cosmo = FLARE.core.default_cosmo()

    c = np.zeros((len(bin_centres['z']), len(bin_centres['log10L'])))

    for i, z in enumerate(bin_centres['z']):
        for j, log10L in enumerate(bin_centres['log10L']):
            z_sample = np.random.uniform(z - dz / 2, z + dz / 2, N_samples)
            log10L_sample = np.random.uniform(log10L - dlog10L / 2, log10L + dlog10L / 2, N_samples)
            f = lum_to_flux(10 ** log10L_sample, cosmo, z_sample)

            N_fract = 0.

            for q in range(N_samples):
                N_fract += 0.5 * (1.0 + cps.erf(stretch * (f[q] - flux_limit)))

            c[i, j] = N_fract / N_samples

    return c.T


def completeness_erf(bin_centres, flux_limit, stretch=1.0, cosmo=False, N_samples=100):
    # calculates a grid for completeness.
    # Should be done to the same size as N, take bin_centres from N.
    # takes flux_limit and background (in nJy), and snr_target
    # stretch is variable.

    dz = bin_centres['z'][1] - bin_centres['z'][0]
    dlog10L = bin_centres['log10L'][1] - bin_centres['log10L'][0]

    if not cosmo: cosmo = FLARE.core.default_cosmo()

    c = np.zeros((len(bin_centres['z']), len(bin_centres['log10L'])))

    for i, z in enumerate(bin_centres['z']):
        for j, log10L in enumerate(bin_centres['log10L']):

            f = np.linspace(lum_to_flux(10 ** (log10L - dlog10L / 2), cosmo, z),
                                lum_to_flux(10 ** (log10L + dlog10L / 2), cosmo, z), N_samples)
            N_fract = 0.

            for fluxx in f:
                N_fract += 0.5 * (1.0 + cps.erf(stretch * (fluxx - flux_limit)))

            c[i, j] = N_fract / N_samples

    return c.T


def completeness_erf_legacy(bin_centres, flux_limit, stretch=1.0, cosmo=False, N_samples=100):
    # calculates a grid for completeness. Due to a design fla.....feature, this only calculates completeness in the bin
    # that gets dissected by the flux cut. Will most likely be removed in a future revision.
    # Should be done to the same size as N, take bin_centres from N.
    # takes flux_limit and background (in nJy), and snr_target
    # stretch is variable.

    dz = bin_centres['z'][1] - bin_centres['z'][0]
    dlog10L = bin_centres['log10L'][1] - bin_centres['log10L'][0]

    if not cosmo: cosmo = FLARE.core.default_cosmo()

    c = np.zeros((len(bin_centres['z']), len(bin_centres['log10L'])))

    for i, z in enumerate(bin_centres['z']):
        for j, log10L in enumerate(bin_centres['log10L']):

            if abs(log10L - np.log10(flux_to_L(flux_limit, cosmo, z))) <= dlog10L:
                # print(abs(log10L-np.log10(flux_to_L(flux_limit, cosmo, z))))
                f = np.linspace(lum_to_flux(10 ** (log10L - dlog10L / 2), cosmo, z),
                                lum_to_flux(10 ** (log10L + dlog10L / 2), cosmo, z), N_samples)
                N_fract = 0.

                for fluxx in f:
                    N_fract += 0.5 * (1.0 + cps.erf(stretch * (fluxx - flux_limit)))
                c[i, j] = N_fract / N_samples

    return c.T

