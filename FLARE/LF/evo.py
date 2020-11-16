# ---
from scipy.stats import linregress
import numpy as np

import scipy.integrate as cp
import scipy.interpolate as cpi
import scipy.special as cps

from mpmath import gammainc

import matplotlib.pyplot as plt

import csv

from FLARE.photom import flux_to_L, lum_to_flux, M_to_lum, lum_to_M
import FLARE.core

import FLARE.LF.lf_parameters as lf_parameters
import FLARE.LF.binned_lf as binned_lf


#
# geo = (4. * np.pi * (100. * 10. * 3.0867 * 10 ** 16) ** 2)  # factor relating the L to M in cm^2
#
#
# def L(M):
#     return 10 ** (-0.4 * (M + 48.6)) * geo
#
#
# def M_to_log10L(M):
#     return -0.4 * (M + 48.6) + np.log10(geo)
#
#
# def M(log10L):
#     return -2.5 * (log10L - np.log10(geo)) - 48.6


def dVc(z, cosmo):
    return cosmo.differential_comoving_volume(z).value


def _integ(x, a):
    return x ** (a) * np.exp(-x)


def _integ2(x, a):
    return 0.4*np.log(10)*10**(-0.4*x* (a+1)) * np.exp(-10**(-0.4*x))


def _integ_dblpow(x, a, b):
    return 1 / (10 ** (x*(a+1)) + 10 ** (x*(b+1)))


def _integ_dblpow2(x, a, b):
    return x**a * (1+x)**(b-a)


def quadfunct(f, a, b, args):
    return cp.quad(f, a, b, args=args)[0]


def quadfunct2(f, x1, x2, a, b):
    args = (a, b)
    return cp.quad(f, x1, x2, args=args)[0]


def trapzfunct(f, bine1, bine2, alpha):
    x = np.array([bine1, bine2])
    y = f(x, alpha)
    return np.trapz(y, x)
# Main part of module:


class evo_base:

    def __init__(self, verbose = False):

        if verbose:
            print(self.model_style)


    def N(self, area=1., cosmo=False, redshift_limits=[8., 15.], log10L_limits=[27.5, 30.], dz=0.05, dlog10L=0.05, per_arcmin=False, return_volumes=False):

        # calculates the number of galaxies in each bin on a grid defined by redshift_limits, log10L_limits, dz, dlog10L
        # and area based on a luminosity function evolution model.

        area_sm = area  # Area in square arcmin
        area_sd = area_sm / 3600.  # Area in square degrees
        area_sr = (np.pi / 180.) ** 2 * area_sd  # Area in steradian

        if not cosmo: cosmo = FLARE.core.default_cosmo()

        # Setting the bin edges as well as centres for later operations
        bin_edges = {'log10L': np.arange(log10L_limits[0], log10L_limits[-1] + dlog10L, dlog10L),
                     'z': np.arange(redshift_limits[0], redshift_limits[-1] + dz, dz)}
        bin_centres = {'log10L': bin_edges['log10L'][:-1] + dlog10L / 2., 'z': bin_edges['z'][:-1] + dz / 2.}

        # Using astropy.cosmology to calculate the volume in each redshift bin
        volumes = np.asarray([cp.quad(dVc, bin_edges['z'][i - 1], bin_edges['z'][i], args=cosmo)[0] for i in
                              range(1, len(bin_edges['z']))])

        params = self.parameters(bin_centres['z'])

        if 'beta' in params.keys():
            alphas = params['alpha']
            Lstars = M_to_lum(params['M*'])
            Mstars = params['M*']
            phistars = 10 ** params['log10phi*']
            betas = params['beta']

            N = phistars[None, :] * np.vectorize(quadfunct2)(_integ_dblpow,
                                                             0.4 * (lum_to_M(
                                                                 10 ** bin_edges['log10L'][1:, None]) - Mstars[None,
                                                                                                        :]),
                                                             0.4 * (lum_to_M(
                                                                 10 ** bin_edges['log10L'][:-1, None]) - Mstars[None,
                                                                                                         :]),
                                                             alphas[None, :], betas[None, :]) * volumes[None,
                                                                                               :] * area_sr

            ''' Left in for possible future implementation
                        #N = phistars[None, :] * np.vectorize(quadfunct2)(_integ_dblpow2,
                        #                                                10**(bin_edges['log10L'][:-1, None] - Lstars[None, :]),
                        #                                                10**(bin_edges['log10L'][1:, None] - Lstars[None, :]),
                        #                                                alphas[None, :], betas[None, :] ) * volumes[None, :] * area_sr
            '''

        else:
            alphas = params['alpha']
            Lstars = M_to_lum(params['M*'])
            phistars = 10 ** params['log10phi*']
            Mstars = params['M*']
            N = phistars[None, :] * np.vectorize(quadfunct)(_integ,
                                                         10 ** bin_edges['log10L'][:-1, None] / Lstars[None, :],
                                                         10 ** bin_edges['log10L'][1:, None] / Lstars[None, :],

                                                         args=(alphas[None, :])) * volumes[None, :] * area_sr
            ''' Left in for possible future implementation
            N = phistars[None, :] * np.vectorize(quadfunct)(_integ2,
                                                         lum_to_M(10**bin_edges['log10L'][1:, None]) - Mstars[None, :],
                                                         lum_to_M(10**bin_edges['log10L'][:-1, None]) - Mstars[None, :],
                                                         args=(alphas[None, :])) * volumes[None, :] * area_sr
            '''

        if per_arcmin:
            N /= area_sm

        if return_volumes:
            return bin_edges, bin_centres, volumes, N

        else:
            return bin_edges, bin_centres, N


    def sample(self, area=1., cosmo=False, redshift_limits=[8., 15.], log10L_limits=[27., 30.], dz=0.05,
                     seed=False):

        # samples the LF evolution model in a given area

        if not cosmo: cosmo = FLARE.core.default_cosmo()

        area_sm = area  # Area in square arcmin
        area_sd = area_sm / 3600.  # Area in square degrees
        area_sr = (np.pi / 180.) ** 2 * area_sd  # Area in steradian

        # Setting the bin edges as well as centres for later operations
        bin_edges = {'z': np.arange(redshift_limits[0], redshift_limits[-1] + dz, dz)}
        bin_centres = {'z': bin_edges['z'][:-1] + dz / 2.}

        # Using astropy.cosmology to calculate the volume in each redshift bin
        volumes = np.asarray([cp.quad(dVc, bin_edges['z'][i - 1], bin_edges['z'][i], args=cosmo)[0] for i in
                              range(1, len(bin_edges['z']))])

        # Initialising the output array
        sample = {'log10L': np.asarray([]), 'z': np.asarray([])}

        if seed: np.random.seed(seed)

        for i in range(len(bin_centres['z'])):
            params = self.parameters(bin_centres['z'][i])

            sp = {}
            sp['alpha'] = params['alpha']
            sp['phi*'] = 10 ** params['log10phi*']
            sp['log10L*'] = np.log10(M_to_lum(params['M*']))

            LF = LF_interpolation(sp)

            samples = LF.sample(volumes[i] * area_sr, log10L_limits[0])
            sample['log10L'] = np.concatenate((sample['log10L'], samples))
            sample['z'] = np.concatenate(
                (sample['z'], np.array([bin_centres['z'][i] + dz * (np.random.random() - 0.5) for item in samples])))

        return sample


    def interpolate_LF(self, bin_centres, N, zs, log10Ls):

        # returns numbers of objects for input luminosities (log10Ls) and redshifts (zs)
        # takes bin centres and N-array from, e.g. FLARE.evo.linear.N()

        lf_interp = cpi.interp2d(bin_centres['z'], bin_centres['log10L'], N)

        # This returns a 2d surface at the moment. The idea was to return a single value for pairs of L, z
        # Can still be used for this if individual pairs are used, but a bit useless as is.

        return lf_interp(zs, log10Ls)


class linear(evo_base):

    def __init__(self, model):
        # lp is a dictionary of the parameters of the linear evolution model

        self.model = model

        self.lp = model.lp
        self.z_ref = model.z_ref

        self.model_style = 'Linear regression LF evolution method.'

        super().__init__()

    def parameters(self, z):
        # use linear evolution model
        # get parameters as a function of z
        # returns a dictionary of parameters
        p = {}
        for param in self.lp:
            p[param] = self.lp[param][0] * (z - self.z_ref) + self.lp[param][1]

        return p

    def parameters_line(self, zr = [6.,13.]):
        p = {}
        for param in self.lp:
            p[param] = [self.lp[param][0] * (zr[0] - self.z_ref) + self.lp[param][1], self.lp[param][0] * (zr[1] - self.z_ref) + self.lp[param][1]]

        return zr, p


class interp(evo_base):

    def __init__(self, model):
        # lp is a dictionary of the parameters of the linear evolution model

        self.model = model

        self.model_style = 'Interpolation LF evolution method.'

        super().__init__()

    def parameters(self, z=8.):
        # use interpolation model
        # get parameters as a function of z
        # returns a dictionary of parameters

        return self.model.interpolate_parameters(z)


class existing_model:

    def __init__(self, verbose = False):


        # convert M* to L* HERE

        self.lp, self.z_ref = self.calculate_linear_evolution_coeffs()
        if verbose:
            print(self.name, self.LF_model)

    def interpolate_parameters(self, z=8.):
        # interpolates parameters as a function of z
        # returns a dictionary of the Schechter function parameters for given redshift(s)

        z_mod = self.redshifts
        alpha_mod = self.alpha
        log10phi_mod = self.phi_star
        log10M_mod = self.M_star
        if self.LF_model == 'Double Power Law':
            beta_mod = self.beta
            p = {'alpha': np.interp(z, z_mod, alpha_mod), 'log10phi*': np.interp(z, z_mod, log10phi_mod),
                 'M*': np.interp(z, z_mod, log10M_mod), 'beta': np.interp(z, z_mod, beta_mod)}

        else:
            p = {'alpha': np.interp(z, z_mod, alpha_mod), 'log10phi*': np.interp(z, z_mod, log10phi_mod),
                 'M*': np.interp(z, z_mod, log10M_mod)}

        return p

    def calculate_linear_evolution_coeffs(self, zr = [6., 15.], z_ref = 6.):
        # Function that calculates the linear evolution coeffs
        # returns a dictionary of linear model coefficients and goodness of fit

        s = (np.array(self.redshifts)>=zr[0])&(np.array(self.redshifts)<=zr[1])

        z_mod = np.array(self.redshifts)[s] - z_ref
        alpha_mod = np.array(self.alpha)[s]
        log10phi_mod = np.array(self.phi_star)[s]
        M_mod = np.array(self.M_star)[s]

        # The output contains full linregress output (0th and 1st element contain the slope and intercept respectively)
        fit_alpha = linregress(z_mod, alpha_mod)
        fit_log10phi = linregress(z_mod, log10phi_mod)
        fit_M = linregress(z_mod, M_mod)

        if self.LF_model == 'Double Power Law':
            beta_mod = np.array(self.beta)[s]
            fit_beta = linregress(z_mod, beta_mod)
            lp = {'alpha': fit_alpha, 'beta': fit_beta, 'log10phi*': fit_log10phi, 'M*': fit_M}

        else:
            lp = {'alpha': fit_alpha, 'log10phi*': fit_log10phi, 'M*': fit_M}

        return lp, z_ref



class LF_interpolation:
    # --- LF interpolation functions for sampling and predictions

    def __init__(self, sp):
        self.sp = sp

    def CulmPhi(self, log10L):
        y = log10L - self.sp['log10L*']
        x = 10 ** y
        alpha = self.sp['alpha']

        gamma = cp.quad(_integ, x, np.inf, args=alpha)[0]
        num = gamma * self.sp['phi*']

        return num

    def CDF(self, log10L_limit, normed=True):
        log10Ls = np.arange(self.sp['log10L*'] + 5., log10L_limit - 0.01, -0.01)

        CDF = np.array([self.CulmPhi(log10L) for log10L in log10Ls])

        if normed: CDF /= CDF[-1]

        return log10Ls, CDF

    def N_exact(self, volume, bin_edges):
        # --- return the exact number of galaxies expected in each bin

        CulmN = np.array([self.CulmPhi(x) for x in bin_edges]) * volume

        return -(CulmN[1:] - CulmN[0:-1])

    def sample(self, volume, log10L_limit):
        L, CDF = self.CDF(log10L_limit, normed=False)

        n = np.random.poisson(volume * CDF[-1])

        nCDF = CDF / CDF[-1]

        log10L_sample = np.interp(np.random.random(n), nCDF, L)

        return log10L_sample

    def bin(self, log10L_sample, bins):
        # --- bins can either be the number of bins or the bin_edges

        N_sample, bin_edges = np.histogram(log10L_sample, bins=bins, normed=False)

        return N_sample


def get_model(model):
    return getattr(lf_parameters, model)()


def get_lf(binned_lf_name):
    return getattr(binned_lf, binned_lf_name)().lf


model_names = ['bluetides', 'Finkelstein_review', 'Finkelstein_obs', 'Bowler20152020', 'Bowler20152020_DPL', 'Bouwens2015', 'Ma2019', 'Mason15', 'Yung2018', 'FLARES', 'FLARES_DPL', 'TNG_A', 'TNG_B', 'TNG_C']

binned_lf_names = ['Atek18', 'Bouwens15', 'Bouwens16', 'Bouwens17', 'Bowler20', 'FLARES', 'Finkelstein', 'Finkelstein15', 'McLeod15', 'Oesch18', 'Stefanon19']


def print_model_parameters(model_list):
    headers = ['Redshift', 'log10phi*', 'log10L*', 'M*', 'alpha', 'beta']
    #headers = ['Model', 'Reference', 'type', 'LF model', 'Redshift', 'log10phi*', 'log10L*', 'alpha', 'beta']

    if type(model_list) == str:
        model_list = [model_list]

    row_format = '{:>12}' * (len(headers))

    for model in model_list:

        m = get_model(model)

        print(row_format.format(*headers))
        if m.LF_model == 'Double Power Law':
            for i, z in enumerate(m.redshifts):
                redshift = f'{m.redshifts[i]:.2f}'
                phistar = f'{m.phi_star[i]:.2f}'
                log10L = f'{np.log10(M_to_lum(m.M_star[i])):.2f}'
                alpha = f'{m.alpha[i]:.2f}'
                beta = f'{m.beta[i]:.2f}'
                Mstar = f'{m.M_star[i]:.2f}'
                print(row_format.format(redshift, phistar, log10L, Mstar, alpha, beta))

        if m.LF_model == 'Schechter':
            for i, z in enumerate(m.redshifts):
                redshift = f'{m.redshifts[i]:.2f}'
                phistar = f'{m.phi_star[i]:.2f}'
                log10L = f'{np.log10(M_to_lum(m.M_star[i])):.2f}'
                Mstar = f'{m.M_star[i]:.2f}'
                alpha = f'{m.alpha[i]:.2f}'
                print(row_format.format(redshift, phistar, log10L, Mstar, alpha, 'N/A'))

    return


def write_model_parameters_csv(model_list, filename='lf_parameters'):
    #headers = ['Redshift', 'log10phi*', 'log10L*', 'M*', 'alpha', 'beta']
    headers = ['Model', 'Reference', 'type', 'LF model', 'Redshift', 'log10phi*', 'log10L*', 'M*', 'alpha', 'beta', 'ADS', 'arXiv']

    if type(model_list) == str:
        model_list = [model_list]

    row_format = '{},' * (len(headers))

    with open(f'{filename}.csv', mode='w', newline='') as output_file:
        output_writer = csv.writer(output_file, delimiter=',')

        output_writer.writerow(headers)

        for model in model_list:

            m = get_model(model)

            #print(row_format.format(*headers))
            if m.LF_model == 'Double Power Law':
                for i, z in enumerate(m.redshifts):
                    redshift = f'{m.redshifts[i]:.4f}'
                    phistar = f'{m.phi_star[i]:.4f}'
                    log10L = f'{np.log10(M_to_lum(m.M_star[i])):.4f}'
                    alpha = f'{m.alpha[i]:.4f}'
                    beta = f'{m.beta[i]:.4f}'
                    Mstar = f'{m.M_star[i]:.4f}'

                    output_writer.writerow([m.name, m.ref, m.type, m.LF_model, redshift, phistar, log10L, Mstar, alpha, beta, m.ads, m.arxiv])


            if m.LF_model == 'Schechter':
                for i, z in enumerate(m.redshifts):
                    redshift = f'{m.redshifts[i]:.4f}'
                    phistar = f'{m.phi_star[i]:.4f}'
                    log10L = f'{np.log10(M_to_lum(m.M_star[i])):.4f}'
                    Mstar = f'{m.M_star[i]:.4f}'
                    alpha = f'{m.alpha[i]:.4f}'
                    output_writer.writerow([m.name, m.ref, m.type, m.LF_model, redshift, phistar, log10L, Mstar, alpha, 'N/A', m.ads, m.arxiv])

    return
