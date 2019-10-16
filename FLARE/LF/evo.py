# ---
from scipy.stats import linregress
import numpy as np

import scipy.integrate as cp
import scipy.interpolate as cpi
import scipy.special as cps

from mpmath import gammainc

import matplotlib.pyplot as plt

from FLARE.photom import flux_to_L, lum_to_flux, M_to_lum
import FLARE.core


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
    return x ** (a - 1) * np.exp(-x)


def quadfunct(f, a, b, args):
    return cp.quad(f, a, b, args=args)[0]

def trapzfunct(f, bine1, bine2, alpha):
    x = np.array([bine1, bine2])
    y = f(x, alpha)
    return np.trapz(y, x)
# Main part of module:


class evo_base:

    def __init__(self):

        print(self.model_style)


    def N(self, area=1., cosmo=False, redshift_limits=[8., 15.], log10L_limits=[27.5, 30.], dz=0.05, dlog10L=0.05, per_arcmin=False):

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
        alphas = params['alpha']
        Lstars = M_to_lum(params['M*'])
        phistars = 10**params['log10phi*']

        N = phistars[None, :] * np.vectorize(quadfunct)(_integ,
                                                         10 ** bin_edges['log10L'][:-1, None] / Lstars[None, :],
                                                         10 ** bin_edges['log10L'][1:, None] / Lstars[None, :],
                                                         args=(alphas[None, :])+1) * volumes[None, :] * area_sr

        if per_arcmin:
            N /= area_sm

        return bin_edges, bin_centres, N


    def N_trapz(self, area=1., cosmo=False, redshift_limits=[8., 15.], log10L_limits=[27.5, 30.], dz=0.05, dlog10L=0.05, per_arcmin=False):

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
        alphas = params['alpha']
        Lstars = M_to_lum(params['M*'])
        phistars = 10**params['log10phi*']


        N = phistars[None, :] * np.vectorize(trapzfunct)(_integ, 10**bin_edges['log10L'][:-1, None]/Lstars[None, :],
                                                        10 ** bin_edges['log10L'][1:, None] / Lstars[None, :],
                                                        alphas[None, :]+1) * volumes[None, :] * area_sr


        if per_arcmin:
            N /= area_sm

        return bin_edges, bin_centres, N


    def N_old(self, area=1., cosmo=False, redshift_limits=[8., 15.], log10L_limits=[27.5, 30.], dz=0.05, dlog10L=0.05,
          flux_min=False):

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

        # Initialising the output array
        N = np.zeros((len(bin_centres['log10L']), len(bin_centres['z'])))

        # Loop calculates LF for each input z (bin centres) and returns the exact numbers expected in each bin
        # (There may be a better option for generating this)

        for i in range(len(bin_centres['z'])):
            params = self.parameters(bin_centres['z'][i])

            sp = {}
            sp['alpha'] = params['alpha']
            sp['phi*'] = 10 ** params['log10phi*']
            sp['log10L*'] = np.log10(M_to_lum(params['M*']))

            LF = LF_interpolation(sp)

            if not flux_min:

                N_ext = LF.N_exact(volumes[i] * area_sr, bin_edges['log10L'])

                for j in range(len(N_ext)):
                    N[j, i] = N_ext[j] / area_sm
            else:
                log10L_min = np.log10(flux_to_L(flux_min, cosmo, bin_centres['z'][i]))
                log10L_edges = bin_edges['log10L'][np.where(bin_edges['log10L'] >= log10L_min - dlog10L)]
                N_ext = LF.N_exact(volumes[i] * area_sr, log10L_edges)

                for j in range(1, len(N_ext) + 1):
                    N[-j, i] = N_ext[len(N_ext) - j] / area_sm

        return bin_edges, bin_centres, N


    def N_mpmath(self, area=1., cosmo=False, redshift_limits=[8., 15.], log10L_limits=[27.5, 30.], dz=0.05, dlog10L=0.05):

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
        alphas = params['alpha']
        Lstars = M_to_lum(params['M*'])
        phistars = 10**params['log10phi*']

        N = phistars[None, :]*np.vectorize(gammainc)(alphas[None, :]+1., 10**bin_edges['log10L'][:-1, None]/Lstars[None, :], 10**bin_edges['log10L'][1:, None]/Lstars[None, :])*volumes[None, :]*area_sr

        N /= area_sm

        return bin_edges, bin_centres, N


    def sample(self, area=1., cosmo=False, redshift_limits=[8., 15.], dlog10L=0.05, dz=0.05, flux_min=1E-8, seed=False):

        # samples the LF evolution model in a given area

        if not cosmo: cosmo = FLARE.core.default_cosmo()

        area_sm = area  # Area in square arcmin
        area_sd = area_sm / 3600.  # Area in square degrees
        area_sr = (np.pi / 180.) ** 2 * area_sd  # Area in steradian

        # Setting the bin edges as well as centres for later operations
        bin_edges = {'z': np.arange(redshift_limits[0], redshift_limits[-1] + dz, dz)}
        bin_centres = {'log10L': bin_edges['log10L'][:-1] + dlog10L / 2., 'z': bin_edges['z'][:-1] + dz / 2.}

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

            samples = LF.sample(volumes[i] * area_sr, np.log10(flux_to_L(flux_min, cosmo, bin_centres['z'][i])))
            sample['log10L'] = np.concatenate((sample['log10L'], samples))
            sample['z'] = np.concatenate((sample['z'], np.asarray([bin_centres['z'][i]] * len(samples))))

        return sample

    def sample_input(self, area=1., cosmo=False, redshift_limits=[8., 15.], log10L_limits=[27., 30.], dz=0.05,
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

    def bin_sample(self, area=1., cosmo=False, redshift_limits=[8., 15.], log10L_limits=[27., 32.], dz=0.05,
                   dlog10L=0.05, flux_min=1E-8, seed=False):

        # bin the sampled LF
        if not cosmo: cosmo = FLARE.core.default_cosmo()
        # samples the LF evolution model in a given area and bins it

        area_sm = area  # Area in square arcmin
        area_sd = area_sm / 3600.  # Area in square degrees
        area_sr = (np.pi / 180.) ** 2 * area_sd  # Area in steradian

        # Setting the bin edges as well as centres for later operations
        bin_edges = {'log10L': np.arange(log10L_limits[0], log10L_limits[-1] + dlog10L, dlog10L),
                     'z': np.arange(redshift_limits[0], redshift_limits[-1] + dz, dz)}
        bin_centres = {'log10L': bin_edges['log10L'][:-1] + dlog10L / 2., 'z': bin_edges['z'][:-1] + dz / 2.}

        # Using astropy.cosmology to calculate the volume in each redshift bin
        volumes = np.asarray([cp.quad(dVc, bin_edges['z'][i - 1], bin_edges['z'][i], args=cosmo)[0] for i in
                              range(1, len(bin_edges['z']))])

        # Initialising the output array
        N_sample = np.zeros((len(bin_centres['log10L']), len(bin_centres['z'])))

        if seed: np.random.seed(seed)

        # Loop calculates LF for each input z (bin centres) and returns the exact numbers expected in each bin
        # (There may be a better option for generating this)
        for i in range(len(bin_centres['z'])):
            params = self.parameters(bin_centres['z'][i])

            sp = {}
            sp['alpha'] = params['alpha']
            sp['phi*'] = 10 ** params['log10phi*']
            sp['log10L*'] = np.log10(M_to_lum(params['M*']))

            LF = LF_interpolation(sp)

            sample = LF.sample(volumes[i] * area_sr, np.log10(flux_to_L(flux_min, cosmo, bin_centres['z'][i])))
            N_binned = LF.bin(sample, bin_edges['log10L'])

            for j in range(len(N_binned)):
                N_sample[j, i] = N_binned[j]

        return bin_edges, bin_centres, N_sample

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

        self.model_style = 'Linear regression LF evolution method.'

        super().__init__()

    def parameters(self, z):
        # use linear evolution model
        # get parameters as a function of z
        # returns a dictionary of parameters
        p = {}
        for param in self.lp:
            p[param] = self.lp[param][0] * z + self.lp[param][1]

        return p

    def parameters_line(self, zr = [6.,13.]):
        # use linear evolution model
        # get parameters as a function of z
        # returns a dictionary of parameters
        p = {}
        for param in self.lp:
            p[param] = [self.lp[param][0] * zr[0] + self.lp[param][1], self.lp[param][0] * zr[1] + self.lp[param][1]]

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

    def __init__(self):
        self.lp = self.calculate_linear_evolution_coeffs()
        print(self.name)

    def interpolate_parameters(self, z=8.):
        # interpolates parameters as a function of z
        # returns a dictionary of the Schechter function parameters for given redshift(s)

        z_mod = self.redshifts
        alpha_mod = self.alpha
        log10phi_mod = self.phi_star
        log10M_mod = self.M_star
        p = {'alpha': np.interp(z, z_mod, alpha_mod), 'log10phi*': np.interp(z, z_mod, log10phi_mod),
             'M*': np.interp(z, z_mod, log10M_mod)}

        return p

    def calculate_linear_evolution_coeffs(self, zr = [7., 15.]):
        # Function that calculates the linear evolution coeffs
        # returns a dictionary of linear model coefficients and goodness of fit

        s = (np.array(self.redshifts)>=zr[0])&(np.array(self.redshifts)<=zr[1])

        z_mod = np.array(self.redshifts)[s]
        alpha_mod = np.array(self.alpha)[s]
        log10phi_mod = np.array(self.phi_star)[s]
        M_mod = np.array(self.M_star)[s]

        # The output contains full linregress output (0th and 1st element contain the slope and intercept respectively)
        fit_alpha = linregress(z_mod, alpha_mod)
        fit_log10phi = linregress(z_mod, log10phi_mod)
        fit_M = linregress(z_mod, M_mod)

        self.lp = {'alpha': fit_alpha, 'log10phi*': fit_log10phi, 'M*': fit_M}

        return self.lp


class bluetides(existing_model):  # --- based on bluetides simulation

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'Bluetides (Wilkins+2017)'
        self.ref = 'Wilkins+2017'
        self.type = 'hydro' 
        self.redshifts = [8.0, 9.0, 10.0, 11.0, 12.0, 13.0]  # array of redshifts
        self.phi_star = [-3.92, -4.2, -4.7, -4.79, -5.09, -5.71]  # array of log10(phi_star) values
        self.M_star = [-20.93, -20.68, -20.69, -20.17, -19.92, -19.91]  # array of M_star values
        self.alpha = [-2.04, -2.1, -2.27, -2.27, -2.35, -2.54]  # array of alpha values

        super().__init__()


class Finkelstein_review(existing_model):
    # --- LF evolution model based on Finkelstein review (2016)

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'Observational review (Finkelstein2016)'
        self.ref = 'Finkelstein2016'
        self.type = 'empirical'
        # self.redshifts = [4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]  # array of redshifts
        # self.phi_star = [-2.99, -3.18, -3.37, -3.56, -3.75, -3.94, -4.13]  # array of log10(phi_star) values
        # self.M_star = [-21.05, -20.92, -20.79, -20.66, -20.52, -20.39, -20.25]  # array of M_star values
        # self.alpha = [-1.69, -1.80, -1.91, -2.02, -2.13, -2.24, -2.35]  # array of alpha values

        # --- cut down to just include observed redshifts
        self.redshifts = [4.0, 5.0, 6.0, 7.0, 8.0]  # array of redshifts
        self.phi_star = [-2.99, -3.18, -3.37, -3.56, -3.75]  # array of log10(phi_star) values
        self.M_star = [-21.05, -20.92, -20.79, -20.66, -20.52]  # array of M_star values
        self.alpha = [-1.69, -1.80, -1.91, -2.02, -2.13]  # array of alpha values

        super().__init__()


class Finkelstein_obs(existing_model):
    # --- LF evolution model based on Finkelstein+ (2015)

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'Observational (Finkelstein+2015)'
        self.ref = 'Finkelstein+2015'
        self.type = 'observed'
        self.redshifts = [4.0, 5.0, 6.0, 7.0, 8.0]  # array of redshifts
        self.phi_star = [np.log10(14.1 * 10 ** -4), np.log10(8.95 * 10 ** -4), np.log10(1.86 * 10 ** -4),
                         np.log10(1.57 * 10 ** -4), np.log10(0.72 * 10 ** -4)]  # array of log10(phi_star) values
        self.M_star = [-20.73, -20.81, -21.13, -21.03, -20.89]  # array of M_star values
        self.alpha = [-1.56, -1.67, -2.02, -2.03, -2.36]  # array of alpha values

        super().__init__()


class Ma2019(existing_model):
    # --- LF evolution model based on Ma et al. (2019) (f_dust = 0.8)

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'FIRE-2 (Ma+2019)'
        self.ref = 'Ma+2019'
        self.type = 'hydro'
        self.redshifts = [5., 6., 7., 8.0, 9.0, 10.0]  # array of redshifts
        self.phi_star = [-3.55, -3.44, -4.09, -3.98, -4.57, -4.74]  # array of log10(phi_star) values
        self.M_star = [-21.77, -21.34, -21.73, -20.97, -21.30, -20.90]  # array of M_star values
        self.alpha = [-1.9, -1.87, -2.05, -2.08, -2.20, -2.31]  # array of alpha values

        super().__init__()


class Mason15(existing_model):
    # --- LF evolution model based on Mason et al. (2015)

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'Semi-empirical (Mason+2015)'
        self.ref = 'Mason+2015'
        self.type = 'semi-empirical'
        self.redshifts = [0., 2., 4., 5., 6., 7., 8., 9., 10., 12., 14., 16.]    # array of redshifts
        self.phi_star = [-2.97, -2.52, -2.93, -3.12, -3.19, -3.48, -4.03, -4.50, -5.12, -5.94, -7.05, -8.25]    # array of phi_star value to interpolate
        self.M_star = [-19.9, -20.3, -21.2, -21.2, -20.9, -21.0, -21.3, -21.2, -21.1, -21.0, -20.9, -20.7]  # array of M_star values
        self.alpha = [-1.68, -1.46, -1.64, -1.75, -1.83, -1.95, -2.10, -2.26, -2.47, -2.74, -3.11, -3.51] # array of alpha values

        super().__init__()


class Yung2018(existing_model):
    # --- LF evolution model based on Yung et al. (2018)

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'Semi-analytical (Yung+2018)'
        self.ref = 'Yung+2018'
        self.type = 'semi-analytical'
        self.redshifts = [4., 5., 6., 7., 8., 9., 10.]    # array of redshifts
        self.phi_star = [np.log10(3.151 * 10 ** -3), np.log10(2.075 * 10 ** -3), np.log10(1.352 * 10 ** -3),
         np.log10(0.818 * 10 ** -3), np.log10(0.306 * 10 ** -3), np.log10(0.133 * 10 ** -3), np.log10(0.053 * 10 ** -3)]  # array of log10(phi_star) values
        self.M_star = [-20.717, -20.774, -20.702, -20.609, -20.660, -20.584, -20.373]  # array of M_star values
        self.alpha = [-1.525, -1.602, -1.672, -1.715, -1.825, -1.879, -1.967] # array of alpha values

        super().__init__()


class LF_interpolation:
    # --- LF interpolation functions for sampling and predictions

    def __init__(self, sp):
        self.sp = sp

    def CulmPhi(self, log10L):
        y = log10L - self.sp['log10L*']
        x = 10 ** y
        alpha = self.sp['alpha']

        gamma = cp.quad(_integ, x, np.inf, args=alpha + 1)[0]
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
    bar.set_label(r'$\rm log_{10}(N/ [arcmin^{-2}])$', rotation=90)

    plt.ylabel(r"$\rm log_{10}(L_{\nu}/ [erg\, s^{-1}\, Hz^{-1}])$")
    plt.xlabel(r"$\rm z$")
    plt.xlim(min(bin_edges['z']), max(bin_edges['z']))
    plt.ylim(min(bin_edges['log10L']), max(bin_edges['log10L']))

    if save_file == False:
        return fig
    else:
        plt.savefig(save_file + '.png', dpi=300)


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
