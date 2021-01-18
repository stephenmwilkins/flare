from scipy.stats import linregress
import numpy as np



class existing_model:

    def __init__(self, verbose = False):

        # convert M* to L* HERE

        self.lp, self.z_ref = self.calculate_linear_evolution_coeffs()
        if verbose:
            print(self.name, self.type, self.LF_model)

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

    def calculate_linear_evolution_coeffs(self, zr=[6., 15.], z_ref=6.):
        # Function that calculates the linear evolution coeffs
        # returns a dictionary of linear model coefficients and goodness of fit

        s = (np.array(self.redshifts) >= zr[0]) & (np.array(self.redshifts) <= zr[1])

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


# Different LF parameters

class bluetides(existing_model):  # --- based on bluetides simulation

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'Bluetides (Wilkins+2017)'
        self.ref = 'Wilkins+2017'
        self.type = 'hydro'
        self.LF_model = 'Schechter'
        self.label = 'BlueTides'
        self.label2 = 'Hydro\ (Bluetides)'

        self.ads = 'https://ui.adsabs.harvard.edu/abs/2017MNRAS.469.2517W/abstract'
        self.arxiv = 'https://arxiv.org/abs/1704.00954'

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

        self.name = 'Observational review (Finkelstein 2016)'
        self.ref = 'Finkelstein+2016'
        self.type = 'empirical'
        self.LF_model = 'Schechter'
        self.label = 'Empirical (F15)'
        self.label2 = 'Empirical\ (F15)'

        self.ads = 'https://ui.adsabs.harvard.edu/abs/2016PASA...33...37F/abstract'
        self.arxiv = 'https://arxiv.org/abs/1511.05558'

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
        self.LF_model = 'Schechter'
        self.label = 'Finkelstein+2015'

        self.ads = 'https://ui.adsabs.harvard.edu/abs/2015ApJ...810...71F/abstract'
        self.arxiv = 'https://arxiv.org/abs/1410.5439'

        self.redshifts = [4.0, 5.0, 6.0, 7.0, 8.0]  # array of redshifts
        self.phi_star = [np.log10(14.1 * 10 ** -4), np.log10(8.95 * 10 ** -4), np.log10(1.86 * 10 ** -4),
                         np.log10(1.57 * 10 ** -4), np.log10(0.72 * 10 ** -4)]  # array of log10(phi_star) values
        self.M_star = [-20.73, -20.81, -21.13, -21.03, -20.89]  # array of M_star values
        self.alpha = [-1.56, -1.67, -2.02, -2.03, -2.36]  # array of alpha values

        super().__init__()


class Bowler20152020(existing_model):
    # --- LF evolution model based on Bowler et al. (2015, 2020)

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'Observational (Bowler+2015, 2020)'
        self.ref = 'Bowler+2015,2020'
        self.type = 'observed'
        self.LF_model = 'Schechter'
        self.label = 'Bowler+2015, 2020'

        self.ads = 'https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.2059B/abstract'
        self.arxiv = 'https://arxiv.org/abs/1911.12832'

        self.redshifts = [5, 6, 7, 8, 9]  # array of redshifts
        self.phi_star = np.log10(
            np.array([6.4e-4, 5.7e-4, 3.7e-4, 1.92e-4, 0.53e-4]))  # array of log10(phi_star) values
        self.M_star = [-21.07, -20.77, -20.56, -20.48, -20.80]  # array of M_star values
        self.alpha = [-1.81, -1.88, -2.09, -2.18, -2.31]  # array of alpha values

        super().__init__()


class Bowler20152020_DPL(existing_model):
    # --- LF evolution model based on Bowler et al. (2015, 2020)

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'Observational (Bowler+2015, 2020)'
        self.ref = 'Bowler+2015,2020'
        self.type = 'observed'
        self.LF_model = 'Double Power Law'
        self.label = 'Bowler+2015, 2020 DPL'

        self.ads = 'https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.2059B/abstract'
        self.arxiv = 'https://arxiv.org/abs/1911.12832'

        self.redshifts = [5, 6, 7, 8, 9]  # array of redshifts
        self.phi_star = np.log10(
            np.array([2.5e-4, 1.9e-4, 2.2e-4, 4.83e-4, 2.85e-4]))  # array of log10(phi_star) values
        self.M_star = [-21.40, -21.20, -20.61, -19.80, -19.67]  # array of M_star values
        self.alpha = [-2.00, -2.10, -2.19, -1.96, -2.10]  # array of alpha values
        self.beta = np.array([-4.8, -5.1, -4.6, -3.98, -3.75])

        super().__init__()


class Bouwens2015(existing_model):
    # --- LF evolution model based on Bouwens et al. (2015)

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'Observational (Bouwens+2015)'
        self.ref = 'Bouwens+2015'
        self.type = 'observed'
        self.LF_model = 'Schechter'
        self.label = 'Bouwens+2015'

        self.ads = 'https://ui.adsabs.harvard.edu/abs/2015ApJ...803...34B/abstract'
        self.arxiv = 'https://arxiv.org/abs/1403.4295'

        self.redshifts = [3.8, 4.9, 5.9, 6.8, 7.9, 10.4]  # array of redshifts
        self.phi_star = np.log10(
            np.array([1.97e-3, 0.74e-3, 0.50e-3, 0.29e-3, 0.21e-3, 0.008e-3]))  # array of log10(phi_star) values
        self.M_star = [-20.88, -21.17, -20.94, -20.87, -20.63, -20.92]  # array of M_star values
        self.alpha = [-1.64, -1.76, -1.87, -2.06, -2.02, -2.27]  # array of alpha values

        super().__init__()


class Ma2019(existing_model):
    # --- LF evolution model based on Ma et al. (2019) (f_dust = 0.8)

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'FIRE-2 (Ma+2019)'
        self.ref = 'Ma+2019'
        self.type = 'hydro'
        self.LF_model = 'Schechter'
        self.label = 'Ma+2019'
        self.label2 = 'Hydro\ (FIRE)'

        self.ads = 'https://ui.adsabs.harvard.edu/abs/2019MNRAS.487.1844M/abstract'
        self.arxiv = 'https://arxiv.org/abs/1902.10152'

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

        self.name = 'analytical (Mason+2015)'
        self.ref = 'Mason+2015'
        self.type = 'analytical'
        self.LF_model = 'Schechter'
        self.label = 'Mason+2015'
        self.label2 = 'Analytical\ (M15)'

        self.ads = 'https://ui.adsabs.harvard.edu/abs/2015ApJ...813...21M/abstract'
        self.arxiv = 'https://arxiv.org/abs/1508.01204'

        self.redshifts = [0., 2., 4., 5., 6., 7., 8., 9., 10., 12., 14., 16.]  # array of redshifts
        self.phi_star = [-2.97, -2.52, -2.93, -3.12, -3.19, -3.48, -4.03, -4.50, -5.12, -5.94, -7.05,
                         -8.25]  # array of phi_star value to interpolate
        self.M_star = [-19.9, -20.3, -21.2, -21.2, -20.9, -21.0, -21.3, -21.2, -21.1, -21.0, -20.9,
                       -20.7]  # array of M_star values
        self.alpha = [-1.68, -1.46, -1.64, -1.75, -1.83, -1.95, -2.10, -2.26, -2.47, -2.74, -3.11,
                      -3.51]  # array of alpha values

        super().__init__()


class Yung2018(existing_model):
    # --- LF evolution model based on Yung et al. (2018)

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'Semi-analytical (Yung+2018)'
        self.ref = 'Yung+2018'
        self.type = 'semi-analytical'
        self.LF_model = 'Schechter'
        self.label = 'Yung+2018'
        self.label2 = 'SAM\ (Y18)'

        self.ads = 'https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.2983Y/abstract'
        self.arxiv = 'https://arxiv.org/abs/1803.09761'

        self.redshifts = [4., 5., 6., 7., 8., 9., 10.]  # array of redshifts
        self.phi_star = [np.log10(3.151 * 10 ** -3), np.log10(2.075 * 10 ** -3), np.log10(1.352 * 10 ** -3),
                         np.log10(0.818 * 10 ** -3), np.log10(0.306 * 10 ** -3), np.log10(0.133 * 10 ** -3),
                         np.log10(0.053 * 10 ** -3)]  # array of log10(phi_star) values
        self.M_star = [-20.717, -20.774, -20.702, -20.609, -20.660, -20.584, -20.373]  # array of M_star values
        self.alpha = [-1.525, -1.602, -1.672, -1.715, -1.825, -1.879, -1.967]  # array of alpha values

        super().__init__()


class FLARES(existing_model):
    # --- LF evolution model based on Vijayan et al. (2020)

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'FLARES (Vijayan+2020)'
        self.ref = 'Vijayan+2020'
        self.type = 'hydro'
        self.LF_model = 'Schechter'
        self.label = 'FLARES'
        self.label2 = 'Hydro\ (FLARES)'

        self.ads = 'https://ui.adsabs.harvard.edu/abs/2020arXiv200806057V/abstract'
        self.arxiv = 'https://arxiv.org/abs/2008.06057'

        self.redshifts = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        self.phi_star = [-3.6305001277308033, -3.875836361821652, -4.892210648014174, -5.254864977740009,
                         -4.95384659793114,
                         -4.877348845014224]
        self.M_star = [-21.79285493467043, -21.59150068090623, -22.255242417520705, -22.027768435403182,
                       -21.34875019363828,
                       -20.55601510467414]
        self.alpha = [-1.9772870410668615, -2.1213054661801465, -2.492156665448407, -2.7226300850912906,
                      -2.722532096381987,
                      -3.148004062716375]

        super().__init__()


class FLARES_old(existing_model):
    # --- LF evolution model based on Vijayan et al. (2020)

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'FLARES (Vijayan+2020)'
        self.ref = 'Vijayan+2020'
        self.type = 'hydro'
        self.LF_model = 'Schechter'
        self.label = 'FLARES'

        self.ads = 'https://ui.adsabs.harvard.edu/abs/2020arXiv200806057V/abstract'
        self.arxiv = 'https://arxiv.org/abs/2008.06057'

        self.redshifts = [5., 6., 7., 8., 9., 10.]  # array of redshifts
        self.phi_star = [-3.674, -3.869, -4.353, -4.379, -4.299, -4.416]  # array of log10(phi_star) values
        self.M_star = [-21.812, -21.484, -21.465, -20.946, -20.458, -20.084]  # array of M_star values
        self.alpha = [-1.987, -2.141, -2.421, -2.584, -2.671, -3.053]  # array of alpha values

        super().__init__()


class FLARES_DPL(existing_model):
    # --- LF evolution model based on Vijayan et al. (2020)
    # --- Double Power Law

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'FLARES DPL (Vijayan+2020)'
        self.ref = 'Vijayan+2020'
        self.type = 'hydro'
        self.LF_model = 'Double Power Law'
        self.label = 'FLARES DPL'

        self.ads = 'https://ui.adsabs.harvard.edu/abs/2020arXiv200806057V/abstract'
        self.arxiv = 'https://arxiv.org/abs/2008.06057'

        self.redshifts = [5., 6., 7., 8., 9., 10.]  # array of redshifts
        self.phi_star = [-3.767, -4.171, -4.957, -5.262, -3.158, -4.027]  # array of log10(phi_star) values
        self.M_star = [-21.704, -21.766, -22.124, -21.823, -19.026, -19.612]  # array of M_star values
        self.alpha = [-2.032, -2.204, -2.529, -2.766, -2.338, -3.082]  # array of alpha values
        self.beta = [-4.437, -5.184, -5.237, -5.173, -3.669, -4.006]  # array of alpha values

        super().__init__()


class FLARES_DPL_OLD(existing_model):
    # --- LF evolution model based on Vijayan et al. (2020)
    # --- Double Power Law

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'FLARES DPL (Vijayan+2020)'
        self.ref = 'Vijayan+2020'
        self.type = 'hydro'
        self.LF_model = 'Double Power Law'
        self.label = 'FLARES DPL'

        self.ads = 'https://ui.adsabs.harvard.edu/abs/2020arXiv200806057V/abstract'
        self.arxiv = 'https://arxiv.org/abs/2008.06057'

        self.redshifts = [5., 6., 7., 8., 9., 10.]  # array of redshifts
        self.phi_star = [-3.771, -4.054, -4.5, -4.605, -3.812, -4.148]  # array of log10(phi_star) values
        self.M_star = [-21.658, -21.446, -21.380, -20.966, -19.712, -19.658]  # array of M_star values
        self.alpha = [-2.034, -2.218, -2.500, -2.674, -2.567, -3.008]  # array of alpha values
        self.beta = [-4.306, -5.194, -5.190, -4.773, -4.467, -4.864]  # array of alpha values

        super().__init__()


class TNG_A(existing_model):
    # --- LF evolution model based on TNG dust model A, Vogelsberger+(2019)

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'Illustris TNG: Model-A (Vogelsberger+2019)'
        self.ref = 'Vogelsberger+2019'
        self.type = 'hydro'
        self.LF_model = 'Schechter'
        self.label = 'Illustris TNG: Model-A'

        self.ads = 'https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.5167V/abstract'
        self.arxiv = 'https://arxiv.org/abs/1904.07238'

        self.redshifts = [5., 6., 7., 8.]  # array of redshifts
        self.phi_star = [-3.244, -3.079, -3.846, -4.445]  # array of log10(phi_star) values
        self.M_star = [-21.17, -20.61, -21.18, -21.38]  # array of M_star values
        self.alpha = [-1.924, -1.876, -2.133, -2.280]  # array of alpha values

        super().__init__()


class TNG_B(existing_model):
    # --- LF evolution model based on TNG dust model B (2019), Vogelsberger+(2019)

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'Illustris TNG: Model-B (Vogelsberger+2019)'
        self.ref = 'Vogelsberger+2019'
        self.type = 'hydro'
        self.LF_model = 'Schechter'
        self.label = 'Illustris TNG: Model-B'

        self.ads = 'https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.5167V/abstract'
        self.arxiv = 'https://arxiv.org/abs/1904.07238'

        self.redshifts = [5., 6., 7., 8.]  # array of redshifts
        self.phi_star = [-3.107, -3.025, -3.418, -4.111]  # array of log10(phi_star) values
        self.M_star = [-20.95, -20.52, -20.58, -20.86]  # array of M_star values
        self.alpha = [-1.884, -1.833, -1.967, -2.216]  # array of alpha values

        super().__init__()


class TNG_C(existing_model):
    # --- LF evolution model based on TNG dust model C (2019), Vogelsberger+(2019)

    def __init__(self):
        # Contains model redshift range (must be increasing) and corresponding LF evolution model parameters
        # Custom models should be created following the same form

        self.name = 'Illustris TNG: Model-C (Vogelsberger+2019)'
        self.ref = 'Vogelsberger+2019'
        self.type = 'hydro'
        self.LF_model = 'Schechter'
        self.label = 'Illustris TNG: Model-C'

        self.ads = 'https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.5167V/abstract'
        self.arxiv = 'https://arxiv.org/abs/1904.07238'

        self.redshifts = [5., 6., 7., 8.]  # array of redshifts
        self.phi_star = [-3.398, -3.608, -4.209, -4.714]  # array of log10(phi_star) values
        self.M_star = [-21.21, -21.31, -21.47, -21.44]  # array of M_star values
        self.alpha = [-1.941, -2.042, -2.279, -2.455]  # array of alpha values

        super().__init__()
