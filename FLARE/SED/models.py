

import numpy as np
import pickle

from . import core
from ..core import *

import copy

from . import dust_curves

class empty: pass



def beta(lam, slope, normalisation, normalisation_wavelength = 1500., include_ISM = True):

    model = core.sed(lam)

    model.lnu = normalisation*(lam/normalisation_wavelength)**(slope + 2.0)

    if include_ISM: model.lnu[lam<912.] = 0.0 # add ISM absorption

    return model



class SED():

    def A(self,l):

        return -2.5*np.log10(np.interp(l, self.total.lam, self.total.lnu)/np.interp(l, self.total_intrinsic.lam, self.total_intrinsic.lnu))

    def A1500(self):

        return self.A(1500.)




class SPS():

    def __init__(self, grid, path_to_SPS_grid = '/data/SPS/nebular/3.0/'):

        self.grid = pickle.load(open(FLARE_dir + path_to_SPS_grid + grid + '/nebular.p','rb'), encoding='latin1')

        self.lam = self.grid['lam']


    def get_Lnu(self, SFZH, SED_p, dust = False, fast = False):

        sed = SED()


        SFZH = np.expand_dims(SFZH, axis=2)

        sed.stellar = core.sed(self.lam)
        sed.stellar.lnu = np.sum(self.grid['stellar'] * SFZH, axis=(0,1))
        sed.stellar.lnu[self.lam<912.] *= SED_p['fesc']

        sed.nebular = core.sed(self.lam)
        sed.nebular.lnu = np.sum(self.grid['nebular'] * SFZH, axis=(0,1))

        sed.total = core.sed(self.lam)
        sed.total.lnu = sed.stellar.lnu + (1.-SED_p['fesc'])*sed.nebular.lnu

        if dust:

            dust_model, dust_model_params = dust

            if not fast:
                sed.stellar_intrinsic = copy.deepcopy(sed.stellar)
                sed.nebular_intrinsic = copy.deepcopy(sed.nebular)
                sed.total_intrinsic = copy.deepcopy(sed.total)

            tau = 10**(SED_p['log10tau_V']) * getattr(dust_curves, dust_model)(params = dust_model_params).tau(self.lam)

            T = np.exp(-tau)

            sed.stellar.lnu *= T
            sed.nebular.lnu *= T
            sed.total.lnu *= T

        return sed

    def get_Q(self, SFZH):

        return np.log10(np.sum(10**self.grid['log10Q'] * SFZH, axis=(0,1)))





class lines():

    def __init__(self, grid, path_to_SPS_grid = FLARE_dir + '/data/SPS/nebular/3.0', dust = False):

        self.grid = pickle.load(open(path_to_SPS_grid + grid + '/lines.p','rb'), encoding='latin1')


    def get_L(self, SFZH, line):

        L = np.sum(10**self.grid[line]['luminosity'] * SFZH, axis=(0,1))

        return L
