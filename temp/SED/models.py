

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


    def get_Lnu(self, sfzh, SED_p, dust = False, fast = False):

        sed = SED()


        SFZH = np.expand_dims(sfzh, axis=2)

        sed.stellar = core.sed(self.lam)
        sed.stellar.lnu = np.sum(self.grid['stellar'] * SFZH, axis=(0,1))
        sed.stellar.lnu[self.lam<912.] *= SED_p['fesc']

        sed.nebular = core.sed(self.lam)
        sed.nebular.lnu = np.sum(self.grid['nebular'] * SFZH, axis=(0,1))

        sed.total = core.sed(self.lam)
        sed.total.lnu = sed.stellar.lnu + (1.-SED_p['fesc'])*sed.nebular.lnu


        if dust:

            dust_model, dust_model_params = dust


            # --- create SED of young and old components

            if dust_model == 'pacman':

                dmp = dust_model_params

                # --- see Ciaran's thesis example for how to use this


                i = np.where(self.grid['log10age']>dmp['log10age_BC'])[0][0]

                sfzh_young = copy.copy(sfzh)
                sfzh_young[i:,:] = 0.0

                sfzh_old = copy.copy(sfzh)
                sfzh_old[:i,:] = 0.0

                SFZH_young = np.expand_dims(sfzh_young, axis=2)
                SFZH_old = np.expand_dims(sfzh_old, axis=2)


                tau_BC = 10**(dmp['tau_ISM_to_BC']*SED_p['log10tau_V']) * dust_curves.simple(params = {'slope': dmp['alpha_BC']}).tau(self.lam)

                tau_ISM = 10**(SED_p['log10tau_V']) * dust_curves.simple(params = {'slope': dmp['alpha_ISM']}).tau(self.lam)

                T_young = np.exp(-tau_BC)*np.exp(-tau_ISM)
                T_old = np.exp(-tau_ISM)


                # --- young stellar light which escapes with no gas/dust reprocessing
                sed.stellar_young_esc = core.sed(self.lam)
                sed.stellar_young_esc.lnu = SED_p['fesc'] * np.sum(self.grid['stellar'] * SFZH_young, axis=(0,1))

                sed.stellar_young_nesc = core.sed(self.lam)
                sed.stellar_young_nesc.lnu = T_young*(1-SED_p['fesc']) * np.sum(self.grid['stellar'] * SFZH_young, axis=(0,1))
                sed.stellar_young_nesc.lnu[self.lam<912.] *= SED_p['fesc']


                sed.nebular_young_nesc = core.sed(self.lam)
                sed.nebular_young_nesc.lnu = T_young*(1-SED_p['fesc']) * np.sum(self.grid['nebular'] * SFZH_young, axis=(0,1))

                sed.total_young_esc = core.sed(self.lam)
                sed.total_young_esc.lnu = sed.stellar_young_esc.lnu # should be same

                sed.total_young = core.sed(self.lam)
                sed.total_young.lnu = sed.total_young_esc.lnu + sed.stellar_young_nesc.lnu + sed.nebular_young_nesc.lnu

                # --- old stellar light which escapes with no gas/dust reprocessing
                sed.stellar_old_esc = core.sed(self.lam)
                sed.stellar_old_esc.lnu = SED_p['fesc'] * np.sum(self.grid['stellar'] * SFZH_old, axis=(0,1))

                sed.stellar_old_nesc = core.sed(self.lam)
                sed.stellar_old_nesc.lnu = T_old*(1-SED_p['fesc']) * np.sum(self.grid['stellar'] * SFZH_old, axis=(0,1))
                sed.stellar_old_nesc.lnu[self.lam<912.] *= SED_p['fesc']

                sed.nebular_old_nesc = core.sed(self.lam)
                sed.nebular_old_nesc.lnu = T_old*(1-SED_p['fesc']) * np.sum(self.grid['nebular'] * SFZH_old, axis=(0,1))

                sed.total_old_esc = core.sed(self.lam)
                sed.total_old_esc.lnu = sed.stellar_young_esc.lnu # should be same

                sed.total_old = core.sed(self.lam)
                sed.total_old.lnu = sed.total_old_esc.lnu + sed.stellar_old_nesc.lnu + sed.nebular_old_nesc.lnu

                # --- total
                sed.total = core.sed(self.lam)
                sed.total.lnu = sed.total_young.lnu + sed.total_old.lnu


            elif dust_model == 'CF00':

                print('not yet implemented')

            else:

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
