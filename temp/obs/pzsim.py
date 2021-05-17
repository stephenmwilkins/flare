import sys
import os
import copy
from types import SimpleNamespace

import h5py

import numpy as np


import FLARE
import FLARE.filters
import FLARE.obs
import FLARE.obs.misc
import FLARE.obs.photometry
import FLARE.obs.plots
import FLARE.SED.models
import FLARE.SED.SFZH as SFZH
import FLARE.surveys
import FLARE.photom
import FLARE.obs.EAZY as eazy

uniform = lambda x: np.random.uniform(low = x[0], high = x[1])

def removekey(d, key):
    r = dict(d)
    del r[key]
    return r

def get_name_shape(name, item):
    shape = ''
    if hasattr(item, 'value'):
        shape = item.shape
    print(name, shape)





class simulation():

    def __init__(self, Filters, prange, cosmo = FLARE.default_cosmo(),  SPS = False, verbose = False, SimType = 'const'):

        self.SimType = SimType

        self.Filters = Filters

        self.prange = prange
        self.cosmo = cosmo

        self.verbose = verbose

        self.SPS = SPS # necessary for SED generating if not \beta model

    def run(self, N, ID = False, OutputFolder = False):

        Sources = self.run_many(N)

        self.hf = self.write_to_HDF5(Sources, ID = ID, OutputFolder = OutputFolder)


    def run_many(self, N):

        return [self.run_single() for i in range(N)]

    def run_single(self):

        s = SimpleNamespace() # output object

        # --- choose parameters

        p = {}

        # --- determine redshift first

        if self.prange['z'][0] == 'uniform':
            p['z'] = uniform(self.prange['z'][1])
        elif self.prange['z'][0] == 'delta':
            p['z'] = self.prange['z'][1]


        # --- set age of Universe constraint on duration

        if 'duration' in self.prange.keys():
            if self.prange['duration'][0] == 'uniform':
                self.prange['duration'][1][1] = self.cosmo.age(p['z']).to('Myr').value

        # --- for isntantanesous set age = age of the Universe

        if self.SimType == 'instantaneous':
            p['log10age'] = self.cosmo.age(p['z']).to('Myr').value


        # --- set other parameters

        for parameter, value in removekey(self.prange, 'z').items():
            if value[0] == 'uniform':
                p[parameter] = uniform(value[1])
            elif value[0] == 'delta':
                p[parameter] = value[1]
            elif value[0] == 'normal':
                p[parameter] = value[1][0] + value[1][1] * np.random.randn()

        if 'duration' in self.prange.keys():
            p['log10_duration'] = np.log10(p['duration']) + 6.



        if self.verbose:
            print('--- Input parameters')
            for k,v in p.items(): print(f'{k}: {v:.2f}')


        s.input = p

        # --- create SED
        s.Fnu= self.create_SED(p)

        return s





    def create_SED(self, p):


        derived = {}


        if self.SimType == 'beta':

            rest_lam = np.arange(0., 5000., 1.)

            self.F = FLARE.filters.add_filters(self.Filters, new_lam = rest_lam * (1. + p['z']))

            sed = FLARE.SED.models.beta(rest_lam, p['beta'], 1.0, normalisation_wavelength = 1500.)

        elif self.SimType == 'const':

            # --- get SFH for a given choice of


            # -- this was the old behaviour where we later rescaled by SNR
            if 'SNR' in p:
                sfzh, sfr = SFZH.constant(self.SPS.grid['log10age'], self.SPS.grid['log10Z'] , {'log10_duration': p['log10_duration'], 'log10Z': p['log10Z'], 'log10M*': 8.0})
            if 'log10M*' in p:
                sfzh, sfr = SFZH.constant(self.SPS.grid['log10age'], self.SPS.grid['log10Z'] , {'log10_duration': p['log10_duration'], 'log10Z': p['log10Z'], 'log10M*': p['log10M*']})

            # --- generate SED for a given choice of parameters

            SED = self.SPS.get_Lnu(sfzh, {'fesc': p['fesc'], 'log10tau_V': p['log10tau_V']}, dust = ('simple', {'slope': -1.0}))

            sed = SED.total
            self.F = FLARE.filters.add_filters(self.Filters, new_lam = sed.lam * (1.+p['z']))

        elif self.SimType == 'instantaneous': # instantaneous burst as close to age of the Universe as possible

            # --- get SFH for a given choice of

            sfzh, sfr = SFZH.instantaneous(self.SPS.grid['log10age'], self.SPS.grid['log10Z'] , {'log10age': p['log10age'], 'log10Z': p['log10Z'], 'log10M*': p['log10M*']})

            # --- generate SED for a given choice of parameters

            SED = self.SPS.get_Lnu(sfzh, {'fesc': 0.0, 'log10tau_V': p['log10tau_V']}, dust = ('simple', {'slope': -1.0}))

            sed = SED.total
            self.F = FLARE.filters.add_filters(self.Filters, new_lam = sed.lam * (1.+p['z']))


        else:

            print('WARNING! Incorrect simulation type set')



        sed.get_fnu(self.cosmo, p['z']) # --- generate observed frame spectrum (necessary to get broad band photometry)
        sed.get_Fnu(self.F) # --- generate broadband photometry

        Fnu = {f: sed.Fnu[f]for f in self.Filters}

        return Fnu










    def write_to_HDF5(self, Sources, ID = False, OutputFolder = False, return_struct = True):

        if not ID:
            ID = np.random.randint(1E9)

        self.ID = ID

        if not OutputFolder: OutputFolder = f'data'

        OutputFile = f'{OutputFolder}/{ID}.h5'


        if self.verbose:
            print()
            print('-'*10, 'Writing to HDF5')
            print(f'Output file: {OutputFile}')

        Sources = np.array(Sources)

        hf = h5py.File(OutputFile, 'w')

        # --- detected flag

        # --- necessary for future collation

        hf.attrs['total'] = len(Sources)

        hf.attrs['filters'] = np.array(self.Filters, dtype='S')

        # --- input parameters

        for k in Sources[0].input.keys():
            data = np.array([source.input[k] for source in Sources])
            hf.create_dataset(f'input/{k}', data=data)


        # --- fluxes

        for k in Sources[0].Fnu.keys():
            data = np.array([source.Fnu[k] for source in Sources])
            hf.create_dataset(f'Fnu/{k}', data=data)

        if return_struct:
            return hf
        else:
            hf.close()
