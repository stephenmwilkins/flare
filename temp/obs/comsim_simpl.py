import sys
import os
import copy
from types import SimpleNamespace

import h5py

class empty(): pass

import numpy as np
import pickle

import FLARE
import FLARE.filters
import FLARE.obs
import FLARE.obs.photometry
import FLARE.SED.models
import FLARE.SED.SFZH as SFZH
import FLARE.surveys
import FLARE.photom
import FLARE.obs.EAZY as eazy
import FLARE.obs.plots

uniform = lambda x: np.random.uniform(low = x[0], high = x[1])

def removekey(d, key):
    r = dict(d)
    del r[key]
    return r







class simple():

    def __init__(self, surveyName, fieldName, prange = False, cosmo = FLARE.default_cosmo(), SPS = False, verbose = False, make_plots = False):

        self.SimType = 'idealised'
        self.surveyName = surveyName
        self.fieldName = fieldName
        self.prange = prange
        self.cosmo = cosmo

        self.verbose = verbose
        self.make_plots = make_plots

        self.Field = FLARE.surveys.surveys[self.surveyName].fields[self.fieldName]
        self.Filters = self.Field.filters
        self.Backgrounds = FLARE.obs.FieldBackgrounds(self.Field, verbose = self.verbose) # --- create Background object

        self.SPS = SPS # necessary for SED generating if not \beta model


    def run_list(self, list, run_EAZY = False):

        Sources = []

        for i in range(len(list['z'])):
            p = {k:v[i] for k,v in list.items()}
            Sources.append(self.run_single(p))

        hf = self.write_to_HDF5(Sources, OutputFolder = f'data/{self.surveyName}/{self.fieldName}/{self.SimType}')

        if self.verbose:

            def get_name_shape(name, item):
                shape = ''
                if hasattr(item, 'value'):
                    shape = item.shape
                print(name, shape)

            hf.visititems(get_name_shape)


        if run_EAZY:

            hf_EAZY = eazy.eazy(ID = self.ID, create_POFZ_FILE = True).run(hf, self.F)

            # --- append EAZY group to original file

            hf_EAZY.copy('EAZY', hf)


        if self.verbose:

            def get_name_shape(name, item):
                shape = ''
                if hasattr(item, 'value'):
                    shape = item.shape
                print(name, shape)

            hf.visititems(get_name_shape)


        hf.close()



    def run(self, N, run_EAZY = False):

        # run for many galaxies, run EAZY, and output as HDF5

        Sources = self.run_many(N)

        hf = self.write_to_HDF5(Sources)

        if self.verbose:

            def get_name_shape(name, item):
                shape = ''
                if hasattr(item, 'value'):
                    shape = item.shape
                print(name, shape)

            hf.visititems(get_name_shape)


        if run_EAZY:

            hf_EAZY = eazy.eazy(ID = self.ID).run(hf, self.F)

            # --- append EAZY group to original file

            hf_EAZY.copy('EAZY', hf)


        if self.verbose:

            def get_name_shape(name, item):
                shape = ''
                if hasattr(item, 'value'):
                    shape = item.shape
                print(name, shape)

            hf.visititems(get_name_shape)


        hf.close()


    def run_many(self, N):

        return [self.run_single(self.get_p()) for i in range(N)]


    def get_p(self):

        # --- This gets the parameters for the model


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

        return p






    def run_single(self, p):

        if self.verbose:
            print('--- Input parameters')
            for k,v in p.items(): print(f'{k}: {v:.2f}')

        s = SimpleNamespace() # output object

        s.InputParameters = p

        # --- create SED
        Fnu, derived = self.create_SED(p)

        s.Fnu = Fnu


        # --- add noise



        # --- check whether detected


        s.derived = derived

        s.detected = detected
        s.ObservedProperties = ObservedProperties



        return s





    def create_SED(self, p):


        derived = {}


        if 'beta' in p:

            rest_lam = np.arange(0., 5000., 1.)

            self.F = FLARE.filters.add_filters(self.Filters, new_lam = rest_lam * (1. + p['z']))

            sed = FLARE.SED.models.beta(rest_lam, p['beta'], 10**p['log10L1500'], normalisation_wavelength = 1500.)

        else:

            # --- get SFH for a given choice of

            sfzh, sfr = SFZH.constant(self.SPS.grid['log10age'], self.SPS.grid['log10Z'] , {'log10_duration': p['log10_duration'], 'log10Z': p['log10Z'], 'log10M*': p['log10M*']})

            # --- generate SED for a given choice of parameters

            SED = self.SPS.get_Lnu(sfzh, {'fesc': p['fesc'], 'log10tau_V': p['log10tau_V']}, dust_model = 'very_simple')

            sed = SED.total

            self.F = FLARE.filters.add_filters(self.Filters, new_lam = sed.lam * (1.+p['z']))

            restF = FLARE.filters.add_filters(self.Filters, new_lam = sed.lam * (1.+p['z']))

            derived['log10L1500'] = np.log10(sed.return_Lnu(FLARE.filters.add_filters(['FAKE.TH.FUV'], new_lam = sed.lam))['FAKE.TH.FUV'])
            derived['log10SFR'] = np.log10(sfr)

            if self.verbose:
                print(rf'log10(L1500/erg/s/Hz)={derived["log10L1500"]:2f}')
                print(rf'log10(SFR/M/yr)={derived["log10SFR"]:2f}')

        sed.get_fnu(self.cosmo, p['z']) # --- generate observed frame spectrum (necessary to get broad band photometry)
        sed.get_Fnu(self.F) # --- generate broadband photometry
        Fnu = {f: [sed.Fnu[f]] for f in self.Filters}

        # --- print out fluxes

        if self.verbose:
            print('-'*5, 'SED')
            for f in self.Filters: print('{0}: {1:.2f}/nJy'.format(f, Fnu[f][0]))

        return Fnu, derived










    def write_to_HDF5(self, Sources, ID = np.random.randint(1E9), OutputFolder = False):


        self.ID = ID

        if not OutputFolder: OutputFolder = f'data/{self.surveyName}/{self.fieldName}/{self.SimType}/individual'

        OutputFile = f'{OutputFolder}/{ID}.h5'


        if self.verbose:
            print()
            print('-'*10, 'Writing to HDF5')
            print(f'Output file: {OutputFile}')

        Sources = np.array(Sources)

        hf = h5py.File(OutputFile, 'w')

        # --- detected flag

        detected = np.array([source.detected for source in Sources])

        hf.create_dataset('detected', data = detected)

        # --- necessary for future collation

        hf.attrs['total'] = len(Sources)
        hf.attrs['detected'] = len(Sources[detected])

        hf.attrs['filters'] = np.array(self.Filters, dtype='S')

        # --- input parameters

        for k in Sources[0].InputParameters.keys():
            data = np.array([source.InputParameters[k] for source in Sources])
            hf.create_dataset(f'input/{k}', data=data)

        # --- derived properties

        for k in Sources[0].derived.keys():
            data = np.array([source.derived[k] for source in Sources])
            hf.create_dataset(f'derived/{k}', data=data)

        # --- fluxes

        for k in Sources[0].Fnu.keys():
            data = np.array([source.Fnu[k][0] for source in Sources])
            hf.create_dataset(f'Fnu/{k}', data=data)


        # --- Model properties

        SizeTypes = ['COG', 'pixel']

        for f in self.Filters:

            # hf.create_dataset(f'model/{f}/photometry/total_flux', data = np.array([source.ModelProperties[f]['photometry']['total'].flux for source in Sources]))

            for SizeType in SizeTypes:
                hf.create_dataset(f'model/{f}/sizes/{SizeType}', data = np.array([source.ModelProperties[f]['sizes'][SizeType].radius for source in Sources]))


        # --- Observed Properties

        PhotTypes = ['circular_kron', 'optimum_aperture', 'small_circular_kron', 'ISO']

        for f in self.Filters:
            for PhotType in PhotTypes:
                hf.create_dataset(f'obs/{f}/photometry/{PhotType}_flux', data = np.array([source.ObservedProperties[f]['photometry'][PhotType].flux for source in Sources[detected]]))
                hf.create_dataset(f'obs/{f}/photometry/{PhotType}_error', data = np.array([source.ObservedProperties[f]['photometry'][PhotType].error for source in Sources[detected]]))


            for SizeType in SizeTypes:
                hf.create_dataset(f'obs/{f}/sizes/{SizeType}', data = np.array([source.ObservedProperties[f]['sizes'][SizeType].radius for source in Sources[detected]]))


        # --- print structure of file with shape of each object

        return hf
