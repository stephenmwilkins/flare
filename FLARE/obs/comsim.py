import sys
import os
import copy
from types import SimpleNamespace

import scipy.special

import h5py

class empty(): pass

import numpy as np
import pickle

from photutils import detect_sources
from photutils import deblend_sources
from photutils import source_properties

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

import SynthObs
from SynthObs.Morph import measure
from SynthObs.Morph import PSF
from SynthObs.Morph import images

uniform = lambda x: np.random.uniform(low = x[0], high = x[1])

def removekey(d, key):
    r = dict(d)
    del r[key]
    return r







class simulation():


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


        p['r_eff'] = 10**p['log10r_eff']

        if self.verbose:
            print('--- Input parameters')
            for k,v in p.items(): print(f'{k}: {v:.2f}')
            print('r_eff/pix:', p['r_eff']*self.cosmo.arcsec_per_kpc_proper(p['z']).value/self.Field.pixel_scale)
            print('r_eff/arcsec:', p['r_eff']*self.cosmo.arcsec_per_kpc_proper(p['z']).value)

        s = SimpleNamespace() # output object

        s.InputParameters = p

        # --- create SED
        Fnu, derived = self.create_SED(p)

        s.Fnu = Fnu
        s.derived = derived

        # --- create the original CutoutImages
        CutoutImages = self.create_OriginalCutoutImages()

        # --- create the ModelImages
        ModelImages = self.create_ModelImages(Fnu, p, CutoutImages)

        # --- Calculate intrinsic properties of the Model
        s.ModelProperties = {filter: FLARE.obs.photometry.measure_model_properties(ModelImages[filter], verbose = self.verbose) for filter in self.Filters}

        # --- add model images to original CutoutImages
        for filter in self.Filters:
            CutoutImages[filter].sci += ModelImages[filter].sci
            CutoutImages[filter].filter = filter

        # --- Make detection images
        DetectionImages = {filter: CutoutImages[filter] for filter in self.DetectionFilters}
        DetectionImage = FLARE.obs.create_stack(DetectionImages)

        # --- detect sources

        detected, DetectionProperties, Mask, ExclusionMask, ObservedProperties = self.detect(DetectionImage, CutoutImages, threshold = self.threshold, npixels = self.npixels)

        s.detected = detected
        s.ObservedProperties = ObservedProperties

        if self.make_plots:

            #FLARE.obs.plots.make_plots({**OriginalCutoutImages}, signficance_plot = True, filename = 'temp/original.png', show = self.test)
            FLARE.obs.plots.make_plots({**ModelImages}, signficance_plot = False, filename = 'temp/model.png', show = self.verbose, fixed_range = True)
            FLARE.obs.plots.make_plots({**CutoutImages}, signficance_plot = True, filename = 'temp/stamps.png', show = self.verbose)

        # if test and detected:

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



    def create_ModelImages(self, Fnu, p, CutoutImages):

        width_arcsec = self.CutoutWidth * self.Field.pixel_scale # size of generated source in "

        ModelImages = images.Sersic(Fnu, p, self.Filters, self.cosmo, p['z'], width_arcsec, pixel_scale = self.Field.pixel_scale, PSFs = self.PSFs)

        for filter in self.Filters:
            ModelImages[filter].sci = ModelImages[filter].data * CutoutImages[filter].nJy_to_es # add a new entry to ModelImages that has units of e/s
            ModelImages[filter].nJy_to_es = CutoutImages[filter].nJy_to_es

        return ModelImages


    def detect(self, DetectionImage, CutoutImages, threshold, npixels):

        SegmentationImage = detect_sources(DetectionImage.sig, threshold, npixels = npixels)

        if type(SegmentationImage) is not type(None):

            print('HERE')

            DeblendedSegmentationImage = deblend_sources(DetectionImage.sig, SegmentationImage, npixels = npixels, nlevels=32)
            AllSourceProperties = source_properties(DetectionImage.sig, DeblendedSegmentationImage)

            Cat = AllSourceProperties.to_table()
            x, y = Cat['xcentroid'].value, Cat['ycentroid'].value
            r = np.sqrt((x - self.CutoutWidth/2)**2 + (y - self.CutoutWidth/2)**2)
            tol = 2. # pixels, impossible for more than one source I think, but at this tolerance the source could be shifted.
            s = r<tol

            if len(x[s])==1:

                detected = True
                idx = np.where(s==True)[0][0]
                SourceProperties = AllSourceProperties[idx]
                DetectionProperties, Mask, ExclusionMask = FLARE.obs.photometry.measure_core_properties(SourceProperties, DetectionImage, DeblendedSegmentationImage, verbose = self.verbose)

                if self.verbose:
                    print()
                    print('-'*10, 'Observed Properties')
                ObservedProperties = {filter: FLARE.obs.photometry.measure_properties(DetectionProperties, CutoutImages[filter], Mask, ExclusionMask, verbose = self.verbose) for filter in self.Filters}

                return detected, DetectionProperties, Mask, ExclusionMask, ObservedProperties

            else:

                detected = False
                if self.verbose: print('**** SOURCES DETECTED BUT NOT IN MIDDLE OF NO SOURCE DETECTED')

                return detected, None, None, None, None
        else:

            detected = False
            if self.verbose: print('**** NO SOURCES DETECTED AT ALL')

            return detected, None, None, None, None







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


class idealised(simulation):


    def __init__(self, surveyName, fieldName, DetectionFilters, prange = False, cosmo = FLARE.default_cosmo(), CutoutWidth = 101, threshold = 2.5, npixels = 5,  SPS = False, verbose = False, make_plots = False):



        self.SimType = 'idealised'
        self.surveyName = surveyName
        self.fieldName = fieldName
        self.DetectionFilters = DetectionFilters
        self.prange = prange
        self.cosmo = cosmo
        self.CutoutWidth = CutoutWidth
        self.threshold = threshold
        self.npixels = npixels

        self.verbose = verbose
        self.make_plots = make_plots

        self.Field = FLARE.surveys.surveys[self.surveyName].fields[self.fieldName]
        self.Filters = self.Field.filters
        self.Backgrounds = FLARE.obs.FieldBackgrounds(self.Field, verbose = self.verbose) # --- create Background object

        # --- define PSF
        TargetPSF = PSF.PSF(self.Filters[-1])
        self.PSFs = {f: TargetPSF for f in self.Filters} # force all to use the PSF of the last filter

        self.SPS = SPS # necessary for SED generating if not \beta model

    def create_OriginalCutoutImages(self):

        return {filter: self.Backgrounds[filter].create_background_image(self.CutoutWidth) for filter in self.Filters}


# class real(simulation):
#
#     def __init__():
#





















class simple():

    def __init__(self, surveyName, fieldName, detection_criteria, prange = False, cosmo = FLARE.default_cosmo(), SPS = False, verbose = False, make_plots = False):

        self.SimType = 'simple'
        self.surveyName = surveyName
        self.fieldName = fieldName
        self.prange = prange
        self.cosmo = cosmo

        self.verbose = verbose
        self.make_plots = make_plots

        self.Field = FLARE.surveys.surveys[self.surveyName].fields[self.fieldName]
        self.Filters = self.Field.filters

        self.SPS = SPS # necessary for SED generating if not \beta model

        self.depth =  {f: FLARE.photom.m_to_flux(self.Field.depths[f]) for f in self.Filters}

        self.detection_filter = detection_criteria[0]
        self.detection_snr = detection_criteria[1]
        self.detection_stretch = detection_criteria[2]

        if verbose:
            for f in self.Filters: print(f'{f} {self.depth[f]:.2f}')


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

            hf_EAZY = eazy.eazy(ID = self.ID, create_POFZ_FILE = False).run_new(hf, self.F, path = lambda f: f'ObsFnu/{f}/')

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

            hf_EAZY = eazy.eazy(ID = self.ID).run_new(hf, self.F, path = lambda f: f'obs/{f}/')

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
        s.derived = derived

        # --- add noise
        s.ObsFnu = {}

        if self.verbose: print('-'*5, 'Observed SED')

        for f in self.Filters:

            s.ObsFnu[f] = SimpleNamespace()
            s.ObsFnu[f].error = self.depth[f]
            s.ObsFnu[f].flux = s.Fnu[f][0] + s.ObsFnu[f].error*np.random.randn()

            if self.verbose:
                print(f'{f}: {s.ObsFnu[f].flux:.2f} {s.ObsFnu[f].flux/s.ObsFnu[f].error:.2f}')


        # --- check whether detected


        snr = s.ObsFnu[self.detection_filter].flux/s.ObsFnu[self.detection_filter].error

        if np.random.random() < 0.5*(1.0+scipy.special.erf((snr-self.detection_snr)/self.detection_stretch)):
            s.detected = True
        else:
            s.detected = False

        if self.verbose:
            print(f'SNR: {snr:.2f} Detected: {s.detected}')


        # s.detected = detected
        # s.ObservedProperties = ObservedProperties

        return s


    def create_SED(self, p):


        derived = {}


        if 'beta' in p:

            rest_lam = np.arange(0., 5000., 1.)

            self.F = FLARE.filters.add_filters(self.Filters, new_lam = rest_lam * (1. + p['z']))

            sed = FLARE.SED.models.beta(rest_lam, p['beta'], 10**p['log10L1500'], normalisation_wavelength = 1500.)

        elif 'log10M*' in p:

            print('WARNING: not yet implemented')

        elif 'log10L1500' in p:

            if 'duration' in p.keys():
                p['log10_duration'] = np.log10(p['duration']) + 6.

            sfzh, sfr = SFZH.constant(self.SPS.grid['log10age'], self.SPS.grid['log10Z'] , {'log10_duration': p['log10_duration'], 'log10Z': p['log10Z'], 'log10M*': 0.0})

            # --- generate SED for a given choice of parameters

            SED = self.SPS.get_Lnu(sfzh, {'fesc': p['fesc'], 'log10tau_V': p['log10tau_V']}, dust_model = 'very_simple')

            sed = SED.total

            log10L1500 = np.log10(sed.return_Lnu(FLARE.filters.add_filters(['FAKE.TH.FUV'], new_lam = sed.lam))['FAKE.TH.FUV'])

            derived['log10M*'] = p['log10L1500'] - log10L1500
            derived['log10SFR'] = np.log10(sfr) + derived['log10M*']

            sed.lnu *= 10**derived['log10M*']

            self.F = FLARE.filters.add_filters(self.Filters, new_lam = sed.lam * (1.+p['z']))


        sed.get_fnu(self.cosmo, p['z']) # --- generate observed frame spectrum (necessary to get broad band photometry)
        sed.get_Fnu(self.F) # --- generate broadband photometry
        Fnu = {f: [sed.Fnu[f]] for f in self.Filters}

        # --- print out fluxes



        if self.verbose:
            print('-'*5, 'Derived Quantities')
            for k,v in derived.items(): print(f'{k} {v:.2f}')
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



        # --- Observed Properties

        for f in self.Filters:
            hf.create_dataset(f'obs/{f}/flux', data = np.array([source.ObsFnu[f].flux for source in Sources])) # --- all sources, not just those detected
            hf.create_dataset(f'obs/{f}/error', data = np.array([source.ObsFnu[f].error for source in Sources]))


        # --- print structure of file with shape of each object

        return hf
