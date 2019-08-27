import sys
import os
import copy 
from types import SimpleNamespace

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
import FLARE.obs.plots
import FLARE.SED.models
import FLARE.surveys
import FLARE.photom
import FLARE.obs.EAZY as eazy

import SynthObs
from SynthObs.Morph import measure 
from SynthObs.Morph import PSF
from SynthObs.Morph import images

uniform = lambda x: np.random.uniform(low = x[0], high = x[1])








class simulation():


    def run(self, N, threshold = 2.5, npixels = 5, run_EAZY = False):
    
        Sources = self.run_many(N, threshold =threshold, npixels = npixels)
    
        hf = self.write_to_HDF5(Sources)
    
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
    
    
    
    def run_many(self, N, threshold = 2.5, npixels = 5):
    
        return [self.run_single(threshold = threshold, npixels = npixels) for i in range(N)]  
     
    def run_single(self, threshold = 2.5, npixels = 5, make_plots = False):
      
        s = SimpleNamespace() # output object
      
        # --- choose parameters

        p = {}
        for parameter, value in self.prange.items():
            if value[0] == 'uniform':
                p[parameter] = uniform(value[1])
            elif value[0] == 'delta':
                p[parameter] = value[1]
                    
       
        p['r_eff'] = 10**p['log10r_eff']
    
        if self.verbose:
            print('--- Input parameters')
            for k,v in p.items(): print(f'{k}: {v:.2f}')
            print('r_eff/pix:', p['r_eff']*self.cosmo.arcsec_per_kpc_proper(p['z']).value/self.Field.pixel_scale)
            print('r_eff/arcsec:', p['r_eff']*self.cosmo.arcsec_per_kpc_proper(p['z']).value)

        
        s.InputParameters = p

        # --- create SED
        Fnu = self.create_SED(p)
        
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

        detected, DetectionProperties, Mask, ExclusionMask, ObservedProperties = self.detect(DetectionImage, CutoutImages, threshold = threshold, npixels = npixels)
    
        s.detected = detected
        s.ObservedProperties = ObservedProperties

#         try:
#         
#             detected, DetectionProperties, Mask, ExclusionMask, ObservedProperties = self.detect(DetectionImage, CutoutImages, threshold = threshold, npixels = npixels)
#     
#             s.detected = detected
#             s.ObservedProperties = ObservedProperties
#     
#         except:
#         
#             if self.verbose: print('**** SOURCE DETECTION FAILED')
#             s.detected = False
        
    
        if make_plots:
        
            #FLARE.obs.plots.make_plots({**OriginalCutoutImages}, signficance_plot = True, filename = 'temp/original.png', show = self.test)
            FLARE.obs.plots.make_plots({**ModelImages}, signficance_plot = False, filename = 'temp/model.png', show = self.verbose, fixed_range = True)
            FLARE.obs.plots.make_plots({**CutoutImages}, signficance_plot = True, filename = 'temp/stamps.png', show = self.verbose)
    
        # if test and detected:
        
        
    
        return s
    
    
        
        

    def create_SED(self, p):

        rest_lam = np.arange(0., 5000., 1.)

        self.F = FLARE.filters.add_filters(self.Filters, new_lam = rest_lam * (1. + p['z'])) 
    
        # -------------------------------------------------------   
        # --- make model SED

        sed = FLARE.SED.models.beta(rest_lam, p['beta'], 10**p['log10L1500'], normalisation_wavelength = 1500.)
        sed.get_fnu(self.cosmo, p['z']) # --- generate observed frame spectrum (necessary to get broad band photometry)
        sed.get_Fnu(self.F) # --- generate broadband photometry
        Fnu = {f: [sed.Fnu[f]] for f in self.Filters}

        if self.verbose: 
            print('-'*5, 'SED')
            for f in self.Filters: print('{0}: {1:.2f}/nJy'.format(f, Fnu[f][0]))

        return Fnu
        
    
    
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
        
        print(self.Filters)
        
        hf.attrs['filters'] = np.array([self.Filters], dtype='S')

        print(hf.attrs['filters'])

        # --- input parameters

        for k in Sources[0].InputParameters.keys():
            data = np.array([source.InputParameters[k] for source in Sources])
            hf.create_dataset(f'input/{k}', data=data)


        # --- Model properties

        SizeTypes = ['COG', 'pixel']

        for f in self.Filters:  

            hf.create_dataset(f'model/{f}/photometry/total_flux', data = np.array([source.ModelProperties[f]['photometry']['total'].flux for source in Sources]))

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

    
    def __init__(self, surveyName, fieldName, DetectionFilters, prange, cosmo = FLARE.default_cosmo(), CutoutWidth = 101, verbose = False):

        
        self.SimType = 'idealised'
        self.surveyName = surveyName
        self.fieldName = fieldName
        self.DetectionFilters = DetectionFilters
        self.prange = prange
        self.cosmo = cosmo
        self.CutoutWidth = CutoutWidth
        self.verbose = verbose
        
        self.Field = FLARE.surveys.surveys[self.surveyName].fields[self.fieldName] # deepest portion of the XDF
        self.Filters = self.Field.filters
        self.Backgrounds = FLARE.obs.FieldBackgrounds(self.Field, verbose = self.verbose) # --- create Background object

        # --- define PSF
        TargetPSF = PSF.HubblePSF(self.Filters[-1])
        self.PSFs = {f: TargetPSF for f in self.Filters} # force all to use the PSF of the last filter


    def create_OriginalCutoutImages(self):
    
        return {filter: self.Backgrounds[filter].create_background_image(self.CutoutWidth) for filter in self.Filters} 





# class real(simulation):
# 
#     def __init__():
# 




    

    
    





  
  
  
  

