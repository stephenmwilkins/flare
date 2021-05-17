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




    def __init__(self, DepthModel, ReferenceFilter, Filters, prange, cosmo = FLARE.default_cosmo(),  SPS = False, verbose = False):

        
        self.DepthModel = DepthModel
        self.depths = FLARE.surveys.depthmodel[self.DepthModel]

        self.ReferenceFilter = ReferenceFilter

        self.Filters = Filters

        self.prange = prange
        self.cosmo = cosmo

        self.verbose = verbose
        
        self.SPS = SPS # necessary for SED generating if not \beta model




    def run_EAZY(self, EAZYFilters, EAZYID):
    
    
        F = FLARE.filters.add_filters(EAZYFilters)
    
        hf_EAZY = eazy.eazy(ID = self.ID).run(self.hf, F, phot_type = '')
    
        # --- append EAZY group to original file
        
        hf_EAZY.copy(f'EAZY', self.hf, name = f'EAZY/{EAZYID}')

        self.hf[f'EAZY/{EAZYID}'].attrs['filters'] = EAZYFilters




    def run(self, N):
    
        Sources = self.run_many(N)
        
        self.hf = self.write_to_HDF5(Sources)
            
    
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
        s.Fnu, s.obs = self.create_SED(p)
        
        
        # --- measure beta
        
        s.derived = {}
        s.derived['beta_int'] = FLARE.obs.misc.measure_beta(p['z'], s.Fnu, self.F)
        s.derived['beta'] = FLARE.obs.misc.measure_beta(p['z'], s.obs.flux, self.F)  # NOT CORRECT, NEED TO RUN PZ FIRST
        
        if self.verbose:
            print(rf"intrinsic \beta: {s.derived['beta_int']:.2f}")
            print(rf"observed \beta: {s.derived['beta']:.2f}")
        
        return s
    
    
        
        

    def create_SED(self, p):


        derived = {}


        if 'beta' in p:

            rest_lam = np.arange(0., 5000., 1.)

            self.F = FLARE.filters.add_filters(self.Filters, new_lam = rest_lam * (1. + p['z'])) 
    
            sed = FLARE.SED.models.beta(rest_lam, p['beta'], 10**p['log10L1500'], normalisation_wavelength = 1500.)

        else:
        
            # --- get SFH for a given choice of 
            
            sfzh, sfr = SFZH.constant(self.SPS.grid['log10age'], self.SPS.grid['log10Z'] , {'log10_duration': p['log10_duration'], 'log10Z': p['log10Z'], 'log10M*': 8.0})

            # --- generate SED for a given choice of parameters

            SED = self.SPS.get_Lnu(sfzh, {'fesc': p['fesc'], 'log10tau_V': p['log10tau_V']}, dust_model = 'very_simple')

            sed = SED.total 

            self.F = FLARE.filters.add_filters(self.Filters, new_lam = sed.lam * (1.+p['z'])) 
            
            restF = FLARE.filters.add_filters(self.Filters, new_lam = sed.lam * (1.+p['z'])) 
            


        sed.get_fnu(self.cosmo, p['z']) # --- generate observed frame spectrum (necessary to get broad band photometry)
        sed.get_Fnu(self.F) # --- generate broadband photometry
        
        Fnu = {f: sed.Fnu[f]/sed.Fnu[self.ReferenceFilter] for f in self.Filters}

        # --- add in noise
        
        
        ref_noise = Fnu[self.ReferenceFilter]/p['SNR']
        
        ref_depth =  FLARE.photom.m_to_flux(self.depths[self.ReferenceFilter])
        
        noise_scaling = {f: FLARE.photom.m_to_flux(self.depths[f])/ref_depth for f in self.Filters}

        obs = SimpleNamespace()

        

        obs.error = {f: ref_noise * noise_scaling[f] for f in self.Filters}

        obs.flux = {f: Fnu[f] + obs.error[f]*np.random.randn() for f in self.Filters}

        obs.flux = {f: obs.flux[f]/obs.flux[self.ReferenceFilter] for f in self.Filters}


        # --- print out fluxes

        if self.verbose: 
            print('-'*5, 'SED')
            for f in self.Filters: print(f'{f}: {Fnu[f]:.2f} | {obs.flux[f]:.2f} {obs.error[f]:.2f} {obs.flux[f]/obs.error[f]:.2f}')




        return Fnu, obs
        
    
    

            
        




    def write_to_HDF5(self, Sources, ID = np.random.randint(1E9), OutputFolder = False, return_struct = True):
          
        self.ID = ID
        
        if not OutputFolder: OutputFolder = f'data/{self.DepthModel}/'

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


        # --- Observed Properties

        for f in self.Filters:  
            hf.create_dataset(f'obs/{f}/flux', data = np.array([source.obs.flux[f] for source in Sources]))
            hf.create_dataset(f'obs/{f}/error', data = np.array([source.obs.error[f] for source in Sources]))
      
      
        if return_struct:
            return hf
        else:
            hf.close()
        



    
    





  
  
  
  

