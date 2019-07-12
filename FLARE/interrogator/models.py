

import sys
import os
import numpy as np
import matplotlib.pyplot as plt 
import pickle

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

class empty(): pass

import FLARE
from ..core import * 
from .. import filters
from ..SED import models
from ..SED import SFZH
from ..SED import IGM






def default_priors(model_definitions):

    # ------------------------- set default priors
        
    prior_def = {}

    # --- redshift

    prior_def['z'] = {'type': 'uniform', 'limits': [0., 15.]} 

    # --- stellar mass

    prior_def['log10M*'] = {'type': 'uniform', 'limits': [5., 13.]}

    if model_definitions['SFZH'] == 'constant_constant':

         # --- constant star formation and metallicity

        prior_def['log10_duration'] = {'type': 'uniform', 'limits': [7.,10.]}
        prior_def['log10Z'] = {'type': 'uniform', 'limits': [-4.,-1.5]}

#    
#         # --- evolving metallicity
#     
#         self.prior_def['log10Z0'] = {'type': 'uniform', 'limits': [-4.,-1.5]}
#         self.prior_def['dlog10Z/dt'] = {'type': 'uniform', 'limits': [-0.01,0.0]}
#         self.prior_def['sigma_Z'] = {'type': 'uniform', 'limits': [0.,1.0]}


    if model_definitions['dust'] == 'very_simple':    

        prior_def['log10tau_V'] = {'type': 'uniform', 'limits': [-2.,1.]}

    if model_definitions['dust'] == 'simple':    

        prior_def['log10tau_V'] = {'type': 'uniform', 'limits': [-2.,1.]}
        prior_def['dust_slope'] = {'type': 'uniform', 'limits': [-2.,-0.5]}

    # --- nebular emission

    prior_def['log10fesc'] = {'type': 'uniform', 'limits': [-2.,0.]}

    return prior_def


def generate_photometry_grid(sps, imf, filters, redshifts = np.arange(0.,20.,0.01)):

    SPS = models.SPS(sps+'/'+imf)
    lam = SPS.lam

    F = FLARE.filters.add_filters(filters)
    
    outdir = FLARE_dir + '/data/SPS/SED_grids/{0}_{1}'.format(sps, imf)
    if not os.path.exists(outdir): os.mkdir(outdir)
    
    np.save('{0}/redshifts.npy'.format(outdir), redshifts)
    np.save('{0}/log10ages.npy'.format(outdir), SPS.grid['log10age'])
    np.save('{0}/log10Zs.npy'.format(outdir), SPS.grid['log10Z'])
    
    grid = {}
    for f in filters:
        grid[f] = {}
        for grid_type in ['stellar', 'nebular']: 
            grid[f][grid_type] = np.zeros((*SPS.grid[grid_type].shape[0:-1], len(redshifts)))
    
    
    for iz, z in enumerate(redshifts):
    
        lamz = lam * (1.+z)
    
        dlamz = lamz[1:] - lamz[:-1]
        dlamz = np.append(dlamz, dlamz[-1])
    
        igm = IGM.madau(lamz, z)
    
        print(z)
    
        for grid_type in ['stellar', 'nebular']: 
                    
            for f in filters:
            
                # --- map filter on to 
    
                T = np.interp(lamz, F[f].l, F[f].t)

                grid[f][grid_type][:,:,iz] = np.sum(igm * SPS.grid[grid_type] * T.T * dlamz, axis=2) / np.trapz(T, lamz)

    for f in filters: pickle.dump(grid[f], open('{0}/{1}.p'.format(outdir, f),'wb'))





class photo():


    def __init__(self, model_definitions, filters = False):
    
    
        self.model_definitions = model_definitions
        self.cosmo = self.model_definitions['cosmo']
        
        self.filters = filters
        nfilters = len(self.filters)
                
        F = FLARE.filters.add_filters(self.filters)
        self.pivwv = np.array([F[f].pivwv() for f in filters])
        
        # --- read in grids
        
        grid_dir = FLARE_dir + '/data/SPS/SED_grids/{0}_{1}'.format(model_definitions['SPS'], model_definitions['IMF'])

        self.redshifts = np.load('{0}/redshifts.npy'.format(grid_dir))
        self.log10ages = np.load('{0}/log10ages.npy'.format(grid_dir))
        self.log10Zs = np.load('{0}/log10Zs.npy'.format(grid_dir))
        
        self.grid = {}
        for grid_type in ['stellar', 'nebular']:
            grids = [pickle.load(open(FLARE_dir + '/data/SPS/SED_grids/{0}_{1}/{2}.p'.format(model_definitions['SPS'], model_definitions['IMF'], f),'rb'))[grid_type] for f in filters]
            self.grid[grid_type] = np.stack(grids, axis=3)
          
        self.prior_def = default_priors(model_definitions)
        self.parameters = list(self.prior_def.keys())


    def p(self, params):
    
        if self.model_definitions['SFZH'] == 'constant_constant': sfzh, sfr = SFZH.constant(self.log10ages, self.log10Zs, params)  # --- generate SFZH

        sfzh = np.expand_dims(sfzh, axis=2)

        z = params['z']
        params['iz'] = (np.abs(self.redshifts - z)).argmin()

        mod = empty()
        mod.properties = {'SFR10': sfr} 
        
        if self.filters: mod.F = self.F(sfzh, params) # calculate broadband photometry
        
        return mod


    def F(self, sfzh, params): # --- calculate broadband photometry

        z = params['z']
        iz =  params['iz']

        stellar = np.sum(self.grid['stellar'][:,:,iz,:] * sfzh, axis=(0,1))
  
        nebular = (1.-10**params['log10fesc'])*np.sum(self.grid['nebular'][:,:,iz,:] * sfzh, axis=(0,1))
  
        total = stellar + nebular
          
        properties = {}
  
        total *= 1E23 * 1E9 * (1.+z) / (4 * np.pi * self.model_definitions['cosmo'].luminosity_distance(z).to('cm').value**2) # nJy
  
        # -- add dust
  
        dust_model = self.model_definitions['dust']
        
        if dust_model:
                       
            tau_V = 10**params['log10tau_V']
                        
            if dust_model == 'simple':
            
                tau = tau_V*(((self.pivwv/(1.+z))/5500.)**SED_p['slope'])
                T = np.exp(-tau) 
                
            if dust_model == 'very_simple':
            
                tau = tau_V*(((self.pivwv/(1.+z))/5500.)**(-1.))
                T = np.exp(-tau) 
                   
            total *= T

  
        return total





class full():


    def __init__(self, model_definitions, filters = False, lines = False):
    
        self.model_definitions = model_definitions
        
        self.cosmo = self.model_definitions['cosmo']
        
        self.filters = filters
        
        self.SPS = models.SPS(model_definitions['SPS']+'/'+model_definitions['IMF'])

        self.F = FLARE.filters.add_filters(self.filters) # ---filter object # THIS IS HORRIBLY INEFFICIENT 

        self.prior_def = default_priors(model_definitions)

        self.parameters = self.prior_def.keys()



    def p(self, params, include_IGM = True):
    
        if self.model_definitions['SFZH'] == 'constant_constant': sfzh, sfr = SFZH.constant(self.SPS.grid['log10age'], self.SPS.grid['log10Z'], params)  # --- generate SFZH

        params['fesc'] = 10**params['log10fesc']

        mod = empty()

        mod.properties = {'SFR10': sfr} 

        if self.filters:
        
            SED = self.SPS.get_Lnu(sfzh, params, dust_model = self.model_definitions['dust'], fast = True) # --- generate intrinsic SED  ** THIS IS THE SPEED BOTTLENECK
            SED.total.get_fnu(self.model_definitions['cosmo'], params['z'], include_IGM = include_IGM)

            for f in self.filters: self.F[f].apply_new_lam(SED.total.lamz)
            
            SED.total.get_Fnu(self.F) # generates Fnu (broad band fluxes)

            mod.F = SED.total.Fnu_array
            
        return mod
            
        
