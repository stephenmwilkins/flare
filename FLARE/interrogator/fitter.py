


import numpy as np
import scipy.stats
import sys
import os
import emcee
from . import models
import pickle


class empty: pass # --- output class

class delta():
    
    def __init__(self, value):     
        self.value = value
        
    def rvs(self):    
        return self.value
        
    def logpdf(self, v):  
    
        if v == self.value:
            return 0.0
        else:
            return -np.inf

    
    
class source():


    def __init__(self, obs, model, verbose = False):
    
    
        self.verbose = verbose
    
        # --- define model
        
        self.model = model
     
        # --- define observations
     
        self.obs = obs 
     
        # --- define model parameters
     
        self.parameters = self.model.parameters
     
        
        # --- copy over the default priors for the model
        
        self.prior_def = {}
        
        for parameter in self.parameters:
        
            self.prior_def[parameter] = self.model.prior_def[parameter]
            
              
              
    def update_priors(self):
    
        # --- update priors

        self.priors = {}

        if self.verbose: 
            print('----------------------------------------------------------------------')
            print('------------------- Priors -------------------------------------------')

        for parameter in self.parameters:
        
            if self.prior_def[parameter]['type'] == 'uniform':
            
                self.priors[parameter] = scipy.stats.uniform(loc = self.prior_def[parameter]['limits'][0], scale = self.prior_def[parameter]['limits'][1] - self.prior_def[parameter]['limits'][0])  
                
                if self.verbose: print(parameter, self.prior_def[parameter]['type'], self.prior_def[parameter]['limits'])
                
                                   
            if self.prior_def[parameter]['type'] == 'delta':
    
                self.priors[parameter] = delta(value = self.prior_def[parameter]['value'])
                
                if self.verbose: print(parameter, self.prior_def[parameter]['type'], self.prior_def[parameter]['value'])


            if self.prior_def[parameter]['type'] == 'norm': 
            
                self.priors[parameter] = scipy.stats.norm(loc = self.prior_def[parameter]['loc'], scale = self.prior_def[parameter]['scale'])

                if self.verbose: print(parameter, self.prior_def[parameter]['type'], self.prior_def[parameter]['loc'], self.prior_def[parameter]['scale'])




    def lnprob(self, params):
        """Log probability function"""

        p = {parameter:params[i] for i,parameter in enumerate(self.parameters)}
      
        lp = np.sum([self.priors[parameter].logpdf(p[parameter]) for parameter in self.parameters])
  
        # --- check to see if posterior propbability is non-zero
  
        if not np.isfinite(lp): 
            return -np.inf , [-99.]
     
        # --- age of Universe constraint
        
        if p['log10_duration']>np.log10(self.model.cosmo.age(p['z']).to('yr').value):
            return -np.inf , [-99.]
            

        # --- generate model observations

        mod = self.model.p(p)
        
        if np.isnan(mod.F).any():
            return -np.inf , [-99.]
        else:  
            lnlike = -0.5*np.sum(((self.obs.fluxes - mod.F)/self.obs.flux_errors)**2 - np.log(self.obs.flux_errors**2))
            return lp + lnlike, [mod.properties['SFR10']]



    def fit(self, nwalkers = 50, nsamples = 1000, burn = 200):
    
        ndim = len(self.parameters)
    
        self.update_priors()
    
        self.ndim = ndim
        self.nwalkers = nwalkers        
        self.nsamples = nsamples
    
        p0 = [ [self.priors[parameter].rvs() for parameter in self.parameters] for i in range(nwalkers)]
        
        self.sampler = emcee.EnsembleSampler(nwalkers, ndim, self.lnprob, args=())
                
        pos, prob, state, blobs = self.sampler.run_mcmc(p0, burn)
          
        self.sampler.reset()
        
        self.sampler.run_mcmc(pos, nsamples)

            
        
    def save(self):
    
        outdir =  'outputs/' + self.obs.ID
    
        if not os.path.exists(outdir): os.system('mkdir '+outdir)  # --- make directory for each object

        
        summary = empty()
    
        summary.model_definitions = self.model.model_definitions # --- 

        summary.model_parameters = self.parameters 
        
        summary.derived_parameters = ['SFR10']  # --- list of derived parameters
                
        summary.prior_def = self.prior_def # --- prior definitions
        
        summary.obs = self.obs  # --- input observations 
        
    
        summary.ndim = self.ndim # --- number of parameters
        summary.nwalkers = self.nwalkers # --- number of walkers
        summary.nsamples = self.nsamples # --- number of samples
        
        
        
        # --- reformat everything just just have a dictionary of samples
        
        samples = {}
        
        chains = self.sampler.chain[:, :, ].reshape((-1, self.ndim))
        
        for ip, p in enumerate(self.parameters): samples[p] = chains[:,ip] 
        
        # --- read blobs
        
        blob_list = ['SFR10']
        
        b = {}     
        for i, p in enumerate(blob_list):
            b[p] = []    
            for k in np.arange(self.nwalkers):
                for j in np.arange(self.nsamples):      
                    b[p].append(self.sampler.blobs[j][k][i])   
            b[p] = np.array(b[p])
        
            samples[p] = b[p]
        
        del b
        
        pickle.dump(samples, open(outdir+'/samples.p','wb'))
        
        
        # ---- determine percentiles
        
        summary.percentiles = {}

        if self.verbose: print('----------------------------------------------------------------------')
        if self.verbose: print('------------------- Results ------------------------------------------')

        for i, p in enumerate(summary.model_parameters + summary.derived_parameters): 
        
            summary.percentiles[p] = np.array([np.percentile(samples[p], x) for x in range(0,101)])
        
            if self.verbose: print('{0} {1:.2f} {2:.2f} {3:.2f}'.format(p, summary.percentiles[p][16], summary.percentiles[p][50], summary.percentiles[p][84]))
                     
        pickle.dump(summary, open(outdir+'/summary.p','wb'))
        
     
        
        
       
