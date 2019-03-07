

import numpy as np
import pickle

from . import core
from ..core import * 

import copy

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

    def __init__(self, grid, path_to_SPS_grid = FLARE_dir + 'data/SPS/nebular/1.0/Z/', dust = False):
    
        self.grid = pickle.load(open(path_to_SPS_grid + grid + '/nebular.p','rb'), encoding='latin1')

        self.lam = self.grid['lam']

        self.dust = dust
      
  
    def get_Lnu(self, SFZH):
    
        sed = SED()
        
        
        SFZH = np.expand_dims(SFZH, axis=2)
         
        sed.stellar = core.sed(self.lam)
        sed.stellar.lnu = np.sum(self.grid['stellar'] * SFZH, axis=(0,1))
        
        sed.nebular = core.sed(self.lam)
        sed.nebular.lnu = np.sum(self.grid['nebular'] * SFZH, axis=(0,1))
        
        sed.total = core.sed(self.lam)
        sed.total.lnu = np.sum((self.grid['stellar']+self.grid['nebular']) * SFZH, axis=(0,1))

        if self.dust:
        
            sed.stellar_intrinsic = copy.deepcopy(sed.stellar)
            sed.nebular_intrinsic = copy.deepcopy(sed.nebular)
            sed.total_intrinsic = copy.deepcopy(sed.total)
                        
            if self.dust['model'] == 'simple':
            
                tau = self.dust['tau_V']*(self.lam/5500.)**self.dust['slope']
                T = np.exp(-tau) 
        
            sed.stellar.lnu *= T
            sed.nebular.lnu *= T
            sed.total.lnu *= T
                
        return sed
        
        
        
        
        