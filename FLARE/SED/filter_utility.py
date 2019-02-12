

import numpy as np

def add_filters(filters, new_lam = False,  filter_path = '.'):

    F = {f: filter(f, new_lam, filter_path) for f in filters}
    
    return F
    
    

class filter():

    def __init__(self, f, new_lam = False, filter_path = '.'):

        self.l, self.t = np.loadtxt(filter_path + '/'+'/'.join(f.split('.'))+'.txt').T 
        
        if isinstance(new_lam, np.ndarray):
        
            self.lam = new_lam
            self.T = np.interp(self.lam, self.l, self.t)
        
        else:
        
            self.lam = l
            self.T = t
            
    def pivwv(self): # -- calculate pivot wavelength using original l, t
    
        return np.sqrt(np.trapz(self.l * self.t , x = self.l)/np.trapz(self.t / self.l, x = self.l))  
        
      


