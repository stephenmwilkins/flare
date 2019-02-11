

import numpy as np




# --- filter helper






def add_filters(filters, new_lam = False,  data_dir = 'data/filters'):

    F = {f: filter(f, new_lam, data_dir) for f in filters}
    
    return F
    
    

class filter():

    def __init__(self, f, new_lam = False, data_dir = 'data/filters'):

        self.l, self.t = np.loadtxt(data_dir + '/'+'/'.join(f.split('.'))+'.txt').T 
        
        if isinstance(new_lam, np.ndarray):
        
            self.lam = new_lam
            self.T = np.interp(self.lam, self.l, self.t)
        
        else:
        
            self.lam = l
            self.T = t
            
    def pivwv(self): # -- calculate pivot wavelength using original l, t
    
        return np.sqrt(np.trapz(self.l * self.t , x = self.l)/np.trapz(self.t / self.l, x = self.l))  
        
      


