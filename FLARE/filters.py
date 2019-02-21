

import numpy as np




WFC3UV_W = ['HST.WFC3.'+f for f in ['f225w','f275w','f336w']]
ACS_W = ['HST.ACS.'+f for f in ['f814w', 'f606w', 'f775w', 'f814w']]
WFC3NIR_W = ['HST.WFC3.'+f for f in ['f105w', 'f125w', 'f140w', 'f160w']]

HST = WFC3UV_W + ACS_W + WFC3NIR_W

NIRCam_W = ['JWST.NIRCAM.'+f for f in ['F070W','F090W','F115W','F150W','F200W','F277W','F356W','F444W']]
NIRCam_M = ['JWST.NIRCAM.'+f for f in ['F140M','F162M','F182M','F210M','F250M','F300M','F360M','F410M','F430M','F460M','F480M']]

 


def add_filters(filters, new_lam = False,  filter_path = '.'):

    F = {f: filter(f, new_lam, filter_path) for f in filters}
    
    F['filters'] = filters
    
    return F
    
    
class filter():

    def __init__(self, f, new_lam = False, filter_path = '.'):

        self.l, self.t = np.loadtxt(filter_path + '/'+'/'.join(f.split('.'))+'.txt', skiprows = 1).T 
        
        if f.split('.')[0] == 'JWST': self.l *= 1E4 # convert from microns to \AA
        
        
        if isinstance(new_lam, np.ndarray):
        
            self.lam = new_lam
            self.T = np.interp(self.lam, self.l, self.t)
        
        else:
        
            self.lam = self.l
            self.T = self.t
    
    
    #http://stsdas.stsci.edu/stsci_python_epydoc/SynphotManual.pdf   #got to page 42 (5.1)
            
    def pivwv(self): # -- calculate pivot wavelength using original l, t
    
        return np.sqrt(np.trapz(self.l * self.t , x = self.l)/np.trapz(self.t / self.l, x = self.l))  
    
    
    def pivT(self): # -- calculate pivot wavelength using original l, t
    
        return np.interp(self.pivwv(), self.l, self.t)
    
    def meanwv(self):
    
        return np.exp(np.trapz(np.log(self.l) * self.t / self.l, x = self.l)/np.trapz(self.t / self.l, x = self.l))
    
    def bandw(self): 
    
        return self.meanwv() * np.sqrt(np.trapz((np.log(self.l/self.meanwv())**2)*self.t/self.l,x=self.l))/np.sqrt(np.trapz(self.t/self.l,x=self.l))

    def fwhm(self):
    
        return np.sqrt(8.*np.log(2))*self.bandw()
        
    def Tpeak(self):
    
        return np.max(self.t)
        
    def rectw(self):
    
        return np.trapz(self.t, x=self.l)/self.Tpeak()