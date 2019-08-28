
from .core import *

class empty: pass

import numpy as np


# --- artificial

FAKE = ['FAKE.FAKE.'+f for f in ['1500','2500','Uth','Bth','Vth','Ith','Zth','Yth','Jth','Hth']] 
TH = ['FAKE.TH.'+f for f in ['FUV','MUV', 'NUV','U','B','V','R','I','Z','Y','J','H','K']] 

# --- Hubble

WFC3UVIS_W = ['HST.WFC3.'+f for f in ['f225w','f275w','f336w']]
WFC3NIR_W = ['HST.WFC3.'+f for f in ['f105w', 'f125w', 'f140w', 'f160w']]
WFC3 = WFC3UVIS_W + WFC3NIR_W

ACS_W = ['HST.ACS.'+f for f in ['f435w', 'f606w', 'f775w', 'f814w', 'f850lp']]
ACS = ACS_W

HST = ACS + WFC3
Hubble = HST

# --- Spitzer

# IRAC = ['Spitzer.IRAC.'+f for f in ['ch1', 'ch2', 'ch3', 'ch4']]
IRAC = ['Spitzer.IRAC.'+f for f in ['ch1', 'ch2']]
Spitzer = IRAC

# --- Euclid

Euclid_VIS = ['Euclid.VIS.'+f for f in ['VIS']]
Euclid_NISP = ['Euclid.NISP.'+f for f in ['Y','J','H']]
Euclid = Euclid_VIS + Euclid_NISP

# --- Subaru

HSC = ['Subaru.HSC.'+f for f in ['g','r','i','z','y']]
Subaru = HSC

# --- Webb

NIRCam_s_W = ['JWST.NIRCAM.'+f for f in ['F070W','F090W','F115W','F150W','F200W']]
NIRCam_s_M = ['JWST.NIRCAM.'+f for f in ['F140M','F162M','F182M','F210M']]
NIRCam_l_W = ['JWST.NIRCAM.'+f for f in ['F277W','F356W','F444W']]
NIRCam_l_M = ['JWST.NIRCAM.'+f for f in ['F250M','F300M','F360M','F410M','F430M','F460M','F480M']]
NIRCam_s = NIRCam_s_W + NIRCam_s_M
NIRCam_l = NIRCam_l_W + NIRCam_l_M 
NIRCam_W = NIRCam_s_W + NIRCam_l_W
NIRCam =  NIRCam_s + NIRCam_l
MIRI = ['JWST.MIRI.'+f for f in ['F560W','F770W','F1000W','F1130W','F1280W','F1500W','F1800W','F2100W','F2550W']]
Webb = NIRCam + MIRI


# --- All *real* filters

all_filters = FAKE + TH + Hubble + Spitzer + Euclid + Subaru + Webb 


# --- filter info

info = {filter:empty() for filter in all_filters}

for filter in all_filters: info[filter].zeropoint = None # photometric zeropoint (e/s -> erg/s/Hz)
for filter in all_filters: info[filter].pixel_scale = None # arcsec
for filter in all_filters: info[filter].FOV = None # field of view of instrument filter is on

info['HST.ACS.f435w'].zeropoint = 25.684
info['HST.ACS.f606w'].zeropoint = 26.505
info['HST.ACS.f775w'].zeropoint = 25.678
info['HST.ACS.f814w'].zeropoint = 25.959
info['HST.ACS.f850lp'].zeropoint = 24.867
info['HST.WFC3.f105w'].zeropoint = 26.269
info['HST.WFC3.f125w'].zeropoint = 26.230
info['HST.WFC3.f140w'].zeropoint = 26.452
info['HST.WFC3.f160w'].zeropoint = 25.946


for filter in NIRCam: info[filter].zeropoint = 26.0 # need to think about this

for filter in IRAC: info[filter].pixel_scale = 1.22
for filter in Euclid_NISP: info[filter].pixel_scale = 0.3
for filter in ACS: info[filter].pixel_scale = 0.05 
for filter in WFC3: info[filter].pixel_scale = 0.13
for filter in NIRCam_s: info[filter].pixel_scale = 0.031
for filter in NIRCam_l: info[filter].pixel_scale = 0.063
for filter in MIRI: info[filter].pixel_scale = 0.11
# for filter in : info[filter].pixel_scale = 


# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------

# below is deprecated in favour of the above

pixel_scale = {}

# --- Spizer

pixel_scale.update({f: 1.22 for f in IRAC}) 

# --- Euclid

pixel_scale.update({f: 0.3 for f in Euclid_NISP}) 

# --- Hubble

pixel_scale.update({f: 0.05 for f in ACS_W}) # this is wrong
pixel_scale.update({f: 0.13 for f in WFC3UVIS_W}) # this is wrong
pixel_scale.update({f: 0.13 for f in WFC3NIR_W})

# --- Webb

pixel_scale.update({f: 0.031 for f in NIRCam_s}) 
pixel_scale.update({f: 0.063 for f in NIRCam_l}) 
pixel_scale.update({f: 0.11 for f in MIRI}) 



zeropoints = {}
zeropoints['HST.ACS.f435w'] = 25.684
zeropoints['HST.ACS.f606w'] = 26.505
zeropoints['HST.ACS.f775w'] = 25.678
zeropoints['HST.ACS.f814w'] = 25.959
zeropoints['HST.ACS.f850lp'] = 24.867

zeropoints['HST.WFC3.f105w'] = 26.269
zeropoints['HST.WFC3.f125w'] = 26.230
zeropoints['HST.WFC3.f140w'] = 26.452
zeropoints['HST.WFC3.f160w'] = 25.946










def add_filters(filters, new_lam = False,  filter_path = FLARE_dir + '/data/filters/'):

    F = {f: filter(f, new_lam, filter_path) for f in filters}
    
    F['filters'] = filters
    
    return F
    
    
class filter():

    def __init__(self, f, new_lam = False, filter_path = FLARE_dir + '/data/filters/'):

        # l, t are the original wavelength and transmission grid

        self.l, self.t = np.loadtxt(filter_path + '/'+'/'.join(f.split('.'))+'.txt', skiprows = 1).T 
        
        if f.split('.')[1] == 'NIRCAM': self.l *= 1E4 # convert from microns to \AA
        
        if isinstance(new_lam, np.ndarray):
        
            self.lam = new_lam
            self.T = np.interp(self.lam, self.l, self.t)
        
        else:
        
            self.lam = self.l
            self.T = self.t
    
        self.info = info[f]
    
    
    def apply_new_lam(self, new_lam):
    
        self.lam = new_lam
        self.T = np.interp(self.lam, self.l, self.t)
    

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
        
        


def create_EAZY_filter_res(F, filter_res_file = 'FILTER.RES'):
    
    o = []

    for f in F['filters']:
        
        o.append('{n} {f}\n'.format(n=len(F[f].l), f=f)) # filter header
    
        for i,l in enumerate(F[f].l):  
            o.append('{i:>5}   {l}   {t}\n'.format(i=i+1, l=l,t=F[f].t[i]))
        
    open(filter_res_file,'w').writelines(o)
        






