
import numpy as np
from astropy.io import fits
from photutils import CircularAperture
from photutils import aperture_photometry


import FLARE
import FLARE.filters






class empty: pass


class image:

    def __init__(self, filename, filter, mask = False, convert_to_nJy = True, pixel_scale = 0.06, verbose = False, sci_suffix = 'sci', wht_suffix = 'wht'):

        """read in sci/wht image and appy mask if specified"""
        
        self.filter = filter
        self.pixel_scale = pixel_scale
        self.mask = mask
        self.sci = fits.getdata('{0}/{1}_{2}.fits'.format(FLARE.FLARE_dir, filename, sci_suffix))
        self.wht = fits.getdata('{0}/{1}_{2}.fits'.format(FLARE.FLARE_dir, filename, wht_suffix))
        
        self.zeropoint = FLARE.filters.zeropoints[filter] # AB magnitude zeropoint
        self.nJy_to_es = 1E-9 * 10**(0.4*(self.zeropoint-8.9)) # conversion from nJy to e/s 
        
        if type(mask) == np.ndarray:
            self.masked = True
        else:
            self.masked = False
        
        if self.masked:
            self.sci = np.ma.masked_array(self.sci, mask = self.mask)
            self.wht = np.ma.masked_array(self.wht, mask = self.mask)
            
        if convert_to_nJy:
            self.sci /= self.nJy_to_es
            self.wht /= self.nJy_to_es
            self.unit = 'nJy'
        else:
            self.unit = 'e/s'
            
    def get_random_locations(self, N):
    
        """get N random locations on the image"""
    
        if self.masked:
            pos = np.random.choice(self.sci.count(), size=N)
            return np.take((~self.sci.mask).nonzero(), pos, axis=1)
        
        
        
    def make_cutout(self, x, y, width):
    
        """make """
    
        cutout = empty()

        return cutout
    
    
    
    def determine_depth(self, N = 10000, aperture_diameter_arcsec = 0.35, sigma = 5.):
    
        """determine depth using random apertures"""
    
        aperture_centres = tuple(self.get_random_locations(N).T)
        apertures = [CircularAperture(aperture_centres, r=r) for r in [(aperture_diameter_arcsec/self.pixel_scale)/2.]] # r in pixels
        phot_table = aperture_photometry(self.sci, apertures) 
        aperture_fluxes = phot_table['aperture_sum_0'].quantity
        negative_aperture_fluxes = aperture_fluxes[aperture_fluxes<0]
        return -np.percentile(negative_aperture_fluxes, 100.-68.3) * sigma

    
    
    
       
    
        
# 
# def create_stack():

