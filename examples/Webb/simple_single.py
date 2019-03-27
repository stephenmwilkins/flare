

# NOTE: this doesn't appear to be working. It does not reproduce numbers from Pandeia. sigh.


import numpy as np
from jwst_backgrounds import jbt

from photutils import source_properties
from photutils import CircularAperture
from photutils import aperture_photometry

from astropy.constants import h

from astropy import units as u
 
import matplotlib.pyplot as plt
 
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))


import FLARE.filters 
import SynthObs.Morph

filters = FLARE.filters.NIRCam_W
filters = ['JWST.NIRCAM.F070W']


F = FLARE.filters.add_filters(filters)

PSF = SynthObs.Morph.webbPSFs(filters, 4.) # creates a dictionary of instances of the webbPSF class

wavelengths = np.array([F[f].pivwv()/1E4 for f in F['filters']])


# --- hopefully this is the same as the benchmark background
 
ra = 261.6833333
dec = -73.3322222
day = 170

# calculate background

bg = jbt.background(ra,dec, wavelengths)


# identify index of day

day_i = np.argwhere(bg.bkg_data['calendar'] == day)[0][0]
   
   
# --- telescope area - precise?
area = 25.*100**2 # cm2 

pixel_s = 0.031 
pixel_l = 0.063 


source_flux_original = 100 # nJy

t_exp = 9974.46 # s


for i,f in enumerate(filters):
    
    
    if f in FLARE.filters.NIRCam_s: pixel_scale = pixel_s
    if f in FLARE.filters.NIRCam_l: pixel_scale = pixel_l
        
    pixel_area = pixel_scale**2 * u.arcsec * u.arcsec 
      
    Tint = np.trapz(F[f].t/(F[f].l*h.value), x=F[f].l)
    
    
    
    # --- aperture
    
    aperture_radius_pix = 2.5 
    aperture_area_pix = np.pi*aperture_radius_pix**2 # in pixels
    
    # --- source
     
    positions = [(PSF[f].PSF.shape[0]/2., PSF[f].PSF.shape[0]/2.)] # centre
    apertures = [CircularAperture(positions, r=r) for r in [aperture_radius_pix]] #r in pixels
    
    phot_table = aperture_photometry(PSF[f].PSF, apertures) 
    frac = phot_table[0][3]/np.sum(PSF[f].PSF) # fraction of source's light in aperture
    print('fraction of source flux in aperture:', frac)
    
    
    source_flux = source_flux_original * frac # nJy
    source_flux /= 1E9 # Jy
    source_flux *= 1E-30 # J/s/cm2/Hz
    
    source_e_s =  source_flux*Tint*area # e-/s in each pixel
    
    print('Total flux:', source_e_s / frac, 'e-/s')
    print('Extracted flux:', source_e_s, 'e-/s')
    
    source_e = source_e_s * t_exp # e-
    
    # --- background
    
    # background_flux = bg.bathtub['total_thiswave'][day_i, i] * 1E6 # Jy/sr 
    
    background_flux = 0.21 * 1E6
    
    print('background flux/MJy:', background_flux/1E6)
    
    background_flux *= 1E-30 # J/s/cm2/Hz 
    
    background_flux_pix = background_flux * pixel_area.to('sr').value # J/s/cm2/Hz per pixel    

    background_e_s_pix = background_flux_pix * Tint * area   # e-/s

    background_e_pix = background_e_s_pix * t_exp  # e-
    
    print('total background in each pixel (e-):', background_e_pix)
    
    background_e_aperture = background_e_pix * aperture_area_pix
    
    print('total background in extraction aperture (e-/s):', background_e_aperture/t_exp)
    
    background_e_aperture_STD = np.sqrt(background_e_aperture)
    
    background_e_s_aperture_STD = background_e_aperture_STD / t_exp
    
    SNR = source_e_s/background_e_s_aperture_STD
    
    print('SNR',SNR)
    
    


