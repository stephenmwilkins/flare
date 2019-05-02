

# NOTE: this doesn't appear to be working. It does not reproduce numbers from Pandeia. sigh.


import numpy as np


from astropy.constants import h

from astropy import units as u
 
import matplotlib.pyplot as plt
 
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))


import FLARE.filters 
# import SynthObs.Morph

filters = FLARE.filters.NIRCam
# filters = ['JWST.NIRCAM.F200W']


detailed_output = False

F = FLARE.filters.add_filters(filters)

# wavelengths = np.array([F[f].pivwv()/1E4 for f in F['filters']])
wavelengths = np.arange(0.5, 5.5, 0.01)

from jwst_backgrounds import jbt

# --- hopefully this is the same as the benchmark background (it is not)
 
ra = 261.6833333
dec = -73.3322222
day = 170

bg = jbt.background(ra,dec, wavelengths) # calculate background
day_i = np.argwhere(bg.bkg_data['calendar'] == day)[0][0] # identify index of day
   
   
# --- telescope area - precise?
area = 25.4*100**2 # cm2 

pixel_s = 0.031 
pixel_l = 0.063 

t_exp = 9974.46 # s


pandeia = {'JWST.NIRCAM.F070W': 1126.0721474889592, 'JWST.NIRCAM.F090W': 1786.7394403286296, 'JWST.NIRCAM.F115W': 1958.3717674718248, 'JWST.NIRCAM.F150W': 2094.999837765901, 'JWST.NIRCAM.F200W': 2061.5719837568827, 'JWST.NIRCAM.F140M': 981.5902707204624, 'JWST.NIRCAM.F162M': 942.9066381746721, 'JWST.NIRCAM.F182M': 1185.001923425075, 'JWST.NIRCAM.F210M': 834.308953687194, 'JWST.NIRCAM.F277W': 5783.452503476042, 'JWST.NIRCAM.F356W': 5598.499153909868, 'JWST.NIRCAM.F444W': 10433.88191793632, 'JWST.NIRCAM.F250M': 1812.3164234617484, 'JWST.NIRCAM.F300M': 2043.4940798266255, 'JWST.NIRCAM.F360M': 2468.310748196025, 'JWST.NIRCAM.F410M': 3277.565422796658, 'JWST.NIRCAM.F430M': 1862.4507469018095, 'JWST.NIRCAM.F460M': 2187.9391545638887, 'JWST.NIRCAM.F480M': 3220.657767138416}

background = {}

for i,f in enumerate(filters):
    
    if detailed_output: print('-'*5, f)
    
    if f in FLARE.filters.NIRCam_s: pixel_scale = pixel_s
    if f in FLARE.filters.NIRCam_l: pixel_scale = pixel_l
        
    pixel_area = pixel_scale**2 * u.arcsec * u.arcsec 
      
    # Tint = np.trapz(F[f].t/(F[f].l*h.value), x=F[f].l)
    
    # --- aperture
    
    aperture_radius_pix = 2.5 # in pixels
    aperture_area_pix = np.pi*aperture_radius_pix**2 # in pixels
     
    # --- background
    
    # background_flux = 0.21 * 1E6 # Jy/sr
    
    # background_flux = bg.bathtub['total_thiswave'][day_i, i] * 1E6 # Jy/sr 
    
    background_flux = bg.bathtub['total_thiswave'][day_i, :] * 1E6 # Jy/sr 
    
    if detailed_output: print('background flux/MJy:', background_flux/1E6)
    
    background_flux *= 1E-30 # J/s/cm2/Hz 
    
    background_flux_pix = background_flux * pixel_area.to('sr').value # J/s/cm2/Hz per pixel

    # background_e_s_pix = background_flux_pix * Tint * area   # e-/s
    
    background_interp = np.interp(F[f].l, wavelengths*1E4, background_flux_pix)
    background_e_s_pix = np.trapz(background_interp*F[f].t/(F[f].l*h.value), x=F[f].l) * area   # e-/s


    background_e_pix = background_e_s_pix * t_exp  # e-
    
    if detailed_output: print('total background in each pixel (e-):', background_e_pix)
    
    background_e_aperture = background_e_pix * aperture_area_pix
    
    if detailed_output: print('total background flux in extraction aperture (e-/s):', background_e_s_pix * aperture_area_pix)
    
    background_e_aperture_STD = np.sqrt(background_e_aperture)
    
    background_e_s_aperture_STD = background_e_aperture_STD / t_exp
    
    if not detailed_output: print(f, background_e_pix, background_e_pix/pandeia[f])
    
    background[f] = background_e_pix


