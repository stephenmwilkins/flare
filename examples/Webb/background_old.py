import numpy as np
from jwst_backgrounds import jbt

from astropy.constants import h

from astropy import units as u
 
import matplotlib.pyplot as plt
 
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))


import FLARE.filters 

filters = FLARE.filters.NIRCam



F = FLARE.filters.add_filters(filters)

wavelengths = np.array([F[f].pivwv()/1E4 for f in F['filters']])

print(wavelengths)
 
# ra = 53.1625
# dec = -27.7914
 

ra = 261.6833333
dec = -73.3322222

day = 170

# jbt.get_background(ra, dec, 4.4, plot_background=True, plot_bathtub=True, write_bathtub=True)



bg = jbt.background(ra,dec, wavelengths)

day_i = np.argwhere(bg.bkg_data['calendar'] == day)[0][0]
print(day_i)

print(dir(bg))
print(bg.bathtub.keys())
print(bg.bkg_data.keys())
print(bg.bkg_data['calendar'])
print(bg.bkg_data['total_bg'].shape)
print(bg.bathtub['total_thiswave'].shape)
   
area = 25.*100**2 # cm2
pixel_s = 0.03 * 0.03 * u.arcsec * u.arcsec # size of pixel in sq. arcsec
pixel_l = 0.06 * 0.06 * u.arcsec * u.arcsec # size of pixel in sq. arcsec
# aperture = np.pi*0.1**2/(0.03**2) # using 0.1" radius aperture to match ETC default
aperture = np.pi*2.5**2 # in pixels



t_exp = 1E4


for i,f in enumerate(filters):
    
    # print(f, np.min(bg.bathtub['total_thiswave'][:, i]))
    
    # B = np.min(bg.bathtub['total_thiswave'][:, i]) * 1E6 * u.Jy / u.sr # Jy/sr # MIN
    B = bg.bathtub['total_thiswave'][day_i, i] * 1E6 # Jy/sr 
    
    if f in FLARE.filters.NIRCam_s: pixel = pixel_s
    if f in FLARE.filters.NIRCam_l: pixel = pixel_l
    
    B_pix = B*pixel.to('sr').value # Jy per pixel
    
    B_pix *= 1E-30 # J/s/cm2/Hz per pixel
     
    
    Tint = np.trapz(F[f].t/(F[f].l*h.value), x=F[f].l)
    
    e = 1.2*B_pix*Tint # e-/s/cm2 in each pixel

    e_bkg_app = e*area*aperture*t_exp # e-
    
    # print(e_bkg_app)
    
    e_bkg_app_std = np.sqrt(e_bkg_app) # e-  

    flux_limit = e_bkg_app_std/(area*t_exp) # e-/s/cm2

    flux_limit /= Tint #J/s/cm2/Hz
    
    flux_limit /= 1E-30 # Jy
    
    flux_limit *= 1E9 # nJy
    
    flux_limit *= 10 # 10\sigma
    
    print(f, flux_limit)
   
   



   
# plt(bg.bathtub(
# 
# print(bg)