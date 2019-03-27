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
filters = ['JWST.NIRCAM.F200W']


F = FLARE.filters.add_filters(filters)


PSF = SynthObs.Morph.webbPSFs(filters, 1.) # creates a dictionary of instances of the webbPSF class

# print(PSF[filters[0]].PSF.shape)
# print(np.sum(PSF[filters[0]].PSF))




#     
# self.profile = {'r_pix': radii, 'r_kpc': radii*self.img.resolution, 'I': np.array([float(phot_table[0][j+3]) for j in range(len(radii))])}
#         

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
   
area = 25.*100**2 # cm2
pixel_s = 0.031 
pixel_l = 0.063 
# aperture = np.pi*0.1**2/(0.03**2) # using 0.1" radius aperture to match ETC default
aperture = np.pi*2.5**2 # in pixels



t_exp = 1E4


for i,f in enumerate(filters):
    
    
    
    positions = [(PSF[f].PSF.shape[0]/2., PSF[f].PSF.shape[0]/2.)] # centre
    apertures = [CircularAperture(positions, r=r) for r in [2.5]] #r in pixels
    
    phot_table = aperture_photometry(PSF[f].PSF, apertures) 
    frac = phot_table[0][3]

    
    # print(f, np.min(bg.bathtub['total_thiswave'][:, i]))
    
    
    # B = 1.2*bg.bathtub['zodi_thiswave'][day_i, i] * 1E6 # Jy/sr 
    B = 1.2*np.min(bg.bathtub['total_thiswave'][:, i]) * 1E6 # * u.Jy / u.sr # Jy/sr # MIN
    
    
    if f in FLARE.filters.NIRCam_s: pixel = pixel_s
    if f in FLARE.filters.NIRCam_l: pixel = pixel_l
    
    
    pixel_area = pixel**2 * u.arcsec * u.arcsec 
    
    
    B_pix = B*pixel_area.to('sr').value # Jy per pixel
    
    B_pix *= 1E-30 # J/s/cm2/Hz per pixel
     
    
    Tint = np.trapz(F[f].t/(F[f].l*h.value/1E4), x=F[f].l/1E4)
    
    e = B_pix*Tint # e-/s/cm2 in each pixel

    e_bkg_app = e*area*aperture*t_exp # e- in the aperture
    
    # print(e_bkg_app)
    
    e_bkg_app_std = np.sqrt(e_bkg_app) # e-  

    flux_limit = e_bkg_app_std/(area*t_exp) # e-/s/cm2

    flux_limit /= Tint #J/s/cm2/Hz
    
    flux_limit /= 1E-30 # Jy
    
    flux_limit *= 1E9 # nJy
    
    flux_limit /= frac # accounts for fraction of object flux expected in the aperture
    
    flux_limit *= 10 # 10\sigma
    
    print(f, B, frac, np.sum(PSF[f].PSF), flux_limit)
   
   



   
# plt(bg.bathtub(
# 
# print(bg)