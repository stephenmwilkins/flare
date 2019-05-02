



import numpy as np
from . import photom

# https://euclid.roe.ac.uk/attachments/download/1698/puwg-high-z-counts.pdf

area_deg2 = 40.
area = 60.*60.*area_deg2 # arcmin2
Euclid_area_sky = area_deg2/(360**2/np.pi) # fraction of sky

print('area of Euclid deep: {0} arcmin2'.format(area))

m_limit_deep = 25.3

f_limit_deep = photom.m_to_flux(m_limit_deep) 

f_limit_deep *= 1.0 # taking into account light outside the aperture




noise_in_aperture = f_limit_deep/10.

background_in_aperture = noise_in_aperture**2
aperture_radius = 3.0
area_of_aperture = np.pi * aperture_radius**2
background_in_pixel = background_in_aperture/area_of_aperture
noise_in_pixel = np.sqrt(background_in_pixel)

print('assumed aperture radius: {0:.2f} pix'.format(aperture_radius))
print('noise in aperture: {0:.2f} nJy'.format(noise_in_aperture))
print('noise in pixel: {0:.2f} nJy'.format(noise_in_pixel))

filters = ['Euclid.NISP.'+f for f in ['Y','J','H']]



