# -------- make Figure

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm

# plt.style.use('simple') # --- makes nicer plots 

class empty(): pass


from photutils import detect_sources
from photutils import source_properties

from photutils import CircularAperture
from photutils import EllipticalAperture
from photutils import aperture_photometry

from photutils import data_properties

import FLARE.filters



  

def measure_properties(idx, cat, detection, segm, verbose = False):
    
    ndim = detection.sci.shape[0]
    
    source = empty()
    
    source.detected = True
    
    source.mask = segm.data==idx+1 # mask just including the object (used for ISO)
                  
    source.other_objects_mask = np.isin(segm.data, [idx+1, 0], invert = True) # produce a mask excluding all other objects 
                  
    x, y = cat['xcentroid'][idx].value, cat['ycentroid'][idx].value
    
    source.x, source.y = x, y
    
    source.area = len(segm.data[segm.data==idx+1])
    source.radius = np.sqrt(source.area/np.pi)
    
    if verbose: print('detection area/pix: {0}'.format(source.area))
    if verbose: print('detection radius/pix: {0:.2f}'.format(source.radius))
    
    source.A = cat['semimajor_axis_sigma'][idx].value
    source.B = cat['semiminor_axis_sigma'][idx].value
    source.theta = cat['orientation'][idx].value
    source.ellipticity = cat['ellipticity'][idx].value
    
    if verbose: print('A: {0:.2f} B:{1:.2f} theta:{2:.2f}'.format(source.A,source.B,source.theta))

    # --- determine Kron radius
    
    X = np.linspace(0, ndim -1, ndim) - x
    Y = np.linspace(0, ndim -1, ndim) - y
    
    XX, YY = np.meshgrid(X,Y)
    
    R = np.sqrt(XX**2 + YY**2)
    
    aper = EllipticalAperture((x,y), 6.0*source.A, 6.0*source.B, source.theta)
    phot_table = aperture_photometry(detection.sci, aper, mask = source.other_objects_mask) 
    f2 = phot_table['aperture_sum'][0]
    phot_table = aperture_photometry(detection.sci * R, aper, mask = source.other_objects_mask) 
    f1 = phot_table['aperture_sum'][0]
    source.kron_radius = np.max([1., f1/f2])
    
    
    if verbose: print('ISO flux: {0:.2f}'.format(np.sum(detection.sci[source.mask])))  
    if verbose: print('flux at 6\sigma: {0:.2f}'.format(f2))
    if verbose: print('Kron radius: {0:.2f}'.format(source.kron_radius))

#     source.photometry = {}
#     source.intrinsic = {}
#     
#     if verbose: 
#         print()
#         print('-'*5, 'detection')
#     source.intrinsic['detection'] = measure_intrinsic(source, img, verbose)
#     source.photometry['detection'] = measure_photometry(source, img, verbose)
    
    
    return source 


 

def measure_intrinsic(s, img, verbose = False):           
    
    o = empty()
    
    data = img.no_PSF
    
    # --- determine centroid
    
    cat = data_properties(data)

    x, y = cat.xcentroid.value, cat.ycentroid.value

    o.xcentroid = x
    o.xcentroid = y
    o.semimajor_axis_sigma = cat.semimajor_axis_sigma.value 
    o.semimajor_axis_sigma = cat.semimajor_axis_sigma.value 
    
    o.a = cat.semimajor_axis_sigma.value  * img.pixel_scale_kpc
    o.b = cat.semiminor_axis_sigma.value  * img.pixel_scale_kpc
    o.ellipticity = cat.ellipticity.value
    o.orientation = cat.orientation.value

    if verbose:
        print('    '+'-'*5, 'intrinsic (no PSF, no noise)')
        print('    x,y: {0:.2f} {1:.2f}'.format(x,y)) # only equal to r_e_major for n=1 
        print('    a: {0:.2f}'.format(cat.semimajor_axis_sigma.value * img.pixel_scale_kpc)) # only equal to r_e_major for n=1 
        print('    b: {0:.2f}'.format(cat.semiminor_axis_sigma.value * img.pixel_scale_kpc)) # only equal to r_e_minor for n=1 
        print('    e: {0:.2f}'.format(cat.ellipticity.value))
        print('    theta: {0:.2f}'.format(cat.orientation.value /np.pi))


    # --- measure COG

    o.apertures = empty()
    o.apertures.radii_pix = np.arange(0.5, 20, 0.5)        

    apertures = [CircularAperture((x,y), r=r) for r in o.apertures.radii_pix] #r in pixels

    phot_table = aperture_photometry(data, apertures) 
    o.apertures.flux = np.array([phot_table['aperture_sum_{0}'.format(i)][0] for i, r in enumerate(o.apertures.radii_pix)])
    
    # --------------------------------------------------------------------------------------------------------
    # --- measure sizes

    o.sizes = empty()

    # --- Curve of growth

    o.sizes.COG = np.interp(0.5, o.apertures.flux/img.flux, o.apertures.radii_pix) * img.pixel_scale_kpc
    
    sortpix = np.array(sorted(data.flatten())[::-1])
    cumsum = np.cumsum(sortpix)/img.flux
    o.sizes.pixel = np.sqrt(len(cumsum[cumsum<0.5])/np.pi)  * img.pixel_scale_kpc
 
    if verbose: 
        print('    COG: {0:.2f}'.format(o.sizes.COG))
        print('    pixel: {0:.2f}'.format(o.sizes.pixel))


    return o




def measure_photometry(s, img, verbose = False):           
    
    photo = empty()
            
    x, y = s.x, s.y  
        
    mask = s.mask
      
        
    # --- true flux

    photo.true = np.sum(img.sci)

    if verbose: print('    '+'-'*5, 'observed (with PSF)')

    if verbose: print('    true flux: {0:.2f}'.format(photo.true))
        

    # --- detection region ISO

    photo.ISO = empty()
    photo.ISO.flux = np.sum(img.sci[mask])/img.nJy_to_es
    photo.ISO.noise = np.sqrt(np.sum(img.noise[mask]**2))/img.nJy_to_es

    if verbose: print('    ISO:', '{0:.2f}'.format(photo.ISO.flux), 'S/N: {0:.2f}'.format(photo.ISO.flux/photo.ISO.noise) )
   
      
    # --- Kron flux (AUTO)
 
    k = 2.5
    kron_correction = 1.16 # this needs to be updated for k = 2.5

    aperture = CircularAperture((x,y), r=k*s.kron_radius) 
    photo.circular_kron = empty()
    phot_table = aperture_photometry(img.sci, aperture, img.noise, mask = s.other_objects_mask) 
    photo.circular_kron.flux = phot_table['aperture_sum'][0]/img.nJy_to_es
    photo.circular_kron.noise = phot_table['aperture_sum_err'][0]/img.nJy_to_es

    if verbose: print('    AUTO:', '{0:.2f}'.format(photo.circular_kron.flux), 'S/N: {0:.2f}'.format(photo.circular_kron.flux/photo.circular_kron.noise) )

    
    # --- small Kron flux (AUTO)
 
    k = 1.
    aperture = CircularAperture((x,y), r=k*s.kron_radius) 
    photo.small_circular_kron = empty()
    phot_table = aperture_photometry(img.sci, aperture, img.noise, mask = s.other_objects_mask) 
    photo.small_circular_kron.flux = phot_table['aperture_sum'][0]/img.nJy_to_es
    photo.small_circular_kron.noise = phot_table['aperture_sum_err'][0]/img.nJy_to_es

    if verbose: print('    small AUTO:', '{0:.2f}'.format(photo.small_circular_kron.flux), 'S/N: {0:.2f}'.format(photo.small_circular_kron.flux/photo.small_circular_kron.noise) )
    
    
    # --- series of apertures (to give COG)

    photo.apertures = empty()
    photo.apertures.radii = np.arange(0.5, 20, 0.5)        

    apertures = [CircularAperture((x,y), r=r) for r in photo.apertures.radii] #r in pixels

    phot_table = aperture_photometry(img.sci, apertures, img.noise, mask = s.other_objects_mask) 
    photo.apertures.flux = np.array([phot_table['aperture_sum_{0}'.format(i)][0] for i, r in enumerate(photo.apertures.radii)])/img.nJy_to_es
    photo.apertures.noise = np.array([phot_table['aperture_sum_err_{0}'.format(i)][0] for i, r in enumerate(photo.apertures.radii)])/img.nJy_to_es 
    
    
    # --- determine optimum aperture

    argopt = np.argmax(photo.apertures.flux/photo.apertures.noise)

    photo.apertures.optimum = empty()
    photo.apertures.optimum.radius = photo.apertures.radii[argopt]
    photo.apertures.optimum.flux = photo.apertures.flux[argopt] 
    photo.apertures.optimum.noise = photo.apertures.noise[argopt]

    if verbose: print('    optimum aperture size: {0:.2f}'.format(photo.apertures.optimum.radius))

# 
#     # --------------------------------------------------------------------------------------------------------
#     # --- measure sizes
# 
#     photo.sizes = empty()
# 
# 
#     # --- Curve of growth
# 
#     photo.sizes.true_COG = np.interp(0.5, photo.apertures.noiseless_flux/photo.true, photo.apertures.radii) * img.pixel_scale_kpc
#     photo.sizes.noiseless_COG = np.interp(0.5, photo.apertures.noiseless_flux/(kron_correction*photo.circular_kron.noiseless_flux), photo.apertures.radii) * img.pixel_scale_kpc
#     photo.sizes.COG = np.interp(0.5, photo.apertures.flux/(kron_correction*photo.circular_kron.flux), photo.apertures.radii) * img.pixel_scale_kpc
# 
#     if verbose: print('    true COG: {0:.2f} noiseless COG: {1:.2f} COG: {2:.2f}'.format( photo.sizes.true_COG, photo.sizes.noiseless_COG, photo.sizes.COG))
# 
#     # --- pixels inside k*r_kron
# 
#     sortpix = np.array(sorted(img.img.flatten())[::-1])
#     cumsum = np.cumsum(sortpix)/photo.true
#     photo.sizes.true_pixel = np.sqrt(len(cumsum[cumsum<0.5])/np.pi) * img.pixel_scale_kpc
# 
#     mask = CircularAperture((x,y), r=k*s.kron_radius).to_mask(method='center') # --- mask all pixels in 2*r_kron 
#  
#     masked_image = mask[0].multiply(img.img)
#     sortpix = np.array(sorted(masked_image.flatten())[::-1])
#     cumsum = np.cumsum(sortpix)/(kron_correction*photo.circular_kron.flux)    
#     photo.sizes.noiseless_pixel = np.sqrt(len(cumsum[cumsum<0.5])/np.pi) * img.pixel_scale_kpc
# 
#     masked_image = mask[0].multiply(img.nimg)
#     sortpix = np.array(sorted(masked_image.flatten())[::-1])
#     cumsum = np.cumsum(sortpix)/(kron_correction*photo.circular_kron.flux)        
#     photo.sizes.pixel = np.sqrt(len(cumsum[cumsum<0.5])/np.pi) * img.pixel_scale_kpc
# 
#     if verbose: print('    true pixel: {0:.2f} noiseless pixel: {1:.2f} pixel: {2:.2f}'.format(photo.sizes.true_pixel, photo.sizes.noiseless_pixel, photo.sizes.pixel))
    
    return photo
