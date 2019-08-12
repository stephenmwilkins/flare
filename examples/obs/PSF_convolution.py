

# --- core modules

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# --- astro-specific

from astropy.io import fits
from astropy.convolution import convolve, convolve_fft
from photutils import create_matching_kernel, TopHatWindow

# --- own installed

from SynthObs.Morph import PSF

# --- own, this package

sys.path.insert(0, os.path.abspath('../../'))

import FLARE
import FLARE.obs
import FLARE.surveys
import FLARE.photom




survey = FLARE.surveys.XDF
field = survey.fields['XDF']
ndim = 21
window = TopHatWindow(1.0)
plot = False


target_filter = field.filters[-1]
psf = PSF.PSF(target_filter) 
pixel_scale = field.pixel_scale
native_pixel_Scale = FLARE.filters.pixel_scale[target_filter]
x = y = np.linspace(-(ndim/2.)*(pixel_scale/native_pixel_Scale), (ndim/2.)*(pixel_scale/native_pixel_Scale), ndim) # in original pixels
target_kernel = psf.f(x,y) 
target_kernel /= np.sum(target_kernel)



for filter in field.filters:

    # --- calculate transfer kernel

    print('-'*10, filter)

    native_pixel_Scale = FLARE.filters.pixel_scale[filter]
    x = y = np.linspace(-(ndim/2.)*(pixel_scale/native_pixel_Scale), (ndim/2.)*(pixel_scale/native_pixel_Scale), ndim) # in original pixels

    psf = PSF.PSF(filter) 
    kernel = psf.f(x,y) 
    kernel /= np.sum(kernel)

    if window:
        transfer_kernel = create_matching_kernel(kernel, target_kernel, window = window)
    else:
        transfer_kernel = create_matching_kernel(kernel, target_kernel)


    np.save('{0}/{1}/{2}_kernel.npy'.format(FLARE.FLARE_dir, survey.datadir, filter.split('.')[-1]), transfer_kernel) # save transfer kernel

    # --- do convolution
    
    hdu = fits.open('{0}/{1}_sci.fits'.format(FLARE.FLARE_dir, field.filename[filter])) # --- open original file
    original_data = hdu[0].data
    hdu[0].data = convolve(original_data, transfer_kernel) # --- convolve with transfer kernel
    hdu.writeto('{0}/{1}_sci_convolved.fits'.format(FLARE.FLARE_dir, field.filename[filter]))
    hdu.close()
    
    
    if plot:
    
        fig, axes = plt.subplots(1,4,figsize = (12,3))
        plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.0, hspace=0.0)

        vmax = 0.3
        axes[0].imshow(target_kernel, vmin = 0, vmax = vmax) 
        axes[1].imshow(kernel, vmin = 0, vmax = vmax) 
        axes[2].imshow(transfer_kernel, vmin = 0, vmax = vmax)
        axes[3].imshow(np.fabs(convolve_fft(kernel,transfer_kernel) - target_kernel)/np.max(target_kernel), vmin = 0, vmax = 1.0)      
            
        plt.show()
        fig.clf()
        # input()   

        fig, axes = plt.subplots(1,2,figsize = (10,5))
        plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.0, hspace=0.0)

        axes[0].imshow(old_data[2800:3200,2800:3200]) 
        axes[1].imshow(hdu[0].data[2800:3200,2800:3200]) 
            
        plt.show()
        fig.clf()
        # input()    
    
    
    

