

import sys

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





    
    
    

def measure_core_properties(SourceProperties, DetectionImage, SegmentationImage, verbose = False):
    
    """Extracts some of the SourceProperties quantities and calculates the Kron radius"""
    
    ndim = DetectionImage.sci.shape[0]
    
    p = {}
    
    label = SourceProperties.label 
    p['label'] = label
    
    
    if verbose: print('-'*5, label)
    
    Mask = SegmentationImage.data == label
    ExclusionMask = np.isin(SegmentationImage.data, [label, 0], invert = True) # produce a mask excluding all other objects 
    

              
    x, y = SourceProperties.xcentroid.value, SourceProperties.ycentroid.value
    
    if verbose: print('x,y: {0} {1}'.format(x,y))
    
    p['x'], p['y'] = x, y
    
    p['area'] = len(SegmentationImage.data[SegmentationImage.data == label])
    p['radius'] = np.sqrt(p['area']/np.pi)
    
    if verbose: print('detection area/pix: {0}'.format(p['area']))
    if verbose: print('detection radius/pix: {0:.2f}'.format(p['radius']))
    
    p['A'] = SourceProperties.semimajor_axis_sigma.value
    p['B'] = SourceProperties.semiminor_axis_sigma.value
    p['theta'] = SourceProperties.orientation.value
    p['ellipticity'] = SourceProperties.ellipticity.value
    
    if verbose: print('A: {0:.2f} B:{1:.2f} theta:{2:.2f}'.format(p['A'],p['B'],p['theta']))

    # --- determine Kron radius
    
    X = np.linspace(0, ndim -1, ndim) - x
    Y = np.linspace(0, ndim -1, ndim) - y
    
    XX, YY = np.meshgrid(X,Y)
    
    R = np.sqrt(XX**2 + YY**2)
    
    aper = EllipticalAperture((x,y), 6.0*p['A'], 6.0*p['B'], p['theta'])
    phot_table = aperture_photometry(DetectionImage.sci, aper, mask = ExclusionMask) 
    f2 = phot_table['aperture_sum'][0]
    phot_table = aperture_photometry(DetectionImage.sci * R, aper, mask = ExclusionMask) 
    f1 = phot_table['aperture_sum'][0]
    p['kron_radius'] = np.max([1., f1/f2])
    
    
    
    # --- AUTO (large Kron)
    
    k = 2.5
    
    aperture = CircularAperture((x,y), r=k*p['kron_radius']) 
    
    phot_table = aperture_photometry(DetectionImage.sci, aperture, DetectionImage.noise, mask = ExclusionMask) 
    
    p['circular_kron'] = empty()
    p['circular_kron'].flux = phot_table['aperture_sum'][0]
    p['circular_kron'].noise = phot_table['aperture_sum_err'][0]
    
    kron_correction = 1.16 # this needs to be updated for k = 2.5
    p['circular_kron'].TOTAL_flux = p['circular_kron'].flux * kron_correction
    p['circular_kron'].TOTAL_noise = p['circular_kron'].noise * kron_correction
    
    
    # --- series of apertures (to give COG)

    p['aperture'] = empty()
    p['aperture'].radii = np.arange(0.5, 20, 0.5) # r in pixels     

    apertures = [CircularAperture((x,y), r=r) for r in p['aperture'].radii] 

    phot_table = aperture_photometry(DetectionImage.sci, apertures, DetectionImage.noise, mask = ExclusionMask) 
    p['aperture'].flux = np.array([phot_table['aperture_sum_{0}'.format(i)][0] for i, r in enumerate(p['aperture'].radii)])
    p['aperture'].noise = np.array([phot_table['aperture_sum_err_{0}'.format(i)][0] for i, r in enumerate(p['aperture'].radii)])
    
    p['aperture'].argopt = np.argmax(p['aperture'].flux/p['aperture'].noise)
    p['aperture'].optimum_radius = p['aperture'].radii[p['aperture'].argopt]
    
    if verbose: print('ISO flux: {0:.2f}'.format(np.sum(DetectionImage.sci[Mask])))  
    if verbose: print('flux at 6\sigma: {0:.2f}'.format(f2))
    if verbose: print('Kron radius: {0:.2f}'.format(p['kron_radius']))
    if verbose: print('optimum aperture radius: {0:.2f}'.format(p['aperture'].optimum_radius))
    
    
    # --- curve-of-growth size
    
    p['sizes'] = {}
    
    p['sizes']['COG'] = empty()
    p['sizes']['COG'].radius = np.interp(0.5, p['aperture'].flux/p['circular_kron'].TOTAL_flux, p['aperture'].radii)
    
    # --- pixels inside k*r_kron
    k = 2.5
    KronMask = CircularAperture((x,y), r=k*p['kron_radius']).to_mask(method='center')# --- mask all pixels in k*r_kron 
    if type(KronMask) is list: KronMask = KronMask[0]
    MaskedImage = DetectionImage.sci * (1-ExclusionMask.astype('float'))
    MaskedImage = KronMask.multiply(MaskedImage) 
    sortpix = np.array(sorted(MaskedImage.flatten())[::-1])
    cumsum = np.cumsum(sortpix)/p['circular_kron'].TOTAL_flux   
    
    p['sizes']['pixel'] = empty()
    p['sizes']['pixel'].radius = np.sqrt(len(cumsum[cumsum<0.5])/np.pi) 
    p['sizes']['pixel'].minflux = np.min(sortpix[cumsum<0.5]) # faintest pixel contributing to the size


    if verbose: print('SIZES: COG: {0:.2f} Pixel: {1:.2f}'.format(p['sizes']['COG'].radius, p['sizes']['pixel'].radius))
    
     
    return p, Mask, ExclusionMask

    


def measure_properties(p, img, Mask, ExclusionMask, verbose = False):           
    
    if verbose: print('    -----')
    
    photo = {}
            
    x, y = p['x'], p['y'] 
        
    # --- detection region ISO

    photo['ISO'] = empty()
    photo['ISO'].flux = np.sum(img.sci[Mask])/img.nJy_to_es
    photo['ISO'].error = np.sqrt(np.sum(img.noise[Mask]**2))/img.nJy_to_es

    if verbose: print('    ISO:', '{0:.2f}'.format(photo['ISO'].flux), 'S/N: {0:.2f}'.format(photo['ISO'].flux/photo['ISO'].error) )
     
    # --- Kron flux (AUTO)
 
    k = 2.5
    
    aperture = CircularAperture((x,y), r=k*p['kron_radius']) 
    
    phot_table = aperture_photometry(img.sci, aperture, img.noise, mask = ExclusionMask) 
    
    photo['circular_kron'] = empty()
    photo['circular_kron'].flux = phot_table['aperture_sum'][0]/img.nJy_to_es
    photo['circular_kron'].error = phot_table['aperture_sum_err'][0]/img.nJy_to_es
    
    kron_correction = 1.16 # this needs to be updated for k = 2.5
    photo['circular_kron'].TOTAL_flux = photo['circular_kron'].flux * kron_correction
    photo['circular_kron'].TOTAL_error = photo['circular_kron'].error * kron_correction

    if verbose: print('    AUTO:', '{0:.2f}'.format(photo['circular_kron'].flux), 'S/N: {0:.2f}'.format(photo['circular_kron'].flux/photo['circular_kron'].error) )

    
    # --- small Kron flux (AUTO)
 
    k = 1.
    aperture = CircularAperture((x,y), r=k*p['kron_radius'])  
    phot_table = aperture_photometry(img.sci, aperture, img.noise, mask = ExclusionMask) 
    photo['small_circular_kron'] = empty()
    photo['small_circular_kron'].flux = phot_table['aperture_sum'][0]/img.nJy_to_es
    photo['small_circular_kron'].error = phot_table['aperture_sum_err'][0]/img.nJy_to_es

    if verbose: print('    small AUTO:', '{0:.2f}'.format(photo['small_circular_kron'].flux), 'S/N: {0:.2f}'.format(photo['small_circular_kron'].flux/photo['small_circular_kron'].error) )
    
    
    # --- series of apertures (to give COG)

    photo['aperture'] = empty()
    photo['aperture'].radii = p['aperture'].radii # r in pixels     
    apertures = [CircularAperture((x,y), r=r) for r in photo['aperture'].radii] 
    phot_table = aperture_photometry(img.sci, apertures, img.noise, mask = ExclusionMask) 
    photo['aperture'].flux = np.array([phot_table['aperture_sum_{0}'.format(i)][0] for i, r in enumerate(photo['aperture'].radii)])/img.nJy_to_es
    photo['aperture'].error = np.array([phot_table['aperture_sum_err_{0}'.format(i)][0] for i, r in enumerate(photo['aperture'].radii)])/img.nJy_to_es 
    photo['aperture'].argopt = np.argmax(photo['aperture'].flux/photo['aperture'].error)
    photo['aperture'].optimum_radius = photo['aperture'].radii[photo['aperture'].argopt]
    
    
    # --- use optimum detection image aperture

    argopt = p['aperture'].argopt
    photo['optimum_aperture'] = empty()
    photo['optimum_aperture'].radius = photo['aperture'].radii[argopt]
    photo['optimum_aperture'].flux = photo['aperture'].flux[argopt] 
    photo['optimum_aperture'].error = photo['aperture'].error[argopt]

    if verbose: print('    optimum aperture:', '{0:.2f}'.format(photo['optimum_aperture'].flux), 'S/N: {0:.2f}'.format(photo['optimum_aperture'].flux/photo['optimum_aperture'].error) )
 
    if photo['circular_kron'].TOTAL_flux > 0.0: 

        x,y = p['x'],p['y']

        sizes = {}
    
        # --- curve-of-growth size
        sizes['COG'] = empty()
        sizes['COG'].radius = np.interp(0.5, photo['aperture'].flux/photo['circular_kron'].TOTAL_flux, photo['aperture'].radii)
    
        if verbose: print('    r_e (COG)/pix: {0:.2f}'.format(sizes['COG'].radius))
    
        # --- pixels inside k*r_kron
        k = 2.5
        KronMask = CircularAperture((x,y), r=k*p['kron_radius']).to_mask(method='center') # --- mask all pixels in k*r_kron 
        if type(KronMask) is list: KronMask = KronMask[0]
        MaskedImage = img.sci * (1-ExclusionMask.astype('float'))
        MaskedImage = KronMask.multiply(MaskedImage) 
        sortpix = np.array(sorted(MaskedImage.flatten())[::-1])/img.nJy_to_es
        
        sortpix = sortpix[sortpix>0.0] # ignore negative pixels
        
        cumsum = np.cumsum(sortpix)/photo['circular_kron'].TOTAL_flux  
    
    
        if len(cumsum[cumsum<0.5]) > 0:
    
            sizes['pixel'] = empty()
            sizes['pixel'].radius = np.sqrt(len(cumsum[cumsum<0.5])/np.pi) 
            sizes['pixel'].minflux = np.min(sortpix[cumsum<0.5]) # faintest pixel contributing to the size
            
            if verbose: print('    r_e (Pixel)/pix: {0:.2f}'.format(sizes['pixel'].radius))

        else:
    
            sizes['pixel'] = False
    
            if verbose: print('    r_e (Pixel)/pix: UNDEFINED')
        
    else:
    
        sizes = False
        
    return {'photometry': photo, 'sizes': sizes}










def measure_model_properties(img, verbose = False):           
    
    if verbose: print('    -----')
    
    photo = {}
    
    photo['total'] = empty()
    photo['total'].flux = np.sum(img.sci)/img.nJy_to_es
    
    if verbose: print('    total flux/nJy: {0:.2f}'.format(photo['total'].flux))
     
    # --- series of apertures (to give COG)


    x = y = img.sci.shape[0] / 2.

    photo['aperture'] = empty()
    photo['aperture'].radii = np.arange(0.5, 20., 0.5) # r in pixels     
    apertures = [CircularAperture((x,y), r=r) for r in photo['aperture'].radii] 
    phot_table = aperture_photometry(img.sci, apertures) 
    photo['aperture'].flux = np.array([phot_table['aperture_sum_{0}'.format(i)][0] for i, r in enumerate(photo['aperture'].radii)])/img.nJy_to_es
    
    sizes = {}
    
    # --- curve-of-growth size
    sizes['COG'] = empty()
    sizes['COG'].radius = np.interp(0.5, photo['aperture'].flux/photo['total'].flux, photo['aperture'].radii)

    if verbose: print('    r_e (COG)/pix: {0:.2f}'.format(sizes['COG'].radius))

    # --- pixels inside 

    sortpix = np.array(sorted(img.sci.flatten())[::-1])/img.nJy_to_es
    sortpix = sortpix[sortpix>0.0] # ignore negative pixels
    cumsum = np.cumsum(sortpix)/photo['total'].flux  
    
    sizes['pixel'] = empty()
    sizes['pixel'].radius = np.sqrt(len(cumsum[cumsum<0.5])/np.pi) 
    sizes['pixel'].minflux = np.min(sortpix[cumsum<0.5]) # faintest pixel contributing to the size
    
    if verbose: print('    r_e (Pixel)/pix: {0:.2f}'.format(sizes['pixel'].radius))

    return {'photometry': photo, 'sizes': sizes}









def COG_plots(Properties, ModelProperties = False, filename = False, show = False):
   
    nfilters = len(Properties)
   
    fig, axes = plt.subplots(1, nfilters, figsize = (3*(nfilters),3), dpi = 200)
    
    plt.subplots_adjust(left=0.025, top=0.85, bottom=0.2, right=0.9, wspace=0.2, hspace=0.0)
          
    for ax, (filter, properties) in zip(axes, Properties.items()): 
    
        ax.set_title(filter, fontsize = 10)
    
        ax.plot(properties['photometry']['aperture'].radii, properties['photometry']['aperture'].flux, c = '0.5', label = 'curve-of-growth')
        ax.axvline(properties['photometry']['aperture'].optimum_radius, color = '0.5', alpha = 0.5)
    
        if ModelProperties is not False:
            
            ax.axhline(ModelProperties[filter]['photometry']['total'].flux, color = '0.5', alpha = 0.5, ls = ':')
            ax.plot(ModelProperties[filter]['photometry']['aperture'].radii, ModelProperties[filter]['photometry']['aperture'].flux, c = '0.5', ls = ':', label = 'true curve-of-growth')
            
    
    
        del properties['photometry']['aperture']
    
        color_idx = np.linspace(0, 1, len(properties['photometry']))
        
        for c_idx, (phot_type, p) in zip(color_idx, properties['photometry'].items()):
    
            ax.axhline(p.flux, label = phot_type, color = cm.viridis(c_idx))
            ax.axhspan(p.flux-p.error, p.flux+p.error, color = cm.viridis(c_idx), alpha=0.5)
            if phot_type == 'optimum_aperture': ax.axvline(p.radius, color = cm.viridis(c_idx), alpha = 0.5)
    
    ax.legend(bbox_to_anchor=(1.1, 1.0), fontsize = 8)
    
    if filename:
        plt.savefig(filename)
    if show:
        plt.show()
    
    plt.close(fig)



def SED_plot(Properties,  ModelProperties = False, FilterInfo = False, phot_type = 'optimum_aperture', filename = False, show = False):
    
        
    # if not FilterInfo: 
   
    fig, ax = plt.subplots(1, 1, figsize = (3,2), dpi = 200)
    plt.subplots_adjust(left=0.2, top=0.85, bottom=0.25, right=0.9, wspace=0.2, hspace=0.0)
    
    color_idx = np.linspace(0, 1, len(Properties))
    
    for c_idx, (filter, properties) in zip(color_idx, Properties.items()): 
    
        pivwv = FilterInfo[filter].pivwv()/1E4
    
        ax.scatter(pivwv, properties['photometry'][phot_type].flux, color = cm.viridis(c_idx))
        ax.plot([pivwv]*2, [properties['photometry'][phot_type].flux - properties['photometry'][phot_type].error, properties['photometry'][phot_type].flux + properties['photometry'][phot_type].error], color = 'k', lw = 1)
    
        if ModelProperties is not False:
        
            ax.scatter(pivwv, ModelProperties[filter]['photometry']['total'].flux, color = cm.viridis(c_idx), alpha = 0.5)
            
    
    
    ax.set_xlabel(r'$\lambda/\mu m$')
    ax.set_ylabel(r'$f_{\nu}/nJy$')
        
    if filename:
        plt.savefig(filename)
    if show:
        plt.show()
    
    plt.close(fig)




def size_plot(img, p, ExclusionMask, threshold = 2.5, signficance_plot = False, filename = False, show = False, add_contours = False):

    
    width = img.sci.shape[0]
    
    fig, ax = plt.subplots(1, 1, figsize = (3,3), dpi = width*2)
    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.01, hspace=0.0)
        
    ax.set_axis_off()

    sig = (img.sci/img.noise)

    ax.imshow(sig, cmap = cm.Greys, vmin = -5.0, vmax = 5.0, origin = 'lower') 
    ax.imshow(np.ma.masked_where(sig <= threshold, sig), cmap = cm.plasma, vmin = threshold, vmax = 100, origin = 'lower') 

    k = 2.5
    
    # --- make mask image including Kron Mask and Exclusion mask
    x = np.linspace(-(width//2), (width//2), width)
    X, Y = np.meshgrid(x, x)
    R2 = X**2 + Y**2
    alpha = np.zeros(img.sci.shape)
    alpha[R2>(k*p['kron_radius'])**2] = 1
    alpha[img.sci<p['sizes']['pixel'].minflux] = 1
    alpha[ExclusionMask] = 1

    # --- make RGBA image which is white everywhere but transparent in regions that are counted
    RGBA = np.ones((width, width, 4))
    RGBA[:,:,3] = alpha

    ax.imshow(RGBA, origin = 'lower', alpha = 0.8)     
    
    # --- add contours around pixels contributing to r_pix (VERY SLOW FOR SOME REASON)
    if add_contours:
        image = alpha[::-1]
        f = lambda x,y: image[int(y),int(x) ]
        g = np.vectorize(f)

        x = np.linspace(0,image.shape[1], image.shape[1]*100)
        y = np.linspace(0,image.shape[0], image.shape[0]*100)
        X, Y= np.meshgrid(x[:-1],y[:-1])
        Z = g(X[:-1],Y[:-1])

        ax.contour(Z[::-1], [0.5], colors='k', linewidths=[2], extent=[0-0.5, x[:-1].max()-0.5,0-0.5, y[:-1].max()-0.5])

    
    # --- add k*r_Kron radius
    kKronRadius = plt.Circle((width//2,width//2), k*p['kron_radius'], alpha = 0.2)
    ax.add_artist(kKronRadius)
    # --- add COG radius
    COGRadius = plt.Circle((width//2,width//2), p['sizes']['COG'].radius, alpha = 1.0, fill = False, lw = 1, color='1')
    ax.add_artist(COGRadius)
    
    # --- add COG radius
    PixelRadius = plt.Circle((width//2,width//2), p['sizes']['pixel'].radius, alpha = 1.0, fill = False, lw = 1, ls = '--', color='1')
    ax.add_artist(PixelRadius)


    if filename:
        plt.savefig(filename)
    if show:
        plt.show()
    
    plt.close(fig)

