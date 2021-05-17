
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from photutils import detect_sources
from photutils import source_properties
from photutils import CircularAperture
from photutils import EllipticalAperture
from photutils import aperture_photometry
from photutils import data_properties

import FLARE.filters

from types import SimpleNamespace



aperture_radii = np.arange(0.5, 50, 0.5) # r in pixels



def measure_core_properties(SourceProperties, DetectionImage, SegmentationImage, verbose = False, measure_apertures = True, save_apertures = False):

    """Extracts some of the SourceProperties quantities and calculates the Kron radius"""

    ndim = DetectionImage.sci.shape[0]

    p = {}

    label = SourceProperties.label
    p['label'] = label


    if verbose:
        print()
        print('-'*10, 'Detection Image Properties')

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

    p['circular_kron'] = SimpleNamespace()
    p['circular_kron'].flux = phot_table['aperture_sum'][0]
    p['circular_kron'].noise = phot_table['aperture_sum_err'][0]

    kron_correction = 1.16 # this needs to be updated for k = 2.5
    p['circular_kron'].TOTAL_flux = p['circular_kron'].flux * kron_correction
    p['circular_kron'].TOTAL_noise = p['circular_kron'].noise * kron_correction


    # --- series of apertures (to give COG)

    if measure_apertures:

        p['aperture'] = SimpleNamespace()
        p['aperture'].radii = aperture_radii

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

        p['sizes']['COG'] = SimpleNamespace()
        p['sizes']['COG'].radius = np.interp(0.5, p['aperture'].flux/p['circular_kron'].TOTAL_flux, p['aperture'].radii)

        # --- pixels inside k*r_kron
        k = 2.5
        KronMask = CircularAperture((x,y), r=k*p['kron_radius']).to_mask(method='center')# --- mask all pixels in k*r_kron
        if type(KronMask) is list: KronMask = KronMask[0]
        MaskedImage = DetectionImage.sci * (1-ExclusionMask.astype('float'))
        MaskedImage = KronMask.multiply(MaskedImage)
        sortpix = np.array(sorted(MaskedImage.flatten())[::-1])
        cumsum = np.cumsum(sortpix)/p['circular_kron'].TOTAL_flux

        p['sizes']['pixel'] = SimpleNamespace()
        p['sizes']['pixel'].radius = np.sqrt(len(cumsum[cumsum<0.5])/np.pi)
        p['sizes']['pixel'].minflux = np.min(sortpix[cumsum<0.5]) # faintest pixel contributing to the size


        if verbose: print('SIZES: COG: {0:.2f} Pixel: {1:.2f}'.format(p['sizes']['COG'].radius, p['sizes']['pixel'].radius))

        if not save_apertures:
            del p['aperture'].radii
            del p['aperture'].flux
            del p['aperture'].noise

    return p, Mask, ExclusionMask




def measure_properties(p, img, Mask, ExclusionMask, verbose = False, save_apertures = False):

    if verbose:
        if hasattr(img, 'filter'):
            print(f'-----{img.filter}')
        else:
            print('-----')

    photo = {}

    x, y = p['x'], p['y']

    # --- detection region ISO

    photo['ISO'] = SimpleNamespace()
    photo['ISO'].flux = np.sum(img.sci[Mask])/img.nJy_to_es
    photo['ISO'].error = np.sqrt(np.sum(img.noise[Mask]**2))/img.nJy_to_es

    if verbose: print('ISO:', '{0:.2f}'.format(photo['ISO'].flux), 'S/N: {0:.2f}'.format(photo['ISO'].flux/photo['ISO'].error) )

    # --- Kron flux (AUTO)

    k = 2.5

    aperture = CircularAperture((x,y), r=k*p['kron_radius'])

    phot_table = aperture_photometry(img.sci, aperture, img.noise, mask = ExclusionMask)

    photo['circular_kron'] = SimpleNamespace()
    photo['circular_kron'].flux = phot_table['aperture_sum'][0]/img.nJy_to_es
    photo['circular_kron'].error = phot_table['aperture_sum_err'][0]/img.nJy_to_es

    kron_correction = 1.16 # this needs to be updated for k = 2.5
    photo['circular_kron'].TOTAL_flux = photo['circular_kron'].flux * kron_correction
    photo['circular_kron'].TOTAL_error = photo['circular_kron'].error * kron_correction

    if verbose: print('AUTO:', '{0:.2f}'.format(photo['circular_kron'].flux), 'S/N: {0:.2f}'.format(photo['circular_kron'].flux/photo['circular_kron'].error) )


    # --- small Kron flux (AUTO)

    k = 1.
    aperture = CircularAperture((x,y), r=k*p['kron_radius'])
    phot_table = aperture_photometry(img.sci, aperture, img.noise, mask = ExclusionMask)
    photo['small_circular_kron'] = SimpleNamespace()
    photo['small_circular_kron'].flux = phot_table['aperture_sum'][0]/img.nJy_to_es
    photo['small_circular_kron'].error = phot_table['aperture_sum_err'][0]/img.nJy_to_es

    if verbose: print('small AUTO:', '{0:.2f}'.format(photo['small_circular_kron'].flux), 'S/N: {0:.2f}'.format(photo['small_circular_kron'].flux/photo['small_circular_kron'].error) )


    # --- series of apertures (to give COG)

    photo['aperture'] = SimpleNamespace()
    photo['aperture'].radii = aperture_radii # r in pixels
    apertures = [CircularAperture((x,y), r=r) for r in photo['aperture'].radii]
    phot_table = aperture_photometry(img.sci, apertures, img.noise, mask = ExclusionMask)
    photo['aperture'].flux = np.array([phot_table['aperture_sum_{0}'.format(i)][0] for i, r in enumerate(photo['aperture'].radii)])/img.nJy_to_es
    photo['aperture'].error = np.array([phot_table['aperture_sum_err_{0}'.format(i)][0] for i, r in enumerate(photo['aperture'].radii)])/img.nJy_to_es
    photo['aperture'].argopt = np.argmax(photo['aperture'].flux/photo['aperture'].error)
    photo['aperture'].optimum_radius = photo['aperture'].radii[photo['aperture'].argopt]


    # --- use optimum detection image aperture

    argopt = p['aperture'].argopt
    photo['optimum_aperture'] = SimpleNamespace()
    photo['optimum_aperture'].radius = photo['aperture'].radii[argopt]
    photo['optimum_aperture'].flux = photo['aperture'].flux[argopt]
    photo['optimum_aperture'].error = photo['aperture'].error[argopt]

    if verbose: print('optimum aperture:', '{0:.2f}'.format(photo['optimum_aperture'].flux), 'S/N: {0:.2f}'.format(photo['optimum_aperture'].flux/photo['optimum_aperture'].error) )

    if photo['circular_kron'].TOTAL_flux > 0.0:

        x,y = p['x'],p['y']

        sizes = {}

        # --- curve-of-growth size
        sizes['COG'] = SimpleNamespace()
        sizes['COG'].radius = np.interp(0.5, photo['aperture'].flux/photo['circular_kron'].TOTAL_flux, photo['aperture'].radii)

        if verbose: print('r_e (COG)/pix: {0:.2f}'.format(sizes['COG'].radius))

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

            sizes['pixel'] = SimpleNamespace()
            sizes['pixel'].radius = np.sqrt(len(cumsum[cumsum<0.5])/np.pi)
            sizes['pixel'].minflux = np.min(sortpix[cumsum<0.5]) # faintest pixel contributing to the size

            if verbose: print('r_e (Pixel)/pix: {0:.2f}'.format(sizes['pixel'].radius))

        else:

            sizes['pixel'] = SimpleNamespace()
            sizes['pixel'].radius = -99
            sizes['pixel'].minflux = -99

            if verbose: print('r_e (Pixel)/pix: UNDEFINED')

    else:

        sizes = {}

        sizes['COG'] = SimpleNamespace()
        sizes['COG'].radius = -99

        sizes['pixel'] = SimpleNamespace()
        sizes['pixel'].radius = -99
        sizes['pixel'].minflux = -99

    if not save_apertures:
        del photo['aperture'].radii
        del photo['aperture'].flux
        del photo['aperture'].error


    return {'photometry': photo, 'sizes': sizes}



def measure_properties_quick(p, img, Mask, ExclusionMask, verbose = False):

    if verbose:
        if hasattr(img, 'filter'):
            print(f'-----{img.filter}')
        else:
            print('-----')

    photo = {}

    x, y = p['x'], p['y']

    # --- Kron flux (AUTO)
    k = 2.5
    aperture = CircularAperture((x,y), r=k*p['kron_radius'])
    phot_table = aperture_photometry(img.sci, aperture, img.noise, mask = ExclusionMask)
    photo['circular_kron'] = SimpleNamespace()
    photo['circular_kron'].flux = phot_table['aperture_sum'][0]/img.nJy_to_es
    photo['circular_kron'].error = phot_table['aperture_sum_err'][0]/img.nJy_to_es

    # --- small Kron flux (AUTO)
    k = 1.
    aperture = CircularAperture((x,y), r=k*p['kron_radius'])
    phot_table = aperture_photometry(img.sci, aperture, img.noise, mask = ExclusionMask)
    photo['small_circular_kron'] = SimpleNamespace()
    photo['small_circular_kron'].flux = phot_table['aperture_sum'][0]/img.nJy_to_es
    photo['small_circular_kron'].error = phot_table['aperture_sum_err'][0]/img.nJy_to_es

    return photo






def measure_model_properties(img, verbose = False, save_apertures = False):

    if verbose:
        if hasattr(img, 'filter'):
            print(f'-----{img.filter}')
        else:
            print('-----')

    photo = {}

    photo['total'] = SimpleNamespace()
    photo['total'].flux = np.sum(img.sci)/img.nJy_to_es

    if verbose: print('    total flux/nJy: {0:.2f}'.format(photo['total'].flux))

    # --- series of apertures (to give COG)


    x = y = img.sci.shape[0] / 2.

    photo['aperture'] = SimpleNamespace()
    photo['aperture'].radii = aperture_radii # r in pixels
    apertures = [CircularAperture((x,y), r=r) for r in photo['aperture'].radii]
    phot_table = aperture_photometry(img.sci, apertures)
    photo['aperture'].flux = np.array([phot_table['aperture_sum_{0}'.format(i)][0] for i, r in enumerate(photo['aperture'].radii)])/img.nJy_to_es

    sizes = {}

    # --- curve-of-growth size
    sizes['COG'] = SimpleNamespace()
    sizes['COG'].radius = np.interp(0.5, photo['aperture'].flux/photo['total'].flux, photo['aperture'].radii)

    if verbose: print('    r_e (COG)/pix: {0:.2f}'.format(sizes['COG'].radius))

    # --- pixels inside

    sortpix = np.array(sorted(img.sci.flatten())[::-1])/img.nJy_to_es
    sortpix = sortpix[sortpix>0.0] # ignore negative pixels
    cumsum = np.cumsum(sortpix)/photo['total'].flux

    if len(cumsum[cumsum<0.5]) > 0: # this shouldn't be needed here but seemingly is.

        sizes['pixel'] = SimpleNamespace()
        sizes['pixel'].radius = np.sqrt(len(cumsum[cumsum<0.5])/np.pi)
        sizes['pixel'].minflux = np.min(sortpix[cumsum<0.5]) # faintest pixel contributing to the size

        if verbose: print('r_e (Pixel)/pix: {0:.2f}'.format(sizes['pixel'].radius))

    else:

        sizes['pixel'] = SimpleNamespace()
        sizes['pixel'].radius = -99
        sizes['pixel'].minflux = -99

        if verbose: print('r_e (Pixel)/pix: UNDEFINED')


    if not save_apertures:
        del photo['aperture'].radii
        del photo['aperture'].flux


    return {'photometry': photo, 'sizes': sizes}
