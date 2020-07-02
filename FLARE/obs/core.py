
import numpy as np

from astropy.io import fits
# from photutils import CircularAperture
# from photutils import aperture_photometry

import FLARE
import FLARE.filters

class empty: pass

def make_cutout(data, x, y, width):

    """extract cut out from arbitrary data"""

    cutout = np.zeros((width, width))

    xmin = x - width // 2
    xmax = x + width // 2
    ymin = y - width // 2
    ymax = y + width // 2

    xstart = 0
    ystart = 0
    xend = width
    yend = width

    if xmin < 0:
        xstart = -xmin
        xmin = 0
    if ymin < 0:
        ystart = -ymin
        ymin = 0
    if xmax > data.shape[0]:
        xend -= xmax - data.shape[0]
        xmax = data.shape[0]
    if ymax > data.shape[1]:
        yend -= ymax - data.shape[1]
        ymax = data.shape[1]

    data = np.array(data)
    cutout[xstart:xend,ystart:yend] = data[xmin:xmax,ymin:ymax]

    return cutout







class image:

    def get_random_location(self):

        """get (single) random location on the image"""

        pos = np.random.choice(self.sci.count())
        return np.take((~self.sci.mask).nonzero(), pos, axis=1)

    def get_random_locations(self, N):

        """get N random locations on the image"""

        pos = np.random.choice(self.sci.count(), size=N)
        return np.take((~self.sci.mask).nonzero(), pos, axis=1)


    def get_area(self):

        """calculated non-masked area in units of arcmin2"""

        return self.sci.count()*self.pixel_scale**2/3600.


    def make_cutout(self, x, y, width):

        """extract cut out"""

        sci = np.zeros((width, width))
        wht = np.zeros((width, width))

        xmin = x - width // 2
        xmax = x + width // 2
        ymin = y - width // 2
        ymax = y + width // 2

        xstart = 0
        ystart = 0
        xend = width
        yend = width

        if xmin < 0:
            xstart = -xmin
            xmin = 0
        if ymin < 0:
            ystart = -ymin
            ymin = 0
        if xmax > self.sci.shape[0]:
            xend -= xmax - self.sci.shape[0]
            xmax = self.sci.shape[0]
        if ymax > self.sci.shape[1]:
            yend -= ymax - self.sci.shape[1]
            ymax = self.sci.shape[1]

#         if self.verbose: print(xmin, xmax, ymin, ymax)
#         if self.verbose: print(xstart, xend)
#         if self.verbose: print(ystart, yend)
#         if self.verbose: print(sci.shape)
#         if self.verbose: print(self.sci[xmin:xmax,ymin:ymax].shape)

        if (width % 2) != 0:
            xmax += 1
            ymax += 1

        sci[xstart:xend,ystart:yend] = self.sci[xmin:xmax,ymin:ymax]
        wht[xstart:xend,ystart:yend] = self.wht[xmin:xmax,ymin:ymax]

        return image_from_arrays(sci, wht, self.pixel_scale, zeropoint = self.zeropoint, nJy_to_es = self.nJy_to_es, verbose = self.verbose)



    def determine_depth(self, N = 10000, aperture_diameter_arcsec = 0.35, sigma = 5.):

        """determine depth using random apertures"""

        aperture_centres = tuple(self.get_random_locations(N).T)
        apertures = [CircularAperture(aperture_centres, r=r) for r in [(aperture_diameter_arcsec/self.pixel_scale)/2.]] # r in pixels
        phot_table = aperture_photometry(self.sci, apertures)
        aperture_fluxes = phot_table['aperture_sum_0'].quantity
        negative_aperture_fluxes = aperture_fluxes[aperture_fluxes<0]
        return -np.percentile(negative_aperture_fluxes, 100.-68.3) * sigma


#     def make_significance_plot(self, threshold = 2.5):
#
#         """make significance plot"""
#
#         sig = (self.sci/self.noise)
#
#         fig, ax = plt.subplots(1, 1, figsize = (4,4))
#
#         ax.set_axis_off()
#         ax.imshow(sig, cmap = cm.Greys, vmin = -5.0, vmax = 5.0, origin = 'lower')
#         ax.imshow(np.ma.masked_where(sig <= threshold, sig), cmap = cm.plasma, vmin = threshold, vmax = 100, origin = 'lower')
#
#         plt.show()



    def write_to_fits(self, filename = 'temp/'):

        sci_hdu = fits.PrimaryHDU(self.sci)
        sci_hdu.writeto(f'{filename}sci.fits')

        wht_hdu = fits.PrimaryHDU(self.wht)
        wht_hdu.writeto(f'{filename}wht.fits')

        rms_hdu = fits.PrimaryHDU(self.noise)
        rms_hdu.writeto(f'{filename}rms.fits')



class image_from_file(image):

    def __init__(self, filename, filter, mask = False, pixel_scale = 0.06, verbose = False, sci_suffix = 'sci', wht_suffix = 'wht'):

        """generate instance of image class from file"""

        if verbose:
            print('-'*40)
            print('filter: {0}'.format(filter))
            print('reading image from: {0}'.format(filename))


        self.verbose = verbose

        self.filter = filter
        self.pixel_scale = pixel_scale

        self.sci = fits.getdata('{0}/{1}_{2}.fits'.format(FLARE.FLARE_dir, filename, sci_suffix))
        self.wht = fits.getdata('{0}/{1}_{2}.fits'.format(FLARE.FLARE_dir, filename, wht_suffix))

        self.zeropoint = FLARE.filters.info[filter].zeropoint # AB magnitude zeropoint
        if verbose: print('zeropoint: {0}'.format(self.zeropoint))

        self.nJy_to_es = 1E-9 * 10**(0.4*(self.zeropoint-8.9)) # conversion from nJy to e/s
        if verbose: print('nJy_to_es: {0}'.format(self.nJy_to_es))


        self.mask = mask

        if type(mask) == np.ndarray:
            self.mask = mask
        else:
            self.mask = (self.wht == 0)

        self.sci = np.ma.masked_array(self.sci, mask = self.mask)
        self.wht = np.ma.masked_array(self.wht, mask = self.mask)

        self.noise = 1./np.sqrt(self.wht)
        self.sig = self.sci/self.noise


class image_from_arrays(image):

    def __init__(self, sci, wht, pixel_scale, zeropoint = False, nJy_to_es = False,  verbose = False):

        """generate instance of image class from cutout"""

        self.verbose = verbose

        self.pixel_scale = pixel_scale
        self.zeropoint = zeropoint # AB magnitude zeropoint
        self.nJy_to_es = nJy_to_es # conversion from nJy to e/s

        self.sci = sci
        self.wht = wht
        self.noise = 1./np.sqrt(self.wht)
        self.sig = self.sci/self.noise


def open_image(field, filter, verbose = False, sci_suffix = 'sci'):

    if field.mask_file:
        mask = fits.getdata('{0}/{1}/{2}'.format(FLARE.FLARE_dir, field.datadir, field.mask_file))
    else:
        mask = False

    return image_from_file(field.filename[filter], filter, mask = mask, verbose = verbose, sci_suffix = sci_suffix)



def open_images(field, filters, verbose = False, sci_suffix = 'sci'):

    return {filter: open_image(field, filter, verbose = verbose, sci_suffix = sci_suffix) for filter in filters}


def create_stack(imgs):

    first_img = next(iter(imgs.values()))

    shape = first_img.sci.shape
    sci = np.zeros(shape)
    wht = np.zeros(shape)

    for filter, img in imgs.items():
        sci += img.sci * img.wht
        wht += img.wht

    sci /= wht

    return image_from_arrays(sci, wht, first_img.pixel_scale)







class Background():

    def __init__(self, filter, pixel_scale, aperture_f_limit, aperture_significance, aperture_radius, verbose = False):

        self.filter = filter
        self.pixel_scale = pixel_scale
        self.zeropoint = FLARE.filters.info[filter].zeropoint # AB magnitude zeropoint # AB magnitude zeropoint
        self.nJy_to_es = 1E-9 * 10**(0.4*(self.zeropoint-8.9)) # conversion from nJy to e/s

        self.aperture = empty()
        self.aperture.flux_limit = aperture_f_limit
        self.aperture.radius = aperture_radius
        self.aperture.significance = aperture_significance
        self.aperture.noise = self.aperture.flux_limit/self.aperture.significance # nJy
        self.aperture.background = self.aperture.noise**2
        self.aperture.area = np.pi * self.aperture.radius**2
        self.aperture.noise_es = self.aperture.noise * self.nJy_to_es # convert from nJy to e/s

        self.pixel = empty()
        self.pixel.background = self.aperture.background/self.aperture.area
        self.pixel.noise = np.sqrt(self.pixel.background) # nJy
        self.pixel.noise_es = self.pixel.noise * self.nJy_to_es # convert from nJy to e/s


        if verbose:
            print('assumed aperture radius: {0:.2f} pix'.format(self.aperture.radius))
            print('noise in aperture: {0:.2f} nJy'.format(self.aperture.noise))
            print('noise in pixel: {0:.2f} nJy'.format(self.pixel.noise))
            print('zeropoint: {0}'.format(self.zeropoint))
            print('nJy_to_es: {0}'.format(self.nJy_to_es))

    def create_background_image(self, CutoutWidth):

        img = empty()
        img.nJy_to_es = self.nJy_to_es
        img.pixel_scale = self.pixel_scale
        img.noise = self.pixel.noise_es * np.ones((CutoutWidth,CutoutWidth))
        img.wht = 1./img.noise**2
        img.sci = self.pixel.noise_es * np.random.randn(*img.noise.shape)


        return img


def FieldBackground(Field, filter, verbose = False):

    aperture_f_limit = FLARE.photom.m_to_flux(Field.depths[filter]) # nJy
    aperture_significance = 5. # S/N in aperture
    aperture_radius = Field.depth_aperture_radius_arcsec/Field.pixel_scale # pix

    return Background(filter, Field.pixel_scale, aperture_f_limit, aperture_significance, aperture_radius, verbose = verbose)


def FieldBackgrounds(Field, verbose = False):

    return {filter:FieldBackground(Field, filter, verbose = verbose) for filter in Field.filters}
