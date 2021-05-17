
import copy

class empty: pass

import flare
import flare.filters
import flare.observatories
import flare.photom


surveys  = {}


class Survey:

    def __init__(self, name):

        self.name = name
        self.fields = {}

    def add_field(self, field_name, data_dir = None, data_reference = None, filters = None, depths = None, depths_mag = None, area = None, mask_file = None, depth_aperture_radius_arcsec = None, depth_aperture_radius_pixel = None, depth_aperture_significance = None, pixel_scale = None, detection_filters = None):

        self.fields[field_name] = Field(field_name, self.name, data_dir = data_dir, data_reference = data_reference, filters = filters, depths = depths, depths_mag = depths_mag, area = area, mask_file = mask_file, depth_aperture_radius_arcsec = depth_aperture_radius_arcsec, depth_aperture_radius_pixel = depth_aperture_radius_pixel, depth_aperture_significance = depth_aperture_significance, pixel_scale = pixel_scale, detection_filters = detection_filters)

        pass



class Field:

    def __init__(self, field_name, survey_name, data_dir = None, data_reference = None, filters = None, depths = None, depths_mag = None, area = None, mask_file = None, depth_aperture_radius_arcsec = None, depth_aperture_radius_pixel = None, depth_aperture_significance = None, pixel_scale = None, detection_filters = None):

        # --- set by field

        self.survey_name = survey_name

        self.name = field_name
        self.data_dir = data_dir
        self.data_reference = data_reference
        self.filters = filters

        if depths:
            self.depths = depths
            self.depths_mag = {f: flare.photom.flux_to_m(m) for f, m in self.depths.items()}
        elif depths_mag:
            self.depths_mag = depths_mag
            self.depths = {f: flare.photom.m_to_flux(m) for f, m in self.depths_mag.items()}
        else:
            print('WARNING! No depths have been set')

        self.area = area
        self.mask_file = mask_file
        self.depth_aperture_radius_arcsec = depth_aperture_radius_arcsec
        self.depth_aperture_radius_pixel = depth_aperture_radius_pixel
        self.depth_aperture_significance = depth_aperture_significance
        self.pixel_scale = pixel_scale
        self.detection_filters = detection_filters









# ------------------------------------------------------------------------------------------
# --- Define a simple Hubble consisting of some simple fields


surveys['SimpleHubble'] = Survey('SimpleHubble')

filters = [f'Hubble.ACS.{f}' for f in ['f435w', 'f606w', 'f775w', 'f814w']] + [f'Hubble.WFC3.{f}' for f in ['f105w', 'f125w', 'f160w']]

# --- define 10nJy field

depths = {f: 10. for f in filters} # nJy
surveys['SimpleHubble'].add_field('10nJy', filters = filters, depths = depths, depth_aperture_radius_arcsec = 0.35/2., depth_aperture_significance = 5, pixel_scale = 0.06)


# ------------------------------------------------------------------------------------------
# ------------- Define the test data created by pysep

surveys['test'] = Survey('test')

# --- add test field

filters = [f'Hubble.ACS.{f}' for f in ['f435w', 'f606w', 'f775w', 'f814w', 'f850lp']] + [f'Hubble.WFC3.{f}' for f in ['f105w', 'f125w', 'f140w', 'f160w']]
data_dir = 'test_data/'
depths_mag = {'Hubble.ACS.f435w': 29.8,'Hubble.ACS.f606w':30.3,'Hubble.ACS.f775w':30.3,'Hubble.ACS.f814w':29.1,'Hubble.ACS.f850lp':29.4,'Hubble.WFC3.f105w':30.1,'Hubble.WFC3.f125w':29.8,'Hubble.WFC3.f140w':29.8,'Hubble.WFC3.f160w':29.8}
depth_aperture_radius_arcsec = 0.35/2.
depth_aperture_significance = 5
pixel_scale = 0.06
detection_filters = ['Hubble.WFC3.f105w','Hubble.WFC3.f125w','Hubble.WFC3.f140w','Hubble.WFC3.f160w']

surveys['test'].add_field('test', data_dir = data_dir, filters = filters, depths_mag = depths_mag, depth_aperture_radius_arcsec = depth_aperture_radius_arcsec, detection_filters = detection_filters, depth_aperture_significance = depth_aperture_significance, pixel_scale = pixel_scale)






# ------------------------------------------------------------------------------------------
# ------------- Define the XDF survey. This consists of two fields, though one is a subset of the other defined by a mask

XDF = Survey('XDF')


# --- entire XDF


# --- deepest sub-field

filters = [f'Hubble.ACS.{f}' for f in ['f435w', 'f606w', 'f775w', 'f814w', 'f850lp']] + [f'Hubble.WFC3.{f}' for f in ['f105w', 'f125w', 'f140w', 'f160w']]
filters = [f'Hubble.ACS.{f}' for f in ['f435w', 'f606w', 'f775w', 'f850lp']] + [f'Hubble.WFC3.{f}' for f in ['f105w', 'f125w', 'f140w', 'f160w']]


data_dir = flare.FLARE_dir + '/data/images/hubble/xdf'
depths_mag = {'Hubble.ACS.f435w': 29.8,'Hubble.ACS.f606w':30.3,'Hubble.ACS.f775w':30.3,'Hubble.ACS.f814w':29.1,'Hubble.ACS.f850lp':29.4,'Hubble.WFC3.f105w':30.1,'Hubble.WFC3.f125w':29.8,'Hubble.WFC3.f140w':29.8,'Hubble.WFC3.f160w':29.8}
depth_aperture_radius_arcsec = 0.35/2.
depth_aperture_significance = 5
pixel_scale = 0.06
mask_file = 'hlsp_xdf_hst_deepest_flag_v1.fits'
detection_filters = ['Hubble.WFC3.f105w','Hubble.WFC3.f125w','Hubble.WFC3.f140w','Hubble.WFC3.f160w']

XDF.add_field('dXDF', data_dir = data_dir, filters = filters, depths_mag = depths_mag, depth_aperture_radius_arcsec = depth_aperture_radius_arcsec, mask_file = mask_file, detection_filters = detection_filters, depth_aperture_significance = depth_aperture_significance, pixel_scale = pixel_scale)

surveys['XDF'] = XDF
