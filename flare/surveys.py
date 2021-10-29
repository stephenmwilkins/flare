
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

    def add_field(self, field_name, data_dir = None, data_reference = None, filters = None, depths = None, depths_mag = None, depth_reference_filter = None, area = None, mask_file = None, depth_aperture_radius_arcsec = None, depth_aperture_radius_pixel = None, depth_aperture_significance = None, pixel_scale = None, detection_filters = None):

        self.fields[field_name] = Field(field_name, self.name, data_dir = data_dir, data_reference = data_reference, filters = filters, depths = depths, depths_mag = depths_mag, depth_reference_filter = depth_reference_filter, area = area, mask_file = mask_file, depth_aperture_radius_arcsec = depth_aperture_radius_arcsec, depth_aperture_radius_pixel = depth_aperture_radius_pixel, depth_aperture_significance = depth_aperture_significance, pixel_scale = pixel_scale, detection_filters = detection_filters)

        pass



class Field:

    def __init__(self, field_name, survey_name, data_dir = None, data_reference = None, filters = None, depths = None, depths_mag = None, depth_reference_filter = None, area = None, mask_file = None, depth_aperture_radius_arcsec = None, depth_aperture_radius_pixel = None, depth_aperture_significance = None, pixel_scale = None, detection_filters = None):

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

        self.depth_reference_filter = depth_reference_filter

        self.area = area # --- arcsec2
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


# --- Euclid

Euclid = Survey('Euclid')

filters = ['Euclid.VIS.VIS'] + [f'Euclid.NISP.{f}' for f in ['Y','J','H']]

Euclid.add_field('Deep', filters = filters, depths_mag = {'Euclid.VIS.VIS': 26, 'Euclid.NISP.Y': 26, 'Euclid.NISP.J': 26, 'Euclid.NISP.H': 26}, area = 40*3600, depth_reference_filter = 'Euclid.NISP.H')
Euclid.add_field('Wide', filters = filters, depths_mag = {'Euclid.VIS.VIS': 24, 'Euclid.NISP.Y': 24, 'Euclid.NISP.J': 24, 'Euclid.NISP.H': 24}, area = 15000*3600, depth_reference_filter = 'Euclid.NISP.H')

surveys['Euclid'] = Euclid



# ------------------------------------------------------------
# --- Webb Cycle 1

# --- JADES
JADES = Survey('JADES')
JADES.add_field('deep', depths_mag = {'Webb.NIRCam.F150W': 30.7, 'Webb.NIRCam.F200W': 30.7, 'Webb.NIRCam.F356W': 30.2}, area = 46, depth_reference_filter = 'Webb.NIRCam.F150W')
JADES.add_field('medium', depths_mag = {'Webb.NIRCam.F150W': 29.7, 'Webb.NIRCam.F200W': 29.8, 'Webb.NIRCam.F356W': 29.4}, area = 144, depth_reference_filter = 'Webb.NIRCam.F150W')
surveys['JADES'] = JADES

# --- PANORAMIC
surveys['PANORAMIC'] = Survey('JADES').add_field('', depths_mag = {'Webb.NIRCam.F115W': 28.1, 'Webb.NIRCam.F150W':28.29, 'Webb.NIRCam.F200W': 28.46, 'Webb.NIRCam.F356W': 28.03, 'Webb.NIRCam.F356W': 28.09, 'Webb.NIRCam.F444W': 27.81}, area = 150*9.1, depth_reference_filter = 'Webb.NIRCam.F150W')

# --- CEERS
surveys['CEERS'] = Survey('CEERS').add_field('', depths_mag = {'Webb.NIRCam.F150W': 28.9, 'Webb.NIRCam.F200W': 28.95, 'Webb.NIRCam.F356W': 28.95}, area = 150*9.1, depth_reference_filter = 'Webb.NIRCam.F150W')

# --- COSMOS-Web
surveys['COSMOS-Web'] = Survey('COSMOS-Web').add_field('', depths_mag = {'Webb.NIRCam.F115W': 27.36, 'Webb.NIRCam.F150W': 27.60, 'Webb.NIRCam.F277W': 28.05, 'Webb.NIRCam.F444W': 27.72}, area = 3600*0.6, depth_reference_filter = 'Webb.NIRCam.F150W')

# --- PRIMER
PRIMER = Survey('PRIMER')
PRIMER.add_field('COSMOS-Shallow', depths_mag = {'Webb.NIRCam.F090W': 28.33 ,'Webb.NIRCam.F115W': 28.61, 'Webb.NIRCam.F150W': 28.81, 'Webb.NIRCam.F200W': 28.89, 'Webb.NIRCam.F277W': 28.85, 'Webb.NIRCam.F356W': 28.79, 'Webb.NIRCam.F410M': 28.04, 'Webb.NIRCam.F444W': 28.32}, area = 144.2, depth_reference_filter = 'Webb.NIRCam.F150W')
PRIMER.add_field('COSMOS-Medium', depths_mag = {'Webb.NIRCam.F090W': 28.57,'Webb.NIRCam.F115W': 28.84, 'Webb.NIRCam.F150W': 29.03, 'Webb.NIRCam.F200W': 29.11, 'Webb.NIRCam.F277W': 29.08, 'Webb.NIRCam.F356W': 29.02, 'Webb.NIRCam.F410M': 28.27, 'Webb.NIRCam.F444W': 28.54}, area = 108.3, depth_reference_filter = 'Webb.NIRCam.F150W')
PRIMER.add_field('COSMOS-Deep', depths_mag = {'Webb.NIRCam.F090W': 28.96,'Webb.NIRCam.F115W': 29.23, 'Webb.NIRCam.F150W': 29.43, 'Webb.NIRCam.F200W': 29.51, 'Webb.NIRCam.F277W': 29.47, 'Webb.NIRCam.F356W': 29.42, 'Webb.NIRCam.F410M': 28.65, 'Webb.NIRCam.F444W': 28.93}, area =  33.4, depth_reference_filter = 'Webb.NIRCam.F150W')
PRIMER.add_field('UDS-Shallow', depths_mag = {'Webb.NIRCam.F090W': 27.92,'Webb.NIRCam.F115W': 28.21, 'Webb.NIRCam.F150W': 28.42, 'Webb.NIRCam.F200W': 28.48, 'Webb.NIRCam.F277W': 28.38, 'Webb.NIRCam.F356W': 28.37, 'Webb.NIRCam.F410M': 27.57, 'Webb.NIRCam.F444W': 27.79}, area = 234.02, depth_reference_filter = 'Webb.NIRCam.F150W')
PRIMER.add_field('UDS-Medium', depths_mag = {'Webb.NIRCam.F090W': 28.32,'Webb.NIRCam.F115W': 28.62, 'Webb.NIRCam.F150W': 28.82, 'Webb.NIRCam.F200W': 28.69, 'Webb.NIRCam.F277W': 28.85, 'Webb.NIRCam.F356W': 28.77, 'Webb.NIRCam.F410M': 27.97, 'Webb.NIRCam.F444W': 28.21}, area = 175.17, depth_reference_filter = 'Webb.NIRCam.F150W')

surveys['PRIMER'] = PRIMER

# --- WDEEP
surveys['WDEEP'] = Survey('COSMOS-Web').add_field('', depths_mag = {'Webb.NIRCam.F115W': 30.90, 'Webb.NIRCam.F150W': 30.62, 'Webb.NIRCam.F200W': 30.62, 'Webb.NIRCam.F277W': 30.72, 'Webb.NIRCam.F356W': 30.70, 'Webb.NIRCam.F444W': 30.56}, area = 2*9.1, depth_reference_filter = 'Webb.NIRCam.F150W')
