
import copy

class empty: pass

import FLARE.filters
import FLARE.observatories
import FLARE.photom


surveys  = {}

class survey():

    def __init__(self, name, data_dir = None, data_reference = None):

        self.name = name
        self.fields = {}
        self.data_reference = data_reference
        self.data_dir = data_dir

    def add_field(self, field_name, filters = [], depths = {}, area = None, mask_file = None, depth_aperture_radius_pixel = None):

        self.fields[field_name] = empty()
        self.fields[field_name].filters = filters
        self.fields[field_name].depths = depths
        self.fields[field_name].area = area
        self.fields[field_name].mask_file = mask_file
        self.fields[field_name].depth_aperture_radius_pixel = depth_aperture_radius_pixel



# ------------------------------------------------------------------------------------------
# ------------- Webb baseline
# --- sensitivities are from https://jwst-docs.stsci.edu/near-infrared-camera/nircam-predicted-performance/nircam-sensitivity


filters = [f'Webb.NIRCam.{f}' for f in ['F090W','F115W','F150W','F200W','F277W','F356W','F444W','F410M']]
depths = {}

SNR = 10
depths['Webb.NIRCam.F090W'] = 15.3/SNR # nJy
depths['Webb.NIRCam.F115W'] = 13.2/SNR # nJy
depths['Webb.NIRCam.F150W'] = 10.6/SNR # nJy
depths['Webb.NIRCam.F200W'] = 9.1/SNR # nJy
depths['Webb.NIRCam.F277W'] = 14.3/SNR # nJy
depths['Webb.NIRCam.F356W'] = 12.1/SNR # nJy
depths['Webb.NIRCam.F444W'] = 23.6/SNR # nJy
depths['Webb.NIRCam.F410M'] = 24.7/SNR # nJy

surveys['Webb10k'] = survey('Webb10k')
surveys['Webb10k'].add_field('base', filters, depths, depth_aperture_radius_pixel = 2.5)



# ------------------------------------------------------------------------------------------
# ------------- Euclid deep baseline


filters = FLARE.observatories.Euclid.filterIDs

depths = {}

SNR = 10

depths = {f: FLARE.photom.m_to_flux(26.)/5. for f in FLARE.observatories.Euclid.NISP.get_filter_IDs()}
depths['Euclid.VIS.VIS'] = FLARE.photom.m_to_flux(26.5)/10.

surveys['Euclid'] = survey('Euclid')
surveys['Euclid'].add_field('deep', filters, depths, depth_aperture_radius_pixel = 2.5)













depthmodel = {}

depthmodel['flat_fnu'] = {f: 0.0 for f in FLARE.filters.Webb + FLARE.filters.HST}

depthmodel['2_2_benchmark'] = {'JWST.NIRCAM.F070W': 27.01, 'JWST.NIRCAM.F090W': 27.43, 'JWST.NIRCAM.F115W': 27.59, 'JWST.NIRCAM.F150W': 27.8, 'JWST.NIRCAM.F200W': 27.95, 'JWST.NIRCAM.F277W': 27.7, 'JWST.NIRCAM.F356W': 27.76, 'JWST.NIRCAM.F444W': 27.52}
depthmodel['2_2_benchmark']['exp_time'] = 612.00

depthmodel['3_2_benchmark'] = {'JWST.NIRCAM.F070W': 27.55, 'JWST.NIRCAM.F090W': 27.93, 'JWST.NIRCAM.F115W': 28.08, 'JWST.NIRCAM.F150W': 28.28, 'JWST.NIRCAM.F200W': 28.45, 'JWST.NIRCAM.F277W': 28.08, 'JWST.NIRCAM.F356W': 28.15, 'JWST.NIRCAM.F444W': 27.9}
depthmodel['3_2_benchmark']['exp_time'] = 1041.47

depthmodel['4_2_benchmark'] = {'JWST.NIRCAM.F070W': 27.84, 'JWST.NIRCAM.F090W': 28.2, 'JWST.NIRCAM.F115W': 28.34, 'JWST.NIRCAM.F150W': 28.53, 'JWST.NIRCAM.F200W': 28.71, 'JWST.NIRCAM.F277W': 28.31, 'JWST.NIRCAM.F356W': 28.36, 'JWST.NIRCAM.F444W': 28.1}
depthmodel['4_2_benchmark']['exp_time'] = 1470.94



depthmodel['5_1_benchmark'] = {'JWST.NIRCAM.F070W': 27.64, 'JWST.NIRCAM.F090W': 27.99, 'JWST.NIRCAM.F115W': 28.13, 'JWST.NIRCAM.F150W': 28.32, 'JWST.NIRCAM.F200W': 28.5, 'JWST.NIRCAM.F277W': 28.06, 'JWST.NIRCAM.F356W': 28.13, 'JWST.NIRCAM.F444W': 27.86, 'JWST.NIRCAM.F410M': 27.54} # 5 groups = 944.84 seconds
depthmodel['5_1_benchmark']['exp_time'] = 944.84

depthmodel['5_2_benchmark'] = {'JWST.NIRCAM.F070W': 28.04, 'JWST.NIRCAM.F090W': 28.37, 'JWST.NIRCAM.F115W': 28.52, 'JWST.NIRCAM.F150W': 28.71, 'JWST.NIRCAM.F200W': 28.89, 'JWST.NIRCAM.F277W': 28.44, 'JWST.NIRCAM.F356W': 28.52, 'JWST.NIRCAM.F444W': 28.25, 'JWST.NIRCAM.F410M': 27.92}
depthmodel['5_2_benchmark']['exp_time'] = 1900.41

# -- similar to CEERS strategy (though CEERS does not have F070W, F090W and has double JWST.NIRCAM.F115W

depthmodel['5_3_benchmark'] = {'JWST.NIRCAM.F070W': 28.26, 'JWST.NIRCAM.F090W': 28.6, 'JWST.NIRCAM.F115W': 28.73, 'JWST.NIRCAM.F150W': 28.94, 'JWST.NIRCAM.F200W': 29.1, 'JWST.NIRCAM.F277W': 28.69, 'JWST.NIRCAM.F356W': 28.73, 'JWST.NIRCAM.F444W': 28.46, 'JWST.NIRCAM.F410M': 28.15}
depthmodel['5_3_benchmark']['exp_time'] = 2855.98

depthmodel['5_4_benchmark'] = {'JWST.NIRCAM.F070W': 28.43, 'JWST.NIRCAM.F090W': 28.76, 'JWST.NIRCAM.F115W': 28.91, 'JWST.NIRCAM.F150W': 29.1, 'JWST.NIRCAM.F200W': 29.28, 'JWST.NIRCAM.F277W': 28.83, 'JWST.NIRCAM.F356W': 28.89, 'JWST.NIRCAM.F444W': 28.62, 'JWST.NIRCAM.F410M': 28.31}
depthmodel['5_4_benchmark']['exp_time'] = 3811.55

depthmodel['5_5_benchmark'] = {'JWST.NIRCAM.F070W': 28.56, 'JWST.NIRCAM.F090W': 28.89, 'JWST.NIRCAM.F115W': 29.03, 'JWST.NIRCAM.F150W': 29.24, 'JWST.NIRCAM.F200W': 29.4, 'JWST.NIRCAM.F277W': 28.97, 'JWST.NIRCAM.F356W': 29.03, 'JWST.NIRCAM.F444W': 28.73, 'JWST.NIRCAM.F410M': 28.43}
depthmodel['5_5_benchmark']['exp_time'] = 4767.13

depthmodel['5_6_benchmark'] = {'JWST.NIRCAM.F070W': 28.64, 'JWST.NIRCAM.F090W': 28.99, 'JWST.NIRCAM.F115W': 29.12, 'JWST.NIRCAM.F150W': 29.31, 'JWST.NIRCAM.F200W': 29.5, 'JWST.NIRCAM.F277W': 29.05, 'JWST.NIRCAM.F356W': 29.12, 'JWST.NIRCAM.F444W': 28.84, 'JWST.NIRCAM.F410M': 28.53}
depthmodel['5_6_benchmark']['exp_time'] = 5712 # not exact

depthmodel['5_10_benchmark'] = {'JWST.NIRCAM.F070W': 29.02, 'JWST.NIRCAM.F090W': 29.31, 'JWST.NIRCAM.F115W': 29.45, 'JWST.NIRCAM.F150W': 29.66, 'JWST.NIRCAM.F200W': 29.85, 'JWST.NIRCAM.F277W': 29.36, 'JWST.NIRCAM.F356W': 29.4, 'JWST.NIRCAM.F444W': 29.12, 'JWST.NIRCAM.F410M': 28.81}
depthmodel['5_10_benchmark']['exp_time'] = 9544.99

depthmodel['5_20_benchmark'] = {'JWST.NIRCAM.F070W': 29.4, 'JWST.NIRCAM.F090W': 29.66, 'JWST.NIRCAM.F115W': 29.78, 'JWST.NIRCAM.F150W': 30.0, 'JWST.NIRCAM.F200W': 30.17, 'JWST.NIRCAM.F277W': 29.66, 'JWST.NIRCAM.F356W': 29.72, 'JWST.NIRCAM.F444W': 29.45, 'JWST.NIRCAM.F410M': 29.19}
depthmodel['5_20_benchmark']['exp_time'] = 19100.71



depthmodel['10_1_benchmark'] = {'JWST.NIRCAM.F070W': 28.13, 'JWST.NIRCAM.F090W': 28.43, 'JWST.NIRCAM.F115W': 28.58, 'JWST.NIRCAM.F150W': 28.76, 'JWST.NIRCAM.F200W': 28.94, 'JWST.NIRCAM.F277W': 28.46, 'JWST.NIRCAM.F356W': 28.52, 'JWST.NIRCAM.F444W': 28.23}
depthmodel['10_1_benchmark']['exp_time'] = 2018.51

depthmodel['47_1_benchmark'] = {'JWST.NIRCAM.F070W': 28.81, 'JWST.NIRCAM.F090W': 29.06, 'JWST.NIRCAM.F115W': 29.2, 'JWST.NIRCAM.F150W': 29.4, 'JWST.NIRCAM.F200W': 29.59, 'JWST.NIRCAM.F277W': 29.06, 'JWST.NIRCAM.F356W': 29.13, 'JWST.NIRCAM.F444W': 28.81}
depthmodel['47_1_benchmark']['exp_time'] = 9963.72



# --- CEERS strategy but using benchmark background
depthmodel['flat_time'] = {'JWST.NIRCAM.F070W': 28.26, 'JWST.NIRCAM.F090W': 28.6, 'JWST.NIRCAM.F115W': 28.73, 'JWST.NIRCAM.F150W': 28.94, 'JWST.NIRCAM.F200W': 29.1, 'JWST.NIRCAM.F277W': 28.69, 'JWST.NIRCAM.F356W': 28.73, 'JWST.NIRCAM.F444W': 28.46, 'JWST.NIRCAM.F410M': 28.15}

# --- CEERS strategy but using benchmark background
depthmodel['CEERS_benchmark'] = {'JWST.NIRCAM.F115W': 29.12, 'JWST.NIRCAM.F150W': 28.94, 'JWST.NIRCAM.F200W': 29.1, 'JWST.NIRCAM.F277W': 28.69, 'JWST.NIRCAM.F356W': 28.73, 'JWST.NIRCAM.F444W': 28.46, 'JWST.NIRCAM.F410M': 28.15}

# --- CEERS strategy but swapping F410M for more F444W
depthmodel['CEERS-F410M_benchmark'] = {'JWST.NIRCAM.F115W': 29.12, 'JWST.NIRCAM.F150W': 28.94, 'JWST.NIRCAM.F200W': 29.1, 'JWST.NIRCAM.F277W': 28.69, 'JWST.NIRCAM.F356W': 28.73, 'JWST.NIRCAM.F444W': 28.84}


depthmodel['5_1_HUDF'] = {'JWST.NIRCAM.F070W': 27.69, 'JWST.NIRCAM.F090W': 28.05, 'JWST.NIRCAM.F115W': 28.23, 'JWST.NIRCAM.F150W': 28.43, 'JWST.NIRCAM.F200W': 28.62, 'JWST.NIRCAM.F277W': 28.26, 'JWST.NIRCAM.F356W': 28.29, 'JWST.NIRCAM.F444W': 27.85}
depthmodel['5_1_HUDF']['exp_time'] = 944.84

depthmodel['5_2_HUDF'] = {'JWST.NIRCAM.F070W': 28.09, 'JWST.NIRCAM.F090W': 28.44, 'JWST.NIRCAM.F115W': 28.62, 'JWST.NIRCAM.F150W': 28.83, 'JWST.NIRCAM.F200W': 29.03, 'JWST.NIRCAM.F277W': 28.64, 'JWST.NIRCAM.F356W': 28.67, 'JWST.NIRCAM.F444W': 28.23}
depthmodel['5_2_HUDF']['exp_time'] = 1900.41

depthmodel['5_3_HUDF'] = {'JWST.NIRCAM.F070W': 28.31, 'JWST.NIRCAM.F090W': 28.67, 'JWST.NIRCAM.F115W': 28.86, 'JWST.NIRCAM.F150W': 29.06, 'JWST.NIRCAM.F200W': 29.24, 'JWST.NIRCAM.F277W': 28.86, 'JWST.NIRCAM.F356W': 28.89, 'JWST.NIRCAM.F444W': 28.44}
depthmodel['5_3_HUDF']['exp_time'] = 2855.98

depthmodel['5_4_HUDF'] = {'JWST.NIRCAM.F070W': 28.48, 'JWST.NIRCAM.F090W': 28.83, 'JWST.NIRCAM.F115W': 29.0, 'JWST.NIRCAM.F150W': 29.2, 'JWST.NIRCAM.F200W': 29.4, 'JWST.NIRCAM.F277W': 29.03, 'JWST.NIRCAM.F356W': 29.06, 'JWST.NIRCAM.F444W': 28.6}
depthmodel['5_4_HUDF']['exp_time'] = 3811.55

depthmodel['5_5_HUDF'] = {'JWST.NIRCAM.F070W': 28.6, 'JWST.NIRCAM.F090W': 28.94, 'JWST.NIRCAM.F115W': 29.13, 'JWST.NIRCAM.F150W': 29.32, 'JWST.NIRCAM.F200W': 29.54, 'JWST.NIRCAM.F277W': 29.17, 'JWST.NIRCAM.F356W': 29.17, 'JWST.NIRCAM.F444W': 28.73}
depthmodel['5_5_HUDF']['exp_time'] = 4767.13

depthmodel['5_10_HUDF'] = {'JWST.NIRCAM.F070W': 29.0, 'JWST.NIRCAM.F090W': 29.33, 'JWST.NIRCAM.F115W': 29.54, 'JWST.NIRCAM.F150W': 29.73, 'JWST.NIRCAM.F200W': 29.88, 'JWST.NIRCAM.F277W': 29.54, 'JWST.NIRCAM.F356W': 29.54, 'JWST.NIRCAM.F444W': 29.11}
depthmodel['5_10_HUDF']['exp_time'] = 9544.99

depthmodel['5_20_HUDF'] = {'JWST.NIRCAM.F070W': 29.38, 'JWST.NIRCAM.F090W': 29.73, 'JWST.NIRCAM.F115W': 29.88, 'JWST.NIRCAM.F150W': 30.06, 'JWST.NIRCAM.F200W': 30.26, 'JWST.NIRCAM.F277W': 29.88, 'JWST.NIRCAM.F356W': 29.96, 'JWST.NIRCAM.F444W': 29.48}
depthmodel['5_20_HUDF']['exp_time'] = 19100.71









# ------------------------------------------------------------------------------------------
# ------------- XDF

surveys['XDF'] = empty()
surveys['XDF'].fields = {}
surveys['XDF'].data_reference = ''
surveys['XDF'].datadir = 'data/images/hubble/xdf'

# --- entire field
surveys['XDF'].fields['XDF'] = empty()
surveys['XDF'].fields['XDF'].datadir = surveys['XDF'].datadir
surveys['XDF'].fields['XDF'].mask_file = False
surveys['XDF'].fields['XDF'].depth_aperture_radius_arcsec = (0.35/2.)
# surveys['XDF'].fields['XDF'].filters = ['HST.ACS.{0}'.format(f) for f in ['f435w', 'f606w', 'f775w', 'f814w', 'f850lp']] + ['HST.WFC3.{0}'.format(f) for f in ['f105w', 'f125w', 'f140w', 'f160w']]
surveys['XDF'].fields['XDF'].filters = ['HST.ACS.{0}'.format(f) for f in ['f435w', 'f606w', 'f775w', 'f850lp']] + ['HST.WFC3.{0}'.format(f) for f in ['f105w', 'f125w', 'f140w', 'f160w']]
surveys['XDF'].fields['XDF'].pixel_scale = 0.06
surveys['XDF'].fields['XDF'].filename = {}
for filter in surveys['XDF'].fields['XDF'].filters:
    obs, inst, f = filter.lower().split('.')
    if inst == 'wfc3': INST = 'wfc3ir'
    if inst == 'acs': INST = 'acswfc'
    surveys['XDF'].fields['XDF'].filename[filter] = surveys['XDF'].datadir + '/hlsp_xdf_hst_{0}-60mas_hudf_{1}_v1'.format(INST, f)
surveys['XDF'].fields['XDF'].literature_depths = False

# --- deepest sub-field
surveys['XDF'].fields['dXDF'] = copy.deepcopy(surveys['XDF'].fields['XDF'])
surveys['XDF'].fields['dXDF'].mask_file = 'hlsp_xdf_hst_deepest_flag_v1.fits'
surveys['XDF'].fields['dXDF'].area = 4.04 # arcmin
surveys['XDF'].fields['dXDF'].literature_depths = {'HST.ACS.f435w': 29.8,'HST.ACS.f606w':30.3,'HST.ACS.f775w':30.3,'HST.ACS.f814w':29.1,'HST.ACS.f850lp':29.4,'HST.WFC3.f105w':30.1,'HST.WFC3.f125w':29.8,'HST.WFC3.f140w':29.8,'HST.WFC3.f160w':29.8}
surveys['XDF'].fields['dXDF'].depths = {'HST.ACS.f435w': 29.8,'HST.ACS.f606w':30.3,'HST.ACS.f775w':30.3,'HST.ACS.f814w':29.1,'HST.ACS.f850lp':29.4,'HST.WFC3.f105w':30.1,'HST.WFC3.f125w':29.8,'HST.WFC3.f140w':29.8,'HST.WFC3.f160w':29.8}



# ------------------------------------------------------------------------------------------
# ------------- HLF







# ------------------------------------------------------------------------------------------
# ------------- HFF






# ------------------------------------------------------------------------------------------
# ------------- CANDELS















# ------------------------------------------------------------------------------------------
# ------------- Webb hypothetical


surveys['WebbHypothetical'] = empty()
surveys['WebbHypothetical'].fields = {}
surveys['WebbHypothetical'].data_reference = ''
surveys['WebbHypothetical'].datadir = ''

# --- entire field
surveys['WebbHypothetical'].fields['deep'] = empty()
surveys['WebbHypothetical'].fields['deep'].datadir = surveys['WebbHypothetical'].datadir
surveys['WebbHypothetical'].fields['deep'].mask_file = False
surveys['WebbHypothetical'].fields['deep'].depth_aperture_radius_arcsec = (0.35/2.)
surveys['WebbHypothetical'].fields['deep'].filters = [f'JWST.NIRCAM.{f}' for f in ['F090W','F115W','F150W','F200W','F277W','F356W','F444W','F410M']]
surveys['WebbHypothetical'].fields['deep'].pixel_scale = 0.031 # assumes drizzled to NIRCam_short
surveys['WebbHypothetical'].fields['deep'].depths = {f:depthmodel['5_10_benchmark'][f] for f in surveys['WebbHypothetical'].fields['deep'].filters}



# ------------------------------------------------------------------------------------------
# ------------- CEERS

surveys['CEERS'] = empty()
surveys['CEERS'].fields = {}
surveys['CEERS'].data_reference = ''
surveys['CEERS'].datadir = ''

# --- entire field
surveys['CEERS'].fields['EGS'] = empty()
surveys['CEERS'].fields['EGS'].datadir = surveys['CEERS'].datadir
surveys['CEERS'].fields['EGS'].mask_file = False
surveys['CEERS'].fields['EGS'].depth_aperture_radius_arcsec = (0.031*2.5/2.)
surveys['CEERS'].fields['EGS'].filters = [f'JWST.NIRCAM.{f}' for f in ['F115W','F150W','F200W','F277W', 'F356W','F444W', 'F410M']]
surveys['CEERS'].fields['EGS'].pixel_scale = 0.031 # assumes drizzled to NIRCam_short
surveys['CEERS'].fields['EGS'].depths = {'JWST.NIRCAM.F115W': 29.15, 'JWST.NIRCAM.F150W': 28.90, 'JWST.NIRCAM.F200W': 28.95, 'JWST.NIRCAM.F356W': 28.95, 'JWST.NIRCAM.F444W': 28.96}
