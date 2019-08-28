
import copy

class empty: pass

import FLARE.filters

surveys = {}



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
surveys['XDF'].fields['XDF'].filters = ['HST.ACS.{0}'.format(f) for f in ['f435w', 'f606w', 'f775w', 'f814w', 'f850lp']] + ['HST.WFC3.{0}'.format(f) for f in ['f105w', 'f125w', 'f140w', 'f160w']]
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
surveys['WebbHypothetical'].fields['deep'].filters = FLARE.filters.NIRCam_W
surveys['WebbHypothetical'].fields['deep'].pixel_scale = 0.031 # assumes drizzled to NIRCam_short
surveys['WebbHypothetical'].fields['deep'].depths = {f:31. for f in surveys['WebbHypothetical'].fields['deep'].filters}





# surveys['CEERS'] = empty()
# surveys['CEERS'].fields = {}
# surveys['CEERS'].data_reference = ''
# surveys['CEERS'].datadir = ''
# 
# # --- entire field
# surveys['CEERS'].fields['EGS'] = empty()
# surveys['CEERS'].fields['EGS'].datadir = surveys['EGS'].datadir
# surveys['CEERS'].fields['EGS'].mask_file = False
# surveys['CEERS'].fields['EGS'].depth_aperture_radius_arcsec = (0.35/2.)
# surveys['CEERS'].fields['EGS'].filters = ['HST.ACS.{0}'.format(f) for f in ['f435w', 'f606w', 'f775w', 'f814w', 'f850lp']] + ['HST.WFC3.{0}'.format(f) for f in ['f105w', 'f125w', 'f140w', 'f160w']]
# surveys['CEERS'].fields['EGS'].pixel_scale = 0.031 # assumes drizzled to NIRCam_short
# surveys['CEERS'].fields['EGS'].depths = {'HST.ACS.f435w': 29.8,'HST.ACS.f606w':30.3,'HST.ACS.f775w':30.3,'HST.ACS.f814w':29.1,'HST.ACS.f850lp':29.4,'HST.WFC3.f105w':30.1,'HST.WFC3.f125w':29.8,'HST.WFC3.f140w':29.8,'HST.WFC3.f160w':29.8}
# 

