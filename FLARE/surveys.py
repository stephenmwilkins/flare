
import copy

class empty: pass



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






