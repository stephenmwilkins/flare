 

try:
    from pandeia.engine.calc_utils import build_default_calc
except:
    print('failed')

from pandeia.engine.calc_utils import build_default_calc
from pandeia.engine.perform_calculation import perform_calculation

import numpy as np

# import FLARE.filters



filter = 'f200w'

aperture = 'sw'
mode = 'sw_imaging'

# aperture = 'lw'
# mode = 'lw_imaging'

norm_flux = 100 # nJy
aperture_size = 0.031*2.5 # arcsec
# aperture_size = 0.063*2.5



config = {}

config['background'] = 'minzodi'
config['background_level'] = 'benchmark'

config['configuration'] = {}
config['configuration']['detector'] = {'nexp': 1, 'ngroup':5, 'nint':1, 'readmode': 'deep8', 'subarray': 'full'}
config['configuration']['instrument'] = {'aperture': aperture, 'instrument': 'nircam', 'mode': mode, 'disperser': None, 'filter': filter}

scene = {}
scene['position'] = {'x_offset': 0., 'y_offset': 0., 'orientation': 0., 'position_parameters': ['x_offset', 'y_offset', 'orientation']}
scene['shape'] = {'geometry': 'point'}
scene['spectrum'] = {'name': 'Flat Source', 'spectrum_parameters': ['sed', 'normalization']}
scene['spectrum']['sed'] = {'sed_type': 'flat', 'unit': 'fnu'}
scene['spectrum']['normalization'] = {'type': 'at_lambda', 'norm_wave': 2.0, 'norm_waveunit': 'um', 'norm_flux': norm_flux, 'norm_fluxunit': 'njy'}
config['scene'] = [scene]

config['strategy'] = {"aperture_size":aperture_size, "background_subtraction": False, "method":"imagingapphot", "target_source":"1", "target_type":"source", "target_xy":[0.0,0.0], "units":"arcsec"}

report = perform_calculation(config)
print('S/N', '{0:.2f}'.format(report['scalar']['sn']))
print('-'*30)
for k in report['scalar'].keys():
    v = report['scalar'][k]
    if isinstance(v, float):
        print(k, '{0:.2f}'.format(v))
   