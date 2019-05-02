 



import numpy as np

import FLARE.filters







run = True

lw_filters = FLARE.filters.NIRCam_l
sw_filters = FLARE.filters.NIRCam_s

filters = sw_filters + lw_filters

radius_pix = 2.5


try:
    from pandeia.engine.calc_utils import build_default_calc
except:
    print('failed')

from pandeia.engine.calc_utils import build_default_calc
from pandeia.engine.perform_calculation import perform_calculation

config = {}

config['background'] = 'minzodi'
config['background_level'] = 'benchmark'

config['configuration'] = {}
config['configuration']['detector'] = {'nexp': 1, 'ngroup':47, 'nint':1, 'readmode': 'deep8', 'subarray': 'full'}
config['configuration']['instrument'] = {'aperture': 'sw', 'instrument': 'nircam', 'mode': 'sw_imaging', 'disperser': None}

scene = {}
scene['position'] = {'x_offset': 0., 'y_offset': 0., 'orientation': 0., 'position_parameters': ['x_offset', 'y_offset', 'orientation']}
scene['shape'] = {'geometry': 'point'}
scene['spectrum'] = {'name': 'Flat Source', 'spectrum_parameters': ['sed', 'normalization']}
scene['spectrum']['sed'] = {'sed_type': 'flat', 'unit': 'fnu'}
scene['spectrum']['normalization'] = {'type': 'at_lambda', 'norm_wave': 2.0, 'norm_waveunit': 'um', 'norm_flux': 10., 'norm_fluxunit': 'njy'}
scene['spectrum']['normalization'] ['norm_flux'] = 100. # we don't actually care. We just want the background calculations.
config['scene'] = [scene]

config['strategy'] = {'background_subtraction': False, 'method':'imagingapphot', 'target_source':'1', 'target_type':'source', 'target_xy':[0.0,0.0], 'units':'arcsec'}



for groups in [47]:

    config['configuration']['detector']['ngroup'] = groups

    background = {}

    for f in filters:

        config['configuration']['instrument']['filter'] = f.split('.')[-1].lower()

        if f in sw_filters:

            config['configuration']['instrument']['aperture'] = 'sw'
            config['configuration']['instrument']['mode'] = 'sw_imaging'
            config['strategy']['aperture_size'] = 0.031*radius_pix

        if f in lw_filters:

            config['configuration']['instrument']['aperture'] = 'lw'
            config['configuration']['instrument']['mode'] = 'lw_imaging'
            config['strategy']['aperture_size'] = 0.063*radius_pix
    
    
        report = perform_calculation(config)
            
        background_total = report['scalar']['background_total'] * report['scalar']['total_exposure_time']
    
        background_total_pix = background_total/(np.pi*radius_pix**2)
    
    
        if f == filters[0]: 
            print('-'*10)
            print('groups: {0}; t_exp: {1:.2f}'.format(groups, report['scalar']['total_exposure_time']))
        
        print(f, '{0:.2f}'.format(background_total), '{0:.2f}'.format(background_total_pix))

        background[f] = background_total_pix

    print(background)