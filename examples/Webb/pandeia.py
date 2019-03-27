# 
# # The following section is only needed if the PYSYN_CDBS environment variable is not already set.
# # The PYSYN_CDBS environment variable should point to the path of the CDBS data files
# import os
# location_of_cdbs = "/Users/stephenwilkins/Dropbox/Research/data/pysynphot"
# os.environ['PYSYN_CDBS'] = location_of_cdbs
# # End section
#   
# # The following section is only needed if the pandeia_refdata environment variable is not already set
# # The pandeia_refdata environment variable should point to the path of the pandeia reference data
# import os
# location_of_pandeia_refdata = "/path/to/pandeia/refdata"
# os.environ['pandeia_refdata'] = location_of_pandeia_refdata
# # End section
#  

try:
    from pandeia.engine.calc_utils import build_default_calc
except:
    print('failed')

from pandeia.engine.calc_utils import build_default_calc
from pandeia.engine.perform_calculation import perform_calculation
  
configuration = build_default_calc('jwst', 'nircam', 'sw_imaging')
scene = {}
scene['position'] = {'x_offset': 0., 'y_offset': 0., 'orientation': 0., 'position_parameters': ['x_offset', 'y_offset', 'orientation']}
scene['shape'] = {'geometry': 'point'}
scene['spectrum'] = {'name': 'Flat Source', 'spectrum_parameters': ['sed', 'normalization']}
scene['spectrum']['sed'] = {'sed_type': 'flat', 'unit': 'fnu'}
scene['spectrum']['normalization'] = {'type': 'at_lambda', 'norm_wave': 2.0, 'norm_waveunit': 'um', 'norm_flux': 1000., 'norm_fluxunit': 'njy'}
configuration['scene'][0] = scene
configuration['configuration']['instrument']['filter'] = 'f200w'
  
report = perform_calculation(configuration)
print(report['scalar']['sn'])
print(report['scalar']['exposure_time'])
print(report['scalar']['extracted_flux'])
print(report['scalar']['background'])