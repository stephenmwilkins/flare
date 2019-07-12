
import sys
import os
import numpy as np

sys.path.insert(0, os.path.abspath('../../'))

import FLARE
import FLARE.filters
from FLARE.interrogator import models




filters = FLARE.filters.HST
# filters = ['HST.WFC3.f105w']


SPS = 'BPASSv2.2.1.binary'
IMF =  'ModSalpeter_300'

# models.generate_photometry_grid(SPS, IMF, filters, redshifts = np.arange(0,5,1))
# models.generate_photometry_grid(SPS, IMF, filters, redshifts = [8.])
models.generate_photometry_grid(SPS, IMF, filters)


