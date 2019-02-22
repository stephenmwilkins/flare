
import numpy as np

import os

try:
    FLARE_dir = os.environ['FLARE']
    print('FLARE = '+FLARE_dir)
except:
    print('WARNING: FLARE environment variable not set')




def default_cosmo():

    from astropy.cosmology import WMAP9 as cosmo

    return cosmo