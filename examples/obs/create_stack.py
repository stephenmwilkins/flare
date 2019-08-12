
import sys
import os

from astropy.io import fits
import numpy as np



sys.path.insert(0, os.path.abspath('../../'))

import FLARE
import FLARE.obs
import FLARE.surveys
import FLARE.photom


field = FLARE.surveys.XDF.fields['dXDF'] # deepest portion of the XDF

filters = ['HST.WFC3.f125w','HST.WFC3.f140w','HST.WFC3.f160w']

imgs = FLARE.obs.open_images(field, filters, verbose = True) # opens image and applies mask. img is actually an instance of the obs.core.image class

x,y = next(iter(imgs.values())).get_random_location()
width = 200

print(x,y)

cutouts = {filter: imgs[filter].make_cutout(x,y, width) for filter in filters} 

stack = FLARE.obs.create_stack(cutouts)

FLARE.obs.make_significance_plots({**cutouts, **{'stack':stack}})