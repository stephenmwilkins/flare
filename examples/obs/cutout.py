
import sys
import os

from astropy.io import fits
import numpy as np

np.random.seed(42)

sys.path.insert(0, os.path.abspath('../../'))

import FLARE
import FLARE.obs
import FLARE.surveys
import FLARE.photom


field = FLARE.surveys.XDF.fields['dXDF'] # deepest portion of the XDF

filter = 'HST.WFC3.f160w'

img = FLARE.obs.open_image(field, filter, verbose = True) # opens image and applies mask. img is actually an instance of the obs.core.image class

x,y = img.get_random_location()

width = 200

print(x,y)

cutout = img.make_cutout(x,y, width) 

cutout.make_significance_plot()

x,y = 10,10
width = 40
cutout2 = cutout.make_cutout(x,y, width) 
cutout2.make_significance_plot()

x,y = 190,190
width = 40

cutout3 = cutout.make_cutout(x,y, width) 
cutout3.make_significance_plot()