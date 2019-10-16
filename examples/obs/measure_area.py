
import sys
import os

from astropy.io import fits
import numpy as np

sys.path.insert(0, os.path.abspath('../../'))

import FLARE
import FLARE.obs
import FLARE.surveys
import FLARE.photom

survey = FLARE.surveys.surveys['XDF']

for fieldID, field in survey.fields.items():

    print('-'*5, fieldID)

    if field.mask_file:
        mask = fits.getdata('{0}/{1}/{2}'.format(FLARE.FLARE_dir, survey.datadir, field.mask_file)) 
    else:
        mask = False

    if mask is not False:

        npix = len(mask[mask==0])

        area_arcsec = npix * field.pixel_scale**2
        area_arcmin = area_arcsec / 3600.


        print(mask)
        print(mask.shape)
        print(npix)
        print(area_arcmin)