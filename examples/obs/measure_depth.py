
import sys
import os

from astropy.io import fits
import numpy as np

sys.path.insert(0, os.path.abspath('../../'))

import FLARE
import FLARE.obs
import FLARE.surveys
import FLARE.photom

survey = FLARE.surveys.XDF

for fieldID, field in survey.fields.items():

    print('-'*5, fieldID)

    if field.mask_file:
        mask = fits.getdata('{0}/{1}/{2}'.format(FLARE.FLARE_dir, survey.datadir, field.mask_file)) 
    else:
        mask = False

    for filter in field.filters:

        img = FLARE.obs.image_from_file(field.filename[filter], filter, mask = mask)

        depth = FLARE.photom.flux_to_m(img.determine_depth())

        if filter == field.filters[0]: print('area/arcmin2: {0:.2f}'.format(img.get_area()))

        if field.literature_depths:
            print('{0}: {1:.1f} | {2:.1f})'.format(filter, depth, depth - field.literature_depths[filter])) # 5\sigma limiting magnitude
        else:
            print('{0}: {1:.1f})'.format(filter, depth)) # 5\sigma limiting magnitude


