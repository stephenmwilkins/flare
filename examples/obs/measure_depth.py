


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

    mask = fits.getdata('{0}/{1}/hlsp_xdf_hst_deepest_flag_v1.fits'.format(FLARE.FLARE_dir, survey.datadir)) 

    for filter in field.filters:

        img = FLARE.obs.image(field.filename[filter], filter, mask = mask)

        depth = FLARE.photom.flux_to_m(img.determine_depth())

        print('{0}: {1:.1f} | {2:.1f})'.format(filter, depth, depth - field.literature_depths[filter])) # 5\sigma limiting magnitude



