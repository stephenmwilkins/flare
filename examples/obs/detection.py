
import sys
import os

from astropy.io import fits
import numpy as np
import pickle

np.random.seed(42)

from photutils import detect_sources
from photutils import deblend_sources
from photutils import source_properties

sys.path.insert(0, os.path.abspath('../../'))

import FLARE
import FLARE.filters
import FLARE.obs
import FLARE.obs.photometry as photometry
import FLARE.surveys
import FLARE.photom

test = True


sci_suffix = 'sci' # original images
sci_suffix = 'sci_convolved' # original images




field = FLARE.surveys.XDF.fields['dXDF'] # deepest portion of the XDF

filters = ['HST.ACS.f850lp','HST.WFC3.f105w','HST.WFC3.f125w','HST.WFC3.f140w','HST.WFC3.f160w']
filters = field.filters

FilterInfo = FLARE.filters.add_filters(filters) # Filter info 

generate_img = False

if generate_img:

    imgs = FLARE.obs.open_images(field, filters, verbose = True, sci_suffix = sci_suffix) # opens image and applies mask. img is actually an instance of the obs.core.image class

    x,y = next(iter(imgs.values())).get_random_location()
    width = 400

    CutoutImages = {filter: imgs[filter].make_cutout(x,y, width) for filter in filters} 

    detection_filters = ['HST.WFC3.f125w','HST.WFC3.f140w','HST.WFC3.f160w']

    DetectionCutoutImages = {filter: CutoutImages[filter] for filter in detection_filters}

    DetectionImage = FLARE.obs.create_stack(DetectionCutoutImages)

    pickle.dump(CutoutImages, open('CutoutImages.p','wb'))
    pickle.dump(DetectionImage, open('DetectionImage.p','wb'))

else:

    CutoutImages = pickle.load(open('CutoutImages.p','rb'))
    DetectionImage = pickle.load(open('DetectionImage.p','rb'))
    


threshold = 2.5
npixels = 5

SegmentationImage = detect_sources(DetectionImage.sig, threshold, npixels = npixels)  
DeblendedSegmentationImage = deblend_sources(DetectionImage.sig, SegmentationImage, npixels = npixels, nlevels=32)
AllSourceProperties = source_properties(DetectionImage.sig, DeblendedSegmentationImage)

# FLARE.obs.make_plots({**cutouts, **{'stack':DetectionImage}, **{'segm': DeblendedSegmentationImage.data}}, signficance_plot = [True, True, True, True, False])

print('Number of objects detected: {0}'.format(len(AllSourceProperties)))








if test: AllSourceProperties = [AllSourceProperties[7]]


for idx, SourceProperties in enumerate(AllSourceProperties):

    output_folder = 'detection/{0}'.format(idx)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    CoreProperties, Mask, ExclusionMask = photometry.measure_core_properties(SourceProperties, DetectionImage, DeblendedSegmentationImage, verbose = True)

    #Photometry = {filter: photometry.measure_photometry(CoreProperties, CutoutImages[filter], Mask, ExclusionMask, verbose = True) for filter in filters}

    # pickle.dump(SourceProperties, open('detection/{0}/SourceProperties.pck'.format(idx),'wb')) # NOTE: very large (~2.5Mb per object)
    pickle.dump(CoreProperties, open(output_folder+'/CoreProperties.pck','wb')) 
    #pickle.dump(Photometry, open(output_folder+'/Photometry.pck','wb')) 

    
    # --- make cutout of single object (COULD AUTOMATE FURTHER)

    width = 50

    x,y = round(CoreProperties['y']),round(CoreProperties['x'])

    object_cutouts = {filter: CutoutImages[filter].make_cutout(x,y, width) for filter in filters} 

    #FLARE.obs.make_plots({**object_cutouts}, signficance_plot = True, filename = output_folder+'/stamps.png', show = test)


    DetectionImage_cutout = DetectionImage.make_cutout(x,y, width)

    DeblendedSegmentationImage_cutout = FLARE.obs.make_cutout(DeblendedSegmentationImage.data, x, y, width)

    ExclusionMask_cutout = FLARE.obs.make_cutout(ExclusionMask, x, y, width).astype('bool')

    #FLARE.obs.make_plots({**{'stack':DetectionImage_cutout}, **{'segm':DeblendedSegmentationImage_cutout}, **{'exclusion_mask': ExclusionMask_cutout}}, signficance_plot = [True, False, False], filename = output_folder+'/detection.png', show = test)

    #photometry.COG_plots(Photometry, filename = output_folder+'/COG.png', show = test)
    
    #photometry.SED_plot(Photometry, FilterInfo = FilterInfo, filename = output_folder+'/SED.png', show = test)
    
    photometry.size_plot(DetectionImage_cutout, CoreProperties, ExclusionMask_cutout, filename = output_folder+'/sizes.png', show = test)
    