


import numpy as np
import pickle
import h5py
from types import SimpleNamespace


from photutils import detect_sources
from photutils import deblend_sources
from photutils import source_properties

import FLARE
import FLARE.filters
import FLARE.obs.photometry as photometry
import FLARE.surveys
import FLARE.photom
import FLARE.obs.EAZY as eazy

import FLARE.obs.plots as plots

class detector():

    def __init__(self, surveyName, fieldName, DetectionFilters, sci_suffix = 'sci_convolved', verbose = False, make_figures = False):


        self.surveyName = surveyName
        self.fieldName = fieldName
        self.DetectionFilters = DetectionFilters
        self.sci_suffix = sci_suffix
        self.verbose = verbose
        self.make_figures = make_figures

        self.Field = FLARE.surveys.surveys[self.surveyName].fields[self.fieldName] # deepest portion of the XDF
        self.Filters = self.Field.filters


    # def go(self):
    #
    #     self.make_detect_image(save_detection_image = True)
    #     self.make_segmentation_image(save_segmentation_image = True)
    #     self.detect_sources()

    def go_from_bkup(self, N = False):

        self.DetectionImage = pickle.load(open('temp/DetectionImage.p','rb'))
        if self.verbose: print('read DetectionImage')
        self.DeblendedSegmentationImage = pickle.load(open('temp/DeblendedSegmentationImage.p','rb'))
        if self.verbose: print('read DeblendedSegmentationImage')
        self.detect_sources()

        # self.Filters = [self.Filters[-1]]
        self.get_images()

        Sources = self.get_all_source_properties(N=N)

        hf = self.save_to_HDF5(Sources)

        self.run_EAZY(hf)

        hf.close()


    def test(self, i=42):

        self.DetectionImage = pickle.load(open('temp/DetectionImage.p','rb'))
        if self.verbose: print('read DetectionImage')
        self.DeblendedSegmentationImage = pickle.load(open('temp/DeblendedSegmentationImage.p','rb'))
        if self.verbose: print('read DeblendedSegmentationImage')
        self.detect_sources()

        # self.Filters = [self.Filters[-1]]
        self.get_images()

        s = self.get_source_properties(i)

        self.make_cutout_figure(s)
        self.make_detection_figure(s)



    def make_cutout_figure(self, s):

        CutoutImages = {filter: self.Images[filter].make_cutout(int(s.DetectionProperties['y']),int(s.DetectionProperties['x']), 50) for filter in self.Filters}
        plots.make_plots({**CutoutImages}, signficance_plot = True, filename = 'temp/test_cutouts.png', imsize = 4)

    def make_detection_figure(self, s):

        DetectionImageCutout = self.DetectionImage.make_cutout(int(s.DetectionProperties['y']),int(s.DetectionProperties['x']), 50)
        DeblendedSegmentationImageCutout = self.DeblendedSegmentationImage.make_cutout(int(s.DetectionProperties['y']),int(s.DetectionProperties['x']), 50)
        plots.make_plots({'DetectionImaget':DetectionImageCutout, 'DeblendedSegmentationImage':DeblendedSegmentationImageCutout, }, signficance_plot = [True, False], filename = 'temp/test_detect_segmentaion.png')




    # def go_from_hdf5(self):
    #
    #     OutputFolder = f'data/{self.surveyName}/{self.fieldName}/'
    #
    #     OutputFile = f'{OutputFolder}/all.h5'
    #
    #     hf = h5py.File(OutputFile, 'r+')
    #
    #     hf = self.run_EAZY(hf)
    #
    #     hf.close()


    def get_images(self):

        self.Images = FLARE.obs.open_images(self.Field, self.Filters, verbose = self.verbose, sci_suffix = self.sci_suffix) # opens image and applies mask. img is actually an instance of the obs.core.image class


    def make_detect_image(self, save_detection_image = False):

        DetectionImages = {filter: self.Images[filter] for filter in self.detection_filters}
        self.DetectionImage = FLARE.obs.create_stack(DetectionImages)

        if save_detection_image: pickle.dump(self.DetectionImage, open('temp/DetectionImage.p','wb'))


    def make_segmentation_image(self, threshold = 2.5, npixels = 5, nlevels = 32, save_segmentation_image = False):

        SegmentationImage = detect_sources(self.DetectionImage.sig, threshold, npixels = npixels)
        self.DeblendedSegmentationImage = deblend_sources(self.DetectionImage.sig, SegmentationImage, npixels = npixels, nlevels = nlevels)

        if save_segmentation_image: pickle.dump(DeblendedSegmentationImage, open('temp/DeblendedSegmentationImage.p','wb'))


    def detect_sources(self):

        self.AllSourceProperties = source_properties(self.DetectionImage.sig, self.DeblendedSegmentationImage)

        if self.verbose:
            print(f'{len(self.AllSourceProperties)} objects detected')



    def get_all_source_properties(self, N = False):

        if not N: N = len(self.AllSourceProperties) # do all sources

        return list(filter(None, [self.get_source_properties(i) for i in range(N)]))


    def get_source_properties(self, i):

        try:

            SourceProperties = self.AllSourceProperties[i]

            s = SimpleNamespace() # output object

            s.label = SourceProperties.label

            if self.verbose: print('-'*5+'{0}'.format(s.label))

            # --- This bit does all the work
            s.DetectionProperties, Mask, ExclusionMask = photometry.measure_core_properties(SourceProperties, self.DetectionImage, self.DeblendedSegmentationImage, measure_apertures = False, verbose = self.verbose)

            s.Photometry = {filter: photometry.measure_properties_quick(s.DetectionProperties, self.Images[filter], Mask, ExclusionMask, verbose = self.verbose) for filter in self.Filters}

            return s

        except:

            return


    def save_to_HDF5(self, Sources):

        OutputFolder = f'data/{self.surveyName}/{self.fieldName}/'

        OutputFile = f'{OutputFolder}/all.h5'

        hf = h5py.File(OutputFile, 'w')

        hf.create_dataset('label', data = np.array([source.label for source in Sources]))


        # --- Core Properties

        for p in ['x', 'y', 'area', 'radius', 'A', 'B', 'theta', 'ellipticity', 'kron_radius']:
            hf.create_dataset(f'detection/{p}', data =  np.array([source.DetectionProperties[p] for source in Sources]))

        hf.create_dataset(f'detection/circular_kron/flux', data =  np.array([source.DetectionProperties['circular_kron'].flux for source in Sources]))
        hf.create_dataset(f'detection/circular_kron/error', data =  np.array([source.DetectionProperties['circular_kron'].noise for source in Sources]))


        # --- Observed Properties


        PhotTypes = ['circular_kron', 'small_circular_kron']

        for f in self.Filters:
            for PhotType in PhotTypes:
                hf.create_dataset(f'photometry/{PhotType}/{f}/flux', data = np.array([source.Photometry[f][PhotType].flux for source in Sources]))
                hf.create_dataset(f'photometry/{PhotType}/{f}/error', data = np.array([source.Photometry[f][PhotType].error for source in Sources]))


        # --- print structure of file with shape of each object

        if self.verbose:
            def get_name_shape(name, item):
                shape = ''
                if hasattr(item, 'value'):
                    shape = item.shape
                print(name, shape)

            hf.visititems(get_name_shape)

        return hf



    def run_EAZY(self, hf):

        F = FLARE.filters.add_filters(self.Filters)

        hf_EAZY = eazy.eazy(ID = 'all', create_POFZ_FILE = True).run_new(hf, F, path = lambda f: f'photometry/small_circular_kron/{f}/')

        # --- append EAZY group to original file

        hf_EAZY.copy(f'EAZY', hf, name = f'EAZY')

        return hf
