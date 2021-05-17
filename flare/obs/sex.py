

import os
import h5py
import numpy as np

import FLARE
import FLARE.filters


def run(img, detect_param_file = 'temp/detect.params'):

    cmd = f'{FLARE.FLARE_dir}/software/sextractor-2.19.5/sex {img}  -c {detect_param_file}'

    print(cmd)
    os.system(cmd)


def run_dual(detect_image, photo_image, detect_param_file = 'temp/detect.params'):

    cmd = f'{FLARE.FLARE_dir}/software/sextractor-2.19.5/sex {detect_image},{photo_image} -c {detect_param_file}'

    print(cmd)
    os.system(cmd)






def create_HDF5_cat(output_file, filters, path_to_sex_cat = 'temp/', close = False, verbose = False):

    hf = h5py.File(output_file, 'w')

    # --- Core Properties

    d = np.loadtxt(path_to_sex_cat+'detection.cat').T

    hf.create_dataset(f'label', data = d[0])

    hf.create_dataset(f'detection/x', data = d[1])
    hf.create_dataset(f'detection/y', data = d[2])
    hf.create_dataset(f'detection/ra', data = d[3])
    hf.create_dataset(f'detection/dec', data = d[4])

    hf.create_dataset(f'detection/circular_kron/flux', data =  d[5])
    hf.create_dataset(f'detection/circular_kron/error', data =  d[6])

    # --- Observed Properties

    PhotTypes = ['circular_kron', 'small_circular_kron']

    for f in filters:
        for PhotType in PhotTypes:

            d = np.loadtxt(f'{path_to_sex_cat}{PhotType}_{f}.cat').T

            nJy_to_es = 1E-9 * 10**(0.4*(FLARE.filters.info[f].zeropoint-8.9))

            hf.create_dataset(f'photometry/{PhotType}/{f}/flux', data = d[5]/nJy_to_es)
            hf.create_dataset(f'photometry/{PhotType}/{f}/error', data = d[6]/nJy_to_es)

    # --- print structure of file with shape of each object

    if verbose:
        def get_name_shape(name, item):
            shape = ''
            if hasattr(item, 'value'):
                shape = item.shape
            print(name, shape)

        hf.visititems(get_name_shape)

    if close:
        hf.close()
    else:
        return hf











class cat_params():


    def __init__(self):

        self.params = ['NUMBER','X_IMAGE','Y_IMAGE','ALPHA_J2000','DELTA_J2000','FLUX_AUTO','FLUXERR_AUTO']

    def write(self, filename = 'temp/cat.params'):

        open(filename,'w').writelines([f'{x}\n' for x in self.params])



class detect_params():

    def __init__(self, project = False):

        self.params = self.default_detect_params()

        if project == 'Hubble2020':

            self.params['DETECT_THRESH'] = 2.5
            self.params['ANALYSIS_THRESH'] = 2.5
            self.params['FILTER'] = 'N'
            self.params['WEIGHT_TYPE'] = 'MAP_RMS'


    def default_detect_params(self):

        #### SExtractor Default parameters

        params = {}

        #-------------------------------- Catalog ------------------------------------

        params['CATALOG_NAME'] = 'temp/output.cat'      # name of the output catalog
        params['CATALOG_TYPE'] = 'ASCII_HEAD'     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
                                        # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC
        params['PARAMETERS_NAME'] = 'temp/cat.params'  # name of the file containing catalog contents

        #------------------------------- Extraction ----------------------------------

        params['DETECT_TYPE'] = 'CCD'            # CCD (linear) or PHOTO (with gamma correction)
        params['DETECT_MINAREA'] = 5            # min. # of pixels above threshold
        params['DETECT_MAXAREA'] = 0             # max. # of pixels above threshold (0=unlimited)
        params['THRESH_TYPE'] = 'RELATIVE'       # threshold type: RELATIVE (in sigmas)
                                        # or ABSOLUTE (in ADUs)
        params['DETECT_THRESH'] = 1.5            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
        params['ANALYSIS_THRESH'] = 1.5            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2

        params['FILTER'] = 'Y'              # apply filter for detection (Y or N)?
        params['FILTER_NAME'] = 'default.conv'   # name of the file containing the filter
        params['FILTER_THRESH'] = ''                 # Threshold[s] for retina filtering

        params['DEBLEND_NTHRESH'] = 32             # Number of deblending sub-thresholds
        params['DEBLEND_MINCONT'] = 0.005          # Minimum contrast parameter for deblending

        params['CLEAN '] =  'Y'              # Clean spurious detections? (Y or N)?
        params['CLEAN_PARAM'] = 1.0            # Cleaning efficiency

        params['MASK_TYPE'] = 'CORRECT'        # type of detection MASKing: can be one of
                                        # NONE, BLANK or CORRECT

        #-------------------------------- WEIGHTing ----------------------------------

        params['WEIGHT_TYPE'] = 'NONE'           # type of WEIGHTing: NONE, BACKGROUND,
                                        # MAP_RMS, MAP_VAR or MAP_WEIGHT
        params['RESCALE_WEIGHTS'] = 'Y'              # Rescale input weights/variances (Y/N)?
        params['WEIGHT_IMAGE'] = 'weight.fits'    # weight-map filename
        params['WEIGHT_GAIN'] = 'Y'              # modulate gain (E/ADU) with weights? (Y/N)
        params['WEIGHT_THRESH'] = ''                # weight threshold[s] for bad pixels

        #-------------------------------- FLAGging -----------------------------------

        params['FLAG_IMAGE'] = 'flag.fits'      # filename for an input FLAG-image
        params['FLAG_TYPE'] = 'OR'             # flag pixel combination: OR, AND, MIN, MAX
                                        # or MOST

        #------------------------------ Photometry -----------------------------------

        params['PHOT_APERTURES'] = 5              # MAG_APER aperture diameter(s) in pixels
        params['PHOT_AUTOPARAMS'] = '2.5, 3.5'       # MAG_AUTO parameters: <Kron_fact>,<min_radius>
        params['PHOT_PETROPARAMS'] = '2.0, 3.5'       # MAG_PETRO parameters: <Petrosian_fact>,
                                        # <min_radius>
        params['PHOT_AUTOAPERS'] = '0.0,0.0'        # <estimation>,<measurement> minimum apertures
                                        # for MAG_AUTO and MAG_PETRO
        params['PHOT_FLUXFRAC'] = 0.5            # flux fraction[s] used for FLUX_RADIUS

        params['SATUR_LEVEL'] = 50000.0        # level (in ADUs) at which arises saturation
        params['SATUR_KEY'] = 'SATURATE'       # keyword for saturation level (in ADUs)

        params['MAG_ZEROPOINT'] = 0.0            # magnitude zero-point
        params['MAG_GAMMA'] = 4.0            # gamma of emulsion (for photographic scans)
        params['GAIN'] = 0.0            # detector gain in e-/ADU
        params['GAIN_KEY'] = 'GAIN'           # keyword for detector gain in e-/ADU
        params['PIXEL_SCALE'] = 1.0            # size of pixel in arcsec (0=use FITS WCS info)

        #------------------------- Star/Galaxy Separation ----------------------------

        params['SEEING_FWHM'] = 1.2            # stellar FWHM in arcsec
        params['STARNNW_NAME'] = 'default.nnw'    # Neural-Network_Weight table filename

        #------------------------------ Background -----------------------------------

        params['BACK_TYPE'] = 'AUTO'           # AUTO or MANUAL
        params['BACK_VALUE'] = 0.0            # Default background value in MANUAL mode
        params['BACK_SIZE'] = 64             # Background mesh: <size> or <width>,<height>
        params['BACK_FILTERSIZE'] = 3              # Background filter: <size> or <width>,<height>

        params['BACKPHOTO_TYPE'] = 'GLOBAL'         # can be GLOBAL or LOCAL
        params['BACKPHOTO_THICK'] = 24             # thickness of the background LOCAL annulus
        params['BACK_FILTTHRESH'] = 0.0            # Threshold above which the background-
                                        # map filter operates

        #------------------------------ Check Image ----------------------------------

        params['CHECKIMAGE_TYPE'] = 'NONE'           # can be NONE, BACKGROUND, BACKGROUND_RMS,
                                        # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
                                        # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,
                                        # or APERTURES
        params['CHECKIMAGE_NAME'] = 'check.fits'     # Filename for the check-image

        #--------------------- Memory (change with caution!) -------------------------

        params['MEMORY_OBJSTACK'] = 3000           # number of objects in stack
        params['MEMORY_PIXSTACK'] = 300000         # number of pixels in stack
        params['MEMORY_BUFSIZE'] = 1024           # number of lines in buffer

        #------------------------------- ASSOCiation ---------------------------------

        params['ASSOC_NAME'] = 'sky.list'       # name of the ASCII file to ASSOCiate
        params['ASSOC_DATA'] = '2,3,4'          # columns of the data to replicate (0=all)
        params['ASSOC_PARAMS'] = '2,3,4'          # columns of xpos,ypos[,mag]
        params['ASSOCCOORD_TYPE'] = 'PIXEL'          # ASSOC coordinates: PIXEL or WORLD
        params['ASSOC_RADIUS'] = '2.0'            # cross-matching radius (pixels)
        params['ASSOC_TYPE'] = 'NEAREST'        # ASSOCiation method: FIRST, NEAREST, MEAN,
                                        # MAG_MEAN, SUM, MAG_SUM, MIN or MAX
        params['ASSOCSELEC_TYPE'] = 'MATCHED'        # ASSOC selection type: ALL, MATCHED or -MATCHED

        #----------------------------- Miscellaneous ---------------------------------

        params['VERBOSE_TYPE'] = 'NORMAL'         # can be QUIET, NORMAL or FULL
        params['HEADER_SUFFIX'] = '.head'          # Filename extension for additional headers
        params['WRITE_XML'] = 'N'              # Write XML file (Y/N)?
        params['XML_NAME'] = 'sex.xml'        # Filename for XML output
        params['XSL_URL'] = 'file:///usr/local/share/sextractor/sextractor.xsl'
                                        # Filename for XSL style-sheet
        params['NTHREADS'] = '1'              # 1 single thread

        params['FITS_UNSIGNED'] = 'N'              # Treat FITS integer values as unsigned (Y/N)?
        params['INTERP_MAXXLAG'] = '16'             # Max. lag along X for 0-weight interpolation
        params['INTERP_MAXYLAG'] = '16'             # Max. lag along Y for 0-weight interpolation
        params['INTERP_TYPE'] = 'ALL'            # Interpolation type: NONE, VAR_ONLY or ALL

        #--------------------------- Experimental Stuff -----------------------------

        params['PSF_NAME'] = 'default.psf'   # File containing the PSF model
        params['PSF_NMAX'] = 1              # Max.number of PSFs fitted simultaneously
        params['PATTERN_TYPE'] = 'RINGS-HARMONIC' # can RINGS-QUADPOLE, RINGS-OCTOPOLE,
                                        # RINGS-HARMONICS or GAUSS-LAGUERRE
        params['SOM_NAME'] = 'default.som'    # File containing Self-Organizing Map weights

        return params


    def write(self, filename = 'temp/detect.params'):

        lines = [f'{k}\t {v}\n' for k,v in self.params.items()]

        open(filename,'w').writelines(lines)
