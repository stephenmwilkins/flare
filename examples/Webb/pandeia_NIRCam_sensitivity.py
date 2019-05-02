 



import numpy as np

import FLARE.filters

import pickle





run = True

lw_filters = FLARE.filters.NIRCam_l
sw_filters = FLARE.filters.NIRCam_s

filters = sw_filters + lw_filters

F = FLARE.filters.add_filters(filters)


if run:

    try:
        from pandeia.engine.calc_utils import build_default_calc
    except:
        print('failed')

    from pandeia.engine.calc_utils import build_default_calc
    from pandeia.engine.perform_calculation import perform_calculation

    config = {}

    config['background'] = 'minzodi'
    config['background_level'] = 'benchmark'

    config['configuration'] = {}
    config['configuration']['detector'] = {'nexp': 1, 'ngroup':47, 'nint':1, 'readmode': 'deep8', 'subarray': 'full'}
    config['configuration']['instrument'] = {'aperture': 'sw', 'instrument': 'nircam', 'mode': 'sw_imaging', 'disperser': None}

    scene = {}
    scene['position'] = {'x_offset': 0., 'y_offset': 0., 'orientation': 0., 'position_parameters': ['x_offset', 'y_offset', 'orientation']}
    scene['shape'] = {'geometry': 'point'}
    scene['spectrum'] = {'name': 'Flat Source', 'spectrum_parameters': ['sed', 'normalization']}
    scene['spectrum']['sed'] = {'sed_type': 'flat', 'unit': 'fnu'}
    scene['spectrum']['normalization'] = {'type': 'at_lambda', 'norm_wave': 2.0, 'norm_waveunit': 'um', 'norm_flux': 10., 'norm_fluxunit': 'njy'}
    config['scene'] = [scene]

    config['strategy'] = {'background_subtraction': False, 'method':'imagingapphot', 'target_source':'1', 'target_type':'source', 'target_xy':[0.0,0.0], 'units':'arcsec'}

    limits = {}

    for f in filters:

        config['configuration']['instrument']['filter'] = f.split('.')[-1].lower()

        if f in sw_filters:
    
            config['configuration']['instrument']['aperture'] = 'sw'
            config['configuration']['instrument']['mode'] = 'sw_imaging'
            config['strategy']['aperture_size'] = 0.031*2.5

        if f in lw_filters:
    
            config['configuration']['instrument']['aperture'] = 'lw'
            config['configuration']['instrument']['mode'] = 'lw_imaging'
            config['strategy']['aperture_size'] = 0.063*2.5
        
        
        # --- iteratively find the flux yielding a S/N=10
        
        target_SN = 10.
        min_flux = 5.
        max_flux = 70.
        tol = 0.01
        
        while min_flux <= max_flux:

            current_flux = (min_flux + max_flux)/2.

            scene['spectrum']['normalization'] ['norm_flux'] = current_flux
            config['scene'] = [scene]
            report = perform_calculation(config)
            
            sn = report['scalar']['sn']
        
            if sn <= target_SN: 
                min_flux = current_flux + tol
            else:
                max_flux = current_flux - tol
              
        limit = current_flux
        
        print(f,'{0:.2f}'.format(limit))
    
        limits[f] = limit
    
    pickle.dump(limits, open('NIRCam_limits.p','wb'))
    
    
    
    
    
if not run:
    
    limits = pickle.load(open('NIRCam_limits.p','rb'))
    
    print(limits)
    
    
    
import matplotlib.pyplot as plt 
from matplotlib import cm  


left  = 0.15
bottom = 0.15  
height = 0.8
width = 0.8

fig = plt.figure(figsize = (6,4))

ax = fig.add_axes((left, bottom, width, height))
    
for i,f in enumerate(F['filters']): 

    pivwv = F[f].pivwv()/1E4

    # c = cm.jet(i/len(F['filters']))

    c = cm.jet((pivwv - 0.7)/4.1)
    
    ax.scatter(pivwv, np.log10(limits[f]), zorder = 3, label = f.split('.')[-1], c=[c]) # --- plot pivot wavelength
    ax.text(pivwv, np.log10(limits[f])+0.02, f.split('.')[-1], ha='center', va='center', fontsize=6)
   

ax.set_xlabel(r'$\rm\lambda/\mu m$')
ax.set_ylabel(r'$\rm\log_{10}(f_{\nu}/nJy)\ with\ S/N=10\ in\ 10ks$')


# plt.legend(fontsize=7)

plt.savefig('NIRCam_sensitivity.pdf')
    
    