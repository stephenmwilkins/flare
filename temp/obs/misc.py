


import numpy as np

from .. import filters




def Casternato(filter, aperture_diameter_arcsec = False, output_pixel_size = 0.06, pix_frac = 0.8, verbose = False): 

    """calculate the noise scaling factor given in Casternato+ (https://iopscience.iop.org/article/10.1086/316851/pdf)"""

    native_pixel_size = filters.info[filter].pixel_scale

    s = output_pixel_size/native_pixel_size

    p = pix_frac

    # --- if no aperture given calculate for a single pixel
    
    if not aperture_diameter_arcsec:
    
        if p<s:
            sqrF_A = 1-p/(3*s)
        else:
            sqrF_A = (s/p)*(1-s/(3*p))

        if verbose: print('scaling for single pixel: {0:.2f}'.format(sqrF_A))    
 
    if aperture_diameter_arcsec:
    
    # -- different apertures 

        l = (aperture_diameter_arcsec/2.)/native_pixel_size # SHOULD THIS BE output_pixel_size

        if l>p:
            F_A = (1-p/(3*l))**2
        else:
            F_A = (l/p)**2 * (1-l/(3*p))**2
        
        sqrF_A = np.sqrt(F_A)
        
        if verbose: print('scaling for {0}" diameter aperture: {1:.2f}'.format(aperture_diameter_arcsec, sqrF_A))    


    return sqrF_A
    
    
    
    
 
 
 
 
 
    
    
    


from scipy.stats import linregress


def measure_beta(z, fluxes, F, verbose = False):

    # --- determine which filters can be used

    fitFilters = []
    for f in F['filters']:
        if min(F[f].l[F[f].t > 0.01]) >= 1216.*(1+z) and max(F[f].l[F[f].t > 0.01]) <= 3000.*(1+z):
            fitFilters.append(f)
    
    if verbose: print(fitFilters)
        
    pivwvs = [F[f].pivwv() for f in fitFilters] 
    fluxes = [fluxes[f] for f in fitFilters]

    if len(fluxes) > 1:
        beta = linregress(np.log10(pivwvs), np.log10(fluxes))[0] - 2.
    else:
        beta = None
    
    if verbose: print(f'{beta:.2f}')
    
    return beta
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

