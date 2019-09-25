import sys
import os

import numpy as np

import FLARE
import FLARE.filters
import FLARE.obs
import FLARE.obs.misc

import FLARE.SED.models
import FLARE.SED.SFZH as SFZH


cosmo = FLARE.default_cosmo()

filters = FLARE.filters.WFC3NIR_W

z = 7.

# --- using \beta model 

beta = -2.43
rest_lam = np.arange(0., 5000., 1.)
F = FLARE.filters.add_filters(filters, new_lam = rest_lam * (1. + z))

sed = FLARE.SED.models.beta(rest_lam, beta, 1E28, normalisation_wavelength = 1500.) # generate rest-frame SED
sed.get_fnu(cosmo, z) # --- generate observed frame SED (necessary to get broad band photometry)
Fnu = sed.return_Fnu(F) # --- generate observed frame broadband photometry
        
FLARE.obs.misc.measure_beta(z, Fnu, F, verbose = True)

# --- using SPS model 



SPS = FLARE.SED.models.SPS('BPASSv2.2.1.binary/Chabrier_300')


F = FLARE.filters.add_filters(filters, new_lam = SPS.lam * (1. + z))

sfzh, sfr = SFZH.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': 8., 'log10Z': -8.0, 'log10M*': 8.0})


for log10tau_V in np.arange(-2., 2.0, 0.5):

    SED = SPS.get_Lnu(sfzh, {'fesc': 0.0, 'log10tau_V': log10tau_V}, dust_model = 'very_simple')

    sed = SED.total 

    sed.get_fnu(cosmo, z) # --- generate observed frame SED (necessary to get broad band photometry)
    Fnu = sed.return_Fnu(F) # --- generate observed frame broadband photometry
 
    beta = FLARE.obs.misc.measure_beta(z, Fnu, F)
    
    print(f'{log10tau_V}: {beta:.2f}')