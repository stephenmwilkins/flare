

import numpy as np
import pickle

from . import core





def beta(lam, slope, normalisation, normalisation_wavelength = 1500., include_ISM = True):

    model = core.sed(lam)
    
    model.lnu = normalisation*(lam/normalisation_wavelength)**(slope + 2.0)
    
    if include_ISM: model.lnu[lam<912.] = 0.0 # add ISM absorption
    
    return model
    






