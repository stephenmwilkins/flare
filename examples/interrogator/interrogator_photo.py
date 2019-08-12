
import time
import numpy as np
import sys
import os
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath('../../'))

class empty: pass # --- output class

import FLARE
import FLARE.photom
import FLARE.filters
from FLARE.interrogator import models
from FLARE.interrogator import fitter
from FLARE.interrogator import plots


z = 8.0


# ID = 'XDF_nodust_fixedz'
# filters = FLARE.filters.XDF
# model_definitions = {'SPS': 'BPASSv2.2.1.binary', 'IMF': 'ModSalpeter_300', 'SFZH': 'constant_constant', 'dust': False, 'cosmo': FLARE.default_cosmo()}
# redshift = z

# ID = 'XDF_fixedz'
# filters = FLARE.filters.XDF
# model_definitions = {'SPS': 'BPASSv2.2.1.binary', 'IMF': 'ModSalpeter_300', 'SFZH': 'constant_constant', 'dust': 'very_simple', 'cosmo': FLARE.default_cosmo()}
# redshift = z

ID = 'XDF'
filters = FLARE.filters.XDF
model_definitions = {'SPS': 'BPASSv2.2.1.binary', 'IMF': 'ModSalpeter_300', 'SFZH': 'constant_constant', 'dust': 'very_simple', 'cosmo': FLARE.default_cosmo()}
redshift = False

# --- generate fake data

parameters = {'log10_duration': 8., 'log10Z': -2., 'log10M*': 8., 'log10fesc': -1., 'log10tau_V': -1.0, 'z': z}






# --- full SED

# print('-'*20)
# print('using full SED generator:')
# 
# model = models.full(model_definitions, filters = filters)
# mod = model.p(parameters) 
# 
# for f,flux in zip(filters, mod.F): print('{0}: {1:.3f}'.format(f, flux, FLARE.photom.flux_to_m(flux)))


print('-'*20)
print('using pre-computed photogrid:')

model = models.photo(model_definitions, filters = filters)
mod = model.p(parameters) 

print('SFR10: {0:.2f}'.format(mod.properties['SFR10']))
for f,flux in zip(filters, mod.F): print('{0}: {1:.3f} {2:.2f}'.format(f, flux, FLARE.photom.flux_to_m(flux)))




print('-'*20)

sigma = 0.1 # nJy

obs = empty()
obs.ID = ID
obs.filters = filters
obs.true_p = parameters
obs.true_properties = mod.properties
obs.true_fluxes = mod.F
obs.fluxes = obs.true_fluxes + sigma * np.random.randn(len(obs.filters))
obs.flux_errors = np.ones(obs.fluxes.shape) * sigma

for f,flux,fe in zip(obs.filters, obs.fluxes, obs.flux_errors):
    print('{0}: {1:.2f} {2:.2f} {3:.2f}'.format(f, flux, fe, flux/fe))


if redshift: model.prior_def['z'] = {'type': 'delta', 'value': z} # --- ASSUME REDSHIFT IS KNOWN

# --- test model call time

# start_time = time.time()
# N = 10
# for i in range(N): model.p(parameters) 
# print('{0:.4f}'.format((time.time() - start_time)/N))



# --- run fitter

run = True
if run:
    start_time = time.time()
    fit = fitter.source(obs, model, verbose = True)
    nwalkers = 50
    nsamples = 2000
    burn = 200
    fit.fit(nwalkers = nwalkers, nsamples = nsamples, burn = burn)
    total_time = time.time() - start_time
    print('total time: {0:.4f}'.format(total_time))
    print('time per sample: {0:.4f}'.format(total_time/(nsamples * nwalkers)))
    fit.save()

plots.plots(ID).simple_triangle()


