


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import FLARE.plt as fplt

import FLARE.LF.evo as evo
import FLARE.LF.lf_parameters as lf_params




# --- Models

fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, figsize = (3,6))
fig.subplots_adjust(left=0.2, bottom=0.05, top=0.85, right=0.95, wspace=0.0, hspace=0.0)


axes[-1].set_xlabel(fplt.ml('z'))
axes[0].set_ylabel(fplt.ml(r'M^{\star}'))
axes[1].set_ylabel(fplt.ml(r'\phi^{\star}'))
axes[2].set_ylabel(fplt.ml(r'\alpha'))



axes[0].set_ylim([-21.95, -19.55])
axes[1].set_ylim([-6.45, -2.05])
axes[2].set_ylim([-2.95, -1.55])

for ax in axes:
    ax.set_xlim([4.5,14.])


models = ['bluetides','Finkelstein_review','Ma2019','Mason15','Yung2018','FLARES']
Nmodels = len(models)

s = 30



for model, c, ms in zip(models, cm.viridis(np.arange(Nmodels)/Nmodels), ['o','^','h','d','*','v','p','s']):

    m = getattr(lf_params, model)()

    if model == 'FLARES':
        alpha = 1.0
    else:
        alpha = 0.4

    axes[0].scatter(m.redshifts, m.M_star, marker = ms, c = [c], alpha = alpha, s=s, label = m.name)
    axes[1].scatter(m.redshifts, m.phi_star, marker = ms, c = [c], alpha = alpha, s=s)
    axes[2].scatter(m.redshifts, m.alpha, marker = ms, c = [c], alpha = alpha, s=s)

    # --- add linear fit

    # linear = evo.linear(m)
    #
    # zr, p = linear.parameters_line()
    #
    # axes[0].plot(zr, p['M*'], c = c, alpha = 0.2)
    # axes[1].plot(zr, p['log10phi*'], c = c, alpha = 0.2)
    # axes[2].plot(zr, p['alpha'], c = c, alpha = 0.2)


axes[0].legend(fontsize = 6, bbox_to_anchor=(0.0, 1.05), loc = 'lower left')

# fig.savefig('figs/LF_parameters_models.pdf')

plt.show()


# Printout of all model parameters in FLARE.LF module
models = evo.model_names

evo.print_model_parameters(models)