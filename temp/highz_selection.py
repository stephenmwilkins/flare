

import numpy as np


def basic_pz(self, include_detected_flag = True, path = ''):

    """Basic photometric redshift based selection"""

    sn_cut = 7
    l68_cut = 6

    # --- apply additional selection criteria to catalogue

    sn = self.o[f'{path}detection/sn']
    l68 = self.o[f'{path}eazy/l68']

    if include_detected_flag:
        detected = self.o[f'{path}detected']
        selected = (detected==True)&(sn>sn_cut)&(l68>l68_cut)
    else:
        selected = (sn>sn_cut)&(l68>l68_cut)

    self.o[f'{path}selected'] = selected

    self.N_selected = len(selected[selected])

    if self.verbose:
        print(f'Number of selected high-z galaxies: {self.N_selected}')


# --- add alternative selections based on colour cuts for example
