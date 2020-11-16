
import numpy as np

geo = (4. * np.pi * (100. * 10. * 3.0867 * 10 ** 16) ** 2)  # factor relating the L to M in cm^2

def M_to_log10L(M):
    return -0.4 * (M + 48.6) + np.log10(geo)

class Schechter():

    def __init__(self, p):

        if 'log10L*' not in p.keys():
            p['log10L*'] = M_to_log10L(p['M*'])

        p['L*'] = 10**p['log10L*']
        p['phi*'] = 10**p['log10phi*']
        self.p = p


    def phi_log10L(self, log10L):

        y = 10**(log10L - self.p['log10L*'])

        # n(log10L) dlog10L
        return self.p['phi*'] * np.log(10.) * (y**(self.p['alpha']+1.)) * np.exp(-y)
