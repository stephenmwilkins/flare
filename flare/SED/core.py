
import numpy as np

from . import IGM

from ..photom import *



class sed():

    def __init__(self, lam, description = False):

        self.description = description

        self.lam = lam # \AA
        self.lnu = np.zeros(self.lam.shape) # luminosity ers/s/Hz
        self.nu = 3E8/(self.lam/1E10) # Hz

#     def get_l(self): # luminosity  erg/s
#
#         nu = physics.constants.c / (self.lam * 1E-10)
#
#         return self.Lnu * nu

    def get_Lnu(self, F): # broad band luminosity/erg/s/Hz

        self.Lnu = {f: np.trapz(self.lnu * F[f].T, self.lam) / np.trapz(F[f].T, self.lam) for f in F['filters']}


    def return_Lnu(self, F): # broad band luminosity/erg/s/Hz

        return {f: np.trapz(self.lnu * F[f].T, self.lam) / np.trapz(F[f].T, self.lam) for f in F['filters']}

    
    def get_fnu(self, cosmo, z, include_IGM = True): # flux nJy, depends on redshift and cosmology

        self.lamz = self.lam * (1. + z)

        self.fnu = 1E23 * 1E9 * self.lnu * (1.+z) / (4 * np.pi * cosmo.luminosity_distance(z).to('cm').value**2) # nJy

        if include_IGM:
            self.fnu *= IGM.madau(self.lamz, z)

    def get_Fnu(self, F): # broad band flux/nJy

        self.Fnu = {f: np.trapz(self.fnu * F[f].T, self.lamz) / np.trapz(F[f].T, self.lamz) for f in F['filters']}

        self.Fnu_array = np.array([self.Fnu[f] for f in F['filters']])

    def return_Fnu(self, F): # broad band flux/nJy

        return {f: np.trapz(self.fnu * F[f].T, self.lamz) / np.trapz(F[f].T, self.lamz) for f in F['filters']}

def rebin(l, f, n): # rebin SED [currently destroys original]

    n_len = int(np.floor(len(l)/n))
    l = l[:n_len*n]
    f = f[:n_len*n]
    nl = np.mean(l.reshape(n_len,n), axis=1)
    nf = np.sum(f.reshape(n_len,n), axis=1)/n

    return nl, nf
