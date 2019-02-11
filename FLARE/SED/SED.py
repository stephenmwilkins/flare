

import numpy as np
import pickle

from . import IGM







class sed():

    def __init__(self, lam, description = False):
    
        self.description = description

        self.lam = lam # \AA
        self.lnu = np.zeros(self.lam.shape) # luminosity ers/s/Hz
        
    def get_l(self): # luminosity  erg/s
      
        nu = physics.constants.c / (self.lam * 1E-6)
        
        return self.Lnu * nu
         
    def get_fnu(self, cosmo, z, include_IGM = True): # flux nJy, depends on redshift and cosmology 

        self.lamz = self.lam * (1. + z)

        self.fnu = 1E23 * 1E9 * self.lnu * (1.+z) / (4 * np.pi * cosmo.luminosity_distance(z).to('cm').value**2) # nJy
        
        if include_IGM:
        
            self.fnu *= IGM.madau(self.lamz, z)
        
    def get_Fnu(self, F): # broad band flux/nJy
                        
        self.Fnu = {f: np.trapz(self.fnu * F[f].T, self.lamz) / np.trapz(F[f].T, self.lamz) for f in F.keys()}
             
    def get_Lnu(self, F): # broad band luminosity/erg/s/Hz
      
        self.Lnu = {f: np.trapz(self.lnu * F[f].T, self.lam) / np.trapz(F[f].T, self.lam) for f in F.keys()}
        





    










# 
# 
# 
# 
# class model():
# 
#     def __init__(self, grid, grid_data = 'data/grids', dust = False):
#     
#         self.grid = pickle.load(open(grid_data + '/' + grid + '/nebular.p','rb'), encoding='latin1')
# 
#         self.lam = self.grid['lam']
# 
#         self.dust = dust
#         
# 
# 
# 
# 
# 
# 
# class generate_SED():
# 
#     
#     def __init__(self, model, Masses, Ages, Metallicities, MetSurfaceDensities, include_intrinsic = True, IGM = False):
# 
# 
#         self.model = model
#         self.lam = self.model.grid['lam'] # convenience
# 
#         print(self.model.dust)
# 
#         self.stellar = sed(self.lam)
#         self.nebular = sed(self.lam)
#         self.total = sed(self.lam)
#     
#         if include_intrinsic:
#     
#             self.intrinsic_stellar = sed(self.lam)
#             self.intrinsic_nebular = sed(self.lam)
#             self.intrinsic_total = sed(self.lam)
#     
# 
#         for Mass, Age, Metallicity, MetalSurfaceDensity in zip(Masses, Ages, Metallicities, MetSurfaceDensities):
# 
# 
#             log10age = np.log10(Age) + 6. # log10(age/yr)
#             log10Z = np.log10(Metallicity) # log10(Z)
#         
# 
#             # --- determine dust attenuation
# 
#             if self.model.dust:
#     
#                 tau_V = (10**self.model.dust['A']) * MetalSurfaceDensity                     
# 
#                 tau = tau_V * (self.lam/5500.)**self.model.dust['slope']
#     
#                 T = np.exp(-tau)
#     
#             else:
#     
#                 T = 1.0
#     
#     
#             # --- determine closest SED grid point 
# 
#             ia = (np.abs(self.model.grid['log10age'] - log10age)).argmin()
#             iZ = (np.abs(self.model.grid['log10Z'] - log10Z)).argmin()
#  
#             self.stellar.lnu += Mass * T * self.model.grid['stellar'][ia, iZ] # erg/s/Hz
#             self.nebular.lnu += Mass * T * self.model.grid['nebular'][ia, iZ] # erg/s/Hz
# 
#             if include_intrinsic:
# 
#                 self.intrinsic_stellar.lnu += Mass * self.model.grid['stellar'][ia, iZ] # erg/s/Hz
#                 self.intrinsic_nebular.lnu += Mass * self.model.grid['nebular'][ia, iZ] # erg/s/Hz
# 
# 
# 
#         self.total.lnu = self.stellar.lnu + self.nebular.lnu # erg/s/Hz
#         
#         if include_intrinsic: self.intrinsic_total.lnu = self.intrinsic_stellar.lnu + self.intrinsic_nebular.lnu # erg/s/Hz
# 
#         
# 
#         
#             
#         
        


        
 