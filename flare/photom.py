
import numpy as np

import flare


geo = 4.*np.pi*(100.*10.*3.0867*10**16)**2 # factor relating the L to M in cm^2



def flux_to_m(flux):

    return -2.5*np.log10(flux/1E9) + 8.9 # -- assumes flux in nJy

def m_to_flux(m):

    return 1E9 * 10**(-0.4*(m - 8.9)) # -- flux returned nJy

def flux_to_L(flux, cosmo, z):

    return flux*(4.*np.pi*cosmo.luminosity_distance(z).to('cm').value**2)/(1E9 * 1E23 * (1.+z))

# same as above but with consistent naming
def flux_to_lum(flux, cosmo, z):
    return flux_to_L(flux, cosmo, z)


def lum_to_flux(lum, cosmo, z):

    return 1E9 * 1E23 * lum * (1.+ z)/(4.*np.pi*cosmo.luminosity_distance(z).to('cm').value**2)


def lum_to_M(lum):

    return -2.5*np.log10(lum/geo)-48.6

def M_to_lum(M):

    return 10**(-0.4*(M+48.6)) * geo


# def M_to_m(M, z, cosmo = flare.default_cosmo()):
#
#     return flux_to_m(lum_to_flux(M_to_lum(M), cosmo, z))


def DM(z, cosmo = flare.default_cosmo()):
    luminosity_distance = cosmo.luminosity_distance(z).to('pc').value
    return 5*np.log10(luminosity_distance/(np.sqrt(1.+z)*10.))


def M_to_m(M, z, cosmo = flare.default_cosmo()):
    return M + DM(z, cosmo = cosmo)

def m_to_M(m, z, cosmo = flare.default_cosmo()):
    return m - DM(z, cosmo = cosmo)




def m_to_lum(m, cosmo, z):

    return flux_to_L(m_to_flux(m), cosmo, z)
