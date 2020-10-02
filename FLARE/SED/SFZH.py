

import numpy as np

from scipy import integrate


class empty: pass


def simple(SPS, p):

    # --- returns the SFZH for a single age/metallicity

    SFZH = np.zeros(SPS.grid['log10Q'].shape)

    iZ = (np.abs(SPS.grid['log10Z'] - p['log10Z'])).argmin()

    a = np.abs(SPS.grid['log10age'] - p['log10age'])

    if np.min(a)==0.0:

        ia = a.argmin()

        SFZH[ia, iZ] = 1.0

    else:

        # --- determine bracketing indicies

        ias = a.argsort()[:2]
        SFZH[ias, iZ] = 1./a[ias]

    SFZH /= np.sum(SFZH)

    return SFZH





def CSFH(SPS, p, redshift=0.0, log10Z=-2.):

    log10ages = SPS.grid['log10age']

    iZ = (np.abs(SPS.grid['log10Z'] - log10Z)).argmin()

    SFZH = np.zeros((len(SPS.grid['log10age']), len(SPS.grid['log10Z'])))

    from astropy.cosmology import WMAP9 as cosmo

    dz = 0.01
    z = np.arange(10., redshift, -dz)



    ages = -(cosmo.age(z).to('yr').value - cosmo.age(redshift).to('yr').value)

    csfh = lambda z, p1, p2, p3, p4: p1 * (1+z)**p2 / (1+((1+z)/p3)**p4 )

    sfrd = csfh(z, *p)

    f = lambda x: np.interp(x, ages, sfrd[::-1])


    start = 0.
    for ia, log10age in enumerate(log10ages[:-1]):

        end = 10**log10age

        # determine the amount of star formation in this bin

        sf = integrate.quad(f, start, end)[0]

        SFZH[ia, iZ] = sf
        # print('{0:.1f}-{1:.1f}: SF={2:.2f}'.format(start/1E6, end/1E6, sf))

        start = end

    SFZH /= np.sum(SFZH)

    SFR = csfh(redshift, *p)

    return SFZH, SFR



def constant(log10ages, log10Zs, p):

    # --- returns the SFZH for a constant SFH and single metallicity. At the moment this just chooses a single metallicity.

    a0 = 10**p['log10_duration']


    SFZH = np.zeros((len(log10ages), len(log10Zs)))

    # --- determine the metallicity bin in which to place all the SF

    iZ = (np.abs(log10Zs - p['log10Z'])).argmin()

    prev = 0.0

    cont = True

    for ia, log10age in enumerate(log10ages[:-1]):

        c = 10**np.mean([log10ages[ia+1],log10age])

        if cont:

            if a0<c:

                w = a0 - prev # --- determine age width of bin
                cont = False

            else:

                w = c - prev # --- determine age width of bin

            SFZH[ia, iZ] = w

        else:

            w = 0

        prev = c

        # print(ia,log10age,c,w)

    SFZH /= np.sum(SFZH)

    SFZH *= 10**p['log10M*']

    SFR = 10**p['log10M*']/10**p['log10_duration']

    return SFZH, SFR



def instantaneous(log10ages, log10Zs, p):

    # --- returns the SFZH for an instantaneous SFH and single metallicity. At the moment this just chooses a single metallicity.

    SFZH = np.zeros((len(log10ages), len(log10Zs)))

    # --- determine the metallicity bin in which to place all the SF

    iZ = (np.abs(log10Zs - p['log10Z'])).argmin()

    ia = (np.abs(log10ages - p['log10age'])).argmin()

    SFZH[ia, iZ] = 10**p['log10M*']

    SFR = 0.0

    return SFZH, SFR
