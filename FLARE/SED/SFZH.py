

import numpy as np


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




def constant(SPS, p):


    # --- returns the SFZH for a constant SFH and single metallicity. At the moment this just chooses a single metallicity.
    
    

    a0 = 10**p['log10_duration']

    log10ages = SPS.grid['log10age']

    SFZH = np.zeros(SPS.grid['log10Q'].shape)

    # --- determine the metallicity bin in which to place all the SF

    iZ = (np.abs(SPS.grid['log10Z'] - p['log10Z'])).argmin()
    
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
    


