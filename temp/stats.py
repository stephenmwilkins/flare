
import numpy as np
import scipy.stats

          
def poisson_confidence_interval(n, p=0.68):
    
    #http://ms.mcmaster.ca/peter/s743/poissonalpha.html
        
    #e.g. p=0.68 for 1 sigma
    
    #agrees with http://hyperphysics.phy-astr.gsu.edu/hbase/math/poifcn.html
        
    # see comments on JavaStat page
    
    #  scipy.stats.chi2.ppf((1.-p)/2.,2*n)/2. also known
    
    if n>0:   
        interval=(scipy.stats.chi2.ppf((1.-p)/2.,2*n)/2.,scipy.stats.chi2.ppf(p+(1.-p)/2.,2*(n+1))/2.)       
    
    else:
        
        #this bit works out the case for n=0
        
        ul=(1.-p)/2.
        
        prev=1.0
        for a in np.arange(0.,5.0,0.001):
        
            cdf=scipy.stats.poisson.cdf(n,a)
        
            if cdf<ul and prev>ul:
                i=a
        
            prev=cdf
        
        interval=(0.,i)
    
    
    return interval   