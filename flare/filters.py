
from .core import *
from . import observatories

class empty: pass

import numpy as np

import os
this_dir, this_filename = os.path.split(__file__)






def add_filters(filters, new_lam = False,  filter_path = f'{this_dir}/data/filters/'):

    F = {f: filter(f, new_lam, filter_path) for f in filters}

    F['filters'] = filters

    return F


class filter():

    def __init__(self, f, new_lam = False, filter_path =  f'{this_dir}/data/filters/'):


        observatory =  f.split('.')[0]
        instrument =  f.split('.')[1]
        fs = f.split('.')[2]

        # l, t are the original wavelength and transmission grid

        self.l, self.t = np.loadtxt(filter_path + '/'+'/'.join(f.split('.'))+'.txt', skiprows = 1).T

        self.t[self.t<0] = 0

        if instrument == 'NIRCAM': self.l *= 1E4 # convert from microns to \AA
        if instrument == 'NIRCam': self.l *= 1E4 # convert from microns to \AA
        if instrument == 'WFI':
            self.l *= 1E4 # convert from microns to \AA
            self.t /= 3. # WFI is effective area

        self.zeropoint = False
        self.nJy_to_es = False

        try:
            zeropoints = observatories.observatories[observatory].instrument[instrument].zeropoints
            if zeropoints:
                self.zeropoint = zeropoints[fs]
                self.nJy_to_es = 1E-9 * 10**(0.4*(self.zeropoint-8.9))

        except:
            self.zeropoint = False



        if isinstance(new_lam, np.ndarray):

            self.lam = new_lam
            self.T = np.interp(self.lam, self.l, self.t)

        else:

            self.lam = self.l
            self.T = self.t

        # self.info = info[f]


    def apply_new_lam(self, new_lam):

        self.lam = new_lam
        self.T = np.interp(self.lam, self.l, self.t)


    #http://stsdas.stsci.edu/stsci_python_epydoc/SynphotManual.pdf   #got to page 42 (5.1)

    def pivwv(self): # -- calculate pivot wavelength using original l, t

        return np.sqrt(np.trapz(self.l * self.t , x = self.l)/np.trapz(self.t / self.l, x = self.l))


    def pivT(self): # -- calculate pivot wavelength using original l, t

        return np.interp(self.pivwv(), self.l, self.t)

    def meanwv(self):

        return np.exp(np.trapz(np.log(self.l) * self.t / self.l, x = self.l)/np.trapz(self.t / self.l, x = self.l))

    def bandw(self):

        return self.meanwv() * np.sqrt(np.trapz((np.log(self.l/self.meanwv())**2)*self.t/self.l,x=self.l))/np.sqrt(np.trapz(self.t/self.l,x=self.l))

    def fwhm(self):

        return np.sqrt(8.*np.log(2))*self.bandw()

    def Tpeak(self):

        return np.max(self.t)

    def rectw(self):

        return np.trapz(self.t, x=self.l)/self.Tpeak()

    def mnmx(self):

        return (self.min(), self.max())

    def max(self):

        return self.l[self.T>1E-2][-1]

    def min(self):

        return self.l[self.T>1E-2][0]




def create_EAZY_filter_res(F, filter_res_file = 'FILTER.RES'):

    o = []

    for f in F['filters']:

        o.append('{n} {f}\n'.format(n=len(F[f].l), f=f)) # filter header

        for i,l in enumerate(F[f].l):
            o.append('{i:>5}   {l}   {t}\n'.format(i=i+1, l=l,t=F[f].t[i]))

    open(filter_res_file,'w').writelines(o)









# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------

# below is deprecated in favour of the observatories sub-module


# --- artificial

FAKE = ['FAKE.FAKE.'+f for f in ['1500','2500','Uth','Bth','Vth','Ith','Zth','Yth','Jth','Hth']]
TH = ['FAKE.TH.'+f for f in ['FUV','MUV', 'NUV','U','B','V','R','I','Z','Y','J','H','K']]

# --- Roman

Roman = ['Roman.Roman.'+f for f in ['F062', 'F087', 'F106', 'F129', 'F146', 'F158', 'F184']]
WFIRST = Roman

# --- Hubble

WFC3UVIS_W = ['Hubble.WFC3.'+f for f in ['f225w','f275w','f336w']]
WFC3NIR_W = ['Hubble.WFC3.'+f for f in ['f105w', 'f125w', 'f140w', 'f160w']]
WFC3 = WFC3UVIS_W + WFC3NIR_W

ACS_W = ['HST.ACS.'+f for f in ['f435w', 'f606w', 'f775w', 'f814w', 'f850lp']]
ACS = ACS_W

HST = ACS + WFC3
Hubble = HST

# --- Spitzer

# IRAC = ['Spitzer.IRAC.'+f for f in ['ch1', 'ch2', 'ch3', 'ch4']]
IRAC = ['Spitzer.IRAC.'+f for f in ['ch1', 'ch2']]
Spitzer = IRAC

# --- Euclid

Euclid_VIS = ['Euclid.VIS.'+f for f in ['VIS']]
Euclid_NISP = ['Euclid.NISP.'+f for f in ['Y','J','H']]
Euclid = Euclid_VIS + Euclid_NISP

# --- Subaru

HSC = ['Subaru.HSC.'+f for f in ['g','r','i','z','y']]
Subaru = HSC

# --- Webb

NIRCam_s_W = ['JWST.NIRCAM.'+f for f in ['F070W','F090W','F115W','F150W','F200W']]
NIRCam_s_M = ['JWST.NIRCAM.'+f for f in ['F140M','F162M','F182M','F210M']]
NIRCam_l_W = ['JWST.NIRCAM.'+f for f in ['F277W','F356W','F444W']]
NIRCam_l_M = ['JWST.NIRCAM.'+f for f in ['F250M','F300M','F360M','F410M','F430M','F460M','F480M']]
NIRCam_s = NIRCam_s_W + NIRCam_s_M
NIRCam_l = NIRCam_l_W + NIRCam_l_M
NIRCam_W = NIRCam_s_W + NIRCam_l_W
NIRCam =  NIRCam_s + NIRCam_l
MIRI = ['JWST.MIRI.'+f for f in ['F560W','F770W','F1000W','F1130W','F1280W','F1500W','F1800W','F2100W','F2550W']]
Webb = NIRCam + MIRI
JWST = Webb

CEERS = ['JWST.NIRCAM.'+f for f in ['F115W','F150W','F200W', ]]

# --- All *real* filters

all_filters = FAKE + TH + Hubble + Spitzer + Euclid + Subaru + Webb + Roman
