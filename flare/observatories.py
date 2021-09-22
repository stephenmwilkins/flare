
import numpy as np


class observatory():

    def __init__(self, name, instruments):
        self.name = name
        self.instruments = instruments
        self.instrument = {inst: instrument(name, inst) for inst in instruments}
        self.filters = []
        self.filterIDs = []

    def update(self):

        for inst in self.instruments:
            self.filters += self.instrument[inst].filters
            self.filterIDs += [f'{self.name}.{inst}.{f}' for f in self.instrument[inst].filters]
            setattr(self, inst, self.instrument[inst])



class instrument():

    def __init__(self, obs, inst):
        self.obs = obs
        self.inst = inst
        self.channel = {}
        self.channels = None
        self.pixel_scale = None
        self.zeropoints = None


    def get_filter_IDs(self, exclude = False):

        if not exclude:
            return ['.'.join([self.obs, self.inst, f]) for f in self.filters]
        else:
            filter_IDs = []
            for f in self.filters:
                if f[-1].upper() not in exclude:
                    filter_IDs.append('.'.join([self.obs, self.inst, f]))
            return filter_IDs

    def get_filters(self):
        filters = []
        for c in self.channels:
            for f in self.channel[c].filters:
                filters.append(f)

        return filters

    def get_zeropoints(self):
        zeropoints = {}
        for c in self.channels:
            for f in self.channel[c].filters:
                zeropoints[f] = self.channel[c].zeropoints[f]

        return zeropoints


class channel: pass



FAKE = observatory('FAKE', ['FAKE', 'TH'])
FAKE.instrument['FAKE'] = instrument('FAKE', 'FAKE')
FAKE.instrument['FAKE'].filters = ['1500','1500n','1600','2000','2000n','2500','2500n','3000n','B','BBa','BBb','D4000a','D4000b','FUV','H','I','J','K','MUV','NUV','R','U','Un','V','Vn','Vw','Y','Z']
FAKE.instrument['FAKE'].pixel_scale = None
FAKE.instrument['FAKE'].zeropoints = None
FAKE.instrument['TH'] = instrument('FAKE', 'TH')
FAKE.instrument['TH'].filters = FAKE.instrument['FAKE'].filters
FAKE.instrument['TH'].pixel_scale = None
FAKE.instrument['TH'].zeropoints = None



Webb = observatory('Webb', ['NIRCam', 'MIRI'])
Webb.instrument['NIRCam'].channels = ['Short', 'Long']
Webb.instrument['NIRCam'].channel['Short'] =  channel()
Webb.instrument['NIRCam'].channel['Short'].filters = ['F070W','F090W','F115W','F150W','F200W','F140M','F162M','F182M','F210M']
Webb.instrument['NIRCam'].channel['Short'].pixel_scale = 0.031
Webb.instrument['NIRCam'].channel['Short'].zeropoints = None
Webb.instrument['NIRCam'].Short = Webb.instrument['NIRCam'].channel['Short']
Webb.instrument['NIRCam'].channel['Long'] =  channel()
Webb.instrument['NIRCam'].channel['Long'].filters = ['F277W','F356W','F444W'] + ['F250M','F300M','F360M','F410M','F430M','F460M','F480M']
Webb.instrument['NIRCam'].channel['Long'].pixel_scale = 0.063
Webb.instrument['NIRCam'].channel['Long'].zeropoints = None
Webb.instrument['NIRCam'].filters = Webb.instrument['NIRCam'].get_filters()
Webb.instrument['MIRI'] = instrument('Webb', 'MIRI')
Webb.instrument['MIRI'].filters = ['F560W','F770W','F1000W','F1130W','F1280W','F1500W','F1800W','F2100W','F2550W']
Webb.instrument['MIRI'].pixel_scale = None
Webb.instrument['MIRI'].zeropoints = None

Roman = observatory('Roman', ['WFI'])
Roman.instrument['WFI'].filters = ['F062', 'F087', 'F106', 'F129', 'F146', 'F158', 'F184']
Roman.instrument['WFI'].pixel_scale = 0.11
Roman.instrument['WFI'].zeropoints = None

Hubble = observatory('Hubble', ['ACS','WFC3'])
Hubble.instrument['ACS'].filters = ['f435w', 'f606w', 'f775w', 'f814w', 'f850lp']
Hubble.instrument['ACS'].zeropoints = {'f435w': 25.684, 'f606w': 26.505, 'f775w': 25.678, 'f814w': 25.959, 'f850lp': 24.867}
Hubble.instrument['ACS'].pixel_scale = 0.05
Hubble.instrument['WFC3'].channels = ['UV','NIR']
Hubble.instrument['WFC3'].channel['UV'] = channel()
Hubble.instrument['WFC3'].channel['UV'].filters = ['f225w','f275w','f336w']
Hubble.instrument['WFC3'].channel['UV'].zeropoints = {'f225w': 0.0, 'f275w': 0.0, 'f336w': 0.0}
Hubble.instrument['WFC3'].channel['UV'].pixel_scale = 0.13 # might be wrong?
Hubble.instrument['WFC3'].channel['NIR'] = channel()
Hubble.instrument['WFC3'].channel['NIR'].filters = ['f105w', 'f125w', 'f140w', 'f160w']
Hubble.instrument['WFC3'].channel['NIR'].zeropoints = {'f105w': 26.269, 'f125w': 26.230, 'f140w': 26.452, 'f160w': 25.946}
Hubble.instrument['WFC3'].channel['NIR'].pixel_scale = 0.13
Hubble.instrument['WFC3'].pixel_scale = 0.13
Hubble.instrument['WFC3'].filters = Hubble.instrument['WFC3'].get_filters()
Hubble.instrument['WFC3'].zeropoints = Hubble.instrument['WFC3'].get_zeropoints()

Euclid = observatory('Euclid', ['VIS','NISP'])
Euclid.instrument['NISP'].filters = ['Y','J','H']
Euclid.instrument['NISP'].pixel_scale = 0.3
Euclid.instrument['NISP'].zeropoints = None
Euclid.instrument['VIS'].filters = ['VIS']
Euclid.instrument['VIS'].pixel_scale = 0.3 # probably wrong
Euclid.instrument['VIS'].zeropoints = None

Spitzer = observatory('Spitzer', ['IRAC'])
Spitzer.instrument['IRAC'].filters = ['ch1', 'ch2']
Spitzer.instrument['IRAC'].pixel_scale = 1.22
Spitzer.instrument['IRAC'].zeropoints = None

VISTA = observatory('VISTA', ['VIRCAM'])
VISTA.instrument['VIRCAM'].filters = ['Z', 'Y', 'J', 'H', 'Ks']
VISTA.instrument['VIRCAM'].pixel_scale = None
VISTA.instrument['VIRCAM'].zeropoints = None






observatory_list = ['Webb', 'Roman', 'Hubble', 'Euclid', 'Spitzer','VISTA','FAKE']
for observatory in observatory_list:
    globals()[observatory].update()

observatories = {observatory: globals()[observatory] for observatory in observatory_list}



# ---


filter_info = {}
for observatory in observatory_list:
    for instrument in observatories[observatory].instruments:
        channels = observatories[observatory].instrument[instrument].channels
        if channels:
            for channel in channels:
                for filter in observatories[observatory].instrument[instrument].channel[channel].filters:
                    f = f'{observatory}.{instrument}.{filter}'
                    filter_info[f] = {}
                    filter_info[f]['pixel_scale'] = observatories[observatory].instrument[instrument].channel[channel].pixel_scale

                    if observatories[observatory].instrument[instrument].channel[channel].zeropoints:
                        filter_info[f]['zeropoint'] = observatories[observatory].instrument[instrument].channel[channel].zeropoints[filter]
                        filter_info[f]['nJy_to_es'] = 1E-9 * 10**(0.4*(filter_info[f]['zeropoint']-8.9))

        else:
            for filter in observatories[observatory].instrument[instrument].filters:
                f = f'{observatory}.{instrument}.{filter}'
                filter_info[f] = {}
                filter_info[f]['pixel_scale'] = observatories[observatory].instrument[instrument].pixel_scale

                if observatories[observatory].instrument[instrument].zeropoints:
                    filter_info[f]['zeropoint'] = observatories[observatory].instrument[instrument].zeropoints[filter]
                    filter_info[f]['nJy_to_es'] = 1E-9 * 10**(0.4*(filter_info[f]['zeropoint']-8.9))
                    # print(f, filter_info[f]['pixel_scale'])




def get_zeropoint(filter):

    observatory, instrument, f = filter.split('.')

    return observatories[observatory].instrument[instrument].zeropoints[f]


def get_nJy_to_es(filter):

    zeropoint = get_zeropoint(filter)

    return 1E-9 * 10**(0.4*(zeropoint-8.9))
