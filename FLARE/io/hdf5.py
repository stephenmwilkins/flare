import os
import sys
import numpy as np
import h5py


def get_name_shape(name, item):
    shape = ''
    if hasattr(item, 'value'):
        shape = item.shape
    print(name, shape)







def merge(filenames, output_filename):


    # ---------------------------------------------
    # --- create a new HDF5 file with with the same structure but total length

    os.system(f'cp {filenames[0]}.h5 {output_filename}.h5')

    hf = h5py.File(f'{output_filename}.h5', 'a') # --- create new


        hf.attrs['total'] = N_total
        hf.attrs['detected'] = N_detected

        def create_new(name, item):
            if name.split('/')[0] == 'obs':
                N = N_detected
            elif name.split('/')[0] == 'EAZY':
                N = N_detected
            else:
                N = N_total
            if hasattr(item, 'value'):
                hf.create_dataset(name, (N,))


        Files[0].visititems(create_new)

        # ---------------------------------------------
        # --- append each entry

        def append(name, item):

            if hasattr(item, 'value'):

                if name.split('/')[0] == 'obs':
                    hf[name][running_n_detected:running_n_detected+n_detected] = item.value
                elif name.split('/')[0] == 'EAZY':
                    hf[name][running_n_detected:running_n_detected+n_detected] = item.value
                else:
                    hf[name][running_n_total:running_n_total+n_total] = item.value

        running_n_total = 0
        running_n_detected = 0

        for i,file in enumerate(Files):

            # print(FileNames[i])

            # file.visititems(get_name_shape)

            n_total = file.attrs['total']
            n_detected = file.attrs['detected']

            file.visititems(append)

            running_n_total += n_total
            running_n_detected += n_detected

        # ---------------------------------------------
        # --- print structure of file with shape of each object



        # hf.visititems(get_name_shape)


        # print(hf['obs/HST.WFC3.f160w/sizes/COG'].value)


        hf.close()

        print('deleting files')
        for FileName in FileNames: os.remove(f'{DataFolder}/individual/{FileName}')
