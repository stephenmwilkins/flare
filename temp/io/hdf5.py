import os
import sys
import numpy as np
import h5py






def get_name_shape(name, item):
    shape = ''
    if hasattr(item, 'value'):
        shape = item.shape
    print(name, shape)







def merge(filenames, output_filename, delete_files = False):


    # ---------------------------------------------
    # --- create a new HDF5 file with with the same structure but longer

    hf = h5py.File(f'{output_filename}.h5', 'w') # --- create new

    N_files = len(filenames)
    N_file = int(filenames[0].split('_')[-2])
    N_total = N_files * N_file



    print(N_total)

    def create_new(name, item):
        if hasattr(item, 'value'):
            hf.create_dataset(name, (N_total,))


    first_hf = h5py.File(f'{filenames[0]}.h5', 'r')

    for key, value in first_hf.attrs.items():
        hf.attrs[key] = value

    first_hf.visititems(create_new)
    first_hf.close()

    # ---------------------------------------------
    # --- append each entry

    def append(name, item):
        if hasattr(item, 'value'):
            hf[name][i*N_file:(i+1)*N_file] = item.value

    for i, filename in enumerate(filenames):
        next_hf = h5py.File(f'{filename}.h5', 'r')
        next_hf.visititems(append)
        next_hf.close()

    hf.close()

    # --- delete files

    if delete_files:
        pass
