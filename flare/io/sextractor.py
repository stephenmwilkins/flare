

import numpy as np

def read(filename):

    d = open(filename, 'r').readlines()

    columns = []
    attributes = []

    for i, l in enumerate(d):
        if l[0]=='#':
            columns.append(l.split()[2])
            if len(l.split())>3:
                attributes.append(' '.join(l.split()[3:]))
            else:
                attributes.append('')

    data = np.loadtxt(filename)


    cat = {}

    for i,col in enumerate(columns):
        cat[col] = data[:, i]

    return cat, attributes
