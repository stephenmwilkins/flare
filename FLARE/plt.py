


import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

fancy = lambda x: r'$\rm '+x.replace(' ','\ ')+'$'
ml = lambda x: r'$\rm '+x+'$'

plt.style.use('simple') # --- makes nicer plots

def simple(size = 3):

    fig = plt.figure(figsize = (size, size))

    left  = 0.15
    height = 0.8
    bottom = 0.15
    width = 0.8

    ax = fig.add_axes((left, bottom, width, height))

    return fig, ax


def simple_fig():

    fig = plt.figure(figsize = (3., 3.))

    left  = 0.15
    height = 0.8
    bottom = 0.15
    width = 0.8

    ax = fig.add_axes((left, bottom, width, height))

    return fig, ax

def simple_fig():

    fig = plt.figure(figsize = (3., 3.))

    left  = 0.15
    height = 0.8
    bottom = 0.15
    width = 0.8

    ax = fig.add_axes((left, bottom, width, height))

    return fig, ax


def hist(set_axis_off = True):



    fig = plt.figure(figsize = (4., 4.))

    left  = 0.15
    height = 0.65
    bottom = 0.15
    width = 0.65

    ax = fig.add_axes((left, bottom, width, height)) # main panel

    axx = fig.add_axes((left, bottom + height, width, 0.15))  # x-axis hist panel

    axy = fig.add_axes((left + width, bottom, 0.15, height))  # y-axis hist panel  # y-axis hist panel

    if set_axis_off:
        axx.axis('off')
        axy.axis('off')

    return fig, ax, axx, axy
