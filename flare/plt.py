


import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

fancy = lambda x: r'$\rm '+x.replace(' ','\ ')+'$'
ml = lambda x: r'$\rm '+x+'$'

rcParams = {}
rcParams['savefig.dpi'] = 300
rcParams['path.simplify'] = True

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'stixsans'
rcParams['text.usetex'] = False
rcParams['font.size'] = 9
rcParams['mathtext.fontset'] = 'stixsans'

rcParams['axes.linewidth'] = 0.5

rcParams['xtick.major.size'] = 3
rcParams['ytick.major.size'] = 3
rcParams['xtick.minor.size'] = 1.5
rcParams['ytick.minor.size'] = 1.5
rcParams['xtick.labelsize'] = 7
rcParams['ytick.labelsize'] = 7
rcParams['ytick.direction'] = 'in'
rcParams['xtick.direction'] = 'in'
rcParams['ytick.minor.visible'] = True
rcParams['xtick.minor.visible'] = True
rcParams['xtick.major.width'] = 0.25
rcParams['ytick.major.width'] = 0.25
rcParams['xtick.minor.width'] = 0.25
rcParams['ytick.minor.width'] = 0.25

rcParams['grid.alpha'] = 0.1
rcParams['grid.color'] = 'k'
rcParams['grid.linestyle'] = '-'
rcParams['grid.linewidth'] = 0.8

rcParams['legend.frameon'] = False
# rcParams[''] =
# rcParams[''] =
# rcParams[''] =


mpl.rcParams.update(rcParams)


# print(mpl.rcParams)



def simple_sm(size = 2.5):

    fig = plt.figure(figsize = (size, size))

    left  = 0.25
    height = 0.7
    bottom = 0.25
    width = 0.7

    ax = fig.add_axes((left, bottom, width, height))

    return fig, ax

def simple(size = 3.5):

    fig = plt.figure(figsize = (size, size))

    left  = 0.15
    height = 0.8
    bottom = 0.15
    width = 0.8

    ax = fig.add_axes((left, bottom, width, height))

    return fig, ax


def simple_fig():

    return simple()



def single_wcbar(base_size = 3.5):

    left  = 0.15
    height = 0.8
    bottom = 0.15
    width = 0.65

    fig = plt.figure(figsize = (base_size, base_size*width/height))

    ax = fig.add_axes((left, bottom, width, height))
    cax = fig.add_axes([left+width, bottom, 0.03, height])

    return fig, ax, cax

def simple_wcbar(base_size = 3.5):
    return single_wcbar(base_size = base_size)


def simple_wcbar_whist(base_size = 3.5):

    left  = 0.15
    height = 0.8
    bottom = 0.15
    width = 0.55
    hwidth = 0.15

    fig = plt.figure(figsize = (base_size, base_size*width/height))

    ax = fig.add_axes((left, bottom, width, height))
    hax = fig.add_axes([left+width, bottom, hwidth, height])
    cax = fig.add_axes([left+width+hwidth, bottom, 0.03, height])

    return fig, ax, cax, hax

def simple_wcbart_whist(base_size = 3.5):

    # --- with colour bar on top

    left  = 0.15
    height = 0.6
    bottom = 0.15
    width = 0.6
    hwidth = 0.15

    fig = plt.figure(figsize = (base_size, base_size*width/height))

    ax = fig.add_axes((left, bottom, width, height))
    hax = fig.add_axes([left+width, bottom, hwidth, height])
    cax = fig.add_axes([left, bottom+height, width, 0.03])

    return fig, ax, cax, hax



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
