


import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

fancy = lambda x: r'$\rm '+x.replace(' ','\ ')+'$'
ml = lambda x: r'$\rm '+x+'$'

plt.style.use('simple') # --- makes nicer plots 


def simple_fig():

    fig = plt.figure(figsize = (3., 3.))

    left  = 0.15
    height = 0.8
    bottom = 0.15 
    width = 0.8

    ax = fig.add_axes((left, bottom, width, height))
    
    return fig, ax