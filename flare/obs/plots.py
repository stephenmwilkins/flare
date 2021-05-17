
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# plt.style.use('simple') # --- makes nicer plots




def make_significance_plot(img, threshold = 2.5, show = False, filename = False, imsize = 1):

    fig = plt.figure(figsize=(imsize, imsize))
    ax = fig.add_axes([0,0,1,1])

    ax.axis('off')

    sig = (img.sci/img.noise)

    ax.imshow(sig, cmap = cm.Greys, vmin = -5.0, vmax = 5.0, origin = 'lower')
    ax.imshow(np.ma.masked_where(sig <= threshold, sig), cmap = cm.plasma, vmin = threshold, vmax = 100, origin = 'lower')

    if filename:
        plt.savefig(filename)
    if show:
        plt.show()

    plt.close(fig)




def make_significance_plots(imgs, threshold = 2.5):

    n = len(imgs)

    fig, axes = plt.subplots(1, n, figsize = (4*n,4))
    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.01, hspace=0.0)

    for ax, (filter, img) in zip(axes, imgs.items()):

        sig = (img.sci/img.noise)

        ax.set_axis_off()
        ax.imshow(sig, cmap = cm.Greys, vmin = -5.0, vmax = 5.0, origin = 'lower')
        ax.imshow(np.ma.masked_where(sig <= threshold, sig), cmap = cm.plasma, vmin = threshold, vmax = 100, origin = 'lower')

    plt.show()
    plt.close(fig)




def make_segm_plot(segm, imsize = 1, filename = False, show = False):

    fig, ax = plt.subplots(1, 1, figsize = (imsize,imsize))

    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.0, hspace=0.0)

    new_cmap = rand_cmap(int(np.max(segm)), type='bright', first_color_black=True, last_color_black=False, verbose=False)

    ax.imshow(segm, cmap = new_cmap, origin = 'lower')

    ax.set_axis_off()

    if filename:
        plt.savefig(filename)
    if show:
        plt.show()

    plt.close(fig)



def make_plots(imgs, threshold = 2.5, signficance_plot = False, filter_label = False, filename = False, show = False, use_vmax = True, fixed_range = False, imsize = 1, frame = True):

    n = len(imgs)

    if show:
        imsize = 4
    else:
        imsize = imsize

    if hasattr(next(iter(imgs.values())), 'sci'):
        fig, axes = plt.subplots(1, n, figsize = (n*imsize,1*imsize), dpi = next(iter(imgs.values())).sci.shape[0])
    else:
        fig, axes = plt.subplots(1, n, figsize = (n*imsize,1*imsize))

    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.0, hspace=0.0)

    if type(signficance_plot) != list: signficance_plot = [signficance_plot]*n

    if hasattr(next(iter(imgs.values())), 'sci'):
        if fixed_range:
            vmax = np.max([np.max(img.sci) for img in imgs.values()])
    else:
        if fixed_range:
            vmax = np.max([np.max(img) for img in imgs.values()])


    for ax, (filter, img), sig_plot in zip(axes, imgs.items(), signficance_plot):

        if frame:
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
        else:
            ax.set_axis_off()

        if filter_label: ax.text(0.5, 0.9, filter, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize = 8, color = '1.0')

        if sig_plot:

            sig = (img.sci/img.noise)

            ax.imshow(sig, cmap = cm.Greys, vmin = -5.0, vmax = 5.0, origin = 'lower')
            ax.imshow(np.ma.masked_where(sig <= threshold, sig), cmap = cm.plasma, vmin = threshold, vmax = 100, origin = 'lower')

        else:

            new_cmap = rand_cmap(np.max(img), type='bright', first_color_black=True, last_color_black=False, verbose=False)

            if fixed_range:
                vmin = 0.0
            else:
                vmin = None
                vmax = None

            if hasattr(img, 'sci'):
                ax.imshow(img.sci, cmap = cm.viridis, origin = 'lower', vmin = vmin, vmax = vmax)
            else:
                ax.imshow(img, cmap = new_cmap, origin = 'lower', vmin = vmin, vmax = vmax) # --- assumes img is just a 2D array




    if filename:
        plt.savefig(filename)
    if show:
        plt.show()

    plt.close(fig)



def rand_cmap(nlabels, type='bright', first_color_black=True, last_color_black=False, verbose=True):
    """
    Creates a random colormap to be used together with matplotlib. Useful for segmentation tasks
    :param nlabels: Number of labels (size of colormap)
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :param first_color_black: Option to use first color as black, True or False
    :param last_color_black: Option to use last color as black, True or False
    :param verbose: Prints the number of labels and shows the colormap. True or False
    :return: colormap for matplotlib
    """
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys



    if type not in ('bright', 'soft'):
        print ('Please choose "bright" or "soft" for type')
        return

    if verbose:
        print('Number of labels: ' + str(nlabels))

    # Generate color map for bright colors, based on hsv
    if type == 'bright':
        randHSVcolors = [(np.random.uniform(low=0.0, high=1),
                          np.random.uniform(low=0.2, high=1),
                          np.random.uniform(low=0.9, high=1)) for i in range(nlabels)]

        # Convert HSV list to RGB
        randRGBcolors = []
        for HSVcolor in randHSVcolors:
            randRGBcolors.append(colorsys.hsv_to_rgb(HSVcolor[0], HSVcolor[1], HSVcolor[2]))

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]

        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Generate soft pastel colors, by limiting the RGB spectrum
    if type == 'soft':
        low = 0.6
        high = 0.95
        randRGBcolors = [(np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high)) for i in range(nlabels)]

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Display colorbar
    if verbose:
        from matplotlib import colors, colorbar
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(15, 0.5))

        bounds = np.linspace(0, nlabels, nlabels + 1)
        norm = colors.BoundaryNorm(bounds, nlabels)

        cb = colorbar.ColorbarBase(ax, cmap=random_colormap, norm=norm, spacing='proportional', ticks=None,
                                   boundaries=bounds, format='%1i', orientation=u'horizontal')

    return random_colormap


def COG_plots(Properties, ModelProperties = False, filename = False, show = False):

    nfilters = len(Properties)

    fig, axes = plt.subplots(1, nfilters, figsize = (3*(nfilters),3), dpi = 200)

    plt.subplots_adjust(left=0.025, top=0.85, bottom=0.2, right=0.9, wspace=0.2, hspace=0.0)

    for ax, (filter, properties) in zip(axes, Properties.items()):

        ax.set_title(filter, fontsize = 10)

        ax.plot(properties['photometry']['aperture'].radii, properties['photometry']['aperture'].flux, c = '0.5', label = 'curve-of-growth')
        ax.axvline(properties['photometry']['aperture'].optimum_radius, color = '0.5', alpha = 0.5)

        if ModelProperties is not False:

            ax.axhline(ModelProperties[filter]['photometry']['total'].flux, color = '0.5', alpha = 0.5, ls = ':')
            ax.plot(ModelProperties[filter]['photometry']['aperture'].radii, ModelProperties[filter]['photometry']['aperture'].flux, c = '0.5', ls = ':', label = 'true curve-of-growth')



        del properties['photometry']['aperture']

        color_idx = np.linspace(0, 1, len(properties['photometry']))

        for c_idx, (phot_type, p) in zip(color_idx, properties['photometry'].items()):

            ax.axhline(p.flux, label = phot_type, color = cm.viridis(c_idx))
            ax.axhspan(p.flux-p.error, p.flux+p.error, color = cm.viridis(c_idx), alpha=0.5)
            if phot_type == 'optimum_aperture': ax.axvline(p.radius, color = cm.viridis(c_idx), alpha = 0.5)

    ax.legend(bbox_to_anchor=(1.1, 1.0), fontsize = 8)

    if filename:
        plt.savefig(filename)
    if show:
        plt.show()

    plt.close(fig)



def SED_plot(Properties,  ModelProperties = False, FilterInfo = False, phot_type = 'optimum_aperture', filename = False, show = False):


    # if not FilterInfo:

    fig, ax = plt.subplots(1, 1, figsize = (3,2), dpi = 200)
    plt.subplots_adjust(left=0.2, top=0.85, bottom=0.25, right=0.9, wspace=0.2, hspace=0.0)

    color_idx = np.linspace(0, 1, len(Properties))

    for c_idx, (filter, properties) in zip(color_idx, Properties.items()):

        pivwv = FilterInfo[filter].pivwv()/1E4

        ax.scatter(pivwv, properties['photometry'][phot_type].flux, color = cm.viridis(c_idx))
        ax.plot([pivwv]*2, [properties['photometry'][phot_type].flux - properties['photometry'][phot_type].error, properties['photometry'][phot_type].flux + properties['photometry'][phot_type].error], color = 'k', lw = 1)

        if ModelProperties is not False:

            ax.scatter(pivwv, ModelProperties[filter]['photometry']['total'].flux, color = cm.viridis(c_idx), alpha = 0.5)



    ax.set_xlabel(r'$\lambda/\mu m$')
    ax.set_ylabel(r'$f_{\nu}/nJy$')

    if filename:
        plt.savefig(filename)
    if show:
        plt.show()

    plt.close(fig)




def size_plot(img, p, ExclusionMask, threshold = 2.5, signficance_plot = False, filename = False, show = False, add_contours = False):


    width = img.sci.shape[0]

    fig, ax = plt.subplots(1, 1, figsize = (3,3), dpi = width*2)
    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.01, hspace=0.0)

    ax.set_axis_off()

    sig = (img.sci/img.noise)

    ax.imshow(sig, cmap = cm.Greys, vmin = -5.0, vmax = 5.0, origin = 'lower')
    ax.imshow(np.ma.masked_where(sig <= threshold, sig), cmap = cm.plasma, vmin = threshold, vmax = 100, origin = 'lower')

    k = 2.5

    # --- make mask image including Kron Mask and Exclusion mask
    x = np.linspace(-(width//2), (width//2), width)
    X, Y = np.meshgrid(x, x)
    R2 = X**2 + Y**2
    alpha = np.zeros(img.sci.shape)
    alpha[R2>(k*p['kron_radius'])**2] = 1
    alpha[img.sci<p['sizes']['pixel'].minflux] = 1
    alpha[ExclusionMask] = 1

    # --- make RGBA image which is white everywhere but transparent in regions that are counted
    RGBA = np.ones((width, width, 4))
    RGBA[:,:,3] = alpha

    ax.imshow(RGBA, origin = 'lower', alpha = 0.8)

    # --- add contours around pixels contributing to r_pix (VERY SLOW FOR SOME REASON)
    if add_contours:
        image = alpha[::-1]
        f = lambda x,y: image[int(y),int(x) ]
        g = np.vectorize(f)

        x = np.linspace(0,image.shape[1], image.shape[1]*100)
        y = np.linspace(0,image.shape[0], image.shape[0]*100)
        X, Y= np.meshgrid(x[:-1],y[:-1])
        Z = g(X[:-1],Y[:-1])

        ax.contour(Z[::-1], [0.5], colors='k', linewidths=[2], extent=[0-0.5, x[:-1].max()-0.5,0-0.5, y[:-1].max()-0.5])


    # --- add k*r_Kron radius
    kKronRadius = plt.Circle((width//2,width//2), k*p['kron_radius'], alpha = 0.2)
    ax.add_artist(kKronRadius)
    # --- add COG radius
    COGRadius = plt.Circle((width//2,width//2), p['sizes']['COG'].radius, alpha = 1.0, fill = False, lw = 1, color='1')
    ax.add_artist(COGRadius)

    # --- add COG radius
    PixelRadius = plt.Circle((width//2,width//2), p['sizes']['pixel'].radius, alpha = 1.0, fill = False, lw = 1, ls = '--', color='1')
    ax.add_artist(PixelRadius)


    if filename:
        plt.savefig(filename)
    if show:
        plt.show()

    plt.close(fig)
