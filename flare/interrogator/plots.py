

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import cm
import pickle
import sys
import os
import decimal
from scipy.stats import kde
import scipy
import scipy.stats

# sys.path.insert(0, os.path.abspath('..'))


this_dir, this_filename = os.path.split(__file__)


parameter_labels = {}



parameter_labels['z'] = r'z'

parameter_labels['log10M*'] = r'\log_{10}(M^{*}/M_{\odot})'

parameter_labels['log10Z'] = r'\log_{10}(Z)'

parameter_labels['log10_duration'] = r'\log_{10}(t/yr)'

parameter_labels['SFR10'] = r'(SFR_{10}/M_{\odot} yr^{-1})'


# --- nebular emission parameters

# parameter_labels['fesc'] = r'f_{esc, LyC}'
parameter_labels['log10fesc'] = r'\log_{10}(f_{esc, LyC})'

# --- for Simple dust model

parameter_labels['log10tau_V'] = r'\log_{10}(\tau_{V})' 









class plots():

    def __init__(self, ID, output_dir = 'outputs/'):
            
        self.outdir = output_dir + '/'+ID
            
        self.ID = ID

        # ---- read in files

        self.summary = pickle.load(open('{0}/summary.p'.format(self.outdir), 'rb')) # --- basic info
        self.samples = pickle.load(open('{0}/samples.p'.format(self.outdir), 'rb')) # --- basic info
    
        
    def P_of_z(self):
    
        # --- quick plot showing P(z)
        
        
        plt.style.use('simple')
        
        fig = plt.figure(figsize = (3,3), frameon = False)

        # ---- gender variation

        left  = 0.15  
        bottom = 0.15
        height = 0.8
        width = 0.8

        ax = fig.add_axes((left, bottom, width, height))
        
        ax.hist(self.samples['z'],bins=np.arange(0.,15.,0.1))

        fig.savefig(self.outdir+'/Pz.pdf')
  
  
    def simple_triangle(self, parameters = False, bins = 20):

        if not parameters: parameters = self.summary.model_parameters 
        
        print(parameters)


        plot_priors = True
    
        n = len(parameters)
    
        # ---- initialise figure

        plt.rcParams['mathtext.fontset'] = 'stixsans'
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.size'] = 7 # perhaps should depend on number of parameters to be plotted

        plt.rcParams['ytick.labelsize'] = 5 # perhaps should depend on number of parameters to be plotted
        plt.rcParams['xtick.labelsize'] = 5 # perhaps should depend on number of parameters to be plotted
    
        plt.rcParams['ytick.direction'] = 'in'    # direction: in, out, or inout
        plt.rcParams['xtick.direction'] = 'in'    # direction: in, out, or inout
    
        plt.rcParams['ytick.minor.visible'] = True
        plt.rcParams['xtick.minor.visible'] = True
    

        fig, axes = plt.subplots(n,n, figsize = (1+n*0.7,1+n*0.7))

        left  = 0.125  # the left side of the subplots of the figure
        right = 0.9    # the right side of the subplots of the figure
        bottom = 0.1   # the bottom of the subplots of the figure
        top = 0.9      # the top of the subplots of the figure
        wspace = 0.0   # the amount of width reserved for blank space between subplots
        hspace = 0.0   # the amount of height reserved for white space between subplots
    
        fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

        # ---- generate histograms
    
        hist = {}
      
        # ----- generate histograms using the prior range
    
        for i, pi in enumerate(parameters): hist[pi] = {}
    
        for i, pi in enumerate(parameters):
        
            for j, pj in enumerate(parameters):
            
                # --- create 1D histogram
                
                
                IQR_fac = 3.

                if i == j:
            
                    IQR = np.percentile(self.samples[pi], 75.) - np.percentile(self.samples[pi], 25.)
            
                    range = [np.median(self.samples[pi]) - IQR*IQR_fac, np.median(self.samples[pi]) + IQR*IQR_fac]
            
                    hist[pi][pj] = np.histogram(self.samples[pi], bins = bins, range = range) 
            
            
                # --- create 2D histogram
            
                else:
                
                
                    IQRi = np.percentile(self.samples[pi], 75.) - np.percentile(self.samples[pi], 25.)
                    IQRj = np.percentile(self.samples[pj], 75.) - np.percentile(self.samples[pj], 25.)
                    
                    rangei = [np.median(self.samples[pi]) - IQRi*IQR_fac, np.median(self.samples[pi]) + IQRi*IQR_fac]
                    rangej = [np.median(self.samples[pj]) - IQRj*IQR_fac, np.median(self.samples[pj]) + IQRj*IQR_fac]
                
                    hist[pi][pj] = np.histogram2d(self.samples[pi], self.samples[pj], bins = bins, range = [rangei, rangej]) 

    










        # ---- loop over parameters

        for i in np.arange(n):
            for j in np.arange(n):
    
                axes[i,j].locator_params(axis = 'x', nbins=4)
                axes[i,j].locator_params(axis = 'y', nbins=7)
      
                pi = parameters[i]
                pj = parameters[j]
                    



                if i != n-1: axes[i,j].set_xticklabels([])
    
                if i == n-1: 
                    axes[i,j].set_xlabel(r'${\rm'+parameter_labels[pj]+'}$')
    
                if j == 0 and i!=0: 
                    axes[i,j].set_ylabel(r'${\rm'+parameter_labels[pi]+'}$')
                    axes[i,j].get_yaxis().set_label_coords(-0.1*n,0.5)
  

                if j < i:
    
                    H, xe, ye = hist[pj][pi]
                    

                    X, Y = np.meshgrid(xe, ye)
  
                    xlims = [xe[0], xe[-1]]
                    ylims = [ye[0], ye[-1]]
    
                    axes[i,j].set_xlim(xlims)
                    axes[i,j].set_ylim(ylims)
    
                    axes[i,j].pcolormesh(X, Y, H.T, cmap = 'plasma',linewidth=0,rasterized=True) #

                    if j != 0: axes[i,j].set_yticklabels([])

 
                    if self.summary.obs.true_p:
                
                        axes[i,j].axhline(self.summary.obs.true_p[pi], c = 'white', alpha = 0.5, lw = 0.5)
                        axes[i,j].axvline(self.summary.obs.true_p[pj], c = 'white', alpha = 0.5, lw = 0.5)
                    

                    




                elif j == i: 


                    med = np.median(self.samples[pi])
                    
                    P68 = [np.percentile(self.samples[pi], x) for x in [16.,84.]]


                    H, edges = hist[pj][pi]
                    
                    axes[i,j].set_xlim([edges[0], edges[-1]])
                    axes[i,j].set_yticklabels([])          
                    
                    # hist, bins, patches = axes[i,j].hist(self.samples[pi], 25, normed=1, lw = 0, facecolor='black', alpha=0.2, range = xlims)

                    axes[i,j].bar((edges[:-1]+edges[1:])/2, H, align='center', width=edges[1]-edges[0], color='k', alpha = 0.2)


                    # --- shade the prior region

                    if plot_priors:
                    
                        if pj in self.summary.prior_def.keys():
                        
                            if self.summary.prior_def[pj]['type'] == 'uniform':
                            
                                axes[i,j].fill_between(self.summary.prior_def[pj]['limits'], [np.max(H)*1.3]*2, [0,0], alpha = 0.05, color= 'black')
                                

                            
        
                    if self.summary.obs.true_p:
                        axes[i,j].axvline(self.summary.obs.true_p[pi], c = 'k', alpha = 0.5, lw = 0.5)
                       
                    axes[i,j].plot(P68, [np.max(H)*1.1]*2, c='k', lw=1)
                    axes[i,j].scatter(med, np.max(H)*1.1, c='k', s=5)  

                else:
    
                    axes[i,j].set_axis_off()
     
#         for axi in axes.flat:
#             axi.xaxis.set_major_locator(plt.MaxNLocator(3))
#             axi.yaxis.set_major_locator(plt.MaxNLocator(3))
 
        fig.savefig(self.outdir+'/triangle.pdf', dpi = 500)

        plt.close(fig)

    

    
    

        