# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle, FancyBboxPatch

import planckStyle as s
from pylab import *
import numpy as np
import GetDistPlots, os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc, font_manager
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, \
     grid, savefig, show

outdir='/home/pettorin/codici/cosmomc_plots/output_getdist/clik/clik10/pca/plots/'
g = s.getSinglePlotter(plot_data='/home/pettorin/codici/cosmomc_plots/output_getdist/clik/clik10/pca/4bins/')
g.settings.setWithSubplotSize(4.0000)

# Load data
wt  = np.loadtxt('/home/pettorin/codici/cosmomc_plots/output_getdist/clik/clik10/pca/pythoncode/4bins/output/weights.txt')
pca = np.loadtxt('/home/pettorin/codici/cosmomc_plots/output_getdist/clik/clik10/pca/pythoncode/4bins/output/w_reconstructed.txt')

#f = plt.figure()
#ax = f.add_subplot(1,1,1)
#col = 'blue'
#z0 = pca[:,0] # lower limit
#z1 = pca[:,2] # upper limit
#w0 = pca[:,3] # mean w
#dw = pca[:,4]
#for j in range(len(z0)-1):
#    llc = (z0[j],w0[j]-2.*dw[j]) # lower left corner
#    dz  = z1[j]-z0[j]
#    dx  = 4.*dw[j]
#    bbox1 = Rectangle(llc,dz,dx,transform=ax.transData,ec=col,fc=col,fill=True,alpha=0.5)
#    ax.add_patch(bbox1)
#xlim((0,2.1))
#ylim((-2,1))
#f.show()

# Create the plot
#for width in [18., 12., 10., 8.8]:
for width in [10.]:
    fig = plt.figure(figsize=(cm2inch(width), cm2inch(width*6/8.)))
    # this should be changed for making a panel of multiple figures
    ax = fig.add_subplot(111)
    leg = ax.legend()
    col='#008ae6'
    z0 = pca[:,0] # lower limit
    z1 = pca[:,2] # upper limit
    w0 = pca[:,3] # mean w
    dw = pca[:,4]
    #
    for j in range(len(z0)-2):
        llc = (z0[j],w0[j]-2.*dw[j]) # lower left corner
        dz  = z1[j]-z0[j]
        dx  = 4.*dw[j]
        bbox1 = Rectangle(llc,dz,dx,transform=ax.transData,ec=col,fc=col,fill=True,alpha=0.8)
        ax.add_patch(bbox1)

    for j in range(len(z0)-2,len(z0)-1):
        llc = (z0[j],w0[j]-2.*dw[j]) # lower left corner
        dz  = z1[j]-z0[j]
        dx  = 4.*dw[j]
        bbox1 = Rectangle(llc,dz,dx,transform=ax.transData,ec=col,fc=col,fill=True,alpha=0.8, label='Planck + BSH')
        ax.add_patch(bbox1)


    # x axis  (plots a line hlines(z,xmin,xmax)
    #plt.hlines(0, 1.8, 7.2,color=(.5,.5,.5),label='_nolegend_')

    # legend
    #    leg = plt.legend(frameon=True)
    #leg = plt.legend(loc='upper right')
    #    leg = plt.legend(loc='lower right') # frameon keyword unknown in my version
    # remove box around legend
    #leg.get_frame().set_edgecolor("white")
    #leg.get_frame().set_alpha(.8)

    # labels
    plt.xlabel(r"$z$"); plt.ylabel(r"$w(z)$")
    ax.yaxis.labelpad = 10*width/17.; ax.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

    # reduce ticks for small figures
    if width < 10:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
    # grid
    #    plt.grid(True, which="major", axis="both")  # problem in my version

    # axes limits
    plt.xlim([0, 2.01]); plt.ylim([-1.8, 0.2]);

    # reduce white space around figure
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

    # set vertical y axis ticklables
    for ticklabel in ax.yaxis.get_ticklabels():
        ticklabel.set_rotation("vertical")

    l = plt.legend(loc=1,prop={'size':10}, frameon = False)
    colors = ["#008AE6", "#E03424"]
    for color,text in zip(colors,l.get_texts()):
        text.set_color(color)

#    g.export(os.path.join(outdir,'wPCA.pdf'))

    plt.tight_layout()
    # save to pdf with right bounding box
    plt.savefig("/home/pettorin/codici/cosmomc_plots/output_getdist/clik/clik10/pca/plots/wPCA.pdf", bbox_inches='tight', pad_inches=0.05)
