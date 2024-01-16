import numpy as np
from numpy import ma
from matplotlib import use
use('pdf')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as mpl
import sys



rdeg = 5.
zmin=-20.
zmax=20.
xmin = -rdeg
xmax = rdeg
ymin = -rdeg
ymax = rdeg


def place_map(myaxes, i, j, fname, title=''):
    """put a healpix map fname at location i, j"""    
    grid = np.loadtxt(fname).transpose()
    vg = grid[::-1]
    mymap = myaxes[i, j].imshow(ma.masked_where(vg < -190., vg), extent=[xmin,xmax,ymin, ymax], vmin = zmin, vmax = zmax)
    mymap.axes.get_xaxis().set_visible(False)
    mymap.axes.get_yaxis().set_visible(False)
    mymap.axes.set_frame_on(False)
    mymap.axes.set_xmargin(0.)
    mymap.axes.set_ymargin(0.)
    if(title != ''):
        mymap.axes.set_title(title)    
    return


def make_plot(weight):
    plt.axis('off')
    fig, myaxes = plt.subplots(nrows=3, ncols=3, figsize=(10, 7.2))
    my_cmap = mpl.cm.get_cmap('jet')
    my_cmap.set_bad('0.5')
    mpl.cm.register_cmap(cmap = my_cmap)
    plt.set_cmap(my_cmap)
    place_map(myaxes, 0,0,'t2zeta_' + weight + '_mean.txt', '$\zeta | T$:     mean   ')
    place_map(myaxes, 0,1,'t2zeta_' + weight + '_fluc1.txt', ' realization #1 ')
    place_map(myaxes, 0,2,'t2zeta_' + weight + '_fluc2.txt', ' realization #2 ')
    place_map(myaxes, 1,0,'e2zeta_' + weight + '_mean.txt', '$\zeta | E$:                 ')
    place_map(myaxes, 1,1,'e2zeta_' + weight + '_fluc1.txt')
    place_map(myaxes, 1,2,'e2zeta_' + weight + '_fluc2.txt')
    place_map(myaxes, 2,0,'te2zeta_' + weight + '_mean.txt','$\zeta | T, E$:                  ')
    place_map(myaxes, 2,1,'te2zeta_' + weight + '_fluc1.txt')
    place_map(myaxes, 2,2,'te2zeta_' + weight + '_fluc2.txt')
    ax1 = fig.add_axes([0.25, 0.065, 0.5, 0.012])
    cmap = mpl.cm.jet
    norm = mpl.colors.Normalize(vmin = zmin, vmax = zmax)
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
    cb1.set_label("$10^5\zeta$ (patch mean removed)")
    plt.savefig('zeta_maps_' + weight + '_highres.pdf', format='pdf')


make_plot("recomb_slice")

make_plot("vis")

