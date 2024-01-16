import numpy as np
from matplotlib import use
use('pdf')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as mpl

plt.axis('off')
fig, myaxes = plt.subplots(nrows=3, ncols=3, figsize=(10, 8.2))

def place_map(i, j, fname, title=''):
    """put a healpix map fname at location i, j"""    
    img = mpimg.imread(fname)
    mymap = myaxes[i, j].imshow(img)
    mymap.axes.get_xaxis().set_visible(False)
    mymap.axes.get_yaxis().set_visible(False)
    mymap.axes.set_frame_on(False)
    mymap.axes.set_xmargin(0.)
    mymap.axes.set_ymargin(0.)
    if(title != ''):
        mymap.axes.set_title(title)    
    return

place_map(0,0,'../zeta0256/t2zeta_vis_mean.png', '$\zeta | T$:     mean   ')
place_map(0,1,'../zeta0256/t2zeta_vis_mean_plus_fluc1.png', ' realization #1 ')
place_map(0,2,'../zeta0256/t2zeta_vis_mean_plus_fluc2.png', ' realization #2 ')
place_map(1,0,'../zeta0256/e2zeta_vis_mean.png', '$\zeta | E$:                 ')
place_map(1,1,'../zeta0256/e2zeta_vis_mean_plus_fluc1.png')
place_map(1,2,'../zeta0256/e2zeta_vis_mean_plus_fluc2.png')
place_map(2,0,'../zeta0256/te2zeta_vis_mean.png','$\zeta | T, E$:                  ')
place_map(2,1,'../zeta0256/te2zeta_vis_mean_plus_fluc1.png')
place_map(2,2,'../zeta0256/te2zeta_vis_mean_plus_fluc2.png')


ax1 = fig.add_axes([0.25, 0.065, 0.5, 0.012])
cmap = mpl.cm.jet
norm = mpl.colors.Normalize(vmin = -40., vmax = 40.)
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                   norm=norm,
                                   orientation='horizontal')
cb1.set_label("$10^5\zeta$")
fig.tight_layout()
plt.savefig('zeta_maps_vis.pdf', format='pdf')
