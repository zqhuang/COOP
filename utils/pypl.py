#!/usr/bin/env python
import math
import numpy as np
import sys
from matplotlib import use
use('pdf')
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import colorConverter
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import scipy.interpolate as minterp
import pylab

cdict = {
    'red'  :  ((0., 1., 1.), (0.03, 0.56, 0.56), (0.2, 0., 0.), (0.3, 0., 0.), (0.4, 0.1, 0.1), (0.5, 0., 0.) , (0.6, 0.7, 0.7), (0.7, 1., 1.), (0.8, 1., 1.), (1., 1., 1.)),
    'green':  ((0., 1., 1.), (0.03, 0.21, 0.21), (0.2, 0., 0.), (0.3, 1., 1.), (0.4, 0.9, 0.9) ,(0.5, 1., 1.),  (0.7, 1., 1.),  (0.8, 0.5, 0.5), (1., 0., 0.)),
    'blue' :  ((0., 1., 1.), (0.03, 0.94, 0.94), (0.2, 1., 1.), (0.3, 1., 1.), (0.4, 0.1, 0.1), (0.5, 0., 0.),  (0.7, 0., 0.), (0.8, 0.2, 0.2),  (1., 0., 0.))
    }
#generate the colormap with 1024 interpolated values
my_cmap = pylab.matplotlib.colors.LinearSegmentedColormap('my_cmap', cdict, 1024)
mpl.cm.register_cmap(cmap = my_cmap)

global_im = None
global global_cmap 
global_cmap = ""
global global_zmin
global_zmin = 1.e31
global global_zmax 
global_zmax = -1.e31

if(len(sys.argv)<3):
    print "pypl.py input_file output_file"
    sys.exit()

def insert_items(dic, line):
    items = line.split('=')
    if(len(items) == 2):
        if(items[0].strip != '' and items[1].strip() != '' ):
            dic[items[0].strip()] = items[1].strip()

def load_dictionary(dic, filename):
    fp = open(filename, 'r')
    for line in fp:
        insert_items(dic, line)
    fp.close()

##load the file to a dictionary
global settings
settings = {}
load_dictionary(settings, sys.argv[1])


def str_of(key, default = "", dic = settings):
    try:
        x = dic[key]
        return x
    except:
        return default

def int_of(key, default = 0, dic = settings):
    try:
        x = int(dic[key])
        return x
    except:
        return default

def float_of(key, default = 0., dic=settings):
    try:
        x = float(dic[key])
        return x
    except:
        return default

def logical_of(key, default=False, dic = settings):
    try:
        x = ( dic[key] in ['T', 't', '1', 'Y', 'y', 'yes', 'YES', 'Yes', 'True', 'TRUE'] )
        return x
    except:
        return False


def int_arr_of(key, default = [], dic=settings):
    try:
        x = [ int(s) for s in dic[key].split() ]
        return x
    except:
        return default


def float_arr_of(key, default = [], dic=settings):
    try:
        x = [ float(s) for s in dic[key].split() ]
        return x
    except:
        return default


nrows = int_of('nrows',1)
ncols = int_of('ncols',1)


def read_int(fp, default = 0):
    s = fp.readline()
    try:
        x = int(s)
        return x
    except:
        return default


def read_float(fp, default = 0.):
    s = fp.readline()
    try:
        x = float(s)
        return x
    except:
        return default

def read_int_arr(fp, default = []):
    s = fp.readline().split()
    try:
        x = [ int(u) for u in s ]
        return x
    except:
        return default

def read_float_arr(fp, default = []):
    s = fp.readline().split()
    try:
        x = [ float(u) for u in s ]
        return x
    except:
        return default


def read_logical(fp, default = False):
    s = fp.readline().split()
    try:
        x = (s in ['T', 't', '1', 'Y', 'y', 'yes', 'YES', 'Yes', 'True', 'TRUE'])  
        return x
    except:
        return default

def read_logical_arr(fp, default = []):
    s = fp.readline().split()
    try:
        x = [ (u in ['T', 't', '1', 'Y', 'y', 'yes', 'YES', 'Yes', 'True', 'TRUE'])  for u in s ]
        return x
    except:
        return default

def read_str(fp, default=""):
    s = fp.readline().strip()
    if(s=="NULL" or s == ""):
        return default
    else:
        return s

invisible= ( 1., 1., 1., 0. )
    

def line_config(lc):
    c = lc.split('_')
    cstr = c[0].split(":")
    genre = cstr[0].strip().lower()
    if(len(cstr) == 5 and genre=='rgba'):
        try:
            rval = float(cstr[1])/255.
            gval = float(cstr[2])/255.
            bval = float(cstr[3])/255.
            alpha = float(cstr[4])/255.
        except:
            rval = 0.
            gval = 0.
            bval = 0.
            alpha = 0.
        color = ( rval, gval, bval, alpha )
    elif(len(cstr) == 4 and genre=='rgb'):
        try:
            rval = float(cstr[1])/255.
            gval = float(cstr[2])/255.
            bval = float(cstr[3])/255.
        except:
            rval = 0.
            gval = 0.
            bval = 0.
        color = ( rval, gval, bval )
    elif(len(cstr) == 2 and genre=="hex"):
        color=colorConverter.to_rgb("#" + cstr[1])
    elif(len(cstr) == 2 and (genre=="gray" or genre=="grey" )):
        color = colorConverter.to_rgb( cstr[1] )
    else:
        if(genre in ['b','blue','g','green','r','red','c','cyan','m','magenta','yellow','y', 'k','black','w','white']):
            color=colorConverter.to_rgb(genre)
        elif(genre == 'invisible'):
            color= invisible
        elif(genre in ['violet', 'v']):
            color=(0.55, 0.22, 0.79)
        elif(genre == 'gold'):
            color= (0.83, 0.63, 0.09)
        elif(genre=='skyblue'):
            color = (0.24, 0.6, 1.)
        elif(genre in [ 'orange', 'o']):
            color= (1., 0.647, 0.)
        elif(genre  == 'slateblue'):
            color= (0.21, 0.45, 0.78)
        elif(genre == 'gray','grey'):
            color=(0.5, 0.5, 0.5)
        elif(genre=="maroon"):
            color=(0.51, 0.02, 0.25)
        elif(genre=='purple', 'p'):
            color = colorConverter.to_rgb("#a020f0")
        elif(genre=='brown'):
            color = colorConverter.to_rgb("#a52a2a")
        elif(genre=='pink'):
            color = colorConverter.to_rgb("#ffc0cb")
        elif(genre=='olive'):
            color = colorConverter.to_rgb("#556b2f")
        elif(genre=='turquoise'):
            color = colorConverter.to_rgb("#40e0d0")
        elif(genre=='lightgray', 'lightgrey'):
            color=(0.75, 0.75, 0.75)
        elif(genre=='darkgray', 'darkgrey'):
            color=(0.25, 0.25, 0.25)
        elif(genre=='lightred'):
            color = (1., 0.1, 0.1)
        elif(genre=='lightblue'):
            color = (0.1, 0.1, 1.)
        elif(genre=='lightgreen'):
            color=(0.1,1.,0.1)
        elif(genre=='royalblue'):
            color=(0.255, 0.412, 1.)
        elif(genre in ['lawngreen', 'grassgreen']):
            color = (0.53, 0.97, 0.10)
        elif(genre=='darkgreen'):
            color = colorConverter.to_rgb("#056405")
        elif(genre=='darkblue'):
            color = colorConverter.to_rgb("#050564")
        elif(genre=='darkred'):
            color = colorConverter.to_rgb("#640505")
        elif(genre=='darkbrown'):
            color = colorConverter.to_rgb("#8b4513")
        elif(genre=='darkcyan'):
            color = colorConverter.to_rgb("#00ced1")
        elif(genre=='darkmagenta'):
            color = colorConverter.to_rgb("#9400d3")
        elif(genre=='springgreen'):
            color = colorConverter.to_rgb("#00ff7f")
        else:
            color = invisible
    if(len(c) > 1):
        linetype = c[1]
    else:
        linetype = 'solid'
    if(len(c)>2):
        try:
            linewidth=float(c[2])
        except:
            linewidth = 1.
    else:
        linewidth = 1. 
    if(len(c)>3 and len(color)<4):
        try:
            alpha = float(c[3])
            color.append(alpha)
        except:
            print "Warning: error in alpha transparency parameter. Ignored."
    return ( color, linetype, linewidth )

def get_line_config(fp):
    lc = read_str(fp)
    return line_config(lc)


def plot_curve(ax, fp):
    n = read_int(fp)
    if( n<=0 or n >=100000):
        print "too many points in the curve"
        sys.exit()
    legend = read_str(fp)
    (color, linetype, linewidth) = get_line_config(fp)
    marker = read_str(fp)
    x = []
    y = []
    for i in range(n):
        xy = read_float_arr(fp)
        x.append(xy[0])
        y.append(xy[1])
    if(marker!=''):
        if(legend != ''):
            ax.plot(x, y, color=color, linestyle = linetype, linewidth = linewidth, label=legend, marker= marker)
        else:
            ax.plot(x, y, color=color, linestyle = linetype, linewidth = linewidth, marker= marker)
    else:
        if(legend != ''):
            ax.plot(x, y, color=color, linestyle = linetype, linewidth = linewidth, label=legend)
        else:
            ax.plot(x, y, color=color, linestyle = linetype, linewidth = linewidth)

def plot_legend_advance(ax, fp):
    (color, linetype, linewidth) = get_line_config(fp)
    xmargin = read_float(fp)
    ymargin = read_float(fp)
    linelength = read_float(fp)
    hskip = read_float(fp)
    vskip = read_float(fp)
    cols = read_int(fp)
    cstr = read_str(fp).upper()
    if(cstr == ''):
        loc=read_float_arr(fp)  #ignroe, only treated in Asymptote
        ax.legend(ncol = cols)
    else:
        if(cstr == 'N' or cstr == '9'):
            loc = 9
        elif(cstr == 'E' or cstr == '7'):
            loc = 4
        elif(cstr == 'S' or cstr == '8'):
            loc = 8
        elif(cstr == 'W' or cstr == '6'):
            loc = 6
        elif(cstr == 'NE' or cstr == '1'):
            loc = 1
        elif(cstr == 'NW' or cstr == '2'):
            loc = 2
        elif(cstr == 'SE' or cstr == '4'):
            loc = 4
        elif(cstr == 'SW' or cstr == '3'):
            loc = 3
        elif(cstr == 'C' or cstr == '10'):
            loc = 10
        else:
            loc = 0
        ax.legend(ncol = cols, loc = loc)

def plot_legend(ax, fp):
    cstr = read_str(fp)
    if( cstr == ''):
        loc = read_float_arr(fp)  #ignore, only treated in Asymptote
        cols = read_int(fp)
        ax.legend(ncol = cols)
    elif(cstr == 'VIRTUAL'):
        legend = read_str(fp)
        (color, linetype, linewidth) = get_line_config(fp)
        if(linewidth >= 4.):
            mypatch = mpatches.Patch(color=color, label = legend)
        else:
            mypatch = mlines.Line2D([], [], color=color, ls = linetype, lw=linewidth, label=legend)
        ax.legend(handles=[mypatch])            
    else:
        cols = read_int(fp)
        if(cstr == 'N' or cstr == '9'):
            loc = 9
        elif(cstr == 'E' or cstr == '7'):
            loc = 4
        elif(cstr == 'S' or cstr == '8'):
            loc = 8
        elif(cstr == 'W' or cstr == '6'):
            loc = 6
        elif(cstr == 'NE' or cstr == '1'):
            loc = 1
        elif(cstr == 'NW' or cstr == '2'):
            loc = 2
        elif(cstr == 'SE' or cstr == '4'):
            loc = 4
        elif(cstr == 'SW' or cstr == '3'):
            loc = 3
        elif(cstr == 'C' or cstr == '10'):
            loc = 10
        else:
            loc = 0
        ax.legend(ncol = cols, loc = loc)

def plot_labels(ax, fp):
    n = read_int(fp)
    if(n <= 0 or n>10000):
        print "too many labels"
        sys.exit()
    (color, linetype, linewidth) = get_line_config(fp)
    for i in range(n):
        xy = read_float_arr(fp)
        label = read_str(fp)
        ax.text(xy[0],xy[1], label, color=color, ha='center')

def plot_rightlabels(ax, fp):
    n = read_int(fp)
    if(n <= 0 or n>10000):
        print "too many labels"
        sys.exit()
    (color, linetype, linewidth) = get_line_config(fp)
    for i in range(n):
        xy = read_float_arr(fp)
        label = read_str(fp)
        ax.text(xy[0], xy[1], label, color=color, ha='left')

def plot_leftlabels(ax, fp):
    n = read_int(fp)
    if(n <= 0 or n>10000):
        print "too many labels"
        sys.exit()
    (color, linetype, linewidth) = get_line_config(fp)
    for i in range(n):
        xy = read_float_arr(fp)
        label = read_str(fp)
        ax.text(xy[0],xy[1], label, color=color, ha='right')

def plot_contour(ax, fp):
    ctype=read_int(fp)
    if(ctype == 1):
        (cf, lsf, lwf) = get_line_config(fp)
        (cb, lsb, lwb) = get_line_config(fp)
        legend = read_str(fp)
        smooth = read_logical(fp)
        nc = read_int(fp)
        for ip in range(nc):
            n = read_int(fp)
            x = []
            y = []
            for i in range(n):
                xy = read_float_arr(fp)
                x.append(xy[0])
                y.append(xy[1])
            x.append( x[0] )
            y.append( y[0] )
            if(ip == 0 and legend != '' and legend != 'NULL'):
                ax.plot(x, y, color=cb, ls = lsb, lw = lwb, label=legend) 
            else:
                ax.plot(x, y, color=cb, ls = lsb, lw = lwb)                
            ax.fill(x, y, color=cf)
    else:
        xr = read_float_arr(fp)
        yr = read_float_arr(fp)
        nc = read_int(fp)
        cvals = read_float_arr(fp)
        fcolor = []
        bcolor = []
        fls = []
        bls = []
        flw = []
        blw = []
        falpha = []
        balpha = []
        for i in range(nc):
            (color, ls, lw) = get_line_config(fp)
            if(len(color)>3):
                falpha.append(color[3])
                fcolor.append([color[0], color[1], color[2]])
            else:
                falpha.append(1.)
                fcolor.append(color)
            fls.append(ls)
            flw.append(lw)
            (color, ls, lw) = get_line_config(fp)
            bcolor.append(color)
            bls.append(ls)
            blw.append(lw)
        smooth = read_logical(fp)
        n = read_int_arr(fp)
        grid = []
        for i in range(n[0]):
            x = read_float_arr(fp)
            grid.append(x)
        grid = np.array(grid).transpose()
        grid = grid[::-1]
        mcvals = cvals
        mcvals.append(grid.max())
        ax.contourf(grid, colors = fcolor, linstyles = fls, extent = [xr[0], xr[1], yr[0], yr[1]], levels = cvals, alphas = falpha)
        ax.contour(grid, colors = bcolor, linewidths = blw, linestypes = bls, extent = [xr[0], xr[1], yr[0], yr[1]], levels = cvals)

def plot_lines(ax, fp):
    n = read_int(fp)
    if(n<=0 or n>1000000):
        print 'too many lines'
        sys.exit()
    (color, ls, lw) = get_line_config(fp)
    for i in range(n):
        xyxy=read_float_arr(fp)
        ax.plot( [xyxy[0], xyxy[2]], [xyxy[1], xyxy[3]], color=color, ls = ls, lw = lw)
    return

def plot_dots(ax, fp):
    n = read_int(fp)
    if(n<0 or n > 100000):
        print "too many dots"
        sys.exit()
    (color, ls, lw) = get_line_config(fp)
    marker = read_str(fp)
    x = []
    y = []
    for i in range(n):
        xy = read_float_arr(fp)
        x.append(xy[0])
        y.append(xy[1])
    if(marker == 'dot' or marker == 'DOT'):
        ax.scatter(x, y, c=color)
    else:
        ax.scatter(x, y, c = color, marker=marker)
    return


def plot_density(ax, fp, global_cmap = global_cmap, global_zmin=global_zmin, global_zmax=global_zmax):
    ctbl = read_str(fp).lower()
    if(ctbl == 'my_cmap'):
        cmap = my_cmap
        plt.set_cmap(my_cmap) 
    else:
        cmap = plt.get_cmap(ctbl)
        plt.set_cmap(ctbl)
    zlabel = read_str(fp)
    xr =read_float_arr(fp)
    yr = read_float_arr(fp)
    zr = read_float_arr(fp)
    irr = read_int(fp)
    if(irr == 0 or irr == -1):  # regular points
        grid=[]
        n = read_int_arr(fp)
        for i in range(n[0]):
            x = read_float_arr(fp)
            grid.append(x)
        grid = np.array(grid).transpose()
        grid = grid[::-1]   ##to be the same as Asymptote
        if(abs(zr[0])<1.e30 and abs(zr[1])<1.e30):
            im = ax.imshow(grid, extent = [xr[0], xr[1], yr[0], yr[1]], vmin=zr[0], vmax = zr[1])
        elif(abs(zr[0])<1.e30):
            im = ax.imshow(grid, extent = [xr[0], xr[1], yr[0], yr[1]], vmin=zr[0])
        elif(abs(zr[1])<1.e30):
            im = ax.imshow(grid, extent = [xr[0], xr[1], yr[0], yr[1]], vmax=zr[1])
        else:
            im = ax.imshow(grid, extent = [xr[0], xr[1], yr[0], yr[1]])
        if(irr == 0):
            if(nrows == 1 and ncols==1):
                plt.colorbar(im)
            else:
                if(global_cmap!= ''):
                    if(ctbl != global_cmap):
                        print 'for single figure you cannot use different color bars!'
                        sys.exit()
                else:
                    global_cmap = ctbl
                if(abs(zr[0])<1.e30 and abs(zr[1])<1.e30):
                    if(abs(global_zmin) < 1.e30 and abs(global_zmax)<1.e30):
                        if(z[0]!= global_zmin or z[1] != global_zmax):
                            print 'for density plots sharing a colorbar you must use the same zmin and zmax'
                            sys.exit()
                    else:
                        global_zmin = zr[0]
                        global_zmax = zr[1]
                else:
                    print 'for density plots sharing a colorbar you must specify zmin and zmax'
                    sys.exit()
    elif(irr == 1 or irr == 2):
        n = read_int(fp)
        x = []
        y = []
        z = []
        for i in range(n):
            coor = read_float_arr(fp)
            x.append(coor[0])
            y.append(coor[1])
            z.append(coor[2])
        ns = int(np.sqrt(n+1.))
        if(abs(xr[0])<1.e30):
            xmin = xr[0]
        else:
            xmin = x.min()
        if(abs(xr[1])<1.e30):
            xmax = xr[1]
        else:
            xmax = x.max()
        if(abs(yr[0])<1.e30):
            ymin = yr[0]
        else:
            ymin = y.min()
        if(abs(yr[1])<1.e30):
            ymax = yr[1]
        else:
            ymax = y.max()
        dx = (xmax-xmin)/ns
        dy = (ymax-ymin)/ns
        xs, ys = np.mgrid[xmin:xmax:dx, ymin:ymax:dy]
        resample = mpl.mlab.griddata(x, y, z, xs, ys)
        if(abs(zr[0])<1.e30 and abs(zr[1])<1.e30):
            im = ax.imshow(resample.T, extent = [xmin, xmax, ymin, ymax], vmin=zr[0], vmax = zr[1])
        elif(abs(zr[0])<1.e30):
            im = ax.imshow(resample.T, extent = [xmin, xmax, ymin, ymax], vmin=zr[0])
        elif(abs(zr[1])<1.e30):
            im = ax.imshow(resample.T, extent = [xmin, xmax, ymin, ymax], vmax = zr[1])
        else:
            im = ax.imshow(resample.T, extent = [xmin, xmax, ymin, ymax])
        if(irr == 1):
            if(nrows == 1 and ncols==1):
                plot.colorbar(im)
            else:
                if(global_cmap!= ''):
                    if(ctbl != global_cmap):
                        print 'for single figure you cannot use different color bars!'
                        sys.exit()
                else:
                    global_cmap = ctbl
                if(abs(zr[0])<1.e30 and abs(zr[1])<1.e30):
                    if(abs(global_zmin) < 1.e30 and abs(global_zmax)<1.e30):
                        if(z[0]!= global_zmin or z[1] != global_zmax):
                            print 'for density plots sharing a colorbar you must use the same zmin and zmax'
                            sys.exit()
                    else:
                        global_zmin = zr[0]
                        global_zmax = zr[1]
                else:
                    print 'for density plots sharing a colorbar you must specify zmin and zmax'
                    sys.exit()
    return (im, global_cmap, global_zmin, global_zmax)

def plot_clip(ax, fp):
    print 'Ignored a clip block. Clipping in pyplot is automatic; use Asymptote for the full support'
    sys.exit()

def plot_errorbars(ax, fp):
    n = read_int(fp)
    if(n <=0 or n > 100000):
        print "too many errorbars to plot"
        sys.exit()
    (color, ls, lw) = get_line_config(fp)
    barsize = read_float(fp) #ignroed; use Asymptote for full support.
    center_symbol = read_str(fp)
    (ccolor, cls, clw) = get_line_config(fp)
    x = []
    y = []
    xm = []
    xp = []
    ym = []
    yp = []
    hasxerr = False
    hasyerr = False
    for i in range(n):
        pts = read_float_arr(fp)
        x.append(pts[0])
        y.append(pts[1])
        xm.append(pts[3])
        xp.append(pts[2])
        ym.append(pts[5])
        yp.append(pts[4])
        if(pts[2] > 0. or pts[3] > 0.):
            hasxerr = True
        if(pts[4] > 0. or pts[5] > 0.):
            hasyerr = True
    xerr = [xm, xp]
    yerr = [ym, yp]
    if(center_symbol == ''):
        if(hasxerr and hasyerr):
            ax.errorbar(x, y, ecolor = color, elinewidth = lw, xerr = xerr, yerr = yerr, fmt = None)
        elif(hasxerr):
            ax.errorbar(x, y, ecolor = color, elinewidth = lw, xerr = xerr, fmt = None)
        elif(hasyerr):
            ax.errorbar(x, y, ecolor = color, elinewidth = lw, yerr = yerr, fmt = None)
        else:
            print "Warning: all errorbars are zero?"
    elif(center_symbol == 'DOT') :
        if(hasxerr and hasyerr):
            ax.errorbar(x, y, ecolor = color, elinewidth = lw,  mfc = ccolor, xerr = xerr, yerr = yerr, fmt = 'o', fillstyle = 'full')
        elif(hasxerr):
            ax.errorbar(x, y, ecolor = color, elinewidth = lw,  mfc = ccolor, xerr = xerr, fmt = 'o', fillstyle = 'full')
        elif(hasyerr):
            ax.errorbar(x, y, ecolor = color, elinewidth = lw,  mfc = ccolor, yerr = yerr, fmt = 'o', fillstyle = 'full')
        else:
            print "Warning: all errorbars are zero?"        
    else:
        if(hasxerr and hasyerr):
            ax.errorbar(x, y, ecolor = color, elinewidth = lw,  mfc = ccolor, xerr = xerr, yerr = yerr, fmt = center_symbol )
        elif(hasxerr):
            ax.errorbar(x, y, ecolor = color, elinewidth = lw,  mfc = ccolor, xerr = xerr, fmt = center_symbol)
        elif(hasyerr):
            ax.errorbar(x, y, ecolor = color, elinewidth = lw,  mfc = ccolor, yerr = yerr, fmt = center_symbol)
        else:
            print "Warning: all errorbars are zero?"

    
def plot_arrows(ax, fp):
    n = read_int(fp)
    if(n <=0 or n > 100000):
        print "too many arrows"
        sys.exit()
    (c, ls, lw) = get_line_config(fp)
    for i in range(n):
        xyxy = read_float_arr(fp)
        ax.plot( x= (xyxy[0], xyxy[2]), y = (xyxy[1], xyxy[3]), color = c, ls= ls, lw = lw, clip_on = False)

def plot_annotate(ax, fp):
    (c, ls, lw) = get_line_config(fp)
    xyxy = read_float_arr(fp)
    label = read_str(fp)
    ax.annotate(s=label, xy = (xyxy[0], xyxy[1]), xytext = (xyxy[2], xyxy[3]), xycoords='data', textcoords='data', annotation_clip = False, arrowprops = dict(color=c, lw = lw, ls=ls))


def plot_extra_axis(ax, fp):
    loc = read_str(fp).lower()
    label = read_str(fp)
    islog = read_logical(fp)
    minmax = read_float_arr(fp)
    if(loc == 'right'):
        extra_axis = ax.twiny()
        if(islog):
            extra_axis.set_yscale("log", nonposy = clip)
        extra_axis.set_ylim(ymin = minmax[0], ymax = minmax[1])
    elif(loc == "top"):
        extra_axis = ax.twinx()
        if(islog):
            extra_axis.set_xscale("log", nonposy = clip)            
        extra_axis.set_xlim(xmin = minmax[0], xmax = minmax[1])

def plot_expand(ax, fp):
    expfac = read_float_arr(fp)
    xmin = ax.xmin - (ax.xmax-ax.xmin)*expfac[0]
    xmax = ax.xmax + (ax.xmax-ax.xmin)*expfac[1]
    ymin = ax.ymin - (ax.ymax-ax.ymin)*expfac[2]
    ymax = ax.ymax + (ax.ymax-ax.ymin)*expfac[3]
    ax.set_xlim(xmin = xmin, xmax = xmax)
    ax.set_ylim(ymin = ymin, ymax = ymax)


def loadfig(ax, filename, global_im = global_im, global_cmap = global_cmap, global_zmin = global_zmin, global_zmax = global_zmax) :
    if(filename == ''):
        return (global_im, global_cmap, global_zmin, global_zmax)
    try:
        fp = open(filename, 'r')
    except:
        print filename + 'does not exist'
        sys.exit()
    fp.readline() # skip the xsize, ysize parameters for single figure
    caption = read_str(fp)
    xlabel = read_str(fp)
    ylabel = read_str(fp)
    islog = read_logical_arr(fp, [False, False, False])
    clip = read_logical(fp, False)
    xbounds = read_float_arr(fp, [0., 1.])
    ybounds = read_float_arr(fp, [0., 1.])
    nblocks = read_int(fp, 0)
    if(xlabel == ''):
        plt.setp( [ ax.get_xticklabels() ], visible = False)
    if(ylabel == ''):
        plt.setp( [ ax.get_yticklabels() ], visible = False) 
    if(caption !=''):
        ax.set_title(caption)
    if(xlabel !=''):
        ax.set_xlabel(xlabel)
    if(ylabel !=''):
        ax.set_ylabel(ylabel)
        for tick in ax.get_yticklabels():
            tick.set_rotation(90)
    if(islog[0]):
        ax.set_xscale("log", nonposx='clip')
    if(islog[1]):
        ax.set_yscale("log", nonposy = 'clip')
    if(abs(xbounds[0]) < 1.e30 and abs(xbounds[1]) < 1.e30):
        ax.set_xlim(xmin = xbounds[0], xmax = xbounds[1])
    elif(abs(xbounds[0]) < 1.e30):
        ax.set_xlim(xmin = xbounds[0])
    elif(abs(xbounds[1]) < 1.e30):
        ax.set_xlim(xmax = xbounds[1])
    if(abs(ybounds[0]) < 1.e30 and abs(ybounds[1]) < 1.e30):
        ax.set_ylim(ymin = ybounds[0], ymax = ybounds[1])
    elif(abs(ybounds[0]) < 1.e30):
        ax.set_ylim(ymin = xbounds[0])
    elif(abs(ybounds[1]) < 1.e30):
        ax.set_ylim(ymax = xbounds[1])
    iblocks = 0
    genre = read_str(fp)
    while(genre !=''):
        if(genre == "CURVE"):
            plot_curve(ax, fp)
        elif(genre == "LEGEND_ADVANCE"):
            plot_legend_advance(ax, fp)
        elif(genre == "LEGEND" or genre=="LEGEND_NOBOX"):
            plot_legend(ax, fp)
        elif(genre=="LABELS"):
            plot_labels(ax, fp)
        elif(genre=="RIGHTLABELS"):
            plot_rightlabels(ax, fp)
        elif(genre=="LEFTLABELS"):
            plot_leftlabels(ax, fp)
        elif(genre=="LINES"):
            plot_lines(ax, fp)
        elif(genre=="DOTS"):
            plot_dots(ax, fp)
        elif(genre=="DENSITY"):
            (global_im, global_cmap, global_zmin, global_zmax) = plot_density(ax, fp)
        elif(genre == 'CONTOUR'):
            plot_contour(ax, fp)
        elif(genre=="CLIP"):
            plot_clip(ax, fp)
        elif(genre=="EXPAND"):
            plot_expand(ax, fp)
        elif(genre=="ERRORBARS"):
            plot_errorbars(ax, fp)
        elif(genre=="ANNOTATE"):
            plot_annotate(ax, fp)
        elif(genre=="ARROWS"):
            plot_arrows(ax, fp)
        elif(genre=="EXTRA_AXIS"):
            plot_extra_axis(ax, fp)
        else:
            print genre + ": unknown block name"
            return (global_im, global_cmap, global_zmin, global_zmax)
        iblocks += 1
        if(iblocks == nblocks):
            fp.close()
            return (global_im, global_cmap, global_zmin, global_zmax)
        genre = read_str(fp)        
    fp.close()
    return (global_im, global_cmap, global_zmin, global_zmax)

hspace = float_of('vertical_space', 0.12)
wspace = float_of('horizontal_space', 0.12)
left_space = float_of('left_space', 0.12)
right_space = float_of('right_space', 0.08)
bottom_space = float_of('bottom_space', 0.12)
top_space = float_of('top_space', 0.08)
minhspace = 0.03 #if less than this, remove labels  
minwspace = 0.05   
if(nrows == 1 and ncols==1):
    fp=open(sys.argv[1], 'r')
    sizes = read_float_arr(fp)
    fp.close()
    fig, ax = plt.subplots(figsize=(sizes[0],sizes[1]))
    loadfig(ax, sys.argv[1])
else:
    fig, axarr = plt.subplots( nrows=nrows, ncols=ncols, figsize=( float_of('xsize', 8.),  float_of('ysize', 6.) ), sharex = logical_of('sharex', True), sharey = logical_of('sharey', True))
    fig.subplots_adjust( left = left_space, right = 1.- right_space, bottom = bottom_space, top = 1.-top_space, wspace=wspace, hspace = hspace)
    if(nrows > 1 and ncols > 1):
        for irow in range(nrows):
            for icol in range(ncols):
                if(hspace < minhspace and irow < nrows-1):
                    plt.setp( axarr[irow][icol].get_xticklabels(), visible=False)
                if(wspace<minwspace and icol>0):
                    plt.setp( axarr[irow][icol].get_yticklabels(), visible=False)
                filename = settings['figure['+str(irow)+','+str(icol)+']'].strip() 
                (global_im, global_cmap, global_zmin, global_zmax) = loadfig(axarr[irow,icol], filename)
    elif(nrows > 1):
        for irow in range(nrows):
            if(hspace < minhspace and irow<nrows-1):
                plt.setp( axarr[irow][0].get_xticklabels(), visible=False)
            filename = settings['figure['+str(irow)+',0]'].strip() 
            (global_im, global_cmap, global_zmin, global_zmax) = loadfig(axarr[irow], filename)
    elif(ncols>1):
        for icol in range(ncols):
            if(wspace<minwspace and icol>0):
                plt.setp( axarr[0][icol].get_yticklabels(), visible=False)
            filename = settings['figure[0,'+str(icol)+']'].strip() 
            (global_im, global_cmap, global_zmin, global_zmax) =loadfig(axarr[icol], filename)
    else:
        filename = settings['figure[0,0]'].strip() 
        loadfig(axarr, filename)


if(global_cmap != ''):
    left = 1. - right_space + 0.03
    bottom = bottom_space+hspace/2
    top = 1.-top_space-hspace/2
    cax = fig.add_axes( [ left,  bottom, min(1.-left,0.03), top-bottom])
    fig.colorbar(mappable=global_im, cax = cax, label = str_of('zlabel'))
plt.savefig(sys.argv[2], format='pdf')
