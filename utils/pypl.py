#! /user/bin/python
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
settings = {}
load_dictionary(settings, sys.argv[1])


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



invisible= ( 1., 1., 1., 1. )
    

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
        elif(genre=='purple'):
            color = colorConverter.to_rgb("#a020f0")
        elif(genre=='brown'):
            color = colorConverter.to_rgb("#a52a2a")
        elif(genre=='pink'):
            color = colorConverter.to_rgb("#ffc0cb")
        elif(genre=='olive'):
            color = colorConverter.to_rgb("#556b2f")
        elif(genre=='turquoise'):
            color = colorConverter.to_rgb("#40e0d0")
        elif(genre=='lightgray'):
            color=(0.75, 0.75, 0.75)
        elif(genre=='darkgray'):
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
    cstr = read_str(fp)
    if(cstr == ''):
        loc=read_float_arr(fp)
    ax.legend()

def plot_legend(ax, fp):
    lc = read_str(fp)
    if( lc == ''):
        loc = read_float_arr(fp)
        cols = read_int(fp)
        ax.legend()
    elif(lc=='VIRTUAL'):
        legend = read_str(fp)
        (color, linetype, linewidth) = get_line_config(fp)
        if(linewidth >= 4.):
            mypatch = mpatches.Patch(color=color, label = legend)
        else:
            mypatch = mlines.Line2D([], [], color=color, ls = linetype, lw=linewidth, label=legend)
        ax.legend(handles=[mypatch])            
    else:
        cols = read_int(fp)
        ax.legend()

def plot_labels(ax, fp):
    n = read_int(fp)
    if(n <= 0 or n>10000):
        print "too many labels"
        sys.exit()
    (color, linetype, linewidth) = get_line_config(fp)
    for i in range(n):
        xy = read_float_arr(fp)
        label = read_str(fp)
        ax.text(xy[0],xy[1], label, color=color)

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
        smooth = read_logical(fp)
        np = read_int(fp)
        for ip in range(np):
            n = read_int(fp)
            x = []
            y = []
            for i in range(n):
                xy = read_float_arr(fp)
                x.append(xy[0])
                y.append(xy[1])
            x.append( x[0] )
            y.append( y[0] )
            ax.plot(x, y, color=cb, ls = lsb, lw = lwb)
            ax.fill(x, y, color=cf)
    else:
        print "type II contour is not supported yet; use Asymptote for the full support"
        sys.exit()
            
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


def plot_density(ax, fp):
    ctbl = read_str(fp)
    zlabel = read_str(fp)
    xr =read_float_arr(fp)
    yr = read_float_arr(fp)
    zr = read_float_arr(fp)
    irr = read_int(fp)
    grid=[]
    if(irr == 0):  # regular points
        n = read_int_arr(fp)
        for i in range(n[0]):
            x = read_float_arr(fp)
            grid.append(x)
        if(abs(zr[0])<1.e30 and abs(zr[1])<1.e30):
            im = ax.imshow(grid, extent = [xr[0], xr[1], yr[0], yr[1]], vmin=zr[0], vmax = zr[1])
        elif(abs(zr[0])<1.e30):
            im = ax.imshow(grid, extent = [xr[0], xr[1], yr[0], yr[1]], vmin=zr[0])
        elif(abs(zr[1])<1.e30):
            im = ax.imshow(grid, extent = [xr[0], xr[1], yr[0], yr[1]], vmax=zr[1])
        else:
            im = ax.imshow(grid, extent = [xr[0], xr[1], yr[0], yr[1]])
        #colorbar is ignored in python multiple panels; use Asymptot for the full support
        return
    else:
        print "irregular density plot is not supported yet; use asymptote for the full support"
        sys.exit()
    return

def plot_clip(ax, fp):
    print 'clip not supported yet; use Asymptote for the full support'
    sys.exit()

def plot_errorbars(ax, fp):
    print 'errorbars not supported yet; use Asymptote for the full support'
    sys.exit()

def plot_arrows(ax, fp):
    n = read_int(fp)
    if(n <=0 or n > 100000):
        print "too many arrows"
        sys.exit()
    (c, ls, lw) = get_line_config(fp)
    for i in range(n):
        xyxy = read_float_arr(fp)
        ax.arrow( xyxy[0], xyxy[1], xyxy[2], xyxy[3], head_width=0.05, head_length = 0.1, fc=c, ec = c)

def plot_extra_axis(ax, fp):
    print 'extra_axis not supported yet; use Asymptote for the full support'
    sys.exit()

def plot_expand(ax, fp):
    expfac = read_float_arr(fp)
    #in python i just ignore this; for full support go for asymptote

def loadfig(ax, filename):
    if(filename == ''):
        return False
    try:
        fp = open(filename, 'r')
    except:
        return False
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
            plot_density(ax, fp)
        elif(genre == 'CONTOUR'):
            plot_contour(ax, fp)
        elif(genre=="CLIP"):
            plot_clip(ax, fp)
        elif(genre=="EXPAND"):
            plot_expand(ax, fp)
        elif(genre=="ERRPRBARS"):
            plot_errorbars(ax, fp)
        elif(genre=="ARROWS"):
            plot_arrows(ax, fp)
        elif(genre=="EXTRA_AXIS"):
            plot_extra_axis(ax, fp)
        else:
            print genre + ": unknown block name"
            return
        iblocks += 1
        if(iblocks == nblocks):
            fp.close()
            return
        genre = read_str(fp)        
    fp.close()
    return True
    
nrows = int_of('nrows',0)
ncols = int_of('ncols',1)

if(nrows == 0):
    fp=open(sys.argv[1], 'r')
    sizes = read_float_arr(fp)
    fp.close()
    fig, ax = plt.subplots(figsize=(sizes[0],sizes[1]))
    loadfig(ax, sys.argv[1])
else:
    fig, axarr = plt.subplots( nrows=nrows, ncols=ncols, figsize=( float_of('xsize', 8.),  float_of('ysize', 6.) ), sharex = logical_of('sharex', True), sharey = logical_of('sharey', True))
    fig.subplots_adjust( left=float_of('left_space', 0.125), right = float_of('right_space', 0.9), bottom = float_of('bottom_space', 0.1), top =float_of('top_space', 0.9), wspace=float_of('horizontal_space', 0.2), hspace = float_of('vertical_space', 0.2))
    if(nrows > 1 and ncols > 1):
        for irow in range(nrows):
            for icol in range(ncols):
                filename = settings['figure['+str(irow)+','+str(icol)+']'].strip() 
                try:
                    loadfig(axarr[irow,icol], filename)
                except:
                    print filename + " cannot be plotted."
    elif(nrows > 1):
        for irow in range(nrows):
            filename = settings['figure['+str(irow)+',0]'].strip() 
            loadfig(axarr[irow], filename)
    elif(ncols>1):
        for icol in range(ncols):
            filename = settings['figure[1,'+str(icol)+']'].strip() 
            try:
                loadfig(axarr[icol], filename)
            except:
                print filename + ' cannot be plooted'
    else:
        filename = settings['figure[0,0]'].strip() 
        try:
            loadfig(axarr, filename)
        except:
            print filename + ' cannot be plooted'


plt.savefig(sys.argv[2], format='pdf')
