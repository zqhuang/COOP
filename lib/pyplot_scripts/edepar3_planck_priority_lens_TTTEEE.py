import planckStyle as s
from pylab import *
import numpy as np
import GetDistPlots, os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc, font_manager
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, \
grid, savefig, show


outdir='/home/pettorin/codici/cosmomc_plots/output_getdist/clik/clik10/ede/edepar3/figures/'
g = s.getSinglePlotter(plot_data='/home/pettorin/codici/cosmomc_plots/output_getdist/clik/clik10/ede/edepar3/plot_data/')
g.settings.setWithSubplotSize(4.0000)
n_groups = 4

# Values for TT
redshift = (1./0.1, 1./0.02, 1./0.005, 1./0.001 )
width = (1,1,1,1)
limit2 = (0.0360, 0.0198, 0.0150, 0.0118) #prior on a1dn1

bar_width = (0.6,3.1,11,52)

# Values for WL + RSD
redshift_RSD_WL = (1./0.1, 1./0.02, 1./0.005, 1./0.001 )
redshift_RSD_WL = map(sum,zip(redshift_RSD_WL,bar_width))   #shift just for plotting
limit2_RSD_WL = (0.0381, 0.0208, 0.0158, 0.0131)


# Values for TTTEEE + BSH
redshift_TTTEEE = (1./0.1, 1./0.02, 1./0.005, 1./0.001 )
redshift_TTTEEE = map(sum,zip(redshift_TTTEEE,bar_width, bar_width))   #shift just for plotting
limit2_TTTEEE = (0.0254, 0.0159, 0.00951, 0.00714)

fig, ax = plt.subplots()

#font = {'family':'sans-serif','size':16}
#plt.rc('font',**font)
#plt.rc('legend',**{'fontsize':14})

xlim((8,1500))
ylim((0,0.04))
opacity = 0.8
error_config = {'ecolor': '0.3'}

rects1 = plt.bar(redshift, limit2, bar_width,
                 alpha=opacity,
                 color='#008AE6',
#                 yerr=std_men,
#                 error_kw=error_config)
#,
                 label='Planck (TT + lowP)+ lensing + BSH')


opacity = 0.5
rects2 = plt.bar(redshift_RSD_WL, limit2_RSD_WL, bar_width,
                 alpha=opacity,
                 color='#E03424',
#                 yerr=std_men,
#                 error_kw=error_config)
#,
                 label='Planck + lensing + WL + RSD')

rects3 = plt.bar(redshift_TTTEEE, limit2_TTTEEE, bar_width,
                 alpha=opacity,
                 color='m',
#                 yerr=std_men,
#                 error_kw=error_config)
#,
                 label='Planck (TT TE EE + lowP) + BSH')


leg = ax.legend()



sizeOfFont = 40
fontProperties = {'family':'sans-serif','sans-serif':['sans-serif'],
    'weight' : 'normal', 'size' : sizeOfFont}
ticks_font = font_manager.FontProperties(family='Helvetica', style='normal',
    size=sizeOfFont, weight='normal', stretch='normal')
rc('text', usetex=True)
rc('font',**fontProperties)


labels=['Planck + lensing + BSH']
plt.xlabel('$1+z_e$', fontsize=22,family='sans-serif')
plt.ylabel('$\Omega_e$', fontsize=22,family='sans-serif')
#plt.setp(axes, xticks=[10, 50, 100, 500, 1000])
plt.xticks([10, 20, 50, 100, 200, 500, 1000])

for tick in ax.xaxis.get_major_ticks():
	tick.label.set_fontsize(20)
for tick in ax.yaxis.get_major_ticks():
	tick.label.set_fontsize(20)	
#plt.title('Scores by group and gender')
#plt.xticks(index + bar_width, ('A', 'B', 'C', 'D'))
#plt.legend()
#g.add_legend('Planck + lensing + BSH')
#g.add_legend(labels, legend_loc='upper right', colored_text=True)
#g.add_legend(labels, legend_loc='upper right', colored_text=False)
l = plt.legend(loc=1,prop={'size':20}, frameon = False)
#l = legend()

colors = ["#008AE6", "#E03424", "m"]
for color,text in zip(colors,l.get_texts()):
    text.set_color(color)

ax.set_xscale('log')

#for text in l.get_texts():
#    text.set_color("#008AE6")

a = gca()
a.set_xticklabels(a.get_xticks(), fontProperties)
a.set_yticklabels(a.get_yticks(), fontProperties)

g.export(os.path.join(outdir,'compare_ede3_pol.pdf'))



plt.tight_layout()
#plt.show()
#pp = PdfPages(outdir,'multipage.pdf')



# a1    mean          standdev        limit1        limit2
#0.001  0.5021975E-02  0.3579554E-02  0.00634       0.0117
#0.005  0.5765400E-02  0.4542086E-02  0.00706       0.0145
#0.02   0.7795844E-02  0.5974882E-02  0.00975        0.0188
#0.1    0.1605600E-01  0.1028572E-01  0.0210        0.

#0.001  0.167444E-02  0.451126E-02  0.350380E-02  0.000000E+00  0.553297E-02  0.000000E+00  0.115441E-01  # v2
#0.01   0.283372E-01  0.257334E-01  0.156940E-01  0.000000E+00  0.330228E-01  0.000000E+00  0.537081E-01   # v1 
#0.03   0.423190E-01  0.433657E-01  0.231759E-01  0.000000E+00  0.548644E-01  0.000000E+00  0.847234E-01  # v4  
#0.1    0.120695E-01  0.480778E-01  0.265192E-01  0.000000E+00  0.100000E+00  0.000000E+00  0.100000E+00  # v5
#0.2    0.674810E-01  0.459228E-01  0.280909E-01  0.000000E+00  0.100000E+00  0.000000E+00  0.100000E+00  # v3
#0.3    0.162019E-01  0.493400E-01  0.287854E-01  0.000000E+00  0.100000E+00  0.000000E+00  0.100000E+00  # v8
