import planckStyle as s
from pylab import *

g = s.getSinglePlotter(plot_data='/home/pettorin/codici/cosmomc_plots/output_getdist/clik/clik10/ede/edepar1/plot_data/')

import GetDistPlots, os
#g=GetDistPlots.GetDistPlotter('/home/pettorin/codici/cosmomc_plots/output_getdist/clik/clik9/wwa/plot_data/')
g.settings.setWithSubplotSize(2.0000)
g.settings.param_names_for_labels = '/home/pettorin/codici/git/cosmomcplanck/chains/paper/ede/edepar1/edepar1_lowTEB_plikTT.paramnames'
g.settings.legend_frac_subplot_margin=0.1

#ranges = [-2.1, 1, 0., 1.]

outdir='/home/pettorin/codici/cosmomc_plots/output_getdist/clik/clik10/ede/edepar1/figures/'
#labels=[s.planckTT,'+lensing',s.planckall, '+lensing','+BAO+HST+JLA' ]

roots = ['','_BAO_JLA_HSTlow', '_WL', '_RSD','_RSD_WL']

roots = ['edepar1_lowTEB_plikTT'+root for root in roots]


labels=['Planck', 'Planck + BSH', 'Planck + WL',  'Planck + RSD', 'Planck + WL + RSD']

g.settings.solid_colors=[('#8CD3F5', '#006FED'), ('#F7BAA6', '#E03424'), ('#D1D1D1', '#A1A1A1'), 'g', 'c', 'm']


# Planck: #00007f, -
# Planck + priority1: #008AE6, .
# Planck + WL: g (green), :
# Planck + RSD: #808080, --
# Planck + WL + RSD: #E03424, : -.
#g.settings.legend_fontsize = 20
#g.settings.lw_contour = 20
g.settings.linewidth=3.0
#g.settings.plot_args=[{'ls':'-','color':'#00007f'},{'marker':'.','color':'#008AE6'},{'ls':':','color':'green'},{'ls':'--','color':'#808080'},{'ls':'-.','color':'#E03424'}, {'ls':'-','color':'#008AE6'}]
g.settings.plot_args=[{'ls':'-','color':'#00007f'},{'ls':'-','color':'#008AE6'},{'ls':'-','color':'green'},{'ls':'-','color':'#808080'},{'ls':'-','color':'#E03424'}]


#g.settings.lmstyle = ['--','-.','-','-.']
#g.settings.lmcolor = ['#009966','#555555', '#FF82AB','#990099']
#g.settings.lineM = ['-k', '--b', '-.g', ':m', ':r', '-c', '-y']
#g.plots_1d(roots,['omegaede','w'],legend_labels=labels,legend_ncol=3, nx = 2)
#g.plots_1d(roots,['omegaede'],legend_labels=labels,legend_ncol=2, nx = 1)
#g.plots_1d(roots,['omegaede'])
#g.add_legend(labels, legend_loc='upper right', colored_text=True, label_order = [3,0,1,2]);

g.plot_1d(roots, 'omegaede')
g.add_legend(labels, legend_loc='upper right', colored_text=True)


g.export(os.path.join(outdir,'compare_edepar1_1D_linear.pdf'))
