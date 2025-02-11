'''
This module is used for the visualization of different properties of the CNC.
'''

from CNCstruc.utils import traj_reader as trj

from CNCstruc.structure import CNC_class as CNC

# import traj_reader as trj
import numpy as np
import matplotlib.pyplot as plt
# import CNC_class as cnc
# import cnc_analysis_utils
import seaborn as sns
import pandas as pd

## Plotting property against time
def plot_xvg (ax , x , y , label  , x_label  , y_label  , color , title=None):
    font_size = 10
    leg_fontsize = 8
    labelpad = 0
    ax.plot(x,y,label = label , linewidth = 1 , color = color )
    ax.set_xlabel(x_label , fontsize = font_size , labelpad=labelpad +1 )
    ax.set_ylabel(y_label, fontsize = font_size, labelpad=labelpad)
    ax.set_title(title , fontsize = font_size ,loc='left',fontweight='bold')
    ax.set_ylim([0,100])
    # ax.set_xlim([0.25,1.2])
    leg = ax.legend(fontsize = leg_fontsize, ncol  = 3 , frameon=True,  fancybox=True , loc = 'upper center' , handletextpad=0.5 , handlelength = 1.0 , columnspacing=1.0 ,labelspacing=0.5 ,bbox_to_anchor=(0.5, 1.2))
    leg.get_frame().set_edgecolor('k')
    leg.get_frame().set_linewidth(0.4)
### for the dihedral angles
def Ramachadron_plot(orig_layer,cent_layer):
    dict = {'Phi_orig' : orig_layer[0], 'Psi_orig' : orig_layer[1], 'Phi_cen' : cent_layer[0], 'Psi_cen' : cent_layer[1]}
    df = pd.DataFrame(data = dict)
    levels=5
    alpha=0.6
    cut=3
    thresh=0.001
    gridsize=200
    bw_adjust=0.8
    bw_method="silverman"
    linewidths=0.3
    contour_alpha = 1
    data_point_edge_colour = "#3c3c3c"
    data_point_colour = "#D4AB2D"
    cmap_fill = 'PuBu'
    cmap_cont = "Reds"


    fig, axs = plt.subplots(1, 2,figsize=(5.0,2.5))
    group_vec = ['orig' , 'cen']
    title_vec = ['b) Origin chains' , 'c) Center chains']

    # fig, ax = plt.subplots(figsize=(6, 6))
    for group_iter , group in enumerate(group_vec):

        sns.scatterplot(
            data=df,
            y="Phi_%s" %  group, x="Psi_%s" %  group,
            color=data_point_colour,
            s=1,
            linewidth=0.1,
            edgecolor=data_point_edge_colour,
            ax=axs[group_iter],
            # zorder=2,
            alpha=1
        )
        
        sns.kdeplot(
            data=df,
            y="Phi_%s" %  group, x="Psi_%s" %  group,
            levels=levels,
            fill=True,
            alpha=alpha,
            cut=cut,
            ax=axs[group_iter],
            thresh=thresh,
            gridsize=gridsize,
            bw_adjust=bw_adjust,
            bw_method=bw_method,
            cmap=cmap_fill
        )
        
        sns.kdeplot(
            data=df,
            y="Phi_%s" %  group, x="Psi_%s" %  group,
            levels=levels,
            fill=False,
            alpha=contour_alpha,
            cut=cut,
            ax=axs[group_iter],
            thresh=thresh,
            gridsize=gridsize,
            bw_adjust=bw_adjust,
            bw_method=bw_method,
            zorder=2 , 
            cmap= cmap_cont , 
            linestyles="--",
            # shade=False,
            # shade_lowest=False,
            linewidths = linewidths 
        )

        title_size = 12
        fontsize = 10
        axs[group_iter].set_ylabel(r'$\phi$ ($^\circ)$' , fontsize = fontsize, labelpad=-5)
        axs[group_iter].set_xlabel(r'$\psi$ ($^\circ)$' , fontsize = fontsize, labelpad=2)
        
        x_lim = [-160,-140]
        y_lim = [-103.00, -80.00]
        axs[group_iter].set_ylim(y_lim[0],y_lim[1])
        axs[group_iter].set_xlim(x_lim[0],x_lim[1])
        axs[group_iter].set_xticks(np.arange(x_lim[0],x_lim[1]+1,10))
        axs[group_iter].set_yticks(np.arange(y_lim[0]+3,y_lim[1]+1,10))
        axs[group_iter].tick_params(axis='both', which='major', labelsize=fontsize)
        axs[group_iter].set_title(title_vec[group_iter], fontsize=title_size,loc='left',fontweight='bold' , pad = 10)

    plt.tight_layout(w_pad = 0.1)
    plt.savefig('phi_psi.png', dpi=300 , bbox_inches = 'tight', pad_inches = 0.01 , transparent=True)
    plt.show()

### for the hydrogen bond populations