#%%
from turtle import color
import traj_reader as trj
import numpy as np
import matplotlib.pyplot as plt
from statistics import mode
import CNC_class as cnc
import cnc_analysis_utils
import seaborn as sns
import pandas as pd
from pathlib import Path



file_dir = './simulation_traj_topol/'
gro_file = file_dir+'solute.gro'
Data = trj.gro_reader(gro_file)

# domain = 'interior'
twist_file_dir = './twist_data/'
# CNC_group = cnc.CNC_analys(Data,domain)
y_ind = 2
resid_vec = [x for x in range(1,11)] ## The number of analyzed residues per side. 12 residues are analyzed, where 6 belongs to one side.


CNC_group_vec = []
twist_info_vec= []
domain_vec= ['interior','exterior']
for domain_iter,domain in enumerate(domain_vec):
    CNC_group_vec.append(cnc.CNC_analys(Data,domain))
    twist_file = twist_file_dir + '%s_chain_twist.xvg' % domain
    twist_info_vec.append(cnc_analysis_utils.twist_analysis (twist_file,CNC_group_vec[-1],resid_vec , domain = domain))
# twist_file = twist_file_dir + '%s_chain_twist.xvg' % domain


# twist_info = cnc_analysis_utils.twist_analysis (twist_file,CNC_group,resid_vec , domain = domain)
#%%
from pathlib import Path
def box_plot(data,layer,chain_iter,domain='interior'):
    output_dir = Path('twist_data/Images/%s/' % (domain))
    output_dir.mkdir(parents=True, exist_ok=True)  # This line creates the directory if it doesn't exist
    output_file = str(output_dir / f'{layer}_ch{chain_iter}_{domain}_twist_boxplot.png')
    fig, ax = plt.subplots(figsize=((0.9, 0.55)),dpi=300)
    pd_data = pd.DataFrame(data,columns = ['Twist'])
    flierprops = dict(marker='x', markersize=0.1, linestyle='none', markeredgecolor='red', markerfacecolor='red')
    boxprops = dict(edgecolor='black')
    meanprops = dict(markerfacecolor='black', markeredgecolor='black')
    whiskerprops = dict(color='black')
    medianprops = dict(color='black')
    capprops = dict(color='black')
    box_color = 'skyblue' if domain=='interior' else 'lightcoral'
    sns.boxplot(data=pd_data,linewidth=0.5,ax=ax,color=box_color ,boxprops=boxprops, showfliers=False , \
         flierprops=flierprops, meanprops=meanprops , medianprops=medianprops, whiskerprops=whiskerprops, capprops=capprops)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([0])  
    ax.set(xticklabels=[]) 
    ax.set(yticklabels=[]) 
    ax.tick_params(axis='y', colors='darkred')
    # ax.tick_params(bottom=False)
    ax.legend([],[], frameon=False)
    ax.set_ylim(-30, 30)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.5)
    plt.savefig(output_file,transparent=True, bbox_inches = 'tight',dpi=300, pad_inches = 0.01)
domain_iter = 0
domain = domain_vec[domain_iter]; CNC_group = CNC_group_vec[domain_iter]; twist_info = twist_info_vec[domain_iter]
combined_twist_data = []
for layer in CNC_group.layer_vec:
    if domain=='interior': 
        chain_number_vec = CNC_group.layers[layer][1:-1] # For interior
    else:
        chain_number_vec = [CNC_group.layers[layer][0], CNC_group.layers[layer][-1]] if len(CNC_group.layers[layer]) > 1\
        else [CNC_group.layers[layer][0]] # for the exterior chains

    for chain_iter , chain_number in enumerate(chain_number_vec):
        twist_data = twist_info[layer][chain_iter]
        # if chain_iter == 0 and layer == CNC_group.layer_vec[0]:
            # box_plot(twist_data,layer,chain_iter,domain)
            # break
        # box_plot(twist_data,layer,chain_iter,domain)
        pd_data = pd.DataFrame(twist_data,columns = ['Twist'])

        combined_twist_data += twist_data
        twist_mean_data = np.mean(twist_data)
        twist_std_data = np.std(twist_data)

        # print('The twist angle the layer %s and chain ch%d is %4.2f ± %4.2f\n' % (layer,chain_iter, twist_mean_data,twist_std_data))
        print('The 1st, 2nd, 3rd and 4th quantiles for twist angle the layer %s and chain ch%d is %4.2f , %4.2f , %4.2f , %4.2f\n' % \
        (layer,chain_iter, np.quantile(twist_data,0.25),np.quantile(twist_data,0.5),np.quantile(twist_data,0.75),np.quantile(twist_data,1)))
# print('The combined twist angle is %4.2f ± %4.2f\n' % (np.mean(combined_twist_data), np.std(combined_twist_data)))
# row_num = int(len(combined_twist_data)/len(twist_info[layer][chain_iter]))
# col_num = int(len(twist_info[layer][chain_iter]))
# new_data = np.reshape(combined_twist_data,(row_num,col_num))
# dt = 0.04
# last_time = 0.04*(len(twist_info[layer][chain_iter]))
# time_vec = np.arange(0,last_time,dt)
# plt.plot(time_vec,np.mean(new_data,axis = 0))
# plt.ylim ([-10,10])

# plt.show()
# %%
plt.rcParams['mathtext.fontset'] = 'stixsans'  # A sans-serif math font
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.weight'] = 'regular'  # Adjust font weight as needed
large_list = twist_info['L3'][0]
fontsize = 12
font_leg = 10
color_vec = ['r','blue','orange','green']
fig, axs = plt.subplots(1,2,figsize=(7.5,3.0))
fig.subplots_adjust(wspace=0.3)
resid_vec_plot = [x+4 for x in resid_vec]

time_resid1 = [x*0.04 for x in range(0,7501)]
x_label_vec =['Time (ns)','residue number']
y_label_vec = [r'Twist angle ($^\circ$)',r'Twist angle ($^\circ$)']
title_vec = ['a)' ,'b)']
xlim_vec = [[-10,310],[4,15]]
y_lim = [-26, 23]
for axs_iter , ax in enumerate(axs):
    if axs_iter == 0:
        ax.plot(time_resid1,large_list[:7501], color = 'orange')

    if axs_iter == 1:
        
        sublists = [large_list[i:i + 7501] for i in range(0, len(large_list), 7501)]
        frame_vec_vec = [0,3000,6000,7500]
        time_vec = [int(x*0.04) for x in frame_vec_vec]
        scatter_data = {i: [sublist[i] for sublist in sublists] for i in frame_vec_vec}
        for time_iter , time in enumerate(time_vec):
            # for resids in resid_vec:
            ax.plot(resid_vec_plot , scatter_data[frame_vec_vec[time_iter]],label = '%d ns' % time , color = color_vec[time_iter])
        leg = ax.legend(ncol=2, fontsize=font_leg, frameon=True,  fancybox=True,handlelength = 1,handletextpad=0.3, columnspacing=0.5 , bbox_to_anchor=(0.70, 1.15), loc='upper center' , framealpha=1)
        leg.get_frame().set_edgecolor('k')
        leg.get_frame().set_linewidth(0.5)
        
        x_lim = [4,14]
        # ax.hlines(0,resid_vec_plot[0],resid_vec_plot[1], colors='k', linestyles='dashed',linewidth=1.0)


        ax.set_xticks(np.arange(x_lim[0],x_lim[1]+1,2))
    ax.set_xlim(xlim_vec[axs_iter])
    ax.hlines(0,xlim_vec[axs_iter][0],xlim_vec[axs_iter][1], colors='k', linestyles='dashed',linewidth=1.0,zorder = 2)
    ax.set_yticks(np.arange(y_lim[0]+6,y_lim[1]+1,10))
    ax.set_xlabel(x_label_vec[axs_iter], fontsize = 12, labelpad=5)
    ax.set_ylabel(y_label_vec[axs_iter], fontsize = 12)
    ax.set_ylim(y_lim)
    ax.set_title(title_vec[axs_iter],loc='left',fontweight='bold')
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    
plt.savefig('twist_data/Images/twist_time_resid.png',transparent=True, bbox_inches = 'tight',dpi=300, pad_inches = 0.01)
    

# plt.tick_params(axis='both', which='major', labelsize=fontsize)
# axs[group_iter].(title_vec[group_iter], fontsize=fontsize,loc='left',fontweight='bold' , pad = 10)

# %%
