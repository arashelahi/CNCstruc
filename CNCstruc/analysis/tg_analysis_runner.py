from turtle import color
import traj_reader as trj
import numpy as np
import matplotlib.pyplot as plt
from statistics import mode
import CNC_class as cnc
import cnc_analysis_utils
from pathlib import Path



def bar_plotting_tg_info(conf_vec , tg_info_left ,tg_info_right,layer,chain_iter,png_file_name , domain='interior'):
    fig, axs = plt.subplots(1,2,figsize=(1.28, 0.55), dpi=300)

    output_dir = Path('./tg_data/Images/%s' % domain)
    output_dir.mkdir(parents=True, exist_ok=True)  # This line creates the directory if it doesn't exist

    output_file = str(output_dir / f'{layer}_ch{chain_iter}_tgs{file_ext}.png')

    # output_dir = './tg_data/Images/'
    font_size = 12
    color_vec=['moccasin','coral','powderblue']
    alpha_val = 0.5 if domain =='exterior' else 1
    for tg_iter , tg_type in enumerate(conf_vec):
        bottom_coord_left = np.sum(tg_info_left[layer][chain_iter][:tg_iter]) if tg_iter > 0 else 0
        bottom_coord_right = np.sum(tg_info_right[layer][chain_iter][:tg_iter]) if tg_iter > 0 else 0
        axs[0].bar('hi', tg_info_left[layer][chain_iter][tg_iter] , edgecolor='black',align='center' , color = color_vec[tg_iter] ,linewidth=0.5 , bottom = bottom_coord_left , alpha = alpha_val)
        axs[1].bar('hi', tg_info_right[layer][chain_iter][tg_iter] , edgecolor='black', align='center' , color = color_vec[tg_iter] ,linewidth=0.5 , bottom = bottom_coord_right , alpha = alpha_val)
    plt.subplots_adjust(wspace = 0.2)
    axs[0].axis('off')
    axs[1].axis('off')
    plt.savefig(output_file,dpi=300, bbox_inches = 'tight', pad_inches = 0 , transparent = True)
    # plt.show()



# func_group = 'Bu'
# file_dir = './simulation_traj_topol/'+func_group+'/'
# gro_file = file_dir+'solute.gro'
# tg_data_direct = './tg_data/Modified/' + func_group + '/'


file_dir = './simulation_traj_topol/'
gro_file = file_dir+'solute.gro'
tg_data_direct = './tg_data/'

# gro_file = 'solute.gro'
domain = 'interior' ## 'interior' or 'exterior

Data = trj.gro_reader(gro_file)
CNC_group = cnc.CNC_analys(Data,domain)
file_ext = '' ## if the interior chains are needed
if domain!='interior':
    file_ext = '_ends' ## if the exterior chains are needed
    # CNC_group.layer_vec = list(CNC_group.layers.keys()) #if the other layers than the interior chains are needed, this line is added.

# tg_data_direct = './tg_data/'


chi_file_right = tg_data_direct + 'right_Chis%s_dist.xvg' % file_ext 
chi_p_file_right = tg_data_direct + 'right_Chi_ps%s_dist.xvg'  % file_ext
chi_file_left = tg_data_direct + 'left_Chis%s_dist.xvg'  % file_ext
chi_p_file_left = tg_data_direct + 'left_Chi_ps%s_dist.xvg'  % file_ext

resid_vec = [x for x in range(1,7)] ## The number of analyzed residues per side. 12 residues are analyzed, where 6 belongs to one side.

tg_info_right = cnc_analysis_utils.tg_analysis (chi_file_right,chi_p_file_right, CNC_group,resid_vec,domain=domain) # for the groups pointing to right
tg_info_left = cnc_analysis_utils.tg_analysis (chi_file_left,chi_p_file_left, CNC_group,resid_vec,domain=domain) # for the groups pointing to the left


conf_vec=['gt' , 'gg', 'tg']

conf_column = np.array([])
for layer in CNC_group.layer_vec:
    if domain=='interior': 
        chain_number_vec = CNC_group.layers[layer][1:-1] # For interior
    else:
        chain_number_vec = [CNC_group.layers[layer][0], CNC_group.layers[layer][-1]] if len(CNC_group.layers[layer]) > 1\
            else [CNC_group.layers[layer][0]] # for the exterior chains
    for chain_iter,chain_number in enumerate(chain_number_vec):
        new_array = np.array([[(tg_info_left[layer][chain_iter][0]+tg_info_right[layer][chain_iter][0])/2],\
                                    [(tg_info_left[layer][chain_iter][1]+tg_info_right[layer][chain_iter][1])/2],\
                                    [(tg_info_left[layer][chain_iter][2]+tg_info_right[layer][chain_iter][2])/2]])
        conf_column = new_array if conf_column.size == 0 else np.hstack((conf_column,new_array))
        # print("In layer %s, and chain %d the left gt/gg/tg percentage is %4.2f/%4.2f/%4.2f" % (layer,chain_iter,tg_info_left[layer][chain_iter][0],tg_info_left[layer][chain_iter][1],tg_info_left[layer][chain_iter][2]))
        # print("In layer %s, and chain %d the left gt/gg/tg percentage is %4.2f/%4.2f/%4.2f" % (layer,chain_iter,tg_info_right[layer][chain_iter][0],tg_info_right[layer][chain_iter][1],tg_info_right[layer][chain_iter][2]))
    
        
        # png_file_name = '%s%s_ch%s_tgs.png' % (file_ext , layer , chain_iter)
        # bar_plotting_tg_info(conf_vec , tg_info_left ,tg_info_right,layer,chain_iter,png_file_name,domain)
print('the average value of gt,gg and tg are %f, %f, %5.2f ± %5.2f' % (np.mean(conf_column[0]),np.mean(conf_column[1]),np.mean(conf_column[2]), np.std(conf_column[2])))