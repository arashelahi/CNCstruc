from CNCstruc.structure import CNC_class as CNC
from CNCstruc.utils import traj_reader as trj
from CNCstruc.analysis import cnc_analysis_utils
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['mathtext.fontset'] = 'stixsans'  # A sans-serif math font
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.weight'] = 'regular'  # Adjust font weight as needed



def plot_xvg (ax , x , y , label  , x_label  , y_label  , color , title=None):
    font_size = 10
    leg_fontsize = 8
    labelpad = 0
    ax.plot(x,y,label = label , linewidth = 1 , color = color )
    ax.set_xlabel(x_label , fontsize = font_size , labelpad=labelpad +1 )
    ax.set_ylabel(y_label, fontsize = font_size, labelpad=labelpad)
    ax.set_title(title , fontsize = font_size ,loc='left',fontweight='bold')
    # ax.set_ylim([0,100])
    ax.set_xlim([0.25,1.2])
    leg = ax.legend(fontsize = leg_fontsize, ncol  = 3 , frameon=True,  fancybox=True , loc = 'upper center' , handletextpad=0.5 , handlelength = 1.0 , columnspacing=1.0 ,labelspacing=0.5 ,bbox_to_anchor=(0.5, 1.2))
    leg.get_frame().set_edgecolor('k')
    leg.get_frame().set_linewidth(0.4)
    
    # plt.show()


def plot_rdf(x,y,label = ''):
    # fig, ax = plt.subplots(1,1,figsize=(3,3))
    ax.plot(x,y,linewidth=1,color='k')
    ax.set_xlabel('r (nm)',fontsize=12)
    ax.set_ylabel('g(r)',fontsize=12)
    ax.set_xlim([0,1.0])
    # plt.show()

FF = 'OPLS_CM5'
color_vec = ['k' , 'r' , 'g' , 'b']
main_directory = '/Users/arashelahi/Research/Papers/CNC_model/Review/'
fig, ax  = plt.subplots(1,1,figsize=(3,3))
FF_directory = main_directory +  FF + '/Results/'
# rdf_list = ['O2_O6_rdf.xvg' , 'O3_O5_rdf.xvg' ,  'O3_O6_rdf.xvg']
rdf_list = ['c1_c1_rdf.xvg']
# legs = ['O2-O6' , 'O3-O5' , 'O6-O3']
legs = ['C1-C1']
for iter , rdf_file in enumerate(rdf_list):
    rdf_x , rdf_y = trj.xvg_reader(FF_directory + rdf_file , [1,2])
    plot_xvg(ax , rdf_x , rdf_y, label = legs[iter] , x_label = 'r (nm)' , y_label = 'g(r)' , color = color_vec[iter])
# rdf_x , rdf_y = trj.xvg_reader(FF_directory + 'rdf.xvg' , [1,2])
# plot_rdf(rdf_x , rdf_y)


### plotting parameters 
# x_label = 'Time (ns)'
# y_label = 'Hbond population (%)'

# # y_label = r'$tg$ population (%)'
# dt = 0.04
# legs = ['OPLS-CM5' , 'OPLS-AA' , 'GLYCAM06', 'CHARMM36' ]
# color_vec = ['k' , 'r' , 'g' , 'b']
# main_directory = '/Users/arashelahi/Research/Papers/CNC_model/Review/'
# # fig, axs = plt.subplots(3 , 1 , figsize=(3.0,7.0)) # for hb
# fig, axs = plt.subplots(2 , 2, figsize=(6,5)) # for tg
# axs = axs.flatten()

# title_vec = ['a) HO2-O6' , 'b) HO3-O5' , 'c) HO6-O3']
# time_series = np.arange(0,dt * 7501 , dt)
# for feature_iter , feature in enumerate(['alcohols']):
#     vals = [85 , 15 , 50 , 30]
#     for FF_iter , val in enumerate(vals):
#         y_vec = val * np.ones(7501)
#         time_series = np.arange(0,dt * 7501 , dt)
#         plot_xvg(axs , time_series , y_vec , label = legs[FF_iter] , x_label = x_label , y_label = y_label , color = color_vec[FF_iter])

# plt.savefig('/Users/arashelahi/Research/Papers/CNC_model/Review/' + 'tgs.png',dpi = 600,bbox_inches='tight')
# plt.show()
# title_vec=[r'$\textbf{a) HO2-O6}$',r'$\textbf{b) HO3-O5}$',r'$\textbf{c) HO6-O3}$']

# FFs = ['OPLS_CM5' , 'OPLS_old' , 'GLYCAM-06' , 'charmm']
# # FFs = ['OPLS_old']

# title_vec = ['a)' , 'b)' , 'c)' , 'd)']
# # # FFs = ['OPLS_CM5']
# for feature_iter , feature in enumerate(['unit_cell']):
#     for FF_iter , FF in enumerate(FFs):
#         iter = -1
#         main_directory = '/Users/arashelahi/Research/Papers/CNC_model/Review/'
#         FF_directory = main_directory +  FF + '/'
#         gro_file = FF_directory + 'solute.gro'
#         Data = trj.gro_reader(gro_file)
#         CNC_group = CNC.CNC_analys(Data)

    
        # cnc_analysis_utils.feature_analysis (CNC_group ,feature, FF_directory)
        ## for tg
        # left_tg_values = np.array([value for layer in CNC_group.layer_vec for value in CNC_group.feature_dict[layer][feature]['left']]) == 'tg'
        # right_tg_values = np.array([value for layer in CNC_group.layer_vec for value in CNC_group.feature_dict[layer][feature]['right']]) == 'tg'
        # left_tg_perc = left_tg_values.reshape(1 , -1).reshape(int(16*6) , 7501)
        # right_tg_perc = right_tg_values.reshape(1 , -1).reshape(int(16*6) , 7501)
        # tot_tg_perc = (np.mean(left_tg_perc , axis = 0) + np.mean(right_tg_perc , axis = 0))/2*100
        # plot_xvg(axs, time_series , tot_tg_perc , label = legs[FF_iter] , x_label = x_label , y_label = y_label , color = color_vec[FF_iter])

        ### for HB ###

        # all_values = [value for layer in CNC_group.layer_vec for value in CNC_group.feature_dict[layer]['H_bonds'][feature]]
        # average_values = np.mean(np.array(all_values),axis = 0)
        ## for unit cell
        
        # plot_xvg(axs[feature_iter], time_series , average_values , label = legs[FF_iter] , x_label = x_label , y_label = y_label , color = color_vec[FF_iter] , title = title_vec[feature_iter])
        # # feature = 'unit_cell'
        # for unit_cell_keys in ['angle']: 
        # for unit_cell_keys in ['dimension','angle']:
        #     if unit_cell_keys == 'dimension':
        #         y_label_suffix = 'dimension (nm)'
        #         unit_cell_prop_list = ['a','b','c']
        #     elif unit_cell_keys == 'angle':
        #         y_label = r'$\gamma$ angle ($^\circ$)'
        #         unit_cell_prop_list = ['gamma']
        #         # for unit_cell_prop in ['a','b','c']:
        #     for unit_cell_prop in unit_cell_prop_list:
        #         if unit_cell_keys == 'dimension':
        #             y_label = unit_cell_prop + ' ' + y_label_suffix
        #         # y_label = unit_cell_prop + ' ' + y_label_suffix
        #         print(unit_cell_prop)
        #         iter+=1
        #         unit_dim_vals = CNC_group.feature_dict[feature][unit_cell_keys][unit_cell_prop][:7501]
        #         time_series = np.arange(0 , dt * len(unit_dim_vals) , dt)
        #         plot_xvg(axs[iter], time_series , unit_dim_vals , label = legs[FF_iter] , x_label = x_label , y_label = y_label , color = color_vec[FF_iter] , title = title_vec[iter])
                # print(f'For FF {FF} The average value of {unit_cell_prop} is {np.average(unit_dim_vals):6.3f} ± {np.std(unit_dim_vals):6.3f}')
# # plt.subplots_adjust( hspace = 0.5)
plt.tight_layout()
plt.savefig(main_directory + 'c1_c1_rdf.png',dpi = 600,bbox_inches='tight')
# plt.show()
