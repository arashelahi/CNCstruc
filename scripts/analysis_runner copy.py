from CNCstruc.structure import CNC_class as CNC
from CNCstruc.utils import traj_reader as trj
from CNCstruc.analysis import cnc_analysis_utils
from CNCstruc.analysis import phi_analysis_runner
import matplotlib.pyplot as plt
import numpy as np

# CNC_form = 'Pristine'
# gro_file = f'./data/input/{CNC_form}/solute.gro'

def plot_xvg (x , y , hb_type = 'O3H_O5'):
    font_size = 12
    plt.plot(x,y,label = hb_type)
    # plt.ylim([0,100])
    plt.xlabel('Time (ps)' , fontsize = font_size)
    plt.ylabel('Number of H-bonds', fontsize = font_size)
    plt.legend()
    # plt.show()


# Inputs 
# FFs = ['OPLS_CM5' , 'OPLS_old' , 'GLYCAM-06' , 'charmm']
FFs = ['OPLS_CM5' ]
main_directory = '/Users/arashelahi/Research/Papers/CNC_model/Review/'
T = ''
FF_directory = main_directory +  FFs[0] + '/' + T 
gro_file = main_directory +  FFs[0] + '/' + 'solute.gro'
Data = trj.gro_reader(gro_file)
CNC_group = CNC.CNC_analys(Data)
first_time = 2250
last_time = 2501

# List of features : ['glycosidic', 'alcohols', 'twist', 'O2H_O6','O3H_O5','O6H_O3' , 'unit_cell']

for feature in [ 'glycosidic']:
    cnc_analysis_utils.feature_analysis (CNC_group ,feature, FF_directory)
    # all_values = np.array([value for layer in CNC_group.layer_vec for value in CNC_group.feature_dict[layer]['H_bonds'][feature]])[:,:last_time]
    # average_values = np.mean(np.array(all_values),axis = 0)
    # time_series = np.arange(0,40 * np.shape(all_values)[1] , 40)
    # plot_xvg(time_series , average_values , feature)
    # for chain in range(len(all_values)):
        
    #     plot_xvg(time_series , average_values])

    # print (f"Average number of H-bonds for {feature} : {np.average( np.array(all_values)):6.2f} ± {np.std( np.array(all_values)):6.2f}")
# plt.show()

# ############# For tg conformations:
# feature = 'alcohols'
# left_tg_values = np.array([value for layer in CNC_group.layer_vec for value in CNC_group.feature_dict[layer][feature]['left']]) == 'gt'
# right_tg_values = np.array([value for layer in CNC_group.layer_vec for value in CNC_group.feature_dict[layer][feature]['right']]) == 'gt'
# max_frame = int(np.shape(left_tg_values)[1]/6)
# left_tg_perc = left_tg_values.reshape(1 , -1).reshape(int(16*6) , max_frame)[:,first_time:last_time]
# right_tg_perc = right_tg_values.reshape(1 , -1).reshape(int(16*6) , max_frame)[:,first_time:last_time]
# total_tg = (np.mean(left_tg_perc , axis = 0) + np.mean(right_tg_perc , axis = 0))/2*100
# print("tg statistics  : %6.3f ± %6.3f"  % (np.average(total_tg) , np.std(total_tg)))
# print("right statistics  : %6.3f " % (np.average(right_tg_values) * 100 ))


############ For Phi or Psi values:
feature = 'glycosidic'
# phi_values = [value for layer in CNC_group.layer_vec for value in CNC_group.feature_dict[layer][feature]['Phi']]
# psi_values = [value for layer in CNC_group.layer_vec for value in CNC_group.feature_dict[layer][feature]['Psi']]
orig_phi_vec = np.array([value for layer in ['L4' , 'L6' , 'L8'] for value in CNC_group.feature_dict[layer][feature]['Phi']]).reshape(1 , -1).reshape(8 , 11 , 7501)
orig_psi_vec = np.array([value for layer in ['L4' , 'L6' , 'L8'] for value in CNC_group.feature_dict[layer][feature]['Psi']]).reshape(1 , -1).reshape(8 , 11 , 7501)
cent_phi_vec = np.array([value for layer in ['L3' , 'L5' , 'L7' , 'L9'] for value in CNC_group.feature_dict[layer][feature]['Phi']]).reshape(1 , -1).reshape(8 , 11 , 7501)
cent_psi_vec = np.array([value for layer in ['L3' , 'L5' , 'L7' , 'L9'] for value in CNC_group.feature_dict[layer][feature]['Psi']]).reshape(1 , -1).reshape(8 , 11 , 7501)

orig_layer = np.mean(np.mean(orig_phi_vec , axis = 1) , axis = 0) , np.mean(np.mean(orig_psi_vec , axis = 1) , axis = 0)
cent_layer = np.mean(np.mean(cent_phi_vec , axis = 1) , axis = 0) , np.mean(np.mean(cent_psi_vec , axis = 1) , axis = 0)
phi_analysis_runner.Ramachadron_plot(orig_layer, cent_layer)

# Phi_average = np.mean(phi_values)
# Psi_average = np.mean(psi_values)
# print("%6.3f ± %6.3f" % (Phi_average , np.std(Phi_average)))
# print("%6.3f ± %6.3f" % (Psi_average , np.std(Psi_average)))




############# For twist calculation:
# feature = 'twist'
# FFs = ['OPLS_CM5' , 'OPLS_old' , 'GLYCAM-06' , 'charmm']

# # for feature in ['twist']:
# for FF_iter , FF in enumerate(FFs):
#     FF_directory = main_directory +  FF + '/'
#     gro_file = FF_directory + 'solute.gro'
#     Data = trj.gro_reader(gro_file)
#     CNC_group = CNC.CNC_analys(Data)
#     cnc_analysis_utils.feature_analysis (CNC_group ,feature, FF_directory)
#     twist_values = np.array([value for layer in CNC_group.layer_vec for value in CNC_group.feature_dict[layer][feature]['twist']])
#     twist_per_chain = twist_values.reshape(1 , -1).reshape(16 , 10 , 7501)
#     ave_per_chain = np.mean(twist_per_chain , axis = 1)
#     # average_values = np.mean(abs(twist_values),axis = 0)
#     # average_values = np.mean(twist_values,axis = 0)
#     # time_series = np.arange(0,40 * len(average_values) , 40)
#     # plot_xvg(time_series , average_values , feature)
#     print("for the %s FF, twist statistics without abs : %6.3f ± %6.3f" % (FF , np.average(ave_per_chain) , np.std(ave_per_chain)))
#     # print("for the %s FF, twist statistics with abs.   : %6.3f ± %6.3f" % (FF , np.average(abs(twist_values)) , np.std(abs(twist_values))))
#     # plt.show()

# ############# For unit cell calculation:
# feature = 'unit_cell'
# for unit_cell_keys in ['dimension' , 'angle']: 
#     # for unit_cell_prop in CNC_group.feature_dict[feature][unit_cell_keys].keys():
#     if unit_cell_keys == 'dimension':
#         unit_cell_prop_list = ['a','b','c']
#     elif unit_cell_keys == 'angle':
#         unit_cell_prop_list = ['gamma']
#         # for unit_cell_prop in ['a','b','c']:
#     for unit_cell_prop in unit_cell_prop_list:
#         unit_dim_vals = CNC_group.feature_dict[feature][unit_cell_keys][unit_cell_prop]
#     # time_series = np.arange(0,40 * len(unit_dim_vals) , 40)
#     # plot_xvg(time_series , unit_dim_vals , unit_cell_prop)
# # # plt.show()

#         print(f'The average value of {unit_cell_prop} is {np.average(unit_dim_vals):6.3f} ± {np.std(unit_dim_vals):6.3f}')
# print(CNC_group.feature_dict['unit_cell'])

# print("right statistics  : %6.3f " % (np.average(right_tg_values) * 100 ))
"""
Expected Output:
    L3: [-2.388, -3.478, 8.061, -3.842, -5.488]
    L4: [4.566, 12.224, 21.719, 18.615, 6.833]
    L5: [2.067, -2.57, 10.756, -7.849, 2.624]
    L6: [1.193, -4.351, 4.355, -1.56, 1.349]
    L7: [-3.994, -6.053, -0.008, -3.299, 4.895]
    L8: [21.397, 7.264, 7.869, 7.193, 17.738]
    L9: [6.835, 5.451, -2.087, 3.839, 6.447]
"""
# %%
