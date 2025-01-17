from CNCstruc.structure import CNC_class as CNC
from CNCstruc.utils import traj_reader as trj
from CNCstruc.analysis import cnc_analysis_utils
import matplotlib.pyplot as plt
import numpy as np

# CNC_form = 'Pristine'
# gro_file = f'./data/input/{CNC_form}/solute.gro'

def plot_xvg (x , y):
    font_size = 12
    plt.plot(x,y)
    
    plt.xlabel('Time (ps)' , fontsize = font_size)
    plt.ylabel('Number of H-bonds', fontsize = font_size)


# Inputs 
FF_directory = '/Users/arashelahi/Research/Papers/CNC_model/Review/GLYCAM-06/'
gro_file = FF_directory + 'solute.gro'
Data = trj.gro_reader(gro_file)
CNC_group = CNC.CNC_analys(Data)

# List of features : ['glycosidic', 'alcohols', 'twist', 'O2H_O6','O3H_O5','O6H_O3' , 'unit_cell']
for feature in ['twist']:
    cnc_analysis_utils.feature_analysis (CNC_group ,feature, FF_directory)
    # all_values = [value for layer in CNC_group.layer_vec for value in CNC_group.feature_dict[layer]['H_bonds'][feature]]
    # print (f"Average number of H-bonds for {feature} : {np.average(all_values):6.3f} ± {np.std(all_values):6.3f}")

# ############# For tg conformations:
# feature = 'alcohols'
# left_tg_values = np.array([value for layer in CNC_group.layer_vec for value in CNC_group.feature_dict[layer][feature]['left']]) == 'tg'
# right_tg_values = np.array([value for layer in CNC_group.layer_vec for value in CNC_group.feature_dict[layer][feature]['right']]) == 'tg'

# print("tg statistics  : %6.3f"  % ((np.average(left_tg_values)/2 + np.average(right_tg_values)/2) * 100))
# print("right statistics  : %6.3f " % (np.average(right_tg_values) * 100 ))


############ For Phi or Psi values:
# feature = 'glycosidic'
# phi_values = [value for layer in CNC_group.layer_vec for value in CNC_group.feature_dict[layer][feature]['Phi']]
# psi_values = [value for layer in CNC_group.layer_vec for value in CNC_group.feature_dict[layer][feature]['Psi']]
# Phi_average = np.mean(phi_values)
# Psi_average = np.mean(psi_values)
# print("%6.3f ± %6.3f" % (Phi_average , np.std(Phi_average)))
# print("%6.3f ± %6.3f" % (Psi_average , np.std(Psi_average)))




############# For twist calculation:
feature = 'twist'
twist_values = np.array([value for layer in CNC_group.layer_vec for value in CNC_group.feature_dict[layer][feature]['twist']])
print("twist statistics : %6.3f ± %6.3f" % (np.average(twist_values) , np.std(twist_values)))

############# For unit cell calculation:
# feature = 'unit_cell'
# for unit_cell_keys in ['dimension','angle']: 
#     for unit_cell_prop in CNC_group.feature_dict[feature][unit_cell_keys].keys():
#         unit_dim_vals = CNC_group.feature_dict[feature][unit_cell_keys][unit_cell_prop]
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
