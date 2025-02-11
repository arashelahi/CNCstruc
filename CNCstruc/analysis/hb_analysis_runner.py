import traj_reader as trj
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import CNC_class as cnc
import os
import cnc_analysis_utils


# func_group = 'Bu'
# filename = './simulation_traj_topol/'+func_group+'/solute.gro'
# HB_dir = './HBs/Modified/' + func_group + '/'

file_dir = './simulation_traj_topol/'
filename = file_dir+'solute.gro'
HB_dir = './HBs/'

Data = trj.gro_reader(filename)
CNC_group = cnc.CNC_analys(Data)


HB_vec=[ 'H6O_O3' , 'H3O_O5' , 'H2O_O6']
Data=pd.DataFrame(columns=HB_vec)
max_per_chain_vec = [12,11,11]
domain = 'interior'
for hb_iter , hb_type in enumerate(HB_vec):
    max_hbnum = max_per_chain_vec[hb_iter]
    hb_info =  cnc_analysis_utils.hb_analysis (CNC_group , hb_type , max_hbnum , HB_dir , domain)
    cent_layer = []
    orig_layer = []
    for layer in CNC_group.layer_vec:
        if domain=='interior': 
            chain_number_vec = CNC_group.layers[layer][1:-1] # For interior
        else:
            chain_number_vec = [CNC_group.layers[layer][0], CNC_group.layers[layer][-1]] if len(CNC_group.layers[layer]) > 1\
            else [CNC_group.layers[layer][0]] # for the exterior chains
        for chain_iter , chain_number in enumerate(chain_number_vec):
            if hb_type == 'H6O_O3':               
                if layer == CNC_group.layer_vec[0] or layer==CNC_group.layer_vec[-1]:
                    continue
                elif chain_iter == len(chain_number_vec)-1:
                    if domain == 'interior':
                        continue
            
            hb_data = hb_info[layer][chain_iter] 
            all_hb_data =  np.mean(hb_data)
            all_hb_data_std = np.std(hb_data)
            if layer in ['L3','L5','L7','L9']:
                cent_layer +=  hb_data
            else:
                orig_layer += hb_data
    cen_hb_data = np.mean(cent_layer)
    cen_hb_data_std = np.std(cent_layer)

    orig_hb_data = np.mean(orig_layer)
    orig_hb_data_std =  np.std(orig_layer)

    # print('The %s occupancy for the center layer is %4.2f ± %4.2f\n' % (hb_type,cen_hb_data,cen_hb_data_std))
    # print('The %s occupancy for the origin layer is %4.2f ± %4.2f\n' % (hb_type,orig_hb_data,orig_hb_data_std))
    print('The overal %s occupancy is %4.2f ± %4.2f\n' % (hb_type,np.mean(cent_layer+orig_layer),np.std(cent_layer+orig_layer)))

# plt.show()


