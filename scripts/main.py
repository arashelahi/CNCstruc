
from CNCstruct.structure import CNC_class as CNC


file_dir = './data/input/Pristine/'
filename = file_dir+'solute.gro'
Data = trj.gro_reader(filename)
domain = 'interior'

# Modif_CNC_layers = {
#     'L1': [18], 'L2': [19, 2], 'L3': [14,24,3], 'L4':  [15,25,31,6],
#     'L5': [11,21,26,32,7], 'L6': [12,22,27,33,36,10], 'L7': [13,23,28,34,8],
#     'L8': [16,29,35,9], 'L9': [17,30,4], 'L10':  [20,5], 'L11': [1]
# }
CNC_group = CNC_analys(Data, domain=domain)
# CNC_group.layer_vec = list(CNC_group.layers.keys()) # 
# CNC_group.layer_vec = list(CNC_group.layers.keys())[2:-2] ## This choice could be made when only iterior chains are needed. 
CNC_group.twist_angles()
# print

ndx_file='twist_chain.ndx'
if os.path.isfile(ndx_file):
    os.remove(ndx_file)

for chi_key in ['chain_twist']:
    # confor_per_side = {}
    # confor_per_side['right'] = [right_data for iter,right_data in enumerate(CNC_group.confors_data[chi_key]) if (iter)%8 in (0,1,2,3)]
    # confor_per_side['left']  = [left_data for iter,left_data in enumerate(CNC_group.confors_data[chi_key]) if (iter)%8 in (4,5,6,7)]
    # for side in ['left','right']:        
    trj.ndx_writer(ndx_file,CNC_group.twist_data[chi_key] ,chi_key+'s')