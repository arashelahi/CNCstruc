from CNCstruc.structure import CNC_class as CNC
# from CNCstruc.structure import CNC_class as CNC
from CNCstruc.utils import traj_reader as trj
from CNCstruc.analysis import Index_making_final as indx_gen
import os 

CNC_form = 'Pristine'
filepath= ''
# filepath = '/Users/arashelahi/Research/Papers/CNC_model/Review/Charmm/'
# filename = filepath + 'reordered_solute.gro'
filename = f'./data/input/{CNC_form}/solute.gro'

Data = trj.gro_reader(filename)
CNC_group = CNC.CNC_analys(Data )
feature = 'H_bonds'
CNC_group.get_indices(feature=feature)
# print('hi')
indx_gen.ndx_making(CNC_group, feature , output_path=filepath)
