from CNCstruc.structure import CNC_class as CNC
from CNCstruc.utils import traj_reader as trj, Indexing as indx_gen
# from CNCstruc.analysis import Indexing as indx_gen

import os 

## initialization of the CNC object
CNC_form = 'Pristine'
filepath= ''
filename = f'./data/input/{CNC_form}/solute.gro'
Data = trj.gro_reader(filename)
CNC_group = CNC.CNC_analys(Data)
## indexing a target feature
feature = 'H_bonds'
CNC_group.get_indices(feature=feature)
indx_gen.ndx_making(CNC_group, feature , output_path=filepath)
