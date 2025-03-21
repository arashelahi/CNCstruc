from CNCstruc.structure import CNC_class as CNC
from CNCstruc.utils import traj_reader as trj, Indexing as indx_gen, surf_functionalize as srf
# from CNCstruc.analysis import Indexing as indx_gen

import os 

## initialization of the CNC object
CNC_form = 'Pristine'
filepath= ''
filename = f'./data/input/{CNC_form}/solute.gro'
Data_raw = trj.gro_reader(filename)
# The information for the gro file and spatial conformation is turned into a class
CNC_group = CNC.CNC_analys(Data_raw)

func = 'ethyl'
srf.material_prep(CNC_group , func)
