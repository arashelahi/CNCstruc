from CNCstruc.structure import CNC_class as CNC
# from CNCstruc.structure import CNC_class as CNC
from CNCstruc.utils import traj_reader as trj
import os

CNC_form = 'Pristine'
filename = f'./data/input/{CNC_form}/solute.gro'

Data = trj.gro_reader(filename)
CNC_group = CNC.CNC_analys(Data)

CNC_group.twist_angles()



ndx_file='./data/output/twist/%s/twist_chain.ndx' % CNC_form
if os.path.isfile(ndx_file):
    os.remove(ndx_file)

for chi_key in CNC_group.twist_data.keys():   
    trj.ndx_writer(ndx_file,CNC_group.twist_data[chi_key] ,chi_key+'s')