# from enum import _EnumNames
import traj_reader as trj
import numpy as np
import pandas as pd
import os


func_group = 'Bu'
dir_path = './unit_cell_data/Modified/' + func_group + '/'
# filename = './simulation_traj_topol/'+func_group+'/solute.gro'

# dir_path='unit_cell_data/'

unit_dim_file=dir_path+'unit_dist.xvg'
unit_dim_vec = ['a','b','c']

angle_vec = ['alpha' , 'beta' , 'gamma']

for unit_iter , unit_dim in enumerate(unit_dim_vec) : 
    time,unit_dim_data=trj.xvg_reader(unit_dim_file,[1,unit_iter+2])
    print ('the %s dim is %8.2f ± %8.2f ' % (unit_dim, 10*np.average(unit_dim_data),10*np.std(unit_dim_data)))
for unit_iter , unit_ang in enumerate(angle_vec) : 
    filename = dir_path + unit_ang+'_ang_dist.xvg'
    time,unit_ang_data=trj.xvg_reader(filename,[1,2])
    print ('the %s angle is %8.2f ± %8.2f ' % (unit_ang, np.average(unit_ang_data),np.std(unit_ang_data)))
