from re import X
import numpy as np
import pandas as pd
from pandas.core.frame import DataFrame
import traj_reader as trj
import os
import argparse
## this function gets the mother molecule, the functional to add, the functional atom replacing the atom in originial molecule, the atoms 
## that should be deleted from the mother molecule.
def surf_modifying(data_molec,data_residue,molec_to_delete,molec_to_replace_vec,resid_to_replace):
    Data_final=data_molec.drop(molec_to_delete)
    for molec_to_replace in molec_to_replace_vec:
        Del_r=data_molec.loc[molec_to_replace,'X':'Z']-data_residue.loc[resid_to_replace,'X':'Z']
        Data_resid_new=data_residue
        Data_resid_new.loc[:,'X':'Z']=Data_resid_new.loc[:,'X':'Z']+Del_r
        Data_resid_new=Data_resid_new.drop([resid_to_replace])
        Data_final=pd.concat([Data_final, Data_resid_new],ignore_index=True)
    Data_final.index=Data_final.index+1
    return Data_final
###################################################################
def gro_prep(data_residue,resid_to_replace,data_molec,Molec_Number,Atom_number):
    for molec_iter in range(1,Molec_Number+1):
        data_chain=data_molec[(molec_iter-1)*Atom_number:(molec_iter)*Atom_number]
        molec_to_replace_new=molec_to_replace+(molec_iter-1)*Atom_number
        molec_to_delete_new =molec_to_delete+(molec_iter-1)*Atom_number
        Data_final_new=surf_modifying(data_chain,data_residue,molec_to_delete_new,molec_to_replace_new,resid_to_replace)
        if molec_iter==1:
            Data_final=Data_final_new
            continue
        Data_final=pd.concat([Data_final, Data_final_new],ignore_index=True)
    Data_to_write=Data_final
    Data_to_write.loc[:,'X':'Z']=0.1*Data_final.loc[:,'X':'Z'] ## from angstrom to nm
    return Data_to_write
###################################################################
# parser = argparse.ArgumentParser(description='generate the coordinate file of the modified CNC')

# # parser.add_argument("-s", "--functionside", help="input the cnc molecule file", type=str)
# # parser.add_argument("-d", "--dos", help="density of substitution", type=str)
# parser.add_argument("-f", "--functional", help="input the functional group", type=str)
# parser.add_argument("-dos", "--dos_list",nargs='+',type=int)

# args = parser.parse_args()
# func_group=args.functional
# dos_vec=args.dos_list

func_group='hexyl'
dos_vec=100

# func_group='hexyl' 
Atom_number=331
Molec_Number=1
resid_to_replace=1                 ## Atom index number in the functional molecule replacing the mother atom
side_vec=['left','right']

# dos_vec=[12.5,25,37.5]
molec_to_replace_all=[18, 56, 94, 132, 170, 208, 246, 284] ## atoms that are replaced from functional to mother ## For Left
Site_numbers=np.array([0,1,2,3,4,5,6,7]) # 100 %
first_del_atom=308 ## the index of the first atom to be deleted

for dos in dos_vec:
    for side in side_vec:
        if side=='left':
            molec_to_replace_all=[39,77,115,153,191,229,267,307] ## atoms that are replaced from functional to mother ## For Left
        if side=='right':
            molec_to_replace_all=[18, 56, 94, 132, 170, 208, 246, 284]
        if dos==12:
            Site_numbers=np.array([2]) # 25 %
        if dos==25:
            Site_numbers=np.array([2,5]) # 25 % 
        if dos==37:
            Site_numbers=np.array([0,2,5]) # 25 % 
        if dos==50:
            Site_numbers=np.array([0,2,5,7]) # 50%
        if dos== 75:
            Site_numbers=np.array([0,2,3,5,6,7]) # 75%
        if dos==100:
            Site_numbers=np.array([0,1,2,3,4,5,6,7]) # 100 %
        molec_name='chain_%s.pdb' % side ;data_molec=trj.pdb_reader(molec_name) ## mother molecule
        resid_name='%s_%s.pdb' % (func_group,side) ;data_residue=trj.pdb_reader(resid_name)  ## functional molecule
        molec_to_replace=np.array([molec_to_replace_all[x] for x in  Site_numbers])
        molec_to_delete=np.sort(np.array((first_del_atom+3*Site_numbers).tolist()+(first_del_atom+3*Site_numbers+1).tolist()+(first_del_atom+3*Site_numbers+2).tolist()))###################################################################
        gro_file='chain_%s_%s_%s.gro' % (side,func_group,str(dos))
        Data_to_write=gro_prep(data_residue,resid_to_replace,data_molec,Molec_Number,Atom_number)
        if os.path.exists(str(gro_file)):
            os.remove(gro_file) 
        box_size='10 10 10'
        res_name='BGL'
        trj.gro_writer(gro_file,Data_to_write,box_size,res_name)

###################################################################
# Data_to_write=gro_prep(data_residue,resid_to_replace,data_molec,Molec_Number,Atom_number)
# if os.path.exists(str(gro_file)):
#     os.remove(gro_file) 
# box_size='10 10 10'
# res_name='BGL'
# trj.gro_writer(gro_file,Data_to_write,box_size,res_name)

