import numpy as np
import pandas as pd
from pandas.core.frame import DataFrame
import traj_reader as trj
import os
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
molec_name='cnc_chain.pdb' 
resid_name='ethyl.pdb'    
data_molec=trj.pdb_reader(molec_name) ## mother molecule 
data_residue=trj.pdb_reader(resid_name)  ## functional molecule
molec_to_replace=np.array([287,48]) ## atoms that are replaced from functional to mother
# molec_to_replace=np.array([287,48,245,119,203,161,55,329])
molec_to_delete=molec_to_replace+3 ## atoms in mother molecule needs to be deleted
resid_to_replace=1                 ## Atom index number in the functional molecule replacing the mother atom
# print(data_molec)
# print(data_residue)
Data_final=surf_modifying(data_molec,data_residue,molec_to_delete,molec_to_replace,resid_to_replace)
Data_to_write=Data_final
Data_to_write.loc[:,'X':'Z']=0.1*Data_final.loc[:,'X':'Z'] ## from angstrom to nm
gro_file='surf_ethyl2.gro'
if os.path.exists(str(gro_file)):
    os.remove(gro_file) 
box_size='10 10 10'
res_name='BGL'
trj.gro_writer(gro_file,Data_to_write,box_size,res_name)

