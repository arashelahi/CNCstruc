import traj_reader as trj
import numpy as np
import pandas as pd
import os
filename='solute.gro'
H2O_O6_index='H2O_O6_index.ndx'
H3O_O5_index='H3O_O5_index.ndx'
H6O_O3_index='H6O_O3_index.ndx'


Layer_info = {}

Layer_info['L1']=[11]
Layer_info['L2']=[12,22]
Layer_info['L3']=[4,13,23]
Layer_info['L4']=[5,14,24,31]
Layer_info['L5']=[1,6,15,25,32]
Layer_info['L6']=[2,7,16,26,33,36]
Layer_info['L7']=[3,8,17,27,34]
Layer_info['L8']=[9,18,28,35]
Layer_info['L9']=[10,19,29]
Layer_info['L10']=[20,30]
Layer_info['L11']=[21]


#### Clipping the data
Atom_num=423
Data=trj.gro_reader(filename)
Data_CNC=Data[Data['residue_name']=='BGLC']
res_num_replace=[x//Atom_num+1 for x in range(len(Data_CNC))]
Data_CNC['chain_number']=res_num_replace
resid_vec=np.arange(5,17,1)
Data_CNC=Data_CNC[Data_CNC['residue_number'].isin(resid_vec)]
layer_vec=['L%d' % x for x in range(3,10)]
chain_num_vec_nested= [Layer_info[x][1:-1]  for x in layer_vec]
chain_num_vec=[item for sublist in chain_num_vec_nested for item in sublist ]
Clipped_Data=Data_CNC[Data_CNC['chain_number'].isin(chain_num_vec)]


######################## Phi analysis
ndx_file='Confors.ndx'
# for ndx_file in conf_file:
if os.path.isfile(ndx_file):
        os.remove(ndx_file)

Data_Phi_ini=Clipped_Data[Clipped_Data['atom_name'].isin(['O4' ,'C4' ,'O5', 'C1'])] ## The dataframe for the phi-containing molecules.
Data_Psi_ini=Clipped_Data[Clipped_Data['atom_name'].isin(['O4' ,'C4' ,'C5', 'C1'])] ## The dataframe for the phi-containing molecules.
Phi_vec=[]
Psi_vec=[]

for layer in layer_vec:

    for chain_number in Layer_info[layer][1:-1]:
        for resid in resid_vec[:-1]:
            Phi_vec.append(Data_Phi_ini[Data_Phi_ini['residue_number']==resid][Data_Phi_ini['chain_number']==chain_number][Data_Phi_ini['atom_name']=='O5']['atom_number'].values.tolist()+\
                        Data_Phi_ini[Data_Phi_ini['residue_number']==resid][Data_Phi_ini['chain_number']==chain_number][Data_Phi_ini['atom_name']=='C1']['atom_number'].values.tolist()+\
                        Data_Phi_ini[Data_Phi_ini['residue_number']==resid+1][Data_Phi_ini['chain_number']==chain_number][Data_Phi_ini['atom_name']=='O4']['atom_number'].values.tolist()+\
                        Data_Phi_ini[Data_Phi_ini['residue_number']==resid+1][Data_Phi_ini['chain_number']==chain_number][Data_Phi_ini['atom_name']=='C4']['atom_number'].values.tolist())
      
            Psi_vec.append(Data_Psi_ini[Data_Psi_ini['residue_number']==resid][Data_Psi_ini['chain_number']==chain_number][Data_Psi_ini['atom_name']=='C1']['atom_number'].values.tolist()+\
                       Data_Psi_ini[Data_Psi_ini['residue_number']==resid+1][Data_Psi_ini['chain_number']==chain_number][Data_Psi_ini['atom_name']=='O4']['atom_number'].values.tolist()+\
                       Data_Psi_ini[Data_Psi_ini['residue_number']==resid+1][Data_Psi_ini['chain_number']==chain_number][Data_Psi_ini['atom_name']=='C4']['atom_number'].values.tolist()+\
                       Data_Psi_ini[Data_Psi_ini['residue_number']==resid+1][Data_Psi_ini['chain_number']==chain_number][Data_Psi_ini['atom_name']=='C5']['atom_number'].values.tolist())
Phi_vec_data=np.array(Phi_vec).flatten()
trj.ndx_writer(ndx_file,Phi_vec_data,'Phis') 
Psi_vec_data=np.array(Psi_vec).flatten()
trj.ndx_writer(ndx_file,Psi_vec_data,'Psis')

######################## Tg analysis

ndx_file='Chis.ndx'
# for ndx_file in conf_file:
if os.path.isfile(ndx_file):
        os.remove(ndx_file)
Data_Chi_ini=Clipped_Data[Clipped_Data['atom_name'].isin(['O5' ,'C5' ,'C6', 'O6'])] ## The dataframe for the phi-containing molecules.
Data_Chi_P_ini=Clipped_Data[Clipped_Data['atom_name'].isin(['C4' ,'C5' ,'C6', 'O6'])] ## The dataframe for the phi-containing molecules.
Chi_vec=[]
Chi_P_vec=[]
for layer in layer_vec:
    for chain_number in Layer_info[layer][1:-1]:
        for resid in resid_vec:
            Chi_vec.append(Data_Chi_ini[Data_Chi_ini['residue_number']==resid][Data_Chi_ini['chain_number']==chain_number][Data_Chi_ini['atom_name']=='O5']['atom_number'].values.tolist()+\
                        Data_Chi_ini[Data_Chi_ini['residue_number']==resid][Data_Chi_ini['chain_number']==chain_number][Data_Chi_ini['atom_name']=='C5']['atom_number'].values.tolist()+\
                        Data_Chi_ini[Data_Chi_ini['residue_number']==resid][Data_Chi_ini['chain_number']==chain_number][Data_Chi_ini['atom_name']=='C6']['atom_number'].values.tolist()+\
                        Data_Chi_ini[Data_Chi_ini['residue_number']==resid][Data_Chi_ini['chain_number']==chain_number][Data_Chi_ini['atom_name']=='O6']['atom_number'].values.tolist())
      
            Chi_P_vec.append(Data_Chi_P_ini[Data_Chi_P_ini['residue_number']==resid][Data_Chi_P_ini['chain_number']==chain_number][Data_Chi_P_ini['atom_name']=='C4']['atom_number'].values.tolist()+\
                       Data_Chi_P_ini[Data_Chi_P_ini['residue_number']==resid][Data_Chi_P_ini['chain_number']==chain_number][Data_Chi_P_ini['atom_name']=='C5']['atom_number'].values.tolist()+\
                       Data_Chi_P_ini[Data_Chi_P_ini['residue_number']==resid][Data_Chi_P_ini['chain_number']==chain_number][Data_Chi_P_ini['atom_name']=='C6']['atom_number'].values.tolist()+\
                       Data_Chi_P_ini[Data_Chi_P_ini['residue_number']==resid][Data_Chi_P_ini['chain_number']==chain_number][Data_Chi_P_ini['atom_name']=='O6']['atom_number'].values.tolist())
Chi_vec_data=np.array(Chi_vec).flatten()
trj.ndx_writer(ndx_file,Chi_vec_data,'Chis' ) 
Chi_P_vec_data=np.array(Chi_P_vec).flatten()
trj.ndx_writer(ndx_file,Chi_P_vec_data,'Chi_Ps')
