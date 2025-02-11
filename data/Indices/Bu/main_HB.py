import traj_reader as trj
import numpy as np
import os 
filename='solute.gro'


Num_of_layers=7

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


#### HB data ####################


H2O_O6=['O2','HO2','O6']
H3O_O5=['O3','HO3','O5']
H6O_O3=['O6','HO6','O3']
H2O_O6_index='H2O_O6_index.ndx'
H3O_O5_index='H3O_O5_index.ndx'
H6O_O3_index='H6O_O3_index.ndx'
File_vec=[H2O_O6_index,H3O_O5_index,H6O_O3_index]
for file in File_vec:
    if os.path.isfile(file):
            os.remove(file)

Data_CNC_H2O_O6=Clipped_Data[Clipped_Data['atom_name'].isin(H2O_O6)]
Data_CNC_H3O_O5=Clipped_Data[Clipped_Data['atom_name'].isin(H3O_O5)]
Data_CNC_H6O_O3=Clipped_Data[Clipped_Data['atom_name'].isin(H6O_O3)]



for iter,group in enumerate(layer_vec):
    Data_CNC_H2O_O6_atoms=Data_CNC_H2O_O6[Data_CNC_H2O_O6['chain_number'].isin(Layer_info[group])]['atom_number']
    Data_CNC_H3O_O5_atoms=Data_CNC_H3O_O5[Data_CNC_H3O_O5['chain_number'].isin(Layer_info[group])]['atom_number']
    Data_CNC_H6O_O3_atoms=Data_CNC_H6O_O3[Data_CNC_H6O_O3['chain_number'].isin(Layer_info[group])]['atom_number']
    trj.ndx_writer(H2O_O6_index,Data_CNC_H2O_O6_atoms,group+'_H2O_O6')
    trj.ndx_writer(H3O_O5_index,Data_CNC_H3O_O5_atoms,group+'_H3O_O5')
    trj.ndx_writer(H6O_O3_index,Data_CNC_H6O_O3_atoms,group+'_H6O_O3')