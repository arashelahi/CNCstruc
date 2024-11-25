import traj_reader as trj
import numpy as np
import pandas as pd
import os
filename='solute.gro'
based_data='../../solute.gro'

Layer_info = {}

Layer_info['L1'] =       [18]
Layer_info['L2'] =      [19,2]
Layer_info['L3'] =    [14,24,3]
Layer_info['L4'] =   [15,25,31,6]
Layer_info['L5'] =  [11,21,26,32,7]
Layer_info['L6'] =[12,22,27,33,36,10]
Layer_info['L7'] =  [13,23,28,34,8]
Layer_info['L8'] =   [16,29,35,9]
Layer_info['L9'] =    [17,30,4]
Layer_info['L10']=      [20,5]
Layer_info['L11']=        [1]

########################################################################
Atom_num=423
resid_vec=trj.gro_reader(based_data)[:Atom_num]['residue_number'].tolist()

# #### Clipping the data
Atom_num=423
Data=trj.gro_reader(filename)
Data_CNC=Data[Data['residue_name']=='CNCM']
Data_CNC['chain_number']=Data_CNC['residue_number']
Data_CNC['residue_number']=resid_vec*len(np.unique(Data_CNC['chain_number']))
resid_vec=np.arange(5,17,1)
Data_CNC=Data_CNC[Data_CNC['residue_number'].isin(resid_vec)]
layer_vec=['L%d' % x for x in range(3,10)]
chain_num_vec_nested= [Layer_info[x][1:-1]  for x in layer_vec]
chain_num_vec=[item for sublist in chain_num_vec_nested for item in sublist ]
Clipped_Data=Data_CNC[Data_CNC['chain_number'].isin(chain_num_vec)][Data_CNC['atom_name']=='C1']
####################### Unit dimensions

######################### c-axis
c_vec=[]
for chain_number in chain_num_vec:
    for resid in resid_vec[:-2]:
        c_vec.append(Clipped_Data[Clipped_Data['residue_number'].isin([resid,resid+2])][Clipped_Data['chain_number']==chain_number]['atom_number'].values.tolist())
c_vec_data=np.array(c_vec).flatten()


######################### b-axis
b_vec=[]
for layer in layer_vec[1:-1]: ## the first two layers contain only one chain and the b calculation might not apply.

    for chain_number_iter,chain_number in  enumerate(Layer_info[layer][1:-2]):
        for resid in resid_vec:
            b_vec.append(Clipped_Data[Clipped_Data['residue_number']==resid]\
                        [Clipped_Data['chain_number'].isin(Layer_info[layer][1:-1][chain_number_iter:chain_number_iter+2])]['atom_number'].values.tolist())
b_vec_data=np.array(b_vec).flatten()
######################### a-axis
a_vec=[]
for iter_layer,layer in enumerate(layer_vec): ## the first two layers contain only one chain and the b calculation might not apply.
    if iter_layer<2:continue
    if len(Layer_info[layer])>len(Layer_info[layer_vec[iter_layer-2]]):
        chain_adim_vec=Layer_info[layer][2:-2]
        top_iter=0
    elif len(Layer_info[layer])==len(Layer_info[layer_vec[iter_layer-2]]):
        chain_adim_vec=Layer_info[layer][1:-1]
    else:
        chain_adim_vec=Layer_info[layer][1:-1]
        top_iter=1
    for chain_number_iter,chain_number in  enumerate(chain_adim_vec):
        other_chain=Layer_info[layer_vec[iter_layer-2]][1:-1][top_iter+chain_number_iter]

        for resid in resid_vec:
            a_vec.append(Clipped_Data[Clipped_Data['residue_number']==resid]\
                        [Clipped_Data['chain_number'].isin([chain_number,other_chain])]['atom_number'].values.tolist())
a_vec_data=np.array(a_vec).flatten()

######################### ANGLES #########################

# #### alpha angle
# alpha_vec=[]
# for layer in layer_vec[1:-1]: ## the first two layers contain only one chain and the b calculation might not apply.

#     for chain_number_iter,chain_number in  enumerate(Layer_info[layer][1:-2]):
#         for resid in resid_vec[:-2]:
#             alpha_vec.append(Clipped_Data[Clipped_Data['residue_number']==resid+2][Clipped_Data['chain_number']==chain_number]['atom_number'].values.tolist()+\
#                 Clipped_Data[Clipped_Data['residue_number']==resid]\
#                         [Clipped_Data['chain_number'].isin(Layer_info[layer][1:-1][chain_number_iter:chain_number_iter+2])]['atom_number'].values.tolist())
# alpha_vec_data=np.array(alpha_vec).flatten()

# #### beta angle

# beta_vec=[]
# for layer in layer_vec[2:]: ## the first two layers contain only one chain and the b calculation might not apply.

#     for chain_number_iter,chain_number in  enumerate(Layer_info[layer][1:-1]):
#         if chain_number-2 in chain_num_vec:
#             other_chain=chain_number-2
#         else:continue
#         for resid in resid_vec[:-2]:
#             beta_vec.append(
#                 Clipped_Data[Clipped_Data['residue_number']==resid]\
#                         [Clipped_Data['chain_number'].isin([chain_number,other_chain])]['atom_number'].values.tolist()+\
#             Clipped_Data[Clipped_Data['residue_number']==resid+2][Clipped_Data['chain_number']==chain_number]['atom_number'].values.tolist())
# beta_vec_data=np.array(beta_vec).flatten()


# #### beta angle

# beta_vec=[]
# for layer in layer_vec[2:]: ## the first two layers contain only one chain and the b calculation might not apply.

#     for chain_number_iter,chain_number in  enumerate(Layer_info[layer][1:-1]):
#         if chain_number-2 in chain_num_vec:
#             other_chain=chain_number-2
#         else:continue
#         for resid in resid_vec[:-2]:
#             beta_vec.append(
#                 Clipped_Data[Clipped_Data['residue_number']==resid]\
#                         [Clipped_Data['chain_number'].isin([chain_number,other_chain])]['atom_number'].values.tolist()+\
#             Clipped_Data[Clipped_Data['residue_number']==resid+2][Clipped_Data['chain_number']==chain_number]['atom_number'].values.tolist())
# beta_vec_data=np.array(beta_vec).flatten()

# #### gamma angle


# gamma_vec=[]
# for layer in layer_vec[2:-1]: ## the first two layers contain only one chain and the b calculation might not apply.

#     for chain_number_iter,chain_number in  enumerate(Layer_info[layer][1:-1]):
#         if chain_number-2 in chain_num_vec:
#             top_chain=chain_number-2
#             if Layer_info[layer][Layer_info[layer].index(chain_number)-1] in chain_num_vec:
#                 left_chain=Layer_info[layer][Layer_info[layer].index(chain_number)-1]
#             else:continue
#         else:continue
#         # print('%d   %d   %d' % (top_chain,chain_number,left_chain))
#         for resid in resid_vec:
#             gamma_vec.append(Clipped_Data[Clipped_Data['residue_number']==resid]\
#                            [Clipped_Data['chain_number'].isin([chain_number,top_chain])]['atom_number'].values.tolist()+\
#                Clipped_Data[Clipped_Data['residue_number']==resid][Clipped_Data['chain_number']==left_chain]['atom_number'].values.tolist())
# gamma_vec_data=np.array(gamma_vec).flatten()




# ######################### ANGLES #########################


# # File_list = ['c_dim.ndx','b_dim.ndx','a_dim.ndx','alpha_ang.ndx','beta_ang.ndx','gamma_ang.ndx']
# # group_list=['c_axis','b_axis','a_axis','alpha_ang','beta_ang','gamma_ang']
# # Data_list=[c_vec_data,b_vec_data,a_vec_data,alpha_vec_data,beta_vec_data,gamma_vec_data]
# # for iter,ndx_file in enumerate(File_list):
# #     if os.path.isfile(ndx_file):
# #         os.remove(ndx_file)
# #     trj.ndx_writer(ndx_file,Data_list[iter],group_list[iter])

Dim_data_vec=[a_vec_data,b_vec_data,c_vec_data]
dim_file=['unit_dim.ndx']
group_dim=['a_axis','b_axis','c_axis']
for ndx_file in dim_file:
    if os.path.isfile(ndx_file):
        os.remove(ndx_file)
    for iter,Dim_data in enumerate(Dim_data_vec):
        trj.ndx_writer(ndx_file,Dim_data,group_dim[iter])


# Ang_data_vec=[alpha_vec_data,beta_vec_data,gamma_vec_data]
# angle_file=['unit_angle.ndx']
# group_ang=['alpha_ang','beta_ang','gamma_ang']

# for ndx_file in angle_file:
#     if os.path.isfile(ndx_file):
#         os.remove(ndx_file)
#     for iter,Ang_data in enumerate(Ang_data_vec):
#         trj.ndx_writer(ndx_file,Ang_data,group_ang[iter])