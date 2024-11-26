'''
This module is used for making the index files for any arbitrary choice of atoms. 

'''

import traj_reader as trj
import CNC_class as cnc
import os
import numpy as np
import matplotlib.pyplot as plt
import pathlib

class index_generator:

    def __init__(self,filename,LAYERS):
        self.data = trj.gro_reader(filename)
        self.CNC_group = cnc.CNC_analys(self.data, layers=LAYERS)
        self.property_functions = {
            'hbond' : self.CNC_group.extract_hb_atoms,
            'tg' : self.CNC_group.confor_angles,
            'phi' : self.CNC_group.confor_angles,
            'unit_cell' : self.CNC_group.unit_cells,
            # 'unit_cell_angle' : self.CNC_group.unit_angles,
            'twist' : self.CNC_group.twist_angles,
        }

        self.index_functions = {
            'hbond' : self.generate_hbond_indices,
            'tg' : self.generate_tg_indices,
            'phi' : self.generate_phi_indices,
            'unit_cell' : self.generate_unit_cell_indices,
            'twist' : self.generate_twist_indices,
        }
    def extract_property_functions(self, index_types):
        """Prepare data selectively based on the required index types."""
        executed_preparations = set()
        for index_type in index_types:
            prep_func = self.property_functions.get(index_type)
            if prep_func and prep_func not in executed_preparations:
                prep_func()
                executed_preparations.add(prep_func)

    # def prepare_data(self, index_types):
    #     """Prepare data selectively based on the required index types."""
    #     executed_preparations = set()
    #     for index_type in index_types:
    #         prep_func = self.preparation_functions.get(index_type)
    #         if prep_func and prep_func not in executed_preparations:
    #             prep_func()
    #             executed_preparations.add(prep_func)


    def generate_indices(self, index_types , file_path):
        """Generate specified types of index files."""
        # self.prepare_data(index_types)
        self.extract_property_functions(index_types)
        for index_type in index_types:
            if index_type in self.index_functions:
                self.index_functions[index_type](file_path)

    # def generate_indices(self, index_type):
    #     self.index_functions[index_type](self.CNC_group)



    ########### INDEX MAKING FOR HBONDS ################
    def generate_hbond_indices(self, file_path):
        ''' 
        This file extracts the hydrogen bond indices for the CNC system and creates the ndx file. 
        The hydrogen bond indices are extracted from the CNC_class object. 
        Each block in the index file contains atoms corresponding to each chain within a layer.
        '''

        hb_key_vec = list(self.CNC_group.hbs.keys())
        hb_type_vec =['H2O_O6','H3O_O5' , 'H6O_O3']

        
        for hb_type , hb_key in zip(hb_type_vec,hb_key_vec):
            block_size = 3*len(self.CNC_group.resid_vec)
            ndx_file = file_path + hb_type+'_index.ndx'
            if os.path.isfile(ndx_file):
                os.remove(ndx_file)
            iter = -1
            for layer in self.CNC_group.layer_vec:
                # chain_number_vec = [self.layers[layer][0], self.layers[layer][-1]] if len(self.layers[layer]) > 1 else [self.layers[layer][0]] # for the exterior chains
                chain_number_vec = self.CNC_group.layers[layer][1:-1] # for the interior chains
                for chain_number_iter,chain_number in enumerate(chain_number_vec):
                    if hb_key=='(O6)H--O3':
                        if layer==self.CNC_group.layer_vec[0] or self.CNC_group==self.CNC_group.layer_vec[-1]:
                            continue
                        elif chain_number_iter==len(chain_number_vec)-1:
                            continue
                    iter+=1
                    # data_hb = _extract_hb_indices(CNC_group,hb_type,hb_key, layer, chain_number_vec, chain_number_iter)
                    data_hb = self.CNC_group.hbs[hb_key][iter*block_size:(iter+1)*block_size]
                    trj.ndx_writer(ndx_file,data_hb,"%s_ch%d_%s" % (layer,chain_number_iter,hb_type)) if data_hb else None


    ########### INDEX MAKING FOR tg conformations ################
    def generate_tg_indices(self , file_path):

        ndx_file = file_path +  'Chis.ndx'
        if os.path.isfile(ndx_file):
            os.remove(ndx_file)
        chi_dict={'left_Chis':[] , 'right_Chis' : [] , 'left_Chi_ps' : [] , 'right_Chi_ps' : [] }
        for confor in ['Chi','Chi_p']:
            chi_dict['left_%ss' % confor] =  np.array([self.CNC_group.confors_data[confor][i+4:i+8]\
                    for i in range(0, len(self.CNC_group.confors_data[confor]), 8)]).flatten()
            chi_dict['right_%ss' % confor] = np.array([self.CNC_group.confors_data[confor][i:i+4]\
                        for i in range(0, len(self.CNC_group.confors_data[confor]), 8)]).flatten()
            trj.ndx_writer(ndx_file,chi_dict['left_%ss' % confor],'left_%ss' % confor)
            trj.ndx_writer(ndx_file,chi_dict['right_%ss' % confor],'right_%ss' % confor)

    ########### INDEX MAKING FOR Phi and Psi conformations ################
    def generate_phi_indices(self , file_path):

        ndx_file = file_path + 'Confors.ndx'
        if os.path.isfile(ndx_file):
            os.remove(ndx_file)
        chi_dict={'Phi':[] , 'Psi' : []  }
        for confor in chi_dict:
            chi_dict[confor] =  self.CNC_group.confors_data[confor]
            trj.ndx_writer(ndx_file,chi_dict[confor], '%ss' % confor)

    ########### INDEX MAKING FOR unit cell dimensions ################
    def generate_unit_cell_indices(self , file_path):
        ndx_file = file_path + 'unit_dim.ndx'
        if os.path.isfile(ndx_file):
            os.remove(ndx_file)
        for dim_elem in self.CNC_group.unit_cell_dims.keys():
            trj.ndx_writer(ndx_file,self.CNC_group.unit_cell_dims[dim_elem], '%s_axis' % dim_elem)
        
        ndx_file = file_path + 'unit_angle.ndx'
        if os.path.isfile(ndx_file):
            os.remove(ndx_file)
        for unit_cell_angle in self.CNC_group.unit_cell_angles.keys():
            trj.ndx_writer(ndx_file,self.CNC_group.unit_cell_angles[unit_cell_angle], '%s_ang' % unit_cell_angle)


    ########### INDEX MAKING FOR twists ################
    def generate_twist_indices(self , file_path):
        ndx_file = file_path + 'twist_chain.ndx'
        if os.path.isfile(ndx_file):
            os.remove(ndx_file)
        for twist_elem in self.CNC_group.twist_data.keys():
            trj.ndx_writer(ndx_file,self.CNC_group.twist_data[twist_elem], '%ss' % twist_elem)


    # ########### INDEX MAKING FOR unit cell angles ################
    def selected_inds(O6_atoms,inds_name):
        ndx_file = '%s.ndx' % inds_name
        if os.path.isfile(ndx_file):
            os.remove(ndx_file)
        trj.ndx_writer(ndx_file,O6_atoms, inds_name)



# if __name__ == "__main__":
file_path = './simulation_traj_topol/Bu/'
filename = file_path+'solute.gro'
Modified_LAYERS = {
    'L1': [18], 'L2': [19, 2], 'L3': [14,24,3], 'L4':  [15,25,31,6],
    'L5': [11,21,26,32,7], 'L6': [12,22,27,33,36,10], 'L7': [13,23,28,34,8],
    'L8': [16,29,35,9], 'L9': [17,30,4], 'L10':  [20,5], 'L11': [1]
}
index_types = ['hbond','twist' , 'tg' , 'unit_cell' , 'phi']
index_generator = index_generator(filename , Modified_LAYERS)
index_generator.generate_indices(index_types , file_path)

