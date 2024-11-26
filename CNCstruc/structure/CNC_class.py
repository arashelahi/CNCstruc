'''
This Moldule is created to draw the crystalline properties of the cellulose nanocrystals. The data contains the atom numbers involing in
the hydrogen bond network, unit celll parameters, the phi and psi angles, and the tg conformations. 

'''
from dataclasses import dataclass
from fileinput import filename
# import traj_reader as trj
import os
import numpy as np
import pandas as pd

class CNC_analys:
    """
        Identification of the atoms involved in the hydrogen bond network, primary alcohol conformations, unit cell parameters, and twist angles.

        Parameters:
        ----------
        data (DataFrame): DataFrame containing the atom numbers, atom names, residue numbers, and atom positions

        Attributes:
            layers (dict): Mapping of layer names to chain numbers.
            atom_num (int): Number of atoms per chain.
            data (DataFrame): DataFrame containing CNC structure data.
            clipped_data (DataFrame): Data clipped based on specific criteria.
            confors (dict): Dictionary storing calculated conformational dihedral angles.
    """
    
    DEFAULT_LAYERS = {
        'L1': [11], 'L2': [12, 22], 'L3': [4, 13, 23], 'L4': [5, 14, 24, 31],
        'L5': [1, 6, 15, 25, 32], 'L6': [2, 7, 16, 26, 33, 36], 'L7': [3, 8, 17, 27, 34],
        'L8': [9, 18, 28, 35], 'L9': [10, 19, 29], 'L10': [20, 30], 'L11': [21]
    }
    # {
    #     'L1': [18], 'L2': [19, 2], 'L3': [14,24,3], 'L4':  [15,25,31,6],
    #     'L5': [11,21,26,32,7], 'L6': [12,22,27,33,36,10], 'L7': [13,23,28,34,8],
    #     'L8': [16,29,35,9], 'L9': [17,30,4], 'L10':  [20,5], 'L11': [1]
    # }
    CONFOR_TYPES={'Phi' : ('O5', 'C1', 'O4', 'C4') , 'Psi' : ('C1', 'O4', 'C4','C5'),
                  'Chi' : ('O5' ,'C5' ,'C6', 'O6') , 'Chi_p' : ('C4' ,'C5' ,'C6', 'O6') 

}
    TWIST_TYPES={'chain_twist' : ('C3', 'C5','C5','C3')  
}
    
    def __init__(self, data, domain='interior',layers=None, residue_numbers = 20): 
        self.layers = layers if layers else self.DEFAULT_LAYERS # The layers are defined based on the spatial location of the chains.
        # self.atom_num = atom_num
        self.residue_numbers = residue_numbers # The number of residues in each CNC chain.
        self.data = data                       # The data containing the atom numbers, atom names, residue numbers, and atom positions.
        self.domain = domain                    # The domain of the CNC structure (interior or exterior)
        self.confors_data = {'Phi': [], 'Psi': [], 'Chi' : [] , 'Chi_p' : []}
        self.twist_data = {'chain_twist' : []}
        self.hbs = {'(O2)H--O6': [], '(O3)H--O5': [], '(O6)H--O3' : []}
        self.unit_cell_dims={'a' : [], 'b' : [], 'c' : []}
        self.unit_cell_angles={'alpha' : [], 'beta' : [], 'gamma' : []}
        self.resid_vec = np.arange(5,17,1)
        # self.layer_vec = ['L%d' % x for x in range(3,10)]
        self.layer_vec = list(self.layers.keys()) if self.domain == 'exterior' else list(self.layers.keys())[2:-2] ## This choice could be made when only iterior chains are needed.
        # self.layer_vec = list(CNC_group.layers.keys())[2:-2]
        self.clipped_data = self._clipping()
    
    def _dataframe_preprocess(self):
        """ This module is used to add the chain number, adjust the residue number,
         and clip the data based on the middle residues and the interior chains. """

        data = self.data.copy()
        # Identify where new residues start ('C1' occurrences)
        residue_starts = data['atom_name'] == 'C1'

        # Calculate residue and chain numbers
        data['residue_number'] = ((residue_starts.cumsum() -1) % self.residue_numbers)+1
        
        data['chain_number'] = ((residue_starts.cumsum() -1) // self.residue_numbers)+1

        return data
        

    def _clipping(self):
        """
        Clips the data based on layer specifications.

        Returns:
            DataFrame: Clipped DataFrame.
        """
        # added_chain_data = self._add_chain_num()
        
        refined_data = self._dataframe_preprocess()
        refined_data = refined_data[refined_data['residue_number'].isin(self.resid_vec)]
        chain_num_vec = [num for layer in self.layer_vec for num in self.layers[layer]]
        self.Clipped_Data=refined_data[refined_data['chain_number'].isin(chain_num_vec)]
        return self.Clipped_Data


    def confor_angles(self):
        """
        Obtain the atom numbers for Chi, and Chi_P angles.
        """

        for ang_type, atom_names in self.CONFOR_TYPES.items():
            data_ang=self.clipped_data[self.clipped_data['atom_name'].isin(atom_names)]
            for layer in self.layer_vec:
                chain_number_vec = self.layers[layer][1:-1] # for the interior chains
                # chain_number_vec = [self.layers[layer][0], self.layers[layer][-1]] if len(self.layers[layer]) > 1 else [self.layers[layer][0]] # for the exterior chains
                for chain_number in chain_number_vec: 
                    resid_range = self.resid_vec[:-1] if ang_type in ['Phi' , 'Psi'] else self.resid_vec
                    for resid in resid_range:
                        self.confors_data[ang_type]+=self._extract_ang_atoms(data_ang, chain_number, resid, ang_type)

    def twist_angles(self):
        """
        Obtain the atom numbers for twist angles.
        """

        for ang_type, atom_names in self.TWIST_TYPES.items():
            data_ang=self.clipped_data[self.clipped_data['atom_name'].isin(atom_names)]
            for layer in self.layer_vec:
                chain_number_vec = self.layers[layer][1:-1] if self.domain == 'interior' else [self.layers[layer][0], self.layers[layer][-1]] if len(self.layers[layer]) > 1 else [self.layers[layer][0]] # for the interior chains

                # chain_number_vec = self.layers[layer][1:-1] # for the interior chains
                # # chain_number_vec = [self.layers[layer][0], self.layers[layer][-1]] if len(self.layers[layer]) > 1 else [self.layers[layer][0]] # for the exterior chains
                for chain_number in chain_number_vec: 
                    resid_range = self.resid_vec[:-2]
                    for resid in resid_range:
                        self.twist_data[ang_type]+=self._extract_twist_atoms(data_ang, chain_number, resid, ang_type)


    def _extract_ang_atoms(self, data, chain_number, resid, angle_type):
        """
        Extracts atom numbers for angle calculation.

        Parameters:
            data (DataFrame): DataFrame containing angle-specific data.
            chain_number (int): Chain number.
            resid (int): Residue number.
            angle_type (str): Type of angle ('Phi', 'Psi', 'Chi', or 'Chi_p').

        Returns:

        Returns:
            list: List of atom numbers.
        """
        atom_names = self.CONFOR_TYPES[angle_type]
        atoms = []

        for atom_name in atom_names:
            resid_offset = 1 if angle_type in ['Phi', 'Psi'] and atom_name in ['O4', 'C4', 'C5'] else 0 ## The Phi and Psi angles contain the atoms from two residues, while Chis do not.
            atoms += data[(data['residue_number'] == resid + resid_offset) & 
                        (data['chain_number'] == chain_number) & 
                        (data['atom_name'] == atom_name)]['atom_number'].values.tolist()
        return atoms

    def _extract_twist_atoms(self, data, chain_number, resid, angle_type):
        """
        Extracts atom numbers for angle calculation.

        Parameters:
            data (DataFrame): DataFrame containing angle-specific data.
            chain_number (int): Chain number.
            resid (int): Residue number.
            angle_type (str): Type of angle ('chain_twist').

        Returns:

        Returns:
            list: List of atom numbers.
        """
        atom_names = self.TWIST_TYPES[angle_type]
        atoms = []

        for atom_iter , atom_name in enumerate(atom_names):
            curr_resid = resid if atom_iter < 2 else resid +2
            atoms += data[(data['residue_number'] == curr_resid) & 
                        (data['chain_number'] == chain_number) & 
                        (data['atom_name'] == atom_name)]['atom_number'].values.tolist()
        return atoms

    def extract_hb_atoms(self):
        """
        Extracts the atom numbers that participate in hydrogen bonds.

        The method iterates over predefined hydrogen bond types and extracts    
        corresponding atom numbers from the clipped data.
        """
        hb_type_vec=[['O2','HO2','O6'],['O3','HO3','O5'],['O6','HO6','O3']]
        for hb_type, hb_key in zip(hb_type_vec, list(self.hbs.keys())):
            data_hb=self.clipped_data[self.clipped_data['atom_name'].isin(hb_type)]
            for layer in self.layer_vec:
                # chain_number_vec = [self.layers[layer][0], self.layers[layer][-1]] if len(self.layers[layer]) > 1 else [self.layers[layer][0]] # for the exterior chains
                chain_number_vec = self.layers[layer][1:-1] # for the interior chains
                for chain_number_iter,chain_number in enumerate(chain_number_vec):
                    if hb_key=='(O6)H--O3':
                        if layer==self.layer_vec[0] or layer==self.layer_vec[-1]:
                            continue
                        elif chain_number_iter==len(chain_number_vec)-1:
                            continue
                        odd_glc_numbs = np.unique(self.clipped_data['residue_number'])[0::2]
                        even_glc_numbs = np.unique(self.clipped_data['residue_number'])[1::2]
                        hb_to_add=data_hb[(data_hb['chain_number'] == chain_number_vec[chain_number_iter]) & 
                                          (data_hb['residue_number'].isin(odd_glc_numbs)) &
                                          (data_hb['atom_name'].isin(['O6','HO6']))]['atom_number'].values.tolist()+\
                                  data_hb[(data_hb['chain_number'] == chain_number_vec[chain_number_iter]) & 
                                          (data_hb['residue_number'].isin(even_glc_numbs)) &
                                          (data_hb['atom_name'].isin(['O3']))]['atom_number'].values.tolist()+ \
                                  data_hb[(data_hb['chain_number'] == chain_number_vec[chain_number_iter+1]) & 
                                          (data_hb['residue_number'].isin(odd_glc_numbs)) &
                                          (data_hb['atom_name'].isin(['O3']))]['atom_number'].values.tolist()+\
                                  data_hb[(data_hb['chain_number'] == chain_number_vec[chain_number_iter+1]) & 
                                          (data_hb['residue_number'].isin(even_glc_numbs)) &
                                          (data_hb['atom_name'].isin(['O6','HO6']))]['atom_number'].values.tolist()
                                #       data_hb[(data['chain_number'] == chain_number_vec[chain_number_iter]) & 
                                #               (data['residue_number'].isin([6,8,10,12,14])]['atom_number'].values.tolist()+
                                #   (data['atom_name'].isin(['O6','HO6']))]['atom_number'].values.tolist()
                                
                    else:
                        hb_to_add=data_hb[(data_hb['chain_number']==chain_number)]['atom_number'].values.tolist()
                    self.hbs[hb_key]+=hb_to_add

    def unit_cells(self):
        self._unit_cell_dim()
        self._unit_cell_ang()

    def _unit_cell_dim(self):
        """
        Calculate the unit cells.
        """
        data_unit_dim = self.clipped_data[self.clipped_data['atom_name']=='C1']
        for dim_type in self.unit_cell_dims.keys():
            """ for the 'a' unit dimension, the first two layers are skipped. Also, the ending on the top-half (<L7) do not have the pair chains for a calculation """
            layer_vec = self.layer_vec[2:] if dim_type == 'a' else self.layer_vec[1:-1] if dim_type == 'b' else self.layer_vec
            for layer in layer_vec:
                chain_vec = self.layers[layer][2:-2] if dim_type == 'a' and layer < 'L7' else self.layers[layer][1:-2] if dim_type == 'b' else self.layers[layer][1:-1]
                for chain_number in chain_vec:
                    resid_range = self.resid_vec[:-2] if dim_type == 'c' else self.resid_vec
                    for resid in resid_range:
                        self.unit_cell_dims[dim_type]+=self._extract_unit_dim_atoms(data_unit_dim, layer, chain_number, resid, dim_type)

    def _extract_unit_dim_atoms(self , data, layer, chain_number, resid, dim_type):

        chain_ind = self.layers[layer][1:-1].index(chain_number)
        layer_ind = self.layer_vec.index(layer)


        if dim_type == 'c':
            atoms = data[(data['residue_number'].isin([resid,resid+2])) & 
                         (data['chain_number']==chain_number)]['atom_number'].values.tolist()

        if dim_type == 'b' : 
            atoms = data[(data['residue_number']==resid) &
                         (data['chain_number'].isin(self.layers[layer][1:-1][chain_ind:chain_ind+2]))]['atom_number'].values.tolist()
        if dim_type == 'a' : 

            top_iter = -1 if layer < 'L7' else 0 if layer == 'L7' else +1

            other_chain=self.layers[self.layer_vec[layer_ind-2]][1:-1][top_iter+chain_ind]
            atoms = data[(data['residue_number']==resid) &
                         (data['chain_number'].isin([chain_number,other_chain]))]['atom_number'].values.tolist()

        return atoms

    def _unit_cell_ang(self):
        """
        Calculate the unit cells.
        """
        data_unit_ang = self.clipped_data[self.clipped_data['atom_name']=='C1']
        for ang_type in self.unit_cell_angles.keys():
            """ for the 'gamma' angle, the first two layers are skipped. Also, the ending on the top-half (<L7) do not have the pair chains for a calculation """
            layer_vec = self.layer_vec[2:-1] if ang_type == 'gamma' else self.layer_vec[1:-1] if ang_type == 'alpha' else self.layer_vec[2:]
            for layer in layer_vec:
                chain_vec = self.layers[layer][2:-2] if ang_type == 'gamma' and layer < 'L7'\
                       else self.layers[layer][2:-1] if ang_type == 'gamma' and layer >= 'L7'\
                       else self.layers[layer][2:-2] if ang_type == 'beta'  and layer < 'L7'\
                       else self.layers[layer][1:-1] if ang_type == 'beta'  and layer > 'L6'\
                       else self.layers[layer][1:-2] #  if ang_type == 'gamma' and layer > 'L6' or ang_type = 'alpha'
                for chain_number in chain_vec:
                    resid_range = self.resid_vec if ang_type == 'gamma' else  self.resid_vec[:-2]
                    for resid in resid_range:
                        self.unit_cell_angles[ang_type]+=self._extract_unit_angle_atoms(data_unit_ang, layer, chain_number, resid, ang_type)


    def _extract_unit_angle_atoms(self,data, layer, chain_number, resid, ang_type):


        chain_ind = self.layers[layer][1:-1].index(chain_number)
        layer_ind = self.layer_vec.index(layer)

        if ang_type == 'alpha':
            atoms = data[(data['residue_number']==resid+2) & 
                         (data['chain_number']==chain_number)]['atom_number'].values.tolist()+\
                    data[(data['residue_number']==resid) & 
                         (data['chain_number']==chain_number)]['atom_number'].values.tolist()+\
                    data[(data['residue_number']==resid) &
                         (data['chain_number']==self.layers[layer][1:-1][chain_ind+1])]['atom_number'].values.tolist()

        if ang_type == 'beta' : 
            top_iter = -1 if layer < 'L7' else 0 if layer == 'L7' else +1

            other_chain = self.layers[self.layer_vec[layer_ind-2]][1:-1][top_iter+chain_ind]

            atoms = data[(data['residue_number']==resid) &
                        (data['chain_number']==other_chain)]['atom_number'].values.tolist()+\
                    data[(data['residue_number'].isin([resid,resid+2]))&
                        (data['chain_number']==chain_number)]['atom_number'].values.tolist()

        if ang_type == 'gamma' : 

            top_iter = -1 if layer < 'L7' else 0 if layer == 'L7' else 1
            # if layer>='L7':
            #     print('here')

            other_chain=self.layers[self.layer_vec[layer_ind-2]][1:-1][top_iter+chain_ind]
            atoms = data[(data['residue_number']==resid) &
                        (data['chain_number'].isin([other_chain,chain_number]))]['atom_number'].values.tolist()+\
                    data[(data['residue_number']==resid) &
                         (data['chain_number']==self.layers[layer][1:-1][chain_ind-1])]['atom_number'].values.tolist()

        return atoms

if __name__ == "__main__":
    file_dir = './simulation_traj_topol/'
    filename = file_dir+'solute.gro'
    Data = trj.gro_reader(filename)
    domain = 'interior'

    # Modif_CNC_layers = {
    #     'L1': [18], 'L2': [19, 2], 'L3': [14,24,3], 'L4':  [15,25,31,6],
    #     'L5': [11,21,26,32,7], 'L6': [12,22,27,33,36,10], 'L7': [13,23,28,34,8],
    #     'L8': [16,29,35,9], 'L9': [17,30,4], 'L10':  [20,5], 'L11': [1]
    # }
    CNC_group = CNC_analys(Data, domain=domain)
    # CNC_group.layer_vec = list(CNC_group.layers.keys()) # 
    # CNC_group.layer_vec = list(CNC_group.layers.keys())[2:-2] ## This choice could be made when only iterior chains are needed. 
    CNC_group.twist_angles()
    # print

    ndx_file='twist_chain.ndx'
    if os.path.isfile(ndx_file):
        os.remove(ndx_file)

    for chi_key in ['chain_twist']:
        # confor_per_side = {}
        # confor_per_side['right'] = [right_data for iter,right_data in enumerate(CNC_group.confors_data[chi_key]) if (iter)%8 in (0,1,2,3)]
        # confor_per_side['left']  = [left_data for iter,left_data in enumerate(CNC_group.confors_data[chi_key]) if (iter)%8 in (4,5,6,7)]
        # for side in ['left','right']:        
        trj.ndx_writer(ndx_file,CNC_group.twist_data[chi_key] ,chi_key+'s')