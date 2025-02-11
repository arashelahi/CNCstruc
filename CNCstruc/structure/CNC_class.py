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
        domain (str): Domain of the CNC structure (interior or exterior)
        layers (dict): Mapping of layer names to chain numbers, where L1 layer is the topmost layer and L11 is the bottommost layer.
        residue_numbers (int): Number of residues in each CNC chain

        Attributes:
        ----------
            DEFAULT_LAYERS (dict): Default layer mapping for interior and exterior chains.
            ATOM_TYPES (dict): The dictionary containing the atom names for different structural properties.
            descriptor (dict): The dictionary to store the atom numbers for different structural properties.
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
    
    ATOM_TYPES = { 'glycosidic' : {'Phi' : ('O5' , 'C1' ,'O4' , 'C4') , 'Psi' : ('C1' , 'O4' ,'C4' , 'C5')} ,
                   'alcohols'   : {'Chi' : ('O5' , 'C5' ,'C6' , 'O6') , 'Chi_p' : ('C4' , 'C5' ,'C6' , 'O6')} ,
                   'twist'      : {'twist' : ('C3' , 'C5' ,'C5' , 'C3')} ,
                   'H_bonds'    : {'O2H_O6' : ('O2','HO2','O6') , 'O3H_O5'  : ('O5','O3','HO3')  , 'O6H_O3' : ('O6','HO6','O3')},
                   'unit_cell'  : {'a'     : ('C1','C1')     , 'b'    : ('C1','C1')      , 'c'    : ('C1','C1'),
                                   'alpha' : ('C1','C1','C1'), 'beta' : ('C1','C1','C1'), 'gamma' : ('C1','C1','C1')}
                 }

    
    def __init__(self, data, domain='interior',layers = None, residue_numbers = 20 , FF = 'default'): 

        self.layers = layers if layers else self.DEFAULT_LAYERS # The layers are defined based on the spatial location of the chains.
        self.residue_numbers = residue_numbers                  # The number of residues in each CNC chain.
        self.data = data                                        # The data containing the atom numbers, atom names, residue numbers, and atom positions.
        self.domain = domain                                    # The domain of the CNC structure (interior or exterior)
        self.descriptor = { 'glycosidic' : {'Phi':  [] , 'Psi'   : []} ,
                            'alcohols'   : {'Chi' : [] , 'Chi_p' : []} ,
                            'twist'      : {'twist' : []} ,
                            'H_bonds'    : {'O2H_O6': [], 'O3H_O5': [], 'O6H_O3' : []} ,
                            'unit_cell'  : {'dimension' : {'a' : [], 'b' : [], 'c' : []} ,
                                            'angle'  : {'alpha' : [], 'beta' : [], 'gamma' : []}}
                                            }                   ## The main descriptor that stores the atom numbers for different structural properties.
        # self.descriptor_values = self.descriptor.copy()
        self.resid_vec = np.arange(5,17,1)
        self.layer_vec = list(self.layers.keys()) if self.domain == 'exterior' else list(self.layers.keys())[2:-2] ## This choice could be made when only iterior chains are needed.
        self.feature_dict = { 
                            k: { 
                            'glycosidic' : {'Phi'   :  []               , 'Psi'  :  []} ,
                            'alcohols'   : {'left'  :  []               ,'right' :  []} ,
                            'twist'      : {'twist' :  []                             } ,
                            'H_bonds'    : {'O2H_O6':  [], 'O3H_O5': [], 'O6H_O3' : []} ,
                            } 
                            for k in self.layer_vec
                            }

        self.ff = FF
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
        This module is used to clip the data based on the middle residues and the interior chains.

        Returns:
            DataFrame: Clipped DataFrame.
        """
        
        # Adding the chain number and residue number to the data
        refined_data = self._dataframe_preprocess()
        refined_data = refined_data[refined_data['residue_number'].isin(self.resid_vec)]
        chain_num_vec = [num for layer in self.layer_vec for num in self.layers[layer]]
        self.Clipped_Data=refined_data[refined_data['chain_number'].isin(chain_num_vec)]

        return self.Clipped_Data
    
    def get_indices(self,feature):
        if feature in ['glycosidic', 'alcohols', 'twist']:
            self.get_dihed_indices(feature)
        elif feature == 'H_bonds':
            self.get_hb_indices(feature)
        elif feature == 'unit_cell':
            self.unit_cell_dim(feature)
            self.unit_cell_ang(feature)
        else:
            raise ValueError('Invalid feature name. Please choose from the following: glycosidic, alcohols, twist, H_bonds, or unit_cell.')

    def lay_chain_resid(self, layer, feature):
        """
        Determines the chain_number_vec and resid_range based on the layer and feature.

        Parameters:
            layer (str): Current layer being processed.
            feature (str): The feature type ('glycosidic', 'alcohols', 'twist', 'a', 'b', 'c', 'alpha', 'beta', 'gamma').

        Returns:
            tuple: chain_number_vec (list), resid_range (list).
        """
        # Initialize the variables
        chain_number_vec = []
        resid_range = []

        # Handle unit cell dimensions and angles
        if feature in ['a', 'b', 'c', 'alpha', 'beta', 'gamma']:
            if feature == 'a':
                if layer < 'L7':
                    if layer in self.layer_vec[2:]:
                        chain_number_vec = self.layers[layer][2:-2]
                else:
                    chain_number_vec = self.layers[layer][1:-1]
                resid_range = self.resid_vec

            elif feature == 'b':
                if layer in self.layer_vec[1:-1]:
                    chain_number_vec = self.layers[layer][1:-2]
                resid_range = self.resid_vec

            elif feature == 'c':
                chain_number_vec = self.layers[layer][1:-1]
                resid_range = self.resid_vec[:-2]

            elif feature == 'alpha':
                if layer in self.layer_vec[1:-1]:
                    chain_number_vec = self.layers[layer][1:-2]
                resid_range = self.resid_vec[:-2]

            elif feature == 'beta':
                if layer in self.layer_vec[2:]:
                    if layer < 'L7':
                        chain_number_vec = self.layers[layer][2:-2]
                    else:
                        chain_number_vec = self.layers[layer][1:-1]
                resid_range = self.resid_vec[:-2]

            elif feature == 'gamma':
                if layer in self.layer_vec[2:-1]:
                    if layer < 'L7':
                        chain_number_vec = self.layers[layer][2:-2]
                    else:
                        chain_number_vec = self.layers[layer][2:-1]
                resid_range = self.resid_vec

        # Handle glycosidic, alcohols, and twist
        elif feature in ['glycosidic', 'alcohols', 'twist']:
            if self.ff == 'Charmm':
                self.resid_vec = self.resid_vec[::-1]
            chain_number_vec = self.layers[layer][1:-1]
            if feature == 'glycosidic':
                resid_range = self.resid_vec[:-1]
            elif feature == 'alcohols':
                resid_range = self.resid_vec
            elif feature == 'twist':
                resid_range = self.resid_vec[:-2]
            # if self.ff == 'Charmm':
            #     resid_range = resid_range[::-1]
        
        elif feature in ['O2H_O6', 'O3H_O5', 'O6H_O3']:
        # Extract chain numbers for interior chains
            
            if feature == 'O6H_O3':
                chain_number_vec = self.layers[layer][1:-2]
                # Skip calculation for first/last layer or if there's only one chain
                if layer == self.layer_vec[0] or layer == self.layer_vec[-1]:
                    chain_number_vec = []  # Skip first and last layers
                resid_range = self.resid_vec  # Use all residues for this hydrogen bond type
            else:
                chain_number_vec = self.layers[layer][1:-1]
                resid_range = self.resid_vec  # Default behavior for other hydrogen bond types
        return chain_number_vec, resid_range


    def get_dihed_indices(self , feature = 'glycosidic'):
        """"
        Module to extract the atom numbers for the structural properties described by dihedral angles.
        """
        for ang_type in self.descriptor[feature].keys():
            for layer in self.layer_vec:
                if chain_number_vec == []:
                    continue
                chain_number_vec , resid_range = self.lay_chain_resid(layer, feature)
                for chain_number in chain_number_vec: 
                    for resid in resid_range:
                        self.descriptor[feature][ang_type] += self._extract_dihed_atom_nums(self.clipped_data, chain_number, resid, ang_type , feature)
    def _extract_dihed_atom_nums(self, data , chain_number, resid, angle_type , feature):
        """
        helper function to extract the atom numbers for the structural properties described by dihedral angles.

        Parameters:
            data (DataFrame): DataFrame containing angle-specific data.
            chain_number (int): Chain number.
            resid (int): Residue number.
            angle_type (str): Name of angle type, e.g., Phi.
            feature (str): Name of the structural feature, e.g., glycosidic dihedral angle.
            resid_offset (int): Residue offset indicating the residue number involved in the dihedral angle.
        Returns:
            list: List of atom numbers.
        """
        atom_names = self.ATOM_TYPES[feature][angle_type]

        # if self.ff == 'Charmm':
        #     atom_names = atom_names[::-1]
        atoms = []

        for atom_iter , atom_name in enumerate(atom_names):
            if feature == 'glycosidic':
                resid_offset = 1 if atom_name in ['O4', 'C4', 'C5'] else 0
            elif feature == 'alcohols':
                resid_offset = 0
            elif feature == 'twist':
                resid_offset = 2 if atom_iter > 1 else 0
            if self.ff == 'Charmm':
                resid_offset = -resid_offset
            # resid_offset = 1 if angle_type in ['Phi', 'Psi'] and atom_name in ['O4', 'C4', 'C5'] else 0 ## The Phi and Psi angles contain the atoms from two residues, while Chis do not.
            atoms += data[(data['residue_number'] == resid + resid_offset) & 
                          (data['chain_number']   == chain_number) & 
                          (data['atom_name']      == atom_name)]['atom_number'].values.tolist()
        return atoms


    def get_hb_indices(self , feature):
        """
         Module to extract the atom numbers for the different types of hydrogen bonds.

        The method iterates over predefined hydrogen bond types and extracts    
        corresponding atom numbers from the clipped data.
        """

        for hb_type in self.ATOM_TYPES[feature].keys():
            for layer in self.layer_vec:
                chain_number_vec , resid_range = self.lay_chain_resid(layer, hb_type)
                if chain_number_vec == []:
                    continue
                for chain_number_iter,chain_number in enumerate(chain_number_vec):
                    for resid in self.resid_vec:
                        self.descriptor['H_bonds'][hb_type] += self._extract_hb_atom_nums(self.clipped_data, layer , chain_number_vec, chain_number_iter , resid, feature, hb_type)
    
    def _extract_hb_atom_nums(self, data, layer, chain_number_vec, chain_number_iter , resid, feature , hbtype):

        odd_glc_numbs = np.unique(self.clipped_data['residue_number'])[0::2]
        even_glc_numbs = np.unique(self.clipped_data['residue_number'])[1::2]

        if self.ff == 'Charmm':
            reverse_copy = odd_glc_numbs.copy()
            odd_glc_numbs = even_glc_numbs.copy()
            even_glc_numbs = reverse_copy
            # even_glc_numbs = odd_glc_numbs.copy()
            # odd_glc_numbs = np.unique(self.clipped_data['residue_number'])[1::2]
            # odd_glc_numbs = np.unique(self.clipped_data['residue_number'])[0::2]

        atom_names = self.ATOM_TYPES[feature][hbtype] 
        atoms = []

        for atom_iter , atom_name in enumerate(atom_names):
            if hbtype == 'O6H_O3': 
                if resid in odd_glc_numbs:
                    if atom_name in ['O6','HO6']:
                        chain_number = chain_number_vec[chain_number_iter]
                    else:
                        # chain_number = chain_number_vec[chain_number_iter+1]
                        chain_number = self.layers[layer][1:][chain_number_iter+1]
                if resid in even_glc_numbs:
                    if atom_name == 'O3':
                        chain_number = chain_number_vec[chain_number_iter]
                    else:
                        # chain_number = chain_number_vec[chain_number_iter+1]
                        chain_number = self.layers[layer][1:][chain_number_iter+1]

            else:
                chain_number = chain_number_vec[chain_number_iter]
            atoms += data[(data['chain_number']   == chain_number) & 
                          (data['residue_number']   == resid) & 
                          (data['atom_name']      == atom_name)]['atom_number'].values.tolist()
        return atoms


    def unit_cell_dim(self , feature):
        """
       This module is used to extract the atom numbers for the unit cell dimensions.
        """
        for dim_type in self.descriptor['unit_cell']['dimension'].keys():
            """ for the 'a' unit dimension, the first two layers are skipped. Also, the ending on the top-half (<L7) do not have the pair chains for a calculation """

            for layer_iter_number , layer in enumerate(self.layer_vec):
                chain_number_vec , resid_range = self.lay_chain_resid(layer, dim_type)
                if chain_number_vec == []:
                    continue   
                for chain_number_iter , chain_number in enumerate(chain_number_vec):

                    for resid in resid_range:
                        self.descriptor['unit_cell']['dimension'][dim_type]+=self._extract_unit_dim_atoms(self.clipped_data, layer_iter_number, layer, chain_number_vec, chain_number_iter, resid, feature, dim_type)

    def _extract_unit_dim_atoms(self , data, layer_iter_number , layer,  chain_number_vec, chain_number_iter , resid, feature , dim_type):
        """
        Helper function to extract the atom numbers for the unit cell dimensions.
        """
        atom_names = self.ATOM_TYPES[feature][dim_type] 
        atoms = []
        chain_number =  chain_number_vec[chain_number_iter]
        resid_offset = 0
        for atom_iter , atom_name in enumerate(atom_names):
            if atom_iter == 1:
                if dim_type == 'c':
                    resid_offset = 2
                elif dim_type == 'b' : 
                    chain_number = self.layers[layer][1:-1][chain_number_iter+1]
                elif dim_type == 'a' : 
                    top_iter = 0 if layer <= 'L7' else +1
                    chain_number = self.layers[self.layer_vec[layer_iter_number-2]][1:-1][top_iter+chain_number_iter]
            atoms += data[(data['chain_number']   == chain_number) & 
                          (data['residue_number'] == resid + resid_offset) & 
                          (data['atom_name']      == atom_name)]['atom_number'].values.tolist()
        return atoms

    def unit_cell_ang(self , feature):
        """
        This module is used to extract the atom numbers for the unit cell angles.
        """
        resid_range = self.resid_vec[:-2]
        for ang_type in self.descriptor['unit_cell']['angle'].keys():
            # if ang_type == 'gamma':
            #     resid_range = self.resid_vec
            """ for the 'gamma' angle, the first two layers are skipped. Also, the ending on the top-half (<L7) do not have the pair chains for a calculation """
            for layer_iter_number , layer in enumerate(self.layer_vec):
                chain_number_vec , resid_range = self.lay_chain_resid(layer, ang_type)
                if chain_number_vec == []:
                    continue   
                for chain_number_iter , chain_number in enumerate (chain_number_vec):

                    for resid in resid_range:
                       self.descriptor['unit_cell']['angle'][ang_type]+=self._extract_unit_angle_atoms(self.clipped_data, layer_iter_number, layer, chain_number_vec, chain_number_iter, resid, feature, ang_type)


    def _extract_unit_angle_atoms(self , data, layer_iter_number , layer,  chain_number_vec, chain_number_iter , resid, feature, ang_type):
        """
        Helper function to extract the atom numbers for the unit cell angles.
        """

        
        atom_names = self.ATOM_TYPES[feature][ang_type] 
        atoms = []
        for atom_iter , atom_name in enumerate(atom_names):
            resid_offset = 0
            chain_number =  chain_number_vec[chain_number_iter]
            if atom_iter == 0 :

                if ang_type == 'alpha':
                    resid_offset = 2

                if ang_type in ['beta' , 'gamma']:
                    if ang_type == 'gamma':
                        top_iter = 0 if layer < 'L7' else +1 if layer == 'L7' else +2
                    elif ang_type == 'beta':
                        top_iter = 0 if layer <= 'L7' else +1
                    chain_number = self.layers[self.layer_vec[layer_iter_number-2]][1:-1][top_iter+chain_number_iter]

            if atom_iter == 2:

                if ang_type =='alpha':
                    chain_number = self.layers[layer][1:-1][chain_number_iter+1]

                if ang_type == 'beta':
                    resid_offset =  2
                
                if ang_type == 'gamma':
                    chain_number = self.layers[layer][1:-1][chain_number_iter] # this is one chain to the left of the current chain examined at atom_iter == 1
            atoms += data[(data['chain_number']   == chain_number) & 
                          (data['residue_number'] == resid + resid_offset) & 
                          (data['atom_name']      == atom_name)]['atom_number'].values.tolist()
        return atoms
