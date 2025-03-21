from CNCstruc.structure import CNC_class as CNC
from CNCstruc.utils import traj_reader as trj, Indexing as indx_gen
# from CNCstruc.analysis import Indexing as indx_gen
import pandas as pd
import os 
import numpy as np
# from rdkit import Chem
# from rdkit.Chem import AllChem


def functionalization(parent_molecule , reference_atom , atoms_to_delete , resid_range , func_group_init):
    '''
    This function converts the reactive groups to the modifying functional group:
    input: 
    parent_molecule : the datafile for the parent molecule
    reference_atom : the reference atom where the functional group and the reactive groups meet
    atoms_to_delete : the atom names which should be deleted from the functional group.
    resid_range : the residues where the functional groups will be applied upon
    func_group_init : the functional group dataframe
    # Outputs:
    the resulting data file
    '''
    for resid in resid_range:
        # The translational movement of the functional group replacement
        relative_position = parent_molecule[(parent_molecule['residue_number'] == resid) & 
                                            (parent_molecule['atom_name'] == reference_atom)].loc[:,'x':'z'].values - func_group_init[func_group_init['atom_name'] == reference_atom].loc[:,'x':'z'].values
        # remove the merging atoms from the functional group
        func_group = func_group_init[func_group_init['atom_name'] != reference_atom]
        # the movement of the functional group
        func_group.loc[:,'x':'z'] = func_group_init.loc[:,'x':'z'] + relative_position
        # making the functional group properties aligned with the parent molecule
        func_group[['residue_number', 'residue_name' , 'chain_number']] = parent_molecule[parent_molecule['residue_number'] == resid][['residue_number', 'residue_name','chain_number']].reset_index(drop=True)
        # delete the atoms which are converted
        data_to_delete = parent_molecule[(parent_molecule['residue_number'] == resid) & parent_molecule['atom_name'].isin(atoms_to_delete)].index
        parent_molecule = parent_molecule.drop(data_to_delete) # delete the atoms of the alcohol groups
        # The reordering is based on atom number. Shift the atom numbers of partciles after the removal as well as the functional group
        parent_molecule.loc[data_to_delete[0] + 1 :,'atom_number'] = np.arange(1 , len(parent_molecule.loc[data_to_delete[0] + 1 :,'atom_number'])+1 )\
                                                                        + parent_molecule.loc[data_to_delete[0] - 1,'atom_number'] + len(func_group)  # then the next atom will be shifted by the number of replacing atoms
        func_group.loc[:,'atom_number'] = np.arange(1 , len(func_group)+1 ) + parent_molecule.loc[data_to_delete[0] - 1 ,'atom_number']
        # Concat the functional groups to the cellulose
        parent_molecule = pd.concat([parent_molecule , func_group])
        # Reorder based on the atom number
        parent_molecule = parent_molecule.sort_values(by = 'atom_number' , ascending = True).reset_index(drop = True)
    return parent_molecule

def alkylation (parent_molecule , func):
    side_vec = ['left' , 'right']
    for side in side_vec:
        CNC_COOH_data = carboxylation (parent_molecule , side)
        # The functional group data and the residues to replace
        func_group_init , resid_range = _data_preparation(CNC_COOH_data , func , side)
        # The atom to be replaced and the atoms to be deleted
        reference_atom = 'C6'
        atoms_to_delete = ('O','H')
        # The carboxylated cellulose
        alkylated_CNC = functionalization(CNC_COOH_data , reference_atom , atoms_to_delete , resid_range , func_group_init)
        trj.gro_writer(f'CNC_COOH_{func}_{side}.gro' , alkylated_CNC)

def _data_preparation(parent_molecule , func , side):
    functional_file = f'./data/input/Groups/{func}_{side}.pdb'
    functional_data = trj.pdb_reader(functional_file)
    if side == 'right':
        resid_range = np.arange(parent_molecule['residue_number'].min() , parent_molecule['residue_number'].max() , 2)
    else:
        resid_range = np.arange(parent_molecule['residue_number'].min() + 1, parent_molecule['residue_number'].max()+1 , 2)
    return functional_data , resid_range
def carboxylation (parent_molecule , side):
    func = 'COOH'
    # side_vec = ['left' , 'right']
    # for side in side_vec:
        # The functional group data and the residues to replace
    func_group_init , resid_range = _data_preparation(parent_molecule , func , side)
    # The atom to be replaced and the atoms to be deleted
    reference_atom = 'C6'
    atoms_to_delete = ('H61', 'H62', 'O6' , 'HO6')
    # The carboxylated cellulose
    carboxylated_CNC = functionalization(parent_molecule , reference_atom , atoms_to_delete , resid_range , func_group_init)
    trj.gro_writer(f'CNC_COOH_{side}.gro' , carboxylated_CNC)
    return carboxylated_CNC


## initialization of the CNC object
CNC_form = 'Pristine'
filepath= ''
filename = f'./data/input/{CNC_form}/solute.gro'
Data_raw = trj.gro_reader(filename)
CNC_group = CNC.CNC_analys(Data_raw)
Data = CNC_group.data  # the chain number is added to the data  
one_chain = Data[Data['chain_number'] == 1]
# CNC_COOH_data  = carboxylation (one_chain , 'left')
## alkylation of the cellulose
alkyl_group = 'ethyl'
CNC_COOH_alkyl = alkylation (one_chain , alkyl_group)
# trj.gro_writer('CNC_COOH.gro' , CNC_COOH_data)

