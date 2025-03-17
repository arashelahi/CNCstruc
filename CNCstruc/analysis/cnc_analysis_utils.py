
from CNCstruc.utils import traj_reader as trj

from CNCstruc.structure import CNC_class as CNC

import numpy as np
import matplotlib.pyplot as plt
from statistics import mode


"""
This modules uses the gromacs generated xvg files, and outputs residue per time-step behavior for each of the properties.
chain naming convention:
For the interior chains, the labeling starts from ch0 from each layer, and iteratively increases.
For the exterior chians, the labels are ch0, and ch1, only
"""

def feature_analysis (CNCclass ,feature, FF_direcotory , domain='interior'):
    """

    """
    filepath = FF_direcotory + 'Results/'
    if feature in ['glycosidic', 'alcohols', 'twist']:
        dihed_analysis (CNCclass , feature , filepath)
    if feature in ['O2H_O6','O3H_O5','O6H_O3']:
        hb_analysis (CNCclass ,  feature , filepath)
    if feature == 'unit_cell':
        CNCclass.feature_dict['unit_cell'] = {
            'dimension': {'a': [], 'b': [], 'c': []},
            'angle': {'alpha': [], 'beta': [], 'gamma': []},
        }
        unit_cell_analysis (CNCclass , feature , filepath)

def dihed_analysis (CNCclass , feature , filepath):
    iter = -1 ## counting the chain numbers.
    for layer in CNCclass.layer_vec:
        chain_number_vec , resid_vec = CNCclass.lay_chain_resid(layer, feature) ## The class instance is needed to get the chain numbers
        if feature == 'alcohols':
            resid_vec = resid_vec[::2]
        if chain_number_vec == []:
            continue
        ind_intervs = len(resid_vec)
    
        for chain_num_iter, chain_number in enumerate(chain_number_vec):
            iter += 1
            for subfeature in CNCclass.feature_dict[layer][feature].keys():
                subfeature_values = []
                for resid_iter , resid in enumerate(resid_vec):
                    y_index = (iter * ind_intervs + resid_iter + 2) + 1
                    
                    subfeature_values += dihed_get_values(feature , subfeature , filepath , y_index)
            
                CNCclass.feature_dict[layer][feature][subfeature].append(subfeature_values)

def dihed_get_values(feature , subfeature , filepath,y_index):
    if feature == 'alcohols':
        subfeature_values = alcohol_get_values(feature , subfeature , filepath , y_index)
    # if feature == 'glycosidic':
    else:
        subfeature_values = glyc_twist_get_values(feature , subfeature , filepath , y_index)
    return subfeature_values

def alcohol_get_values(feature , subfeature , filepath , y_index , trans_ang = 150):
    """"
    Outputs the confor_vec, containing tg conformatino over time per residue.
    """

    # if feature 
    chi_file = filepath + f'Chi_{subfeature}_dist.xvg'
    chi_p_file = filepath +  f'Chi_p_{subfeature}_dist.xvg'

    for ang_file in [chi_file,chi_p_file]:
        _ , angle_vec = trj.xvg_reader(ang_file,[1,y_index])
        if ang_file == chi_file: confor_vec = [''] * len(angle_vec)
        for ang_iter , angle_val in enumerate(angle_vec):
            confor_vec[ang_iter] +='t' if abs(angle_val) >= trans_ang else 'g'
            
    return confor_vec
def glyc_twist_get_values(feature , subfeature , filepath , y_index): 
    """"
    Outputs the confor_vec, containing tg conformatino over time per residue.
    """
    
    glyc_file = filepath + f'{subfeature}_dist.xvg'
    _ , glyc_val_vec = trj.xvg_reader(glyc_file,[1,y_index])
    return glyc_val_vec


def _read_hbs(hb_filename,max_hbnum):
    """"
    Outputs the confor_vec, containing tg conformatino over time per residue.
    """
    indice=[1,2]
    _,hb_data=trj.xvg_reader(hb_filename,indice)
    return (np.array(hb_data)/max_hbnum*100).tolist()
def hb_analysis (CNCclass , hb_type  , HB_dir, max_hbnum = 11 , domain='interior' , file_ext = ''):
    """
    Analyzes the hb conformation from xvg files.

    Parameters:
        CNCclass (cnc.CNC_analys): Instance of CNC_analys class.
        hbtype (str): type of hb
        max_hbnum (int): Max number of hb for each chain 
        
    Returns:
        dict: the vector of hb occupancy for each chain
    """
    if hb_type == 'O6H_O3':
        max_hbnum=12 
    for layer in CNCclass.layer_vec:
        
        chain_number_vec , resid_vec = CNCclass.lay_chain_resid(layer, hb_type) ## The class instance is needed to get the chain numbers
        if chain_number_vec == []:
            continue

        for chain_iter , chain_number in enumerate(chain_number_vec):
            filename = HB_dir + '%s_ch%s_%s%s_hbnum.xvg' % (layer , chain_iter, hb_type , file_ext)
            hb_data = _read_hbs(filename,max_hbnum)
            CNCclass.feature_dict[layer]['H_bonds'][hb_type].append(hb_data)

def unit_cell_analysis (CNCclass , feature , filepath):

    for subfeature in CNCclass.feature_dict[feature].keys():
        for unit_cell_prop in CNCclass.feature_dict[feature][subfeature]:
            if subfeature == 'dimension':
                unit_cell_file = filepath + 'unit_dist.xvg'
                if unit_cell_prop == 'a':
                    y_index = 2
                elif unit_cell_prop == 'b':
                    y_index = 3
                elif unit_cell_prop == 'c':
                    y_index = 4

            if subfeature == 'angle':
                unit_cell_file = filepath + f'{unit_cell_prop}_dist.xvg'
                y_index = 2
            _ , subfeature_values = trj.xvg_reader(unit_cell_file,[1,y_index])
            # _get_unit_cell_values(feature , subfeature , unit_cell_file , y_index)
            CNCclass.feature_dict[feature][subfeature][unit_cell_prop] = subfeature_values
