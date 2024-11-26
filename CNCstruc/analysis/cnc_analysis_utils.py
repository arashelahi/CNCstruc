
from CNCstruc.utils import traj_reader as trj
from CNCstruc.structure import CNC_class as cnc
import numpy as np
import matplotlib.pyplot as plt
from statistics import mode


"""
This modules uses the gromacs generated xvg files, and outputs residue per time-step behavior for each of the properties.


chain naming convention:
For the interior chains, the labeling starts from ch0 from each layer, and iteratively increases.
For the exterior chians, the labels are ch0, and ch1, only

"""





def calculate_conformation(chi_file,chi_p_file,y_index, trans_ang):
    """"
    Outputs the confor_vec, containing tg conformatino over time per residue.
    """
    # y_index = resid - 2
    for ang_file in [chi_file,chi_p_file]:
        _ , angle_vec = trj.xvg_reader(ang_file,[1,y_index])
        if ang_file == chi_file: confor_vec = [''] * len(angle_vec)
        for ang_iter , angle_val in enumerate(angle_vec):
            # if  ang_iter == 4874:
            #     print('here')
            confor_vec[ang_iter] +='t' if abs(angle_val) >= trans_ang else 'g'
            
    return confor_vec
def tg_analysis (chi_file,chi_p_file, CNCclass,resid_vec,trans_ang = 150,domain='interior'):
    """
    Analyzes the tg conformation from xvg files.

    Parameters:
        chi_file (str): File path for chi angle data.
        chi_p_file (str): File path for chi prime angle data.
        CNCclass (cnc.CNC_analys): Instance of CNC_analys class.
        resid_vec (list): List of residue indices.
        trans_ang (int): Threshold angle for conformation determination.
        
    Returns:
        dict: Average percentage of tg, gt, gg conformations for chain within each layer.
    """
    if domain!='interior':
        file_ext = '_ends' ## if the exterior chains are needed
    # if domain=='interior':
    #     CNCclass.layer_vec = CNCclass.layer_vec[2:-2]
    tg_info={} 
    iter = -1 ## counting the chain numbers.
    ind_intervs = len(resid_vec)
    for layer in CNCclass.layer_vec:
        tg_info[layer]=[]
        if domain=='interior': 
            chain_number_vec = CNCclass.layers[layer][1:-1] # For interior
        else:
            chain_number_vec = [CNCclass.layers[layer][0], CNCclass.layers[layer][-1]] if len(CNCclass.layers[layer]) > 1\
             else [CNCclass.layers[layer][0]] # for the exterior chains

        for chain_number in chain_number_vec:
            iter += 1
            percent_vec=[]
            for resid in resid_vec:
                y_index = iter * ind_intervs + resid + 2
                # y_index = resid - 2
                confor_vec = calculate_conformation(chi_file,chi_p_file, y_index,trans_ang)
                percent_vec.append([(np.sum(np.array(confor_vec) == confor_key) / len(confor_vec)) * 100 for confor_key in ('gt' , 'gg', 'tg')]) ## an array with rows of resids and cols of tg,gt or ggs

            tg_info[layer].append(np.mean(percent_vec,axis=0))
    return tg_info


def _read_hbs(hb_filename,max_hbnum):
    """"
    Outputs the confor_vec, containing tg conformatino over time per residue.
    """
    indice=[1,2]
    _,hb_data=trj.xvg_reader(hb_filename,indice)
    return (np.array(hb_data)/max_hbnum*100).tolist()
def hb_analysis (CNCclass , hb_type , max_hbnum , HB_dir , domain='interior'):
    """
    Analyzes the hb conformation from xvg files.

    Parameters:
        CNCclass (cnc.CNC_analys): Instance of CNC_analys class.
        hbtype (str): type of hb
        max_hbnum (int): Max number of hb for each chain 
        
    Returns:
        dict: the vector of hb occupancy for each chain
    """
    if domain!='interior':
        file_ext = '_exterior' ## if the exterior chains are needed
    if domain=='interior':
        file_ext = '' ## if the interior chains are needed
        CNCclass.layer_vec = list(CNCclass.layers.keys())[2:-2]

    
    hb_info = {}
    for layer in CNCclass.layer_vec:

        hb_info[layer]=[]
        if domain=='interior': 
            chain_number_vec = CNCclass.layers[layer][1:-1] # For interior
        else:
            chain_number_vec = [CNCclass.layers[layer][0], CNCclass.layers[layer][-1]] if len(CNCclass.layers[layer]) > 1\
            else [CNCclass.layers[layer][0]] # for the exterior chains

        for chain_iter , chain_number in enumerate(chain_number_vec):
            # for HB_iter,HB in enumerate(HB_vec):
            if hb_type=='H6O_O3':
                
                if layer==CNCclass.layer_vec[0] or layer==CNCclass.layer_vec[-1]:
                    continue
                elif chain_iter==len(chain_number_vec)-1:
                    # if domain=='interior':
                    continue

           

            filename = HB_dir + '%s_ch%s_%s%s_hbnum.xvg' % (layer , chain_iter, hb_type , file_ext)
            hb_data = _read_hbs(filename,max_hbnum)
            hb_info[layer].append(hb_data)
    return hb_info

def calculate_phi_psi(phi_file,psi_file,y_index):
    """"
    Outputs the confor_vec, containing tg conformatino over time per residue.
    """
    # y_index = resid - 2
    _ , phi_vec = trj.xvg_reader(phi_file,[1,y_index])
    _ , psi_vec = trj.xvg_reader(psi_file,[1,y_index])
    return [phi_vec[:2600],psi_vec[:2600]]
def phi_psi_analysis (phi_file,psi_file,CNCclass,resid_vec,domain='interior'):
    """
    Analyzes the phi and psi conformation from xvg files.

    Parameters:
        phi_file (str): File path for phi angle data.
        psi_file (str): File path for psi prime angle data.
        CNCclass (cnc.CNC_analys): Instance of CNC_analys class.
        resid_vec (list): List of residue indices.
        
    Returns:
        dict: The vectors of phi and psi angles within each layer.
    """
    if domain =='exterior':
        file_ext = '_exterior' ## if the exterior chains are needed
    if domain=='interior':
        file_ext = '' ## if the interior chains are needed
        # CNCclass.layer_vec = CNCclass.layer_vec[2:-2]
    phi_info={} 
    iter = -1 ## counting the chain numbers.
    ind_intervs = len(resid_vec)
    for layer in CNCclass.layer_vec:
        phi_info[layer]= []
        if domain=='interior': 
            chain_number_vec = CNCclass.layers[layer][1:-1] # For interior
        else:
            chain_number_vec = [CNCclass.layers[layer][0], CNCclass.layers[layer][-1]] if len(CNCclass.layers[layer]) > 1\
             else [CNCclass.layers[layer][0]] # for the exterior chains

        for chain_num_iter, chain_number in enumerate(chain_number_vec):
            resid_data_vec = np.empty((2,0))
            iter += 1
            for resid in resid_vec:
                y_index = iter * ind_intervs + resid + 2
                # y_index = resid - 2
                new_phi_psi = calculate_phi_psi(phi_file,psi_file,y_index)
                resid_data_vec = np.append(resid_data_vec,new_phi_psi,axis = 1)
            phi_info[layer].append(resid_data_vec)
    return phi_info

def calculate_twist(twist_file,y_index):
    """"
    Outputs the confor_vec, containing tg conformatino over time per residue.
    """
    # y_index = resid - 2
    _ , twist_vec = trj.xvg_reader(twist_file,[1,y_index])
    return twist_vec
    # return [abs(x) for x in twist_vec]
def twist_analysis (twist_file,CNCclass,resid_vec,domain='interior'):
    """
    Analyzes the twist conformation from xvg files.

    Parameters:
        twist (str): File path for twist angle data.
        CNCclass (cnc.CNC_analys): Instance of CNC_analys class.
        resid_vec (list): List of residue indices.
        
    Returns:
        dict: The vectors of phi and psi angles within each layer.
    """
    # if domain=='interior':
    #     CNCclass.layer_vec = CNCclass.layer_vec[2:-2]
    twist_info={} 
    iter = -1 ## counting the chain numbers.
    ind_intervs = len(resid_vec)
    for layer in CNCclass.layer_vec:
        twist_info[layer]= []
        if domain=='interior': 
            chain_number_vec = CNCclass.layers[layer][1:-1] # For interior
        else:
            chain_number_vec = [CNCclass.layers[layer][0], CNCclass.layers[layer][-1]] if len(CNCclass.layers[layer]) > 1\
             else [CNCclass.layers[layer][0]] # for the exterior chains

        for chain_num_iter, chain_number in enumerate(chain_number_vec):
            resid_data_vec = []
            iter += 1
            for resid in resid_vec:
                y_index = iter * ind_intervs + resid + 2
                new_twist = calculate_twist(twist_file,y_index)
                resid_data_vec += new_twist
            twist_info[layer].append(resid_data_vec)
    return twist_info