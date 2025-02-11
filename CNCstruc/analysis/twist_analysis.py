from cProfile import label
from turtle import color
import traj_reader as trj
import numpy as np
import matplotlib.pyplot as plt
from statistics import mode
import CNC_class as cnc
import os
os.chdir(os.path.dirname(__file__))

def calculate_conformation(twist_pref,y_index):
    """"
    Outputs the confor_vec, containing tg conformatino over time per residue.
    """
    twist_file = twist_pref+'.xvg'
    _ , angle_vec = trj.xvg_reader(twist_file,[1,y_index])
    return angle_vec
def twist_analysis (twist_pref, CNCclass,resid_vec,domain='interior'):
    """
    Analyzes the twist from xvg files.

    Parameters:
        twist_file (str): File path for twist angle data.
        CNCclass (cnc.CNC_analys): Instance of CNC_analys class.
        resid_vec (list): List of residue indices.
        ind_intervs : the intervals betweens the index numbers of the same residue numbers but on consequent chains.
        
    Returns:
        dict: Average percentage of tg, gt, gg conformations for each layer.
    """
    twist_info={} 
    iter = -1 ## counting the chain numbers.
    # ind_intervs = len(resid_vec)
    file_open = open(twist_pref+'_output.dat','w')
    for layer in CNCclass.layer_vec:
        file_open.write('\n%s\t' % layer)
        twist_info[layer]=[]
        if domain=='interior': 
            chain_number_vec = CNCclass.layers[layer][1:-1] # For interior
        else:
            chain_number_vec = [CNCclass.layers[layer][0], CNCclass.layers[layer][-1]] if len(CNCclass.layers[layer]) > 1\
             else [CNCclass.layers[layer][0]] # for the exterior chains
        for chain_iter,chain_number in enumerate(chain_number_vec):
            iter += 1
            for resid_it, resid in enumerate(resid_vec):
                y_index = iter * len(resid_vec) + resid_it + 3
                twist_data = calculate_conformation(twist_pref, y_index)                              
                if resid_it == 0:
                    twist_vec = np.array(twist_data)
                    continue
                twist_vec += np.array(twist_data)
                
            twist_vec_chain = twist_vec/(resid_it+1)
            all_twist_data = twist_vec_chain if iter==0 else all_twist_data + twist_vec_chain
            mean_twist = np.mean(twist_vec_chain)
            std_twist = np.std(twist_vec_chain)
            file_open.write("%4.2f ± %4.2f\t" % (mean_twist,std_twist))
            print("In layer %s, and chain %d the phi is %4.2f ± %4.2f" % (layer,chain_iter,mean_twist,std_twist))
    file_open.close()
    return all_twist_data/(iter+1)
            

    
if __name__ == "__main__":
    gro_file = 'solute.gro'
    Data = trj.gro_reader(gro_file)
    CNC_group = cnc.CNC_analys(Data)
    # twist_pref_exterior = 'exterior_chain_twist'
    # twist_file = twist_pref_exterior+'.xvg'
    # CNC_group.layer_vec = list(CNC_group.layers.keys()) #if the other layers than the interior chains are needed, this line is added.    
    # twist_data_exter = twist_analysis (twist_pref_exterior, CNC_group,CNC_group.resid_vec[:-2],'exterior')


    dt = 40
    # time_vec = np.arange(0,len(twist_data_exter)*dt,dt)
    # plt.plot(time_vec, twist_data_exter,linewidth = 2, label = 'exterior',color='b')


    twist_pref_interior = 'interior_chain_twist'
    twist_file = twist_pref_interior+'.xvg'
    CNC_group.layer_vec = list(CNC_group.layers.keys())[2:-2] # for interior
    twist_data_interior = twist_analysis (twist_pref_interior, CNC_group,CNC_group.resid_vec[:-2],'interior')
    time_vec = np.arange(0,len(twist_data_interior)*dt,dt)

    # print(twist_data_interior)

    plt.plot(time_vec/1000, abs(twist_data_interior),linewidth = 2 ,color='r')
    # plt.legend()
    plt.xlabel('time (ns)', fontsize = 12)
    plt.ylabel(r'twist ($circ$)', fontsize = 12)


    plt.show()