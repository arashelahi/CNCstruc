import traj_reader as trj
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import pandas as pd
import CNC_class as cnc
import os
import cnc_analysis_utils
from pathlib import Path

def create_violin_plot(hb_df , HB_dir , layer , chain_iter , hb_vec , domain='interior'):
    # hb_type = inter_hb_vec[0]
    file_ext = '' if domain=='interior' else '_exterior'
    output_dir = Path(HB_dir+'Images/%s/%s' % (domain,'intra' if hb_vec[0] != "H6O_O3" else 'inter'))
    output_dir.mkdir(parents=True, exist_ok=True)  # This line creates the directory if it doesn't exist
    
    
    if hb_vec[0] == "H6O_O3":
        fig, ax = plt.subplots(figsize=((0.78, 0.40)),dpi=300)
        output_file = str(output_dir / f'{layer}_ch{chain_iter}_inter{file_ext}.png')
        hb_df['Occupancy'] = pd.to_numeric(hb_df['Occupancy']) # Convert 'Occupancy' column to numeric
        violins = sns.violinplot(data=hb_df, y='Occupancy', x="Chain", hue='HB type',
                            linewidth=0.5, inner=None, fill=True, 
                            bw=0.3, palette=["powderblue"] , ax=ax)
    
    else:
        fig, ax = plt.subplots(figsize=((1.3, 0.7)),dpi=300)
        output_file = str(output_dir / f'{layer}_ch{chain_iter}_intra{file_ext}.png')
        hb_df['Occupancy'] = pd.to_numeric(hb_df['Occupancy']) # Convert 'Occupancy' column to numeric
        violins = sns.violinplot(data=hb_df, y='Occupancy', x="Chain", hue='HB type',
                            linewidth=0.5, inner=None, fill=True, 
                            bw=0.3 , palette=["red", "royalblue"] , ax=ax, split=True)
        
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.set_yticks([])  
    ax.set(xticklabels=[]) 
    ax.tick_params(left=False, bottom=False)
    ax.legend([],[], frameon=False)
    ax.set_ylim(-20, 120)
    ax.axis('off')
    for violin in violins.collections:
        if isinstance(violin, matplotlib.collections.PolyCollection):
            violin.set_edgecolor('black')
            # violin.set_facecolor('powderblue')  # Keep the fill color if desired
            violin.set_linewidth(0.5)
    plt.savefig(output_file, dpi=300 , bbox_inches = 'tight', pad_inches = 0.01 , transparent=True)
    # plt.close()
    
    


def intra_hb_pandas_maker(layer_number_data , chain_iter,layer,intra_hb_vec,HB_dir, domain='interior'):

    indice = [1 , 2]
    intra_hb_occup = []
    intra_hb_type_data = []
    chain_name_data = []
    max_hbnum = 12
    file_ext = '' if domain=='interior' else '_exterior'
    for hb_iter , hb_type in enumerate(intra_hb_vec):
        filename = HB_dir + '%s_ch%s_%s%s_hbnum.xvg' % (layer , chain_iter, hb_type,file_ext)
        _ , hb_data = trj.xvg_reader(filename,indice)
        intra_hb_type_data = intra_hb_type_data + [hb_type]*len(hb_data)
        intra_hb_occup = intra_hb_occup + hb_data
        layer_number_data = layer_number_data + [layer]*len(hb_data)
        chain_name_data = chain_name_data + ['ch%s' % chain_iter]*len(hb_data)
    hb_df = pd.DataFrame([layer_number_data, chain_name_data, intra_hb_type_data, np.array(intra_hb_occup)/max_hbnum*100]).T
    hb_df.columns = ["Layer", "Chain", "HB type", "Occupancy"]
    create_violin_plot(hb_df , HB_dir , layer , chain_iter , intra_hb_vec , domain)
    return hb_df
def inter_hb_pandas_maker(layer_number_data , chain_iter,layer,inter_hb_vec,HB_dir, domain='interior'):
    xvg_file_number = chain_iter
    if layer == CNC_group.layer_vec[0] or layer==CNC_group.layer_vec[-1]:
        return
    if domain=='interior':
        if chain_iter == len(chain_number_vec)-1:
            return
    elif domain=='exterior':
        if chain_iter == 1:
            if len(CNC_group.layers[layer]) <= 2:
                return
            else:
                xvg_file_number = len(CNC_group.layers[layer])-2
                # chain_iter = len(CNC_group.layers[layer])-2
    indice = [1 , 2]
    inter_hb_occup = []
    inter_hb_type_data = []
    chain_name_data = []
    max_hbnum = 11
    file_ext = '' if domain=='interior' else '_exterior'
    for hb_iter , hb_type in enumerate(inter_hb_vec):
        filename = HB_dir + '%s_ch%s_%s%s_hbnum.xvg' % (layer , xvg_file_number, hb_type,file_ext)
        _ , hb_data = trj.xvg_reader(filename,indice)
        inter_hb_type_data = inter_hb_type_data + [hb_type]*len(hb_data)
        inter_hb_occup = inter_hb_occup + hb_data
        layer_number_data = layer_number_data + [layer]*len(hb_data)
        chain_name_data = chain_name_data + ['ch%s' % chain_iter]*len(hb_data)
    hb_df = pd.DataFrame([layer_number_data, chain_name_data, inter_hb_type_data, (np.array(inter_hb_occup)/max_hbnum*100).tolist()]).T
    hb_df.columns = ["Layer", "Chain", "HB type", "Occupancy"]
    create_violin_plot(hb_df , HB_dir , layer , chain_iter , inter_hb_vec , domain)
    return hb_df


filename = './simulation_traj_topol/solute.gro'
domain_vec = ['interior','exterior']
domain = 'exterior'
HB_dir = './HBs/'
intra_hb_vec=['H2O_O6','H3O_O5']
inter_hb_vec=[ 'H6O_O3' ]
indice = [1 , 2]

for domain in domain_vec:
    Data = trj.gro_reader(filename)
    CNC_group = cnc.CNC_analys(Data,domain)


    for layer in CNC_group.layer_vec:
        layer_number_data = []    
        if domain=='interior': 
            chain_number_vec = CNC_group.layers[layer][1:-1] # For interior
        else:
            chain_number_vec = [CNC_group.layers[layer][0], CNC_group.layers[layer][-1]] if len(CNC_group.layers[layer]) > 1\
            else [CNC_group.layers[layer][0]] # for the exterior chains
        for chain_iter , chain_number in enumerate(chain_number_vec):
            intra_hb_df = intra_hb_pandas_maker(layer_number_data , chain_iter,layer,intra_hb_vec,HB_dir,domain)
        # inter_hb_df = inter_hb_pandas_maker(layer_number_data , chain_iter,layer,inter_hb_vec,HB_dir,domain)
        # create_violin_plot(inter_hb_df , HB_dir , layer , chain_iter , inter_hb_vec)