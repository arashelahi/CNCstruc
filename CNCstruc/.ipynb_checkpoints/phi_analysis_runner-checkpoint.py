#%%
import traj_reader as trj
import numpy as np
import matplotlib.pyplot as plt
from statistics import mode
import CNC_class as cnc
import cnc_analysis_utils
import seaborn as sns
import pandas as pd


gro_file = 'solute.gro'
Data = trj.gro_reader(gro_file)
CNC_group = cnc.CNC_analys(Data)
domain = 'interior'


phi_file = 'Phis_dist.xvg'
psi_file = 'Psis_dist.xvg'
y_ind = 2
resid_vec = [x for x in range(1,12)] ## The number of analyzed residues per side. 12 residues are analyzed, where 6 belongs to one side.

phi_psi_info = cnc_analysis_utils.phi_psi_analysis (phi_file,psi_file,CNC_group,resid_vec)
cent_layer = np.empty((2,0))
orig_layer = np.empty((2,0))
for layer in CNC_group.layer_vec:

        if domain=='interior': 
            chain_number_vec = CNC_group.layers[layer][1:-1] # For interior
        else:
            chain_number_vec = [CNC_group.layers[layer][0], CNC_group.layers[layer][-1]] if len(CNC_group.layers[layer]) > 1\
            else [CNC_group.layers[layer][0]] # for the exterior chains

        for chain_iter , chain_number in enumerate(chain_number_vec):
            phi_psi_data = phi_psi_info[layer][chain_iter]
            n_row = int(len(resid_vec))
            n_col = int(len(phi_psi_data[0])/len(resid_vec))
            phi_psi_chain = np.reshape(phi_psi_data,(2,n_row,n_col))

            if layer in ['L3','L5','L7','L9']:
                cent_layer = np.append(cent_layer ,np.mean(phi_psi_chain,axis=1),axis = 1)
            else:
                orig_layer = np.append(orig_layer ,np.mean(phi_psi_chain,axis=1),axis = 1)

cen_phi_psi_data = np.mean(cent_layer,axis = 1)
cen_phi_psi_data_std = np.std(cent_layer,axis = 1)

orig_phi_psi_data = np.mean(orig_layer,axis = 1)
orig_phi_psi_data_std =  np.std(orig_layer,axis = 1)

print('The phi angle in center and origin chains are %4.2f ± %4.2f and %4.2f ± %4.2f, respectively.\n' % (cen_phi_psi_data[0],cen_phi_psi_data_std[0],orig_phi_psi_data[0],orig_phi_psi_data_std[0]))
print('The psi angle in center and origin chains are %4.2f ± %4.2f and %4.2f ± %4.2f, respectively.\n' % (cen_phi_psi_data[1],cen_phi_psi_data_std[1],orig_phi_psi_data[1],orig_phi_psi_data_std[1]))

##
# df = pd.read_csv("Phi_Psis.csv")

#%%
dict = {'Phi' : cent_layer[0], 'Psi' : cent_layer[1]}
df = pd.DataFrame(data = dict)


fig, ax = plt.subplots(figsize=(6, 6))

levels=5
alpha=0.6
cut=3
thresh=0.001
gridsize=200
bw_adjust=0.8
bw_method="silverman"
linewidths=0.3

sns.kdeplot(
    data=df,
    y="Phi", x="Psi",
    levels=levels,
    fill=True,
    alpha=alpha,
    cut=cut,
    ax=ax,
    thresh=thresh,
    gridsize=gridsize,
    bw_adjust=bw_adjust,
    bw_method=bw_method
)

sns.kdeplot(
    data=df,
    y="Phi", x="Psi",
    levels=levels,
    fill=False,
    alpha=alpha,
    cut=cut,
    ax=ax,
    thresh=thresh,
    gridsize=gridsize,
    bw_adjust=bw_adjust,
    bw_method=bw_method,
    linewidths=linewidths
)
sns.scatterplot(
    data=df,
    y="Phi", x="Psi",
    color="gold",
    s=1,
    linewidth=0.1,
    edgecolor="black",
    ax=ax,
)

# ax.set_ylabel(r'$\phi (^\circ)$')
# ax.set_xlabel(r'$\psi (^\circ)$')

x_lim = [-160,-140]
y_lim = [-105.00, -80.00]
ax.set_ylim(y_lim[0],y_lim[1])
ax.set_xlim(x_lim[0],x_lim[1])
ax.set_xticks(np.arange(x_lim[0],x_lim[1]+1,10))
ax.set_yticks(np.arange(y_lim[0],y_lim[1]+1,10))
fig.savefig("phi_psi_seaborn.png", dpi=1200)#, transparent=True)
plt.show()
#%%

# fig.savefig("phi_psi_seaborn.png", dpi=1200)#, transparent=True)
# %%
