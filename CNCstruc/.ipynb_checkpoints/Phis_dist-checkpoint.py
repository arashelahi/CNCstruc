import traj_reader as trj
import numpy as np
import matplotlib.pyplot as plt
from statistics import mode
import pandas as pd
filename='solute.gro'
plt.rcParams['mathtext.fontset'] = 'stixsans'  # A sans-serif math font
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.weight'] = 'regular'  # Adjust font weight as needed
################
Phi_file='Phis_dist.xvg'
Psi_file='Psis_dist.xvg'
time,Phi_vec_data=trj.xvg_reader(Phi_file,[1,2])
time,Psi_vec_data=trj.xvg_reader(Psi_file,[1,2])
Dict_phis={'Phi':Phi_vec_data[:2500],'Psi':Psi_vec_data[:2500]}
Phis_to_write=pd.DataFrame(Dict_phis)
Phis_to_write.to_csv('Phi_Psis.csv',index=False)
