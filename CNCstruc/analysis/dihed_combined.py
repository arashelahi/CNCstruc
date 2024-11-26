import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

plt.rcParams['mathtext.fontset'] = 'stixsans'  # A sans-serif math font
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.weight'] = 'regular'  # Adjust font weight as needed



def RB_ang(ang,cn_vec):
    pot_RB = 0
    for cn_iter,cn_item in enumerate(cn_vec):\
        pot_RB+=cn_item*(np.cos(ang*np.pi/180))**cn_iter
    return pot_RB

def plot_ang_energy(axs,ang_vec,Chi_vec,Chi_p_vec,tit_str):
    text_size = 10
    color_vec = ['green','blue','green','black','orange','purple']
    Fontsize = 12
    Chi_RB = RB_ang(ang_vec,Chi_vec)
    Chi_p_RB = RB_ang(ang_vec,Chi_p_vec)
    axs.plot(ang_vec,Chi_RB,label=r' $\chi$',color = color_vec[0])
    axs.plot(ang_vec,Chi_p_RB,label=r"$\chi$'",color = color_vec[1])
    leg=axs.legend(ncol=1, fontsize=Fontsize, frameon=True,  fancybox=False,handletextpad=0.3, columnspacing=1.0)
    leg.get_frame().set_edgecolor('k')
    leg.get_frame().set_linewidth(0.8)
    conv_vec=['t','g','t']
    region_names = ['trans','gauche','trans']
    region_split = [-174,0,174]
    for iter_x,x in enumerate(region_split):
        text = axs.text(x,12.50, region_names[iter_x], ha='center', va='center', color='black', fontsize=text_size,style='italic')
    axs.set_xlabel(r'degree ($^\circ$)',fontsize=Fontsize)
    axs.set_ylabel(' Potential Energy (kcal/mol)',fontsize=Fontsize)
    axs.tick_params(axis='both', which='major', labelsize=Fontsize)
    axs.set_xlim(-198,198)
    axs.set_xticks(np.arange(-150,151,100))
    axs.set_yticks(np.arange(0,23.0,5))
    axs.set_title(tit_str, loc='left' ,fontweight='bold',pad=7)
    # plt.show()

def two_D_dih_energy(axs,ang_vec,Chi_vec,Chi_p_vec,tit_str):
    chi_grid, chi_p_grid = np.meshgrid(ang_vec, ang_vec)
    energy_grid = RB_ang(chi_grid,Chi_vec) + RB_ang(chi_p_grid,Chi_p_vec)
    Fontsize=12
    imshow_plt = axs.imshow(energy_grid, extent=(-180, 180, -180, 180), origin='lower', cmap='jet')
    cbar = plt.colorbar(imshow_plt ,label='Energy (kcal/mol)',fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=Fontsize)
    imprp_val = [50.0, 70.0]
    access_color = 'white'
    chi_pract_vec = np.arange(-180,181,1)
    for imprp_iter,imprp in enumerate(imprp_val):
        for half_split in [180,-180]:
            chi_p_prac_vec= chi_pract_vec + half_split- imprp
            axs.plot(chi_pract_vec,chi_p_prac_vec,'--',color=access_color,zorder=2)
    axs.axhline(y=-150, color='black', linestyle='--')
    axs.axhline(y=150, color='black', linestyle='--')
    axs.axvline(x=-150, color='black', linestyle='--')
    axs.axvline(x=150, color='black', linestyle='--')
    # # Add text to each cell
    conv_vec=['t','g','t']
    region_split = [-165,0,165]
    for iter_x,x in enumerate(region_split):
        for iter_y,y in enumerate(region_split):
            text = axs.text(x, y, conv_vec[iter_x]+conv_vec[iter_y], ha='center', va='center', color='forestgreen', fontsize=12,style='italic')
    axs.set_xlabel(r' $\chi$ ($^\circ$)',fontsize=Fontsize)
    axs.set_ylabel(r" $\chi$' ($^\circ$)",fontsize=Fontsize)
    axs.tick_params(axis='both', which='major', labelsize=Fontsize)
    axs.set_xlim(-180,180)
    axs.set_ylim(-180,180)
    axs.set_xticks(np.arange(-150,151,100))
    axs.set_yticks(np.arange(-150,151,100))
    # axs.set_title(r'\textbf{a)}', loc='left')
    # axs.set_title(r'\textbf{a)}', loc='left')
    axs.set_title(tit_str, loc='left' ,fontweight='bold',pad=7)
    


    # axs.set_subtitles('a)')
    # plt.tight_layout(w_pad = 0.1)
    # plt.show()

Chi_vec=[9.03534    ,   -9.03534     ,   0.00     , 0.000    ,   0.00000      ,  0.00000] # Chi
Chi_p_vec=[2.87441    ,    0.58158   ,    2.09200      , -5.54799     ,   0.00000 ,       0.00000] # Chi_p
ang_vec=np.arange(-180,181,1)

fig, axs = plt.subplots(2, 1 ,figsize=(3.5,7.5))
subtitl_vec=['b)','d)']
# subtitl_vec=['','']
plot_ang_energy(axs[0],ang_vec,Chi_vec,Chi_p_vec,subtitl_vec[0])
two_D_dih_energy(axs[1],ang_vec,Chi_vec,Chi_p_vec,subtitl_vec[1])

plt.savefig('dih_energy.png', dpi=300, bbox_inches='tight')
plt.show()