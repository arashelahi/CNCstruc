from CNCstruc.structure import CNC_class as CNC
from CNCstruc.utils import traj_reader as trj
from CNCstruc.analysis import cnc_analysis_utils

CNC_form = 'Pristine'
gro_file = f'./data/input/{CNC_form}/solute.gro'


Data = trj.gro_reader(gro_file)
CNC_group = CNC.CNC_analys(Data)

twist_file = './data/output/twist/%s/interior_chain_twist.xvg' % CNC_form

resid_vec = [x for x in range(1,11)] ## The number of analyzed residues per side. 12 residues are analyzed, where 6 belongs to one side.

twist_info = cnc_analysis_utils.twist_analysis(twist_file , CNC_group , resid_vec)
print(twist_info)

"""
Expected Output:
    L3: [-2.388, -3.478, 8.061, -3.842, -5.488]
    L4: [4.566, 12.224, 21.719, 18.615, 6.833]
    L5: [2.067, -2.57, 10.756, -7.849, 2.624]
    L6: [1.193, -4.351, 4.355, -1.56, 1.349]
    L7: [-3.994, -6.053, -0.008, -3.299, 4.895]
    L8: [21.397, 7.264, 7.869, 7.193, 17.738]
    L9: [6.835, 5.451, -2.087, 3.839, 6.447]
"""