# CNCstruc
This package enables the calculation of the different structural properties of cellulose nanocrystals, including primary alcohol conformation, unit cell parameters, twist angles, etc. Please cite the following publication if you are using this software:

Elahi, A.; Yan, X.; Chaudhuri, S. Enhancing the OPLS-AA force field for cellulose Iβ: structural stability and surface functionalization capability with the CM5 charge model. Carbohydr. Polym. 2025, 360, 123572 

This paper also massively uses OPLSCM5 Click [here](quora.com/profile/Ashish-Kulkarni-100).  software developed through this research to create GROMACS-compatible topology files. Here is the link to OPLSCM5 software. 


## Applications
The current version of CNCstruc software is only compatible with GROMACS inputs and outputs, and it has two main applications:
1. It reads the coordinate file (pdb or gro), and extracts the atoms involved for quantification of different properties, e.g., hydrogen bonds, alcohol conformations, etc., and creates GROMACS-compatible index files.  
<div align="center">
<img width="700" alt="Fig1" src="https://github.com/user-attachments/assets/4ca30c20-cdd1-4051-afe0-232b4aa13695">
</div>  
  <br>
  <br>
2. Given the index files, and GROMACS output files, this package is used for analysis of the hydrogen bond occupancy, population of <i>tg</i> conformations, etc.

<div align="center">
<img width="500" alt="Fig2" src="https://github.com/user-attachments/assets/cc14625b-5a0b-442d-a329-d326033a05ed">
</div>

## Installation
Cloning the repository and creating a conda environment.
``` 
git clone https://github.com/arashelahi/CNCstruc.git
cd CNCstruc
conda create -n CNCstruc
```
Installing the prerequisite packages
```
pip install numpy
pip install scikit-learn
pip install pandas
pip install matplotlib
```

Installing the package
```
pip install -e .
```

## Using the Model
### For 1st Application : Making the indices
Loading the structure readers
```
from CNCstruc.structure import CNC_class as CNC
from CNCstruc.utils import traj_reader as trj
```
create the DataFrame of coordinates
```
CNC_form = 'Pristine'
filename = './data/input/{CNC_form}/solute.gro'
Data = trj.gro_reader(filename)
CNC_group = CNC.CNC_analys(Data)
```
Make index files for atoms involved in twist rate calculations

```
CNC_group.twist_angles()
ndx_file='twist_chain.ndx' % CNC_form
for chi_key in CNC_group.twist_data.keys():   
    trj.ndx_writer(ndx_file,CNC_group.twist_data[chi_key] ,chi_key+'s')
```

### For 2nd Application : Analysis of the Structure
Loading the structure readers and analysis module
```
from CNCstruc.structure import CNC_class as CNC
from CNCstruc.utils import traj_reader as trj
from CNCstruc.analysis import cnc_analysis_utils
```
create the DataFrame of coordinates
```
CNC_form = 'Pristine'
gro_file = f'./data/input/{CNC_form}/solute.gro'
Data = trj.gro_reader(gro_file)
CNC_group = CNC.CNC_analys(Data)
```
Read the twist calculations outputed from GROMACS (xvg file) and analyze the quantifications for each layer within the CNC:
```
twist_file = './data/output/twist/%s/interior_chain_twist.xvg' % CNC_form
twist_info = cnc_analysis_utils.twist_analysis(twist_file , CNC_group)
print twist_info
```
```
# Expected Output:
L3: [-2.388, -3.478, 8.061, -3.842, -5.488]
L4: [4.566, 12.224, 21.719, 18.615, 6.833]
L5: [2.067, -2.57, 10.756, -7.849, 2.624]
L6: [1.193, -4.351, 4.355, -1.56, 1.349]
L7: [-3.994, -6.053, -0.008, -3.299, 4.895]
L8: [21.397, 7.264, 7.869, 7.193, 17.738]
L9: [6.835, 5.451, -2.087, 3.839, 6.447]
```
<div align="center">
<img width="500" alt="Fig3" src="https://github.com/user-attachments/assets/938416c8-326e-44aa-848b-c7dca8a1a9ff">
</div>  

