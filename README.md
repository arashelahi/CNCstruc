# CNCstruc
This package enables the calculation of the different structural properties of cellulose nanocrystals, including primary alcohol conformation, unit cell parameters, twist angles, etc.

## Reference
Please cite the following paper if you use the CNCstruc software:

Elahi, A.; Yan, X; Chaudhuri, S. Improving Molecular Dynamics Simulations of Cellulose I𝛽: A Cm5-Opls Model for Crystalline Integrity and Surface Functionalization (November 08, 2024). Available at http://dx.doi.org/10.2139/ssrn.5018578

## Applications
The current version of CNCstruc software is only compatible with GROMACS inputs and outputs, and it has two main applications:
1. It reads the coordinate file (pdb or gro), and extracts the atoms involved for quantification of different properties, e.g., hydrogen bonds, alcohol conformations, etc., and creates GROMACS-compatible index files.  
<div align="center">
<img width="700" alt="Fig1" src="https://github.com/user-attachments/assets/4ca30c20-cdd1-4051-afe0-232b4aa13695">
</div>  
  <br>
  <br>
2. Given the index files, and GROMACS output files, this package is used for analysis of the hydrogen bond occupancy, population of _tg_ conformations, etc.

<div align="center">
<img width="500" alt="Fig2" src="https://github.com/user-attachments/assets/cc14625b-5a0b-442d-a329-d326033a05ed">
</div>

# Installation
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

# Using the Model
## For 1st Application : Making the indices
Loading the structure readers
```
from CNCstruc.structure import CNC_class as CNC
from CNCstruc.utils import traj_reader as trj
```
create the DataFrame of coordinates
```
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

