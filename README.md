# CNCstruc
This package enables the calculation of the different structural properties of cellulose nanocrystals, including primary alcohol conformation, unit cell parameters, twist angles, etc.

## Reference
Please cite the following paper if you use the CNCstruc software:

Elahi, A.; Bidault, X.; Chaudhuri, S. Temperature-transferable coarse-grained model for poly(propylene oxide) to study thermo-responsive behavior of triblock copolymers. J. Phys. Chem. B 2022, 126, 292– 307,  DOI: 10.1021/acs.jpcb.1c06318
Elahi, A.; Yan, X; Chaudhuri, S. Improving Molecular Dynamics Simulations of Cellulose I𝛽: A Cm5-Opls Model for Crystalline Integrity and Surface Functionalization (November 08, 2024). Available at http://dx.doi.org/10.2139/ssrn.5018578

## Applications
The current version of CNCstruc software is only compatible with GROMACS inputs and outputs, and it has two main applications:
1. It reads the coordinate file (pdb or gro), and extracts the atoms involved for quantification of different properties, e.g., hydrogen bonds, alcohol conformations, etc., and creates GROMACS-compatible index files.
![Fig1](https://github.com/user-attachments/assets/4ca30c20-cdd1-4051-afe0-232b4aa13695)

2. Given the index files, and GROMACS output files, this package is used for analysis of the hydrogen bond occupancy, population of _tg_ conformations, etc.

![FigS2](https://github.com/user-attachments/assets/cc14625b-5a0b-442d-a329-d326033a05ed)


