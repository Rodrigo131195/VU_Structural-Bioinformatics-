# VU_Structural-Bioinformatics-

Description:
This repository contains Python scripts developed for an assignment in the field of Structural Bioinformatics. The assignment involves analyzing protein structures to calculate phi-psi angles and amino acid propensities for burial in the protein core. The project aims to implement various algorithms and techniques for extracting structural insights from protein data.

Project Structure:

readDSSP.py: This Python script reads DSSP files containing protein structure data to calculate amino acid propensities for burial in the protein core. It utilizes the unfolded surface accessibility data and DSSP files to perform the analysis. The script generates an output file containing the calculated propensities. Data not provided as is over 25MB.

readPDB.py: This Python script parses PDB files to calculate phi and psi angles and assign secondary structure types to each residue in the protein structure. It utilizes vector operations and geometric calculations to determine the dihedral angles. The script outputs a file containing the calculated angles and secondary structure assignments.

Usage:

Ensure that the required input files (DSSP files for readDSSP.py and PDB files for readPDB.py) are placed in the appropriate directories.
Execute the Python scripts by providing the necessary command-line arguments.
For readDSSP.py, use:

python3 readDSSP.py Data/
For readPDB.py, use:

python3 readPDB.py pdb_filename.txt
Dependencies:

Python 3.x
No additional Python packages are required. The scripts utilize standard libraries such as sys, os, and math.
Note:

Please follow the instructions provided within each script for proper execution and usage.
Ensure that the input files are correctly specified and accessible to the scripts.
The scripts are developed as part of a master's program assignment, focusing on fundamental concepts in structural bioinformatics.
