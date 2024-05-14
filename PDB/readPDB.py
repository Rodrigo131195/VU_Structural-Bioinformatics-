#!/usr/bin/python

"""
Structural Bioinformatics assignment 1 - Phi-psi angles

When you finish the assignment, the script should create an
output file containing the phi and psi angles and secondary
structure assignment of each residue. Please ONLY modify the
code in the three indicated blocks and do NOT use additional
python packages. Use spaces instead of tabs. Do not round
down in any calculation step. The commented #print() lines
can be used to test separate functions, but make sure that
all of them are commented when submtting to CodeGrade.

To run, make 'PDB' your working directory and use:

> python3 readPDB.py pdb_filename.txt
"""

# Packages
from sys import argv
import os
from math import sqrt, atan2, degrees



# Vector functions that we need to calculate the angles
def dot_product(v1, v2):
    """ Calculate the dot product of two vectors """
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
#print(dot_product([1, 2, 3], [1, 3, 2]))

def cross_product(v1, v2):
    """ Calculate the cross product of two vectors """
    i = v1[1]*v2[2] - v1[2]*v2[1]
    j = v1[2]*v2[0] - v1[0]*v2[2]
    k = v1[0]*v2[1] - v1[1]*v2[0]
    return [i,j,k]
#print(cross_product([1, 2, 3], [1, 3, 2]))

def magnitude(v):
    """ Calculate the size of a vector """
    return sqrt(v[0]**2 + v[1]**2 + v[2]**2)
#print(magnitude([1, 2, 2]))

# PDB file parser
def readPDB(PDB_file):
    """ Reads a PDB file and stores the atom
    coordinates and amino acid types of the protein """
    # open the file
    f = open(PDB_file, 'r')

    # dictionaries to store the output
    # pdb atom coordinates:
    #     pdbcoord[chain][residue_number][atom_type] = coordinates
    pdbcoord = {}
    # residue type per chain and residue number (i.e. store the sequence)
    #     pdbseq[chain][resnum] = restype
    pdbseq = {}

    # parse each line in the file
    for line in f:
        # remove whitespace at the end of the line
        line = line.strip()
        # only parse the lines containing atom coordinates
        if line[:4] == 'ATOM':
            # ATOM type (e.g. C-alpha)
            atom_type = line[12:16].strip()
            # AMINO ACID type (e.g. alanine)
            aa_type = line[17:20].strip()
            # residue number
            res_num = int(line[22:26])
            # Protein chain
            chain = line[21]
            # coordinates
            xcoord = float(line[30:38])
            ycoord = float(line[38:46])
            zcoord = float(line[46:54])

            # if chain does not exists create new entry
            if not chain in pdbcoord:
                pdbcoord[chain] = {}
                pdbseq[chain] = {}
            # if resnum does not exists create new entry
            if not res_num in pdbcoord[chain]:
                pdbcoord[chain][res_num] = {}

            # store coordinates as a vector
            pdbcoord[chain][res_num][atom_type] = [xcoord,ycoord,zcoord]
            # store sequence
            pdbseq[chain][res_num] = aa_type

    # close file
    f.close()

    # return dictionaries
    return pdbcoord, pdbseq
#print(readPDB('Data\\1TIM.pdb'))

### THE FOLLOWING THREE FUNCTIONS ARE THE ONES YOU NEED
### TO EDIT FOR THE ASSIGNMENT. ONLY EDIT THE INDICATED
### BLOCKS
def calculateDihedral(a1, a2, a3, a4):

    """ Calculates the normal vector of the planes
    defined by four atom coordinates """
    ### START CODING HERE
    u, v, w, g = [],[],[], []


    for i in range(3):
        u.append(a2[i] - a1[i])
        v.append(a3[i] - a2[i])
        w.append(a3[i] - a2[i])
        g.append(a4[i] - a3[i])

    n1 = cross_product(u , v)
    n2 = cross_product(w , g)

    # Compute the sine and cosine of the dihedral angle
    if dot_product(cross_product(n1, n2), u) < 0:
        sin_value = -magnitude(cross_product(n1, n2)) / (
                    magnitude(n1) * magnitude(n2))
    else:
        sin_value = magnitude(cross_product(n1, n2)) / (
                    magnitude(n1) * magnitude(n2))
    cos_value = dot_product(n1, n2) / (
                magnitude(n1) * magnitude(n2))

    # Calculate the dihedral angle in degrees
    dihedral = degrees(atan2(sin_value, cos_value))

    # Replace the line above with your own code

    return dihedral

    ### END CODING HERE
#print(calculateDihedral([1, 9, 2], [3, 2, 1], [2, 4, 7], [8, 2, 5]))

def assign_ss(phi, psi):
    """ Assign a secondary structure type based on the phi
    and psi angles of a residue """
    ### START CODING HERE
    ### START CODING HERE
    # for code checking purposes use the terms "loop", "alpha" or "beta"
    if (((phi > -160 and phi < 0) or (phi > 30 and phi < 90))) and (psi > -100 and psi < 100):
        secondary_structure = "alpha"
    elif (phi > -180 and phi < -15) or (phi > 160 and phi < 180):
        secondary_structure = "beta"
    else:
        secondary_structure = "loop"
    ### END CODING HERE
    return secondary_structure


# print(assign_ss(-50, -60))
def print_phi_psi(pdbcoord, pdbseq, outfile):
    """ given the PDB coordinates, calculate the dihedral
    angles of all the residues, assign secondary structure
    types and write them into an output file """
    f = open(outfile, 'w')

    # get the chains
    list_chains = sorted(pdbcoord.keys())

    for chain in list_chains:
        # get the sorted residue numbers
        list_residue_numbers = sorted(pdbcoord[chain].keys())
        for res_num in list_residue_numbers:
            try:
                # Get coordinates for the atoms
                a1 = pdbcoord[chain][res_num]["CA"]
                a2 = pdbcoord[chain][res_num]["N"]
                a3 = pdbcoord[chain][res_num]["C"]

                # Check if the previous and next residues exist
                if res_num - 1 in pdbcoord[chain] and res_num + 1 in pdbcoord[chain]:
                    a4 = pdbcoord[chain][res_num - 1]["C"]
                    a5 = pdbcoord[chain][res_num + 1]["N"]

                    # Calculate phi and psi angles
                    phi = calculateDihedral(a4, a2, a1, a3)
                    psi = calculateDihedral(a2, a1, a3, a5)

                    # Assign secondary structure
                    ss = assign_ss(phi, psi)
                else:
                    phi, psi, ss = None, None, None

                    # set N phi angle to None
                    if res_num == min(list_residue_numbers):
                        phi = None
                        # Calculate the psi angle
                        if res_num + 1 in pdbcoord[chain]:
                            a5 = pdbcoord[chain][res_num + 1]["N"]
                            psi = calculateDihedral(a2, a1, a3, a5)

                    elif res_num == max(list_residue_numbers):
                        psi = None
                        # Calculate the phi angle
                        if res_num - 1 in pdbcoord[chain]:
                            a4 = pdbcoord[chain][res_num - 1]["C"]
                            phi = calculateDihedral(a4, a2, a1, a3)

            except KeyError:
                # If any required atom is missing, print a warning message
                print('WARNING: KeyError occurred for residue', chain, res_num)
                # Set phi, psi, and secondary structure to None
                phi, psi, ss = None, None, None

            # get amino acid
            aa_type = pdbseq[chain][res_num]
            # write into output file
            print(chain, res_num, aa_type, phi, psi, ss, file=f)
    f.close()
    print('written:', outfile)



def main():
    # input PDB file
    f_in = argv[1]
    f_out = 'Output/phi_psi.txt'
    # read PDB file
    pdbcoord, pdbseq = readPDB(f_in)
    print_phi_psi(pdbcoord, pdbseq, f_out)
    # for testing
    # for i in ['1TIM', '3PG8']:
    #     f_in = 'student/{}.pdb'.format(i)
    #     print(f_in)
    #     f_out = 'student/output/phi_psi_{}.txt'.format(i)
    #
    #     # read PDB file
    #     pdbcoord, pdbseq = readPDB(f_in)
    #     print_phi_psi(pdbcoord, pdbseq, f_out)

if __name__ == '__main__':
    main()
