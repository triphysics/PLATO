"""
==============================================================================
                              POSCAR_to_CONTROL.py
==============================================================================

Author: Tribhuwan Pandey
Description:
------------
This Python script is designed to convert VASP POSCAR and BORN files into the
CONTROL file format, which is used as input for the ShengBTE code. ShengBTE is
a software package for calculating lattice thermal conductivity using the
Boltzmann transport equation (BTE) approach.

The script provides the following functionalities:

1. Conversion of BORN File:
   - Reads a BORN file containing the dielectric tensor and Born effective charges.
   - Converts the data into the ALMA format, which includes the dielectric tensor
     and Born charges for each atom.
   - Writes the ALMA format data to the CONTROL file.
   - The code for transforming BEC is motivated by "translate_born.py" code from ALMABTE project
   - For the oginial code please see
     https://almabte.bitbucket.io/downloads/translate-born.py

2. Conversion of POSCAR File:
   - Reads a POSCAR file containing structural information, such as lattice vectors,
     elements, and atomic coordinates.
   - Converts the data into the CONTROL file format, which includes lattice information,
     element definitions, atomic types, and atomic positions.
   - Appends the ALMA format data (from the BORN file conversion) to the CONTROL file.

Usage:
------
To use this script, follow these steps:

1. Ensure you have the following input files in the same directory as this script:
   - BORN: The BORN file containing dielectric tensor and Born effective charges.
   - POSCAR: The POSCAR file containing structural information.

2. Run the script with the desired options:
"""

import sys
import argparse
import numpy as np
import scipy as sp
import scipy.linalg as la
import phonopy

def read_born(filename):
    """
    Read the dielectric tensor and a set of Born effective charges from
    a BORN file.

    Args:
        filename (str): The name of the BORN file.

    Returns:
        Tuple: A tuple containing the dielectric tensor (epsilon) and a list of Born tensors (born).
    """
    born = []
    with open(filename, "r") as f:
        next(f)
        fields = [float(i) for i in next(f).split()]
        epsilon = np.reshape(fields, (3, 3))
        for l in f:
            fields = [float(i) for i in l.split()]
            born.append(np.reshape(fields, (3, 3)))
    return (epsilon, born)

def convert_to_ALMA(epsilon, born, cartesian, ops, classes, uclasses, natoms, outf):
    """
    Convert from Phonopy to ALMA format and return the output as a string.

    Args:
        epsilon (np.ndarray): The dielectric tensor.
        born (list): List of Born tensors.
        cartesian (list): List of Cartesian coordinates for symmetry operations.
        ops (list): List of symmetry operations.
        classes (list): List of atom classes.
        uclasses (list): List of unique atom classes.
        natoms (int): Number of atoms.
        outf: Output file stream.

    Returns:
        str: ALMA format data as a string.
    """
    alma_data = []
    for i in range(3):
        alma_data.append("epsilon(:,{}) = {:>10.6f} {:>10.6f} {:>10.6f}".format(i+1, epsilon[i][0], epsilon[i][1], epsilon[i][2]))
    dborn = {uclasses[i]: born[i] for i in range(len(uclasses))}
    for atom_index in range(1, natoms+1):
        newborn = cartesian[ops[atom_index - 1]].dot(
            dborn[classes[atom_index - 1]].dot(cartesian[ops[atom_index - 1]].T))
        for j in range(3):
            alma_data.append("born(:,{},{}) = {:>10.6f} {:>10.6f} {:>10.6f}".format(j+1, atom_index, newborn[j][0], newborn[j][1], newborn[j][2]))

    return '\n'.join(alma_data)

def parse_arguments():
    """
    Check the validity of command-line arguments and turn them into
    useful information.
    """
    def positive_float(argument):
        """
        Raise an error if argument is not a positive floating point
        number.
        """
        value=float(argument)
        if argument<=0:
            raise argparse.ArgumentTypeError(
                "{} is not a positive floating point number".format(argument))
        return value
    parser = argparse.ArgumentParser(
        description = "Convert BORN files between the Phonopy and ALMA formats. "
        "When no output file is specified, print the result to the standard output.",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", help = "input file name",
                        default = "BORN")
    parser.add_argument("-p", "--poscar",
                        help = "name of a file containing structural information",
                        default = "POSCAR")
    parser.add_argument("-t", "--tolerance", type = positive_float,
                        help = "tolerance for symmetry search",
                        default = 1e-5)
    parser.add_argument("-o", "--output",
                        help = "output filename",
                        default = None)
    parser.add_argument("--loto", action="store_true", help="Write epsilon and born charges to the control file.")
    parser.add_argument("--noiso", action="store_true", help="switch off phonon-isotope scattering")
    parser.add_argument("--rta", action="store_true", help="use relaxation time appoximation only")
    parser.add_argument("--scell", type=int, nargs=3, help="Supercell dimensions as an array [x, y, z]")   
    parser.add_argument("--ngrid", type=int, nargs=3, help="Integration grid for thermal conductivity calculations [x, y, z]")
    parser.add_argument("--tgrid", action="store_true", help="use temperature grid for thermal conductivity calculations")
    parser.add_argument("--tmin", type=float, default="100", help="lowest temperature for thermal conductivity calculations")
    parser.add_argument("--tmax", type=float, default="1000", help="highest temperature for thermal conductivity calculations")
    parser.add_argument("--tstep", type=float, default="100", help="step of temperature for thermal conductivity calculations")
    parser.add_argument("--sb", type=float, default="0.1", help="scale broadening parameter")

#    parser.add_argument("--scell", type=int, nargs=3, default=[5, 5, 2],
#                        help="Supercell dimensions as an array [x, y, z]")
#    parser.add_argument("--ngrid", type=int, nargs=3, default=[5, 5, 2], 
#                        help="Supercell dimensions as an array [x, y, z]")
    return parser.parse_args()

def prepare_loto(args):
    outf = sys.stdout if args.output is None else open(args.output, "w")

    epsilon, born = read_born(args.input)
    nborn = len(born)

    cell = phonopy.interface.vasp.read_vasp(args.poscar)
    natoms = cell.get_number_of_atoms()

    symmetry = phonopy.structure.symmetry.Symmetry(cell, symprec=args.tolerance)
    indep = symmetry.get_independent_atoms()
    nindep = len(indep)

    if nindep == natoms:
        sys.exit("No conversion seems necessary")

    format = "ALMA" if nborn == natoms else "Phonopy"

    classes = symmetry.get_map_atoms()
    ops = symmetry.get_map_operations()
    matrices = symmetry.get_symmetry_operations()["rotations"]
    lattvec = cell.get_cell().transpose()
    cartesian = [lattvec.dot(m.dot(la.inv(lattvec))).T for m in matrices]
    uclasses = sorted(list(set(classes)))

    if format == "Phonopy":
        alma_data = convert_to_ALMA(epsilon, born, cartesian, ops, classes, uclasses, natoms, outf)
        return alma_data
    else:
        sys.exit("Cannot convert from ALMA to Phonopy in this context")

# Function to read POSCAR file and extract data
def read_POSCAR(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    lattice_factor = float(lines[1].strip())/10
    lattice_vectors = [list(map(float, line.split())) for line in lines[2:5]]
    elements = lines[5].split()
    num_atoms = list(map(int, lines[6].split()))
    atom_coords = [list(map(float, line.split())) for line in lines[8:]]

    return lattice_factor, lattice_vectors, elements, num_atoms, atom_coords

# Function to write CONTROL file in the desired format
def write_CONTROL(filename, lattice_factor, lattice_vectors, elements, num_atoms, atom_coords, alma_data):
    with open(filename, 'w') as file:
        file.write("&allocations\n")
        file.write(f"        nelements={len(elements)},\n")
        file.write(f"        natoms={sum(num_atoms)},\n")
        if args.ngrid:
            formatted_ngrid = ' '.join(map(str, args.ngrid))
            file.write(f"  ngrid(:) = {formatted_ngrid}\n")
        else:
            sys.exit("please define integration grid for thermal conductivity calculations")
        file.write("&end\n")
        file.write("&crystal\n")
        file.write(f"   lfactor= {lattice_factor:12.18f}\n")
        for i in range(3):
            formatted_vector = ' '.join([f'{coord:12.18f}' for coord in lattice_vectors[i]])
            file.write(f"   lattvec(:,{i+1})= {formatted_vector}\n")
        elements_str = " ".join([f'"{element}"' for element in elements])
        file.write(f"  elements={elements_str}\n")

#        file.write(f"  elements={','.join(['"' + e + '"' for e in elements])}\n")   
#        file.write(f"  elements=\"{' '.join(elements)}\"\n")

        # Calculate and write the types field
        types = []
        for i, num in enumerate(num_atoms):
            types.extend([i + 1] * num)  # Repeat the atom index by the number of atoms
        file.write(f"  types={' '.join(map(str, types))},\n")

        # Write positions with the specified format
        for i in range(sum(num_atoms)):
            formatted_coords = ' '.join([f'{coord:12.18f}' for coord in atom_coords[i]])
            file.write(f"  positions(:,{i+1}) = {formatted_coords}\n")

        # Write ALMA data if --loto flag is specified
        if args.loto:
            file.write(alma_data+"\n")

        if args.scell:
            formatted_scell = ' '.join(map(str, args.scell))
            file.write(f"  scell(:) = {formatted_scell}\n")    
        else:
            sys.exit("please define supercell size used in phonon calculations")
        file.write("&end\n")
        file.write("&parameters\n")
        if args.tgrid:
            file.write("T_min=" + str(args.tmin) + "\n")
            file.write("T_step=" + str(args.tstep) + "\n")
            file.write("T_max=" + str(args.tmax) + "\n")
        else:    
            file.write("        T=300.\n")
        if args.sb:
            file.write("scalebroad=" + str(args.sb) + "\n")
        else:
            file.write("        scalebroad=0.1\n")
        file.write("&end\n")
        file.write("&flags\n")
        if args.rta:
            file.write("        convergence=.FALSE.\n")
        if args.noiso:
            file.write("        isotopes=.FALSE.\n")
        if args.loto:    
            file.write("        nonanalytic=.TRUE.\n")
        else:
            file.write("        nonanalytic=.FALSE.\n")
        
        file.write("        nanowires=.FALSE.\n")
        file.write("&end\n")
        

# Main program
if __name__ == "__main__":
    poscar_filename = "POSCAR"
    control_filename = "CONTROL"

    lattice_factor, lattice_vectors, elements, num_atoms, atom_coords = read_POSCAR(poscar_filename)
    args = parse_arguments()
    alma_data = prepare_loto(args)
    write_CONTROL(control_filename, lattice_factor, lattice_vectors, elements, num_atoms, atom_coords, alma_data)

