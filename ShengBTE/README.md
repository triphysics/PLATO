# VASP POSCAR and BORN to ShengBTE CONTROL File Converter

This Python script is designed to convert VASP POSCAR and BORN files into the CONTROL file format, which is used as input for the ShengBTE code. ShengBTE is a software package for calculating lattice thermal conductivity using the Boltzmann transport equation (BTE) approach.

## Features

- Conversion of BORN File:
  - Reads a BORN file containing the dielectric tensor and Born effective charges.
  - Converts the data into the ALMA format, which includes the dielectric tensor and Born charges for each atom.
  - Writes the ALMA format data to the CONTROL file.

- Conversion of POSCAR File:
  - Reads a POSCAR file containing structural information, such as lattice vectors, elements, and atomic coordinates.
  - Converts the data into the CONTROL file format, which includes lattice information, element definitions, atomic types, and atomic positions.
  - Appends the ALMA format data (from the BORN file conversion) to the CONTROL file.

## Usage

To use this script, follow these steps:

1. Ensure you have the following input files in the same directory as this script:
   - `BORN`: The BORN file containing the dielectric tensor and Born effective charges.
   - `POSCAR`: The POSCAR file containing structural information.

2. Run the script with the desired options:
   ```shell
   python POSCAR_to_CONTROL.py -i BORN -p POSCAR --loto --scell 5 5 2 --ngrid 15 15 6 --noiso --rta
   ```
