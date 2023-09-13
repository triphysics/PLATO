"""
Interatomic Force Constants Extraction Script

This script is designed to extract interatomic force constants with a smaller cutoff from a previous set of calculations
that were performed with a larger cutoff, all at a fixed supercell size. It accomplishes this by comparing the input
POSCAR files for the smaller cutoff calculations with those of the larger cutoff calculations and copying relevant
VASP output files (vasprun.xml) to new directories for further analysis.

Usage:
    - Place this script in the directory where the smaller cutoff calculations were performed.
    - Ensure that the corresponding larger cutoff calculations are in a parent directory (e.g., "../").
    - Run the script to identify matching calculations and create new directories for the smaller cutoff results.

"""

import os
import shutil
import glob

# Use glob to find files matching the pattern
file_list = glob.glob(os.path.join("3RD.POSCAR.*"))
file_list_large = glob.glob(os.path.join("../3RD.POSCAR.*"))
# Count the matching files
f_snn = len(file_list)
f_lnn = len(file_list_large)
print(f_snn, f_lnn)
  
for i in sorted(os.listdir(".")):
    if i.startswith("3RD.POSCAR."):
        s = i.split(".")[2]

        if s.isnumeric() and 0 < int(s) < f_snn+1:
            aa = False
            d = s

            for j in range(1, f_lnn+1):
                with open(i, "r") as file1, open(f"../3RD.POSCAR.{j:03d}", "r") as file2:
                    if file1.read() == file2.read():
                        print(s, j)
                        os.mkdir(d)
                        shutil.copy(f"../{j:03d}/vasprun.xml", d)
                        aa = True
                        break

            if not aa:
                print(s)

