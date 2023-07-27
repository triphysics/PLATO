1. modol_thermal_conductivity_analysis.py
# Thermal Conductivity Data Analysis Script
This script is designed to process thermal conductivity data from a tdep code, specifically from two HDF5 files:
outfile.cumulative_kappa.hdf5: Contains cumulative thermal conductivity data.
outfile.grid_thermal_conductivity.hdf5: Contains spectral thermal conductivity data.

The current version of code performs the following tasks:
1. Loads and processes the cumulative thermal conductivity data from outfile.cumulative_kappa.hdf5.
2. Plots the cumulative thermal conductivity vs. frequency in a separate plot (Part A).
3. Loads and processes the spectral thermal conductivity data from outfile.grid_thermal_conductivity.hdf5.
4. Plots the spectral thermal conductivity vs. frequency for each direction (sk_xx, sk_yy, sk_zz) in a separate plot (Part B).
5. Optionally, dumps the cumulative thermal conductivity data in a file named cumulative_kappa_tensor.dat
   and the spectral thermal conductivity data in a file named Spectral_kappa_tensor.dat for gnuplot.

## Usage

Run the script using the following command:

```bash
python thermal_conductivity_analysis.py [--gnuplot] [--cumulative_xrange XMIN XMAX] [--cumulative_yrange YMIN YMAX]
                                        [--spectral_xrange XMIN XMAX] [--spectral_yrange YMIN YMAX] [--gridon]
 
