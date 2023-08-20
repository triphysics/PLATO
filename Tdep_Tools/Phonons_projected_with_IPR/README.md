# 3. TDEP Phonon Dispersion Plotter with Inverse Participation Ratio (IPR)

This script reads a hdf5 file containing phonon band structure data from TDEP, calculates the Inverse Participation Ratio (IPR) for each phonon mode, and generates a plot with the phonon dispersion overlaid by IPR values.

## Usage
```shell
    python IPR_projected_phonon_tdep.py [--f FILENAME] [--cmap COLORMAP] [--alpha ALPHA] [--gnuplot] [--format FORMAT]
                                         [--width WIDTH] [--height HEIGHT] [--ymax YMAX]
