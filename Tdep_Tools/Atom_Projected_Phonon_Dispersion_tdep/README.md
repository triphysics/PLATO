# Atom-Projected Phonon Dispersion Plotter

This Python script reads phonon dispersion data from a phonopy output file "outfile.dispersion_relations.hdf5" and plots the atom-projected phonon dispersion using the corresponding eigenvectors. The script provides options to generate scatter plots and colormap (cmap) plots for visualizing the atom contributions to the phonon modes.

## Prerequisites

1. Ensure that the "outfile.dispersion_relations.hdf5" file is present in the working directory, containing the phonon dispersion data as generated by the TDEP package.
2. This version best works for simple diatomic structures.

## Features

1. Load phonon dispersion data from the "outfile.dispersion_relations.hdf5" file.
2. Extract eigenvectors representing atom-projected contributions.
3. Normalize amplitudes of eigenvectors for visualization.
4. Plot atom-projected phonon dispersion using scatter plots or colormap plots.

## Usage

1. Install the required dependencies: `numpy`, `matplotlib`.
2. Place the "outfile.dispersion_relations.hdf5" file in the working directory. MgO outfile.dispersion_relations.hdf5 file is supplied with the code.
3. --gnuplot will dump the data for plotting in gnuplot.
4. Run the script using the following command:

```shell
python atom_projected_phonons_tdep.py [--scatter] [--cmap COLORMAP] [--gnuplot] [--atom1_label Mg] [--atom2_label O]

