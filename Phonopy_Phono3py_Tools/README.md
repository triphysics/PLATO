# 1. Atom-Projected Phonon Dispersion Plotter

This Python script reads phonon dispersion data from a phonopy output file "band.yaml" and plots the atom-projected phonon dispersion using the corresponding eigenvectors. The script provides options to generat scatter plots and colormap (cmap) plots for visualizing the atom contributions to the phonon modes.

## Features
1. Load phonon dispersion data from the "band.yaml" file.
2. Extract eigenvectors representing atom-projected contributions.
3. Normalize amplitudes of eigenvectors for visualization.
4. Plot atom-projected phonon dispersion using scatter plots or colormap plots.

## Usage
1. Place the "band.yaml" file in the working directory. YbS band.yaml file is supplied with the code.
2. Run the code using the following command:

```shell
python atom_projected_phonons.py [--scatter] [--cmap COLORMAP] [--gnuplot] [--atom1_label LABEL] [--atom2_label LABEL]
```

# 2.  Group Velocity Projected Phonon Dispersion Plotter

This Python script reads phonon dispersion data from a phonopy output file named "band.yaml" and creates a plot of the average phonon group velocities projected onto the phonon dispersion curve. 

## Features
1. Load phonon dispersion data from the "band.yaml" file.
2. Calculate the norm of phonon group velocity and plot them on phonon dispersion/

## Usage
1. Place the "band.yaml" file in the working directory. MgO band_gv.yaml file is supplied with the code.
2. Run the code using the following command:

```shell
python group_velocity_projected_phonons.py --f band_gv.yaml
```
