"""
Atom-Projected Phonon Dispersion Plotter

This Python script reads the phonon dispersion data from a phonopy output file "band.yaml"
and plots the atom-projected phonon dispersion using the corresponding eigen vectors.

The main function of this code is to visualize the phonon dispersion of a crystalline
material by plotting its band structure. Additionally, it shows the contribution of each
atom to the phonon modes using the eigen vectors obtained from the "band.yaml" file.

Prerequisites:
    1. Ensure that the "band.yaml" file is present in the working directory, containing the
       phonon dispersion data as generated by the phonopy package.
    2. This version works best for diatomic cells
    
Functionality:
    1. Load phonon dispersion data from "band.yaml" file.
    2. Extract the eigen vectors representing the atom-projected contributions.
    3. Calculate and normalize the amplitudes of the eigen vectors for visualization.
    4. Plot the atom-projected phonon dispersion for each band and atom using scatter plots.
    5. Add color bar to represent the magnitude of atom contributions.
    6. Set appropriate x-axis labels for k-points using provided labels or empty strings.
    7. Add vertical lines at the k-point positions to indicate distinct Brillouin zones.
    
Usage:
    Simply run the script in a directory containing the "band.yaml" file to generate the plot.

Note:
    This code assumes that the "band.yaml" file contains the required phonon dispersion data
    in the format generated by the phonopy package. Make sure to have the necessary
    dependencies installed before running the code.

To do:
    Extend the code for any structure.


Author: Tribhuwan Pandey
"""

import yaml
import numpy as np
import matplotlib.pyplot as plt
import copy
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm
import argparse

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

def load_phonon_data(filename):

    """
    Load phonon data from a YAML file.

    Parameters:
        filename (str): The name of the YAML file containing the phonon data.

    Returns:
        tuple: A tuple containing:
            - D1 (numpy.ndarray): Array of distances for each k-point.
            - F1 (numpy.ndarray): Array of phonon frequencies for each k-point and band.
            - normalized_amplitudes (numpy.ndarray): Array of normalized phonon amplitudes for each k-point, band, and atom.
            - num_bnd (int): Number of bands.
            - kpoint_labels (list): List of labels for each k-point.
            - eigvec (list): List of eigenvectors for each k-point and band.
            - Bcell (numpy.ndarray): Reciprocal lattice vectors.
            - Q1 (numpy.ndarray): Array of k-point positions.
            - B1: Number of bandsegments.
    """

    with open(filename, 'r') as file:
        data = yaml.load(file, Loader=Loader)

    num_kpt = len(data['phonon'])
    num_bnd = len(data['phonon'][0]['band'])
    num_atm = len(data['phonon'][0]['band'][0]['eigenvector'])

    phonon_freq = np.zeros([num_kpt, num_bnd])
    phonon_amplitudes = np.zeros([num_kpt, num_bnd, num_atm])
    normalized_amplitudes = np.zeros([num_kpt, num_bnd, num_atm])
    band_distance = [data['phonon'][ikpt]['distance'] for ikpt in range(num_kpt)]

    labels = []
    for j, v in enumerate(data['phonon']):
        if 'label' in v:
            labels.append(v['label'])
        else:
            labels.append(None)

    for ikpt, kpt_data in enumerate(data['phonon']):
        for ibnd, band_data in enumerate(kpt_data['band']):
            phonon_freq[ikpt, ibnd] = band_data['frequency']
            for iatm, root_eigv in enumerate(band_data['eigenvector']):
                vec = np.array(root_eigv).T
                amplitude = (np.sum(vec.real ** 2 + vec.imag ** 2))
                phonon_amplitudes[ikpt, ibnd, iatm] = amplitude
    normalized_amplitudes = phonon_amplitudes
    for iatm in range(num_atm):
        if iatm == 0:
            normalized_amplitudes[:, :, iatm] = normalized_amplitudes[:, :, iatm]/np.max(normalized_amplitudes[:, :, iatm]) * -0.5
        else:
            normalized_amplitudes[:, :, iatm] = normalized_amplitudes[:, :, iatm]/np.max(normalized_amplitudes[:, :, iatm]) * 0.5

    Bcell = np.array(data['lattice'])
    D1 = np.array(band_distance)
    F1 = phonon_freq
    Q1 = np.array([kpt_data['q-position'] for kpt_data in data['phonon']])
    L1 = []

    if all(x is None for x in labels):
        if 'labels' in data:
            ss = np.array(data['labels'])
            labels = list(ss[0])
            for ii, f in enumerate(ss[:-1,1] == ss[1:,0]):
                if not f:
                    labels[-1] += r'|' + ss[ii+1, 0]
                labels.append(ss[ii+1, 1])
        else:
            labels = []


    eigvec = []
    for j, v in enumerate(data['phonon']):
        if 'eigenvector' in v['band'][0]:
            eigvec_data = [np.array(f['eigenvector']) for f in v['band']]
            eigvec.append(eigvec_data)

    kpoint_labels = [label if label is not None else '' for label in labels]
    return D1, F1, normalized_amplitudes, num_kpt, num_bnd, kpoint_labels, eigvec, Bcell, Q1, data['segment_nqpoint'],  data['natom'], phonon_amplitudes

def plot_cmap_phonon_band(cmap_name, norm):
    """
    Plot the cmap phonon band structure for a specific atom.

    Parameters:
        cmap_name (str): Name of the color map.
        norm: Color normalization object.
    """

    for ibnd in range(num_bnd):
        plt.scatter(band_distance1, phonon_freq1[:, ibnd], 
                    s=10 * abs(normalized_amplitudes_all[:, ibnd]), 
                    c=normalized_amplitudes_all[:, ibnd], 
                    cmap=cmap_name, alpha=0.2, norm=norm)

def plot_scatter_phonon_band(atom_colors):

    """
    Plot the projected phonon band structure for a specific atom.

    Parameters:
        atom_index (int): Index of the atom to project (0-indexed).
        color (str): Color of the scatter plot.
        label (str): Label for the scatter plot.
    """
    for ibnd in range(num_bnd):
        for iatm in range(num_atm):
            plt.scatter(band_distance, phonon_freq[:, ibnd],
                    s=10 * abs(phonon_amplitudes[:, ibnd, iatm]),
                    c=atom_colors[iatm], alpha=0.2)

def add_legend(atom_labels):
    """
    Plot the legend based on atom labels and colors.

    Parameters:
        atom_labels (list): List of labels for different atoms.
        atom_colors (list): List of colors for different atoms.
    """
    for iatm in range(num_atm):
        plt.scatter([], [], label=atom_labels[iatm], color=atom_colors[iatm])
        plt.legend()
        
def write_data_to_file(output_filename_prefix):
    """
    Write distance, frequency, and normalized amplitude of each atom to separate files.

    Parameters:
        output_filename_prefix (str): The prefix for the output file names.
    """
    for iatm in range(num_atm):
        atom_output_filename = f"{output_filename_prefix}_atom_{iatm}.dat"
        with open(atom_output_filename, 'w') as file:
            # Write header
            file.write("Distance (Angstrom) Frequency (THz) Normalized Amplitude (0-0.5)\n")

            # Write data for each k-point, band, and atom
            for ikpt in range(len(band_distance)):
                for ibnd in range(num_bnd):
                    distance = band_distance[ikpt]
                    frequency = phonon_freq[ikpt, ibnd]
                    amplitude = normalized_amplitudes[ikpt, ibnd, iatm]
                    file.write(f"{distance:.6f} {frequency:.6f} {amplitude:.6f}\n")
    
filename = "band.yaml"
band_distance, phonon_freq, normalized_amplitudes, num_kpt, num_bnd, x_labels, eigvec_data, Bcell, Q1, B1, num_atm, phonon_amplitudes = load_phonon_data(filename)

# Prepare data for color map. Gather normalized_amplitudes of all atoms to be used in cmap plot. Pretty sure there a 
# better way of doing this ( I am dumb :( )

normalized_amplitudes_all =  []
for iatm in range(num_atm):
    normalized_amplitudes_all =np.append(normalized_amplitudes_all, normalized_amplitudes[:,:,iatm])

normalized_amplitudes_all=np.reshape(normalized_amplitudes_all, (num_atm*num_kpt,num_bnd))

band_distance1 = np.tile(band_distance, num_atm)
phonon_freq1 = np.tile(phonon_freq, (num_atm, 1))

# Some fancy settings

norm = cm.colors.Normalize(vmax=np.max(normalized_amplitudes_all[:,:]), vmin=np.min(normalized_amplitudes_all[:,:]))
parser = argparse.ArgumentParser(description="Atom-Projected Phonon Dispersion Plotter")
parser.add_argument("--scatter", action="store_true", help="Use scatter plot mode")
parser.add_argument("--cmap", default="cool", type=str, help="Use cmap plot mode with specified colormap (default: cool)")
parser.add_argument("--gnuplot", action="store_true", help="Write output file")
parser.add_argument("--atom1_label", default="atom1", type=str, help="Label for the bottom of the color bar (default: atom1)")
parser.add_argument("--atom2_label", default="atom2", type=str, help="Label for the top of the color bar (default: atom2)")
args = parser.parse_args()

if args.scatter:
    scatter_mode = True
else:
    scatter_mode = False

if scatter_mode:
    # Define colors for different atoms (modify as needed)
    atom_colors = ['cyan', 'm', 'red', 'blue', 'green', 'orange', 'purple', 'pink']
    atom_labels = [args.atom1_label, args.atom2_label]
    output_filename = "Phonon_atom_projected_scatter.pdf"
else:
    atom_colors = None
    output_filename = "Phonon_atom_projected_cmap.pdf"

fig, ax = plt.subplots(figsize=(8, 5))

if not scatter_mode:
    cmap_name = args.cmap
    norm = cm.colors.Normalize(vmax=0.5, vmin=-0.5)
    plot_cmap_phonon_band(cmap_name, norm)
    cbar = plt.colorbar()
    cbar.set_label('Atomic Contribution', rotation=270, labelpad=15)
    cbar_ticks = [-0.5, 0, 0.5]
    cbar_ticklabels = [args.atom1_label, "", args.atom2_label]
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels(cbar_ticklabels)
else:
    plot_scatter_phonon_band(atom_colors)
    add_legend( atom_labels)

if args.gnuplot:
    output_filename_prefix = "atom_projected_phonon_data"
    write_data_to_file(output_filename_prefix)

ax.axhline(y=0, color='black', linestyle='--', linewidth=1)  # Zero line
ax.set_ylabel('Phonon Frequency (THz)')
ax.set_title('Projected Phonon Band Structure')

# Calculate midpoints between k-points
x_klabels = (band_distance[np.r_[[0], np.cumsum(B1)-1]])
ax.set_xlim(x_klabels[0], x_klabels[-1])

# Set the x-axis ticks and labels
plt.xticks(x_klabels, x_labels)

x_klabels=x_klabels[1:-1]
for distance in x_klabels[:]:
    plt.axvline(x=distance, ls='--', color='gray', alpha=0.8, lw=0.5)


# Write the data to a dat file
output_filename_prefix = "atom_projected_phonon_data"
write_data_to_file(output_filename_prefix)
plt.tight_layout()
plt.savefig(output_filename)
plt.show()

