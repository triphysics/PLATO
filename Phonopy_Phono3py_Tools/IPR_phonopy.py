"""
Phonopy Phonon Dispersion Plotter with Inverse Participation Ratio (IPR)

This script reads a YAML file containing phonon band structure data from Phonopy, calculates the Inverse Participation Ratio (IPR) for each phonon mode, and generates a plot with the phonon dispersion overlaid by IPR values.

Usage:
    python IPR_phonopy.py [--f FILENAME] [--cmap COLORMAP] [--alpha ALPHA] [--gnuplot] [--format FORMAT]
                                         [--width WIDTH] [--height HEIGHT] [--ymax YMAX]

    Optional arguments:
    --f FILENAME:           Path to the YAML file containing the phonon band structure data (default: 'band.yaml').
    --cmap COLORMAP:        Name of the colormap to be used for coloring the phonon modes by IPR (default: 'cool').
    --alpha ALPHA:          Transparency value for plot markers (default: 0.6).
    --gnuplot:              Include this flag to write data to an output file.
    --format FORMAT:        Output file format for the plot (default: 'png').
    --width WIDTH:          Width of the figure (default: 6).
    --height HEIGHT:        Height of the figure (default: 4).
    --ymax YMAX:            Maximum y-axis limit in THz.

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

def read_band_yaml(filename):

    """
    Load phonon data from a YAML file.

    Parameters:
        filename (str): The name of the YAML file containing the phonon data.

    Returns:
        tuple: A tuple containing:
            - D1 (numpy.ndarray): Array of distances for each k-point.
            - F1 (numpy.ndarray): Array of phonon frequencies for each k-point and band.
            - num_bnd (int): Number of bands.
            - kpoint_labels (list): List of labels for each k-point.
            - Q1 (numpy.ndarray): Array of k-point positions.
            - B1: Number of bandsegments.
            - gv: is phonon group velocity in each band at each q-point
    """

    with open(filename, 'r') as file:
        data = yaml.load(file, Loader=Loader)

    num_kpt = len(data['phonon'])
    num_bnd = len(data['phonon'][0]['band'])
    num_atm = data['natom']

    phonon_freq = np.zeros([num_kpt, num_bnd])
    ipr=np.zeros([num_kpt, num_bnd])
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
            e_add = 0.0
            count = 0.0
            for iatm, root_eigv in enumerate(band_data['eigenvector']):
                vec = np.array(root_eigv).T
                e_sum = (np.sum(vec.real ** 2 + vec.imag ** 2))**2
                e_add += e_sum
                count += 1
                
            if count != num_atm:
                raise ValueError('\nSomething is not right. Please check your band.yaml file!!!')

            else:
                p = num_atm * e_add
                ipr1=(1/p)

            ipr[ikpt, ibnd]=ipr1

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


    kpoint_labels = [label if label is not None else '' for label in labels]
    return D1, F1, num_kpt, num_bnd, kpoint_labels, Bcell, Q1, data['segment_nqpoint'],  data['natom'], ipr

def plot_cmap_phonon_band(cmap_name, norm, alpha,s):
    """
    Plot the cmap phonon band structure for a specific atom.

    Parameters:
        cmap_name (str): Name of the color map.
        norm: Color normalization object.
    """

    for ibnd in range(num_bnd):
        plt.scatter(band_distance, phonon_freq[:, ibnd], 
                    s=s, 
                    c=ipr[:, ibnd], 
                    cmap=cmap_name, alpha=alpha, norm=norm)

def write_data_to_file(output_filename_prefix):
    """
    Write distance, frequency, and normalized amplitude of each atom to separate files.

    Parameters:
        output_filename_prefix (str): The prefix for the output file names.
    """
    with open(output_filename_prefix, 'w') as file:
        # Write header
        file.write("Distance (Angstrom) Frequency (THz) Inverse PR\n")

        # Write data for each k-point, band
        for ikpt in range(len(band_distance)):
            for ibnd in range(num_bnd):
                distance = band_distance[ikpt]
                frequency = phonon_freq[ikpt, ibnd]
                ipr1 = ipr[ikpt, ibnd]
                file.write(f"{distance:.6f} {frequency:.6f} {ipr1:.6f}\n")


# Some fancy settings
parser = argparse.ArgumentParser(description="Phonon Dispersion Plotter")

parser.add_argument("--f", default="band.yaml", type=str, help="File containing phonon frequncies and group velocities. (default: band.yaml) ")
parser.add_argument("--cmap", default="cool", type=str, help="Use cmap plot mode with specified colormap (default: cool)")
parser.add_argument("--alpha", default=0.6, type=float, help="Use transparency to improve plots. Useful when large number of bands are present (default: 0.6)")
parser.add_argument("--gnuplot", action="store_true", help="Write output file")
parser.add_argument("--format", default="png", type=str, help="Output file format (default: png)")
parser.add_argument("--width", default=6, type=float, help="Width of the figure (default: 6)")
parser.add_argument("--height", default=4, type=float, help="Height of the figure (default: 4)")
parser.add_argument("--ymax", type=float, help="Maximum y-axis limit in THz.")
parser.add_argument("--s", default=1, type=float, help="size of symbol in scatter plot. (default: 1)")

args = parser.parse_args()
filename=args.f

band_distance, phonon_freq, num_kpt, num_bnd, x_labels, Bcell, Q1, B1, num_atm, ipr = read_band_yaml(filename)

output_filename = "Phonon_dispersion_with_IPR."+ args.format

fig, ax = plt.subplots(figsize=(args.width, args.height))

s=args.s
cmap_name = args.cmap
alpha =args.alpha
norm = cm.colors.Normalize(vmax=np.max(ipr), vmin=np.min(ipr))
plot_cmap_phonon_band(cmap_name, norm, alpha, s)
cbar = plt.colorbar()
cbar.set_label('IPR', rotation=270, labelpad=35, fontsize=18)
cbar.ax.tick_params(labelsize=18)

if args.ymax is not None:
    ax.set_ylim(0, args.ymax)

if args.gnuplot:
    output_filename_prefix = "Phonon_dispersion_with_IPR.dat"
    write_data_to_file(output_filename_prefix)

ax.axhline(y=0, color='black', linestyle='--', linewidth=1)  # Zero line
ax.set_ylabel('Frequency (THz)', fontsize=18)


# Calculate midpoints between k-points
x_klabels = (band_distance[np.r_[[0], np.cumsum(B1)-1]])
ax.set_xlim(x_klabels[0], x_klabels[-1])

# Set the x-axis ticks and labels
plt.xticks(x_klabels, x_labels)

ax.tick_params(axis="x", labelsize=18, direction="in", top=False, length=8, which="major", width=1)
ax.tick_params(axis="x", which='minor', length=4, direction="in", top=False, width=1)
ax.tick_params(axis="y", labelsize=18, direction="in", right=True, length=8, which="major", width=1)
ax.tick_params(axis="y",which='minor', length=4, direction="in", right=True, width=1)

x_klabels=x_klabels[1:-1]
for distance in x_klabels[:]:
    plt.axvline(x=distance, ls='--', color='gray', alpha=0.8, lw=0.5)


# Write the data to a dat file
plt.tight_layout()
plt.savefig(output_filename)
plt.show()

