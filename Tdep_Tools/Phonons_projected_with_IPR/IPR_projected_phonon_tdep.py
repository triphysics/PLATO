"""
Phonopy Phonon Dispersion Plotter with Inverse Participation Ratio (IPR)

This script reads a hdf5 file containing phonon band structure data from TDEP, calculates the Inverse Participation Ratio (IPR) for each phonon mode, and generates a plot with the phonon dispersion overlaid by IPR values.

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
import h5py

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

filename="outfile.dispersion_relations.hdf5"

def read_band_hdf5(filename):
    data = h5py.File(filename, "r")
    x = np.array(data["q_values"])
    ys = np.array(data["frequencies"])
    q_ticks = np.array(data["q_ticks"])
    q_ticklabels = data.attrs["q_tick_labels"].decode().split()
    num_kpt=len(x)
    num_bnd = data['frequencies'].shape[1]
    #print(num_bnd)
    num_atm = int(num_bnd/3)
    ipr=np.zeros([num_kpt, num_bnd])
    vec_R=data["eigenvectors_re"]
    vec_I=data["eigenvectors_im"]
    for ikpt in range(num_kpt):
        for nbnd in range(num_bnd):
            e_add = 0.0
            count = 0.0
            vec11=vec_R[ikpt,nbnd, :]
            vec22=vec_I[ikpt,nbnd, :]
            vec_T1 = np.split(vec11, num_atm, axis=0) # plug x,y,z together for each atom
            vec_T2 = np.split(vec22, num_atm, axis=0) # plug x,y,z together for each atom
            for iatm in range(num_atm):
                vec1=vec_T1[iatm]
                vec2=vec_T2[iatm]
                e_sum = (np.sum(vec1 ** 2 + vec2 ** 2))**2
                e_add += e_sum
                count += 1

            if count != num_atm:
                raise ValueError('\nSomething is not right. Please check your band.yaml file!!!')

            else:
                p = num_atm * e_add
                ipr1=(1/p)

            ipr[ikpt, nbnd]=ipr1

    return x, ys, q_ticks, q_ticklabels, num_bnd, num_atm,  num_kpt, ipr


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




band_distance, phonon_freq, x_ticks, x_klabels, num_bnd, num_atm,  num_kpt, ipr =  read_band_hdf5(filename)

print(x_klabels)

parser = argparse.ArgumentParser(description="Phonon Dispersion Plotter with IPR")
parser.add_argument("--f", default="outfile.dispersion_relations.hdf5", type=str, help="File containing phonon frequncies and group velocities. (default: band.yaml) ")
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

x_klabels = [w.replace('G','Γ') for w in x_klabels]

ax.axhline(y=0, color='black', linestyle='--', linewidth=1)  # Zero line
ax.set_ylabel('Phonon Frequency (THz)', fontsize=18)
ax.set_xlim(x_ticks[0], x_ticks[-1])

ax.tick_params(axis="x", labelsize=18, direction="in", top=False, length=8, which="major", width=1)
ax.tick_params(axis="x", which='minor', length=4, direction="in", top=False, width=1)
ax.tick_params(axis="y", labelsize=18, direction="in", right=True, length=8, which="major", width=1)
ax.tick_params(axis="y",which='minor', length=4, direction="in", right=True, width=1)

plt.xticks(x_ticks, x_klabels)


for distance in x_ticks[:]:
    plt.axvline(x=distance, ls='--', color='gray', alpha=0.8, lw=0.5)

plt.tight_layout()
plt.savefig(output_filename)
plt.show()
