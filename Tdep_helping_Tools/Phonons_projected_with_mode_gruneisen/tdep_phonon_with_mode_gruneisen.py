"""
Group velocity projected phonon dispersion plotter

This Python script reads the phonon dispersion data from a TDEP output file "outfile.dispersion_relations.hdf5"
and project the mode_gruneisen on phonon dispersion.

Prerequisites:
    1. Ensure that the "outfile.dispersion_relations.hdf5" file is present in the working directory, containing the
       phonon dispersion data as generated by TDEP.
    2. Make sure to run tdep phonon calculation with --gruneisen  option. Note this also requires 
       infile.forceconstant_thirdorder file.
 
Usage:
    Simply run the script in a directory containing the "outfile.dispersion_relations.hdf5" file to generate the plot.
```
python tdep_phonon_with_gv.py --cmap jet
```
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
    gp=np.array(data["mode_gruneisen_parameters"])

    return x, ys, q_ticks, q_ticklabels, num_bnd, num_atm,  num_kpt, gp

def plot_cmap_phonon_band(cmap_name, norm, s, alpha):
    """
    Plot the cmap phonon band structure for a specific atom.

    Parameters:
        cmap_name (str): Name of the color map.
        norm: Color normalization object.
    """

    for ibnd in range(num_bnd):
        plt.scatter(band_distance, phonon_freq[:, ibnd],
                    s=s*abs(gp[:, ibnd]),
                    c=gp[:, ibnd],
                    cmap=cmap_name, alpha=alpha, norm=norm)

def write_data_to_file(output_filename_prefix):
    """
    Write distance, frequency, and normalized amplitude of each atom to separate files.

    Parameters:
        output_filename_prefix (str): The prefix for the output file names.
    """
    with open(output_filename_prefix, 'w') as file:
        # Write header
        file.write("Distance (Angstrom) Frequency (THz)  mode_gruneisen_parameters\n")

        # Write data for each k-point, band, and atom
        for ikpt in range(len(band_distance)):
            for ibnd in range(num_bnd):
                distance = band_distance[ikpt]
                frequency = phonon_freq[ikpt, ibnd]
                avg_gv = gp[ikpt, ibnd]
                file.write(f"{distance:.6f} {frequency:.6f} {avg_gv:.6f}\n")

# Some fancy settings
parser = argparse.ArgumentParser(description="Phonon Dispersion Plotter")

parser.add_argument("--f", default="outfile.dispersion_relations.hdf5", type=str, help="File containing phonon frequncies and group velocities. (default: band.yaml) ")
parser.add_argument("--cmap", default="cool", type=str, help="Use cmap plot mode with specified colormap (default: cool)")
parser.add_argument("--alpha", default=0.6, type=float, help="Use transparency to improve plots. Useful when large number of bands are present (default: 0.6)")
parser.add_argument("--gnuplot", action="store_true", help="Write output file")
parser.add_argument("--format", default="png", type=str, help="Output file format (default: png)")
parser.add_argument("--width", default=6, type=float, help="Width of the figure (default: 6)")
parser.add_argument("--height", default=4, type=float, help="Height of the figure (default: 4)")
parser.add_argument("--ymax", type=float, help="Maximum y-axis limit in THz.")
parser.add_argument("--s", default=10, type=float, help="size of the symbols in scatter plots(default: 5)")
args = parser.parse_args()
filename=args.f

band_distance, phonon_freq, x_ticks, x_klabels, num_bnd, num_atm,  num_kpt,gp =  read_band_hdf5(filename)

THz_A_to_Km_s=0.001

output_filename = "Phonon_dispersion_with_mode_gruneisen_parameters."+ args.format

fig, ax = plt.subplots(figsize=(args.width, args.height))

cmap_name = args.cmap
s=args.s
alpha =args.alpha
norm = cm.colors.Normalize(vmax=np.max(gp), vmin=np.min(gp))
plot_cmap_phonon_band(cmap_name, norm, s, alpha)
cbar = plt.colorbar()
cbar.set_label('$\gamma_{qj}$', rotation=0, labelpad=15, fontsize=18)
cbar.ax.tick_params(labelsize=18)


if args.ymax is not None:
    ax.set_ylim(0, args.ymax)

if args.gnuplot:
    output_filename_prefix = "Phonon_dispersion_with_mode_gruneisen_parameters.dat"
    write_data_to_file(output_filename_prefix)

ax.axhline(y=0, color='black', linestyle='--', linewidth=1)  # Zero line
ax.set_ylabel('Frequency (THz)', fontsize=18)


x_klabels = [w.replace('G','Γ') for w in x_klabels]
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
