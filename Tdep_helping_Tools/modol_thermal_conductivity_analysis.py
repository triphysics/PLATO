"""
Thermal Conductivity Data Analysis Script
-----------------------------------------
Author: Tribhuwan Pandey

This script is designed to process thermal conductivity data from a tdep code, specifically from two HDF5 files:
1. outfile.cumulative_kappa.hdf5: Contains cumulative thermal conductivity data.
2. outfile.grid_thermal_conductivity.hdf5: Contains spectral thermal conductivity data.

The script performs the following tasks:
1. Loads and processes the cumulative thermal conductivity data from outfile.cumulative_kappa.hdf5.
2. Plots the cumulative thermal conductivity vs. frequency in a separate plot (Part A).
3. Loads and processes the spectral thermal conductivity data from outfile.grid_thermal_conductivity.hdf5.
4. Plots the spectral thermal conductivity vs. frequency for each direction (sk_xx, sk_yy, sk_zz) in a separate plot (Part B).
5. Optionally, dumps the cumulative thermal conductivity data in a file named cumulative_kappa_tensor.dat
   and the spectral thermal conductivity data in a file named Spectral_kappa_tensor.dat for gnuplot.

Usage:
python Modol_thermal_conductivity_analysis.py [--gnuplot] [--cumulative_xrange XMIN XMAX] [--cumulative_yrange YMIN YMAX]
                                        [--spectral_xrange XMIN XMAX] [--spectral_yrange YMIN YMAX]

Optional arguments:
--gnuplot              Dump data to files for gnuplot.
--cumulative_xrange    Set the x-axis range for the cumulative plot (Part A).
--cumulative_yrange    Set the y-axis range for the cumulative plot (Part A).
--spectral_xrange      Set the x-axis range for the spectral plot (Part B).
--spectral_yrange      Set the y-axis range for the spectral plot (Part B).

"""
import argparse
import numpy as np
import h5py
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'Times New Roman'

# Set the font size for labels and title
plt.rcParams['font.size'] = 14

# Define functions load_and_process_cumulative_kappa_data(), load_and_process_modal_kappa_data(),

def load_and_process_cumulative_kappa_data(filepath):
    """
    Load and process cumulative kappa data from an HDF5 file.

    Parameters:
        filepath (str): Path to the HDF5 file containing cumulative kappa data.

    Returns:
        numpy.ndarray: Frequency data.
        numpy.ndarray: Spectral kappa data along the diagonal (sk1, sk2, sk3).
    """
    with h5py.File(filepath, "r") as f:
        xx = f['temperature_1']
        freq = xx['frequency_axis'][:]
        sk1 = xx['spectral_kappa_vs_frequency_per_direction'][0, 0, :]
        sk2 = xx['spectral_kappa_vs_frequency_per_direction'][1, 1, :]
        sk3 = xx['spectral_kappa_vs_frequency_per_direction'][2, 2, :]


    return freq, sk1, sk2, sk3

def load_and_process_modal_kappa_data(filepath):
    """
    Load and process modal kappa data from an HDF5 file.

    Parameters:
        filepath (str): Path to the HDF5 file containing modal kappa data.

    Returns:
        numpy.ndarray: Frequency data in THz.
        numpy.ndarray: Modal kappa data (mode_kxx, mode_kyy, mode_kzz).
        int: Number of modes (qp).
    """
    with h5py.File(filepath, "r") as f:
        fp_list = f['frequencies'].shape
        qp = fp_list[0]
        toTHz = 1 / 1E12 / 2 / np.pi
        freq_n = f['frequencies'][:, :] * toTHz
        mode_kxx = f['thermal_conductivity'][:, :, 0, 0]
        mode_kyy = f['thermal_conductivity'][:, :, 1, 1]
        mode_kzz = f['thermal_conductivity'][:, :, 2, 2]


    return freq_n, mode_kxx, mode_kyy, mode_kzz, qp


# plot_cumulative_kappa_data(), plot_spectral_kappa_data(), dump_data_for_gnuplot() here

def plot_cumulative_kappa_data(df, xrange=None, yrange=None, grid=False):
    # Plotting the cumulative data
    plt.plot(df['Freq'], df['cumulative_kxx'], label=r'$\kappa_{xx}$')
    plt.plot(df['Freq'], df['cumulative_kyy'], label=r'$\kappa_{yy}$')
    plt.plot(df['Freq'], df['cumulative_kzz'], label=r'$\kappa_{zz}$')
    plt.xlabel('Frequency (THz)')
    plt.ylabel('Cumulative Thermal Conductivity (W/m-K)')
    plt.legend()
    if grid:
        plt.grid(True)
    else:
        plt.grid(False)
    plt.title('Cumulative Thermal Conductivity vs. Frequency')
    if xrange:
        plt.xlim(xrange[0], xrange[1])
    if yrange:
        plt.ylim(yrange[0], yrange[1])

def plot_spectral_kappa_data(freq, sk1, sk2, sk3, xrange=None, yrange=None, grid=False):
    # Plotting the spectral kappa data
    plt.plot(freq, sk1, label=r'$\kappa_{xx}$')
    plt.plot(freq, sk2, label=r'$\kappa_{yy}$')
    plt.plot(freq, sk3, label=r'$\kappa_{zz}$')
    plt.xlabel('Frequency (THz)')
    plt.ylabel('Spectral Thermal Conductivity (W/m-K)')
    plt.legend()
    if grid:
        plt.grid(True)
    else:
        plt.grid(False)
    plt.title('Spectral Thermal Conductivity vs. Frequency')
    if xrange:
        plt.xlim(xrange[0], xrange[1])
    if yrange:
        plt.ylim(yrange[0], yrange[1])

def dump_data_for_gnuplot(df):
    # Dump data to files for gnuplot
    df.to_csv('cumulative_kappa_tensor.dat', index=False, sep=' ', header=True)

def main():
    parser = argparse.ArgumentParser(description="Process cumulative and modal kappa data and optionally plot or dump data in files.")
    parser.add_argument("--gnuplot", action="store_true", help="Dump data to files for gnuplot.")
    parser.add_argument("--cumulative_xrange", nargs=2, type=float, help="Set the x-axis range for the cumulative plot (Part A).")
    parser.add_argument("--cumulative_yrange", nargs=2, type=float, help="Set the y-axis range for the cumulative plot (Part A).")
    parser.add_argument("--spectral_xrange", nargs=2, type=float, help="Set the x-axis range for the spectral plot (Part B).")
    parser.add_argument("--spectral_yrange", nargs=2, type=float, help="Set the y-axis range for the spectral plot (Part B).")
    parser.add_argument("--gridon", action="store_true", help="Enable grid in the plots.")
    args = parser.parse_args()

    # Load and process cumulative kappa data
    cumulative_kappa_file = "outfile.cumulative_kappa.hdf5"
    freq, sk1, sk2, sk3 = load_and_process_cumulative_kappa_data(cumulative_kappa_file)

    # Save spectral_kappa_tensor.dat
    data_sk1 = np.vstack([freq, sk1, sk2, sk3])
    np.savetxt("Spectral_kappa_tensor.dat", data_sk1.T, fmt="%.8e", header="F(Thz)          sk_xx       sk_yy         sk_zz")

    # Load and process modal kappa data
    modal_kappa_file = "outfile.grid_thermal_conductivity.hdf5"
    freq_n, mode_kxx, mode_kyy, mode_kzz, qp = load_and_process_modal_kappa_data(modal_kappa_file)

    # Create DataFrame for modal kappa
    data_mode_k = {
        'Freq': freq_n.ravel(),
        'mode_kxx': mode_kxx.ravel(),
        'mode_kyy': mode_kyy.ravel(),
        'mode_kzz': mode_kzz.ravel()
    }
    df = pd.DataFrame(data_mode_k)
    df = df.sort_values(by=['Freq'])

    # Make cumulation and save to cumulative_kappa_tensor.dat
    df['cumulative_kxx'] = np.cumsum(df['mode_kxx']) / qp
    df['cumulative_kyy'] = np.cumsum(df['mode_kyy']) / qp
    df['cumulative_kzz'] = np.cumsum(df['mode_kzz']) / qp

    if args.gnuplot:
        # Dump data to files for gnuplot
        dump_data_for_gnuplot(df)

    # Plot both cumulative and spectral kappa data in one figure
    plt.figure(figsize=(12, 6))

    plt.subplot(1, 2, 1)
    plot_cumulative_kappa_data(df, args.cumulative_xrange, args.cumulative_yrange, grid=args.gridon)

    plt.subplot(1, 2, 2)
    plot_spectral_kappa_data(freq, sk1, sk2, sk3, args.spectral_xrange, args.spectral_yrange, grid=args.gridon)

    plt.tight_layout()

    # Save the plot as a PDF file
    plt.savefig("modol_thermal_conductivity_plots.pdf")

    plt.show()


if __name__ == "__main__":
    main()
    
