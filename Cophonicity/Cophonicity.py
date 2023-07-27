"""
Phonon State Moments Calculator
-------------------------------

Author: Tribhuwan Pandey

This code, calculates the moments of the phonon states based on the phonon density of states (DOS).
The phonon DOS data is loaded from a file named "outfile.phonon_dos", which contains information
about the phonon frequencies and their corresponding weights (densities of states) for different elements.

The moments are calculated by raising each phonon frequency to a specified moment order and then
multiplying it by its corresponding weight. The moments are then normalized by dividing them by the
total weight within a given frequency range.

The script calculates three different moments:
1. Overall moment: The moment calculated using all phonon frequencies in the specified range.
2. Moment for Mg: The moment calculated using the phonon frequencies and weights specific to the element Mg.
3. Moment for O: The moment calculated using the phonon frequencies and weights specific to the element O.

Functionality:
- Load phonon DOS data from a file.
- Calculate moments and normalization factors for the overall case, Mg, and O.
- Print the calculated moments for each case.

Usage:
    python phonon_state_moments.py

"""

import numpy as np
import sys


def load_phonon_dos(filename):
    """
    Load phonon density of states data from a file.

    Parameters:
        filename (str): Name of the file containing phonon DOS data.

    Returns:
        freq (numpy.ndarray): Array of phonon frequencies.
        gw (numpy.ndarray): Array of weights (DOS) corresponding to each phonon frequency.
        gw_Mg (numpy.ndarray): Array of weights (DOS) corresponding to the phonon frequencies of the element Mg.
        gw_O (numpy.ndarray): Array of weights (DOS) corresponding to the phonon frequencies of the element O.

    Raises:
        ValueError: If there is an issue loading data from the file.
    """
    try:
        data = np.loadtxt(filename)
        return data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    except OSError:
        raise ValueError("Failed to load data from the file.")

def calculate_moments(freq, weights, moment_order, fmin, fmax):
    """
    Calculate moments of the phonon states within a specified frequency range.

    Parameters:
        freq (numpy.ndarray): Array of phonon frequencies.
        weights (numpy.ndarray): Array of weights (DOS) corresponding to each phonon frequency.
        moment_order (int): Order of the moment calculation.
        fmin (float): Minimum frequency for the calculation range.
        fmax (float): Maximum frequency for the calculation range.

    Returns:
        float: The calculated moment divided by the normalization factor.
    """
    valid_indices = np.where((freq > fmin) & (freq < fmax))[0]
    moments = np.sum(freq[valid_indices] ** moment_order * weights[valid_indices])
    norm = np.sum(weights[valid_indices])
    return moments / norm if norm != 0 else 0.0

if __name__ == "__main__":
    if len(sys.argv) > 1 and (sys.argv[1] == "-h" or sys.argv[1] == "--help"):
        print(__doc__)  # Print the help message
    else:
        filename = "outfile.phonon_dos"
        freq, gw, gw_Mg, gw_O = load_phonon_dos(filename)

        moment_order = 1
        fmin, fmax = 0.01, np.amax(freq)

        moment_overall = calculate_moments(freq, gw, moment_order, fmin, fmax)
        moment_Mg = calculate_moments(freq, gw_Mg, moment_order, fmin, fmax)
        moment_O = calculate_moments(freq, gw_O, moment_order, fmin, fmax)

        print("Overall moment:", moment_overall)
        print("Moment for Mg:", moment_Mg)
        print("Moment for O:", moment_O)

