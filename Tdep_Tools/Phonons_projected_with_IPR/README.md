# 3. TDEP Phonon Dispersion Plotter with Inverse Participation Ratio (IPR)

This script reads a hdf5 file containing phonon band structure data from TDEP, calculates the Inverse Participation Ratio (IPR) for each phonon mode, and generates a plot with the phonon dispersion overlaid by IPR values.

### IPR is defined as:

$$
\frac{1}{PR_{\boldsymbol{q}j}} = N_{\kappa}  \left(\sum_{\kappa}^{N_{\kappa}} \frac{|\boldsymbol{e}(\kappa;\boldsymbol{q}j)|^{2}}{M_{\kappa}}\right)^{2}
$$


where *N* is the total number of atoms in primitive cell. and $e_{i \alpha, \lambda}$ is the αth eigenvector component of eigenmode λ for the ith atom. Each eigenmode λ is specified by the wave-vector **k** and branch index j.  Note that the eigenvectors are already normalized with respect to mass, therefore mass term is ignored in the code.

## Usage
```shell
    python IPR_projected_phonon_tdep.py [--f FILENAME] [--cmap COLORMAP] [--alpha ALPHA] [--gnuplot] [--format FORMAT]
                                         [--width WIDTH] [--height HEIGHT] [--ymax YMAX]
```
