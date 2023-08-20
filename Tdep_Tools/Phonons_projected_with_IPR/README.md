# 3. TDEP Phonon Dispersion Plotter with Inverse Participation Ratio (IPR)

This script reads a hdf5 file containing phonon band structure data from TDEP, calculates the Inverse Participation Ratio (IPR) for each phonon mode, and generates a plot with the phonon dispersion overlaid by IPR values.

### Note the IPR is defined as:

$$
\frac{1}{PR_{\boldsymbol{q}j}} = N_{\kappa}  \left(\sum_{\kappa}^{N_{\kappa}} \frac{|\boldsymbol{e}(\kappa;\boldsymbol{q}j)|^{2}}{M_{\kappa}}\right)^{2}
$$


$$
P_{\lambda}^{-1}=N \sum_{i}\left(\sum_{\alpha} e_{i \alpha, \lambda}^{*} e_{i \alpha, \lambda}\right)^{2}
$$

where *N* is the total number of atoms in primitive cell. and $e_{i \alpha, \lambda}$ is the αth eigenvector component of eigenmode λ for the ith atom. Each eigenmode λ is specified by the wave-vector **k** and branch index s. 


where *N* is the total number of atoms in primitive cell. and $e_{k; \boldsymbol{q}\lambda}$ is the $j^{th}$ eigenvector component of eigenmode $\lambda$ for the $k^{th}$ atom. Each eigenmode λ is specified by the wave-vector **k** and branch index j. 

## Usage
```shell
    python IPR_projected_phonon_tdep.py [--f FILENAME] [--cmap COLORMAP] [--alpha ALPHA] [--gnuplot] [--format FORMAT]
                                         [--width WIDTH] [--height HEIGHT] [--ymax YMAX]
