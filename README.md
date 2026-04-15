# Spin Reconstruction and Lagrangian Remapping Code

## Overview
### spincorr_manga.f90
This code is used for spin reconstruction and spin parameter calculation.

### idsp.f90: 
This code is used for calculating inverse displacement.

---

## System Requirements
- Operating system: Linux (tested on Ubuntu 20.04)
- Fortran compiler: ifort
- FFT library (e.g. FFTW)
- No non-standard hardware is required

---

## Installation (Not required)
No installation is required. The code can be compiled and executed directly.

---

## Demo
### A typical workflow to run the full pipeline is as follows:

1. Generate particle data at redshift z = 0 from an N-body simulation using cosmological initial conditions (e.g. using CUBE2 or other publicly available N-body codes).

2. Run the Lagrangian remapping code to construct the inverse displacement field:
   make idsp.x
   ./idsp.x

3. Provide the input for the spin reconstruction code:
   - Target group positions at z=0
   - The inverse displacement field
   - The initial conditions used in the simulation

4. Run the spin reconstruction code:
   make spincorr.x
   ./spincorr.x

### The code will:
- Map Eulerian positions back to Lagrangian space
- Compute the smoothed tidal field from the initial conditions
- Reconstruct the spin
- Calculate the spin parameter

### Expected output
Reconstructed spin vectors and corresponding spin parameter values at different smoothing scales.

### Runtime 
Depend on the grid size and simulation resolution, typically ranging from minutes to hours on a standard desktop computer.

---

## Instructions For Use
### spincorr_manga.f90
#### Parameter setting
- boxsize = 500 Mpc/h
- grid number = 500
- R = 1.1 * r
- r = 0.2 - 15 Mpc/h

#### Pipeline
- Read manga group positions from 'manga_cat_group.bin' and rotate to ELUCID coordinate system;
- Using ELUCID weight array 'Weight5hR5_0500.bin' to judge which groups in the ELUCID reconstruction volumn;
- Read inverse displacement field from '0.000_xtoq_1.bin' and remapping groups back to Lagrangian space;
- Read ELUCID IC 'cxyz_251_500_500.bin' and convert to potential field;
- Apply Gaussian smooth to potential field and derive smoothed tidal field;
- Calculate reconstructed spin and spin paramter at different smoothing scales.

#### Subroutine 
- spinfield: tide field, reconstructed spin and spin paramter caculation.
- gaussian_fourier_filter: smooth potential field at scale r.

### idsp.f90: 
#### Pipeline
- Read CUBE2 constrained simulation output, which using the ELUCID IC 'cxyz_251_500_500.bin';
- Derive particle q and x position;
- Calculate displacement field x-q;
- Get inverse displacement field at grid.

---

## License
This code is released under the MIT License.
