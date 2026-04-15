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
- Compile and execute spin reconstruction code using: make spincorr.x ; ./spincorr.x
- Compile and execute Lagrangian remapping code using: make idsp.x ; ./idsp.x

---

## Instructions for use
### spincorr_manga.f90
#### parameter setting
- boxsize = 500 Mpc/h
- grid number = 500
- R = 1.1 * r
- r = 0.2 - 15 Mpc/h

#### pipeline
- Read manga group positions from 'manga_cat_group.bin' and rotate to ELUCID coordinate system;
- Using ELUCID weight array 'Weight5hR5_0500.bin' to judge which groups in the ELUCID reconstruction volumn;
- Read inverse displacement field from '0.000_xtoq_1.bin' and remapping groups back to Lagrangian space;
- Read ELUCID IC 'cxyz_251_500_500.bin' and convert to potential field;
- Apply Gaussian smooth to potential field and derive smoothed tidal field;
- Calculate reconstructed spin and spin paramter at different smoothing scales.

#### subroutine 
- spinfield: tide field, reconstructed spin and spin paramter caculation.
- gaussian_fourier_filter: smooth potential field at scale r.

### idsp.f90: 
#### pipeline
- Read CUBE2 constrained simulation output, which using the ELUCID IC 'cxyz_251_500_500.bin';
- Derive particle q and x position;
- Calculate displacement field x-q;
- Get inverse displacement field at grid. 
