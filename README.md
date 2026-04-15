# Spin Reconstruction and Spin Parameter Calculation Code
# spincorr_manga.f90

## Overview
This code is used for spin reconstruction and spin parameter calculation in the associated manuscript. 

---

## System Requirements

- Operating system: Linux (tested on Ubuntu 20.04)
- Fortran compiler: ifort
- No non-standard hardware is required

---

## Installation

---


## Instructions for use
### parameter setting
- boxsize = 500 Mpc/h
- grid number = 500
- R = 1.1 * r
- r = 0.2 - 15 Mpc/h

### pipeline
- read manga group position and rotate to ELUCID coordinate system
- read inverse displacement field and remapping group back to Lagrangian space
- read ELUCID IC and convert to potential field
- apply Gaussian smooth to potential field and derive smoothed tidal field
- calculate reconstructed spin and spin paramter at different smoothing scales

### subroutine 
- spinfield: tide field, reconstructed spin and spin paramter caculation
- gaussian_fourier_filter: smooth potential field at scale r

# Lagrangian Remapping Code
# idsp.f90: 
### pipeline
- read CUBE2 constrained simulation information
- derive particle q and x position
- calculate displacement field x-q
- get inverse displacement field at grid 
