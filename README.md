# spincorr_manga.f90: 
## spin reconstruction and spin parameter calculation code
### parameter setting
- boxsize = 500 Mpc/h
- grid number = 500
- R = 1.1*r
- r = 0.2-15 Mpc/h

### pipeline
- read manga group position and rotate to ELUCID coordinate system
- read inverse displacement field and remapping group back to Lagrangian space
- read ELUCID IC and convert to potential field
- apply Gaussian smooth to potential field and derive smoothed tidal field
- calculate reconstructed spin and spin paramter at different smoothing scales

### subroutine 
- spinfield: tide field, reconstructed spin and spin paramter caculation
- gaussian_fourier_filter: smooth potential field at scale r

# idsp.f90: 
## Lagrangian remapping code
### pipeline
- read CUBE2 constrained simulation information
- derive particle q and x position
- calculate displacement field x-q
- get inverse displacement field at grid 
