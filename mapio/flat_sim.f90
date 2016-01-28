program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  character(LEN=*),parameter::mapdir = "act16/"
  COOP_INT,parameter::lmin = 200
  COOP_INT,parameter::lmax = 2500
  type(coop_fits_image_cea)::cf, noise
  COOP_REAL::fwhm_arcmin = 0.5
  COOP_REAL::Cls(lmin:lmax), ells(lmin:lmax), Cls2(lmin:lmax),junk
  COOP_INT::l, il
  type(coop_file)::fp
  call coop_random_init()
  call fp%open("planck14best_lensedCls.dat")
  do l=2, lmin-1
     read(fp%unit, *) il, junk
  enddo
  do l=lmin, lmax
     ells(l) = l
     read(fp%unit, *) il, Cls(l)
     if(il.ne.l) stop "error in cl file"
     Cls(l) = coop_2pi/l/(l+1.d0)*Cls(l)*coop_Gaussian_filter(fwhm_arcmin = fwhm_arcmin, l=l)**2
  enddo
  call fp%close()
  call cf%open(mapdir//"deep56_coadd_I.fits")
  call cf%simulate(lmin=lmin, lmax=lmax, Cls=Cls)
  call noise%open(mapdir//"deep56_array_2_noise_sim_1_I.fits")
  call cf%write(mapdir//"sim_1_I.fits")
  cf%image = cf%image + noise%image
  call cf%write(mapdir//"sim_with_noise_1_I.fits")

end program test
