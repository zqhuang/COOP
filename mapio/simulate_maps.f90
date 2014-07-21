program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools

  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, imask, polmask
  integer l, m, il
  type(coop_file)::fp
  COOP_REAL::beam_fwhm = 5.*coop_SI_arcmin
  COOP_REAL sigma
  call coop_random_init()
  call map%init(nside = 1024, nmaps=3, spin = (/ 0, 2, 2 /))
  call map%map2alm()
  sigma = coop_sigma_by_fwhm * beam_fwhm
  call fp%open("cls.dat", "r")
  print*, "lmax = ",map%lmax
  do l=2, map%lmax
     read(fp%unit, *) il, map%cl(l, coop_healpix_index_TT), map%cl(l, coop_healpix_index_EE), map%cl(l, coop_healpix_index_BB), map%cl(l, coop_healpix_index_TE)
     map%cl(l, :) = map%cl(l, :)*exp(-l*(l+1.d0)*sigma**2)
  enddo
  print*, map%cl(2:10, 1)
  call map%simulate()
  call map%map2alm()
  print*, map%cl(2:10, 1)
  call map%write("sim2/sim2_iqu_nside1024.fits")

end program test
