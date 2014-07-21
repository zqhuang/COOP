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
  type(coop_healpix_maps)::map, imask, polmask, cutmap
  type(coop_file)::fp
  integer l, il, i
  COOP_REAL::beam_fwhm = 5.*coop_SI_arcmin
  COOP_REAL sigma
  COOP_INT, parameter::nside = 256
  call coop_random_init()
  call map%init(nside = nside, nmaps=1, spin = (/ 0 /) ) ! 2, 2 /))

  map%map = 0.
  call map%map2alm()

  map%cl(0:1,:) = 0.d0
  sigma = coop_sigma_by_fwhm * beam_fwhm
  call fp%open("cls.dat", "r")
  do l=2, map%lmax
     read(fp%unit, *) il, map%cl(l, coop_healpix_index_TT) !, map%cl(l, coop_healpix_index_EE), map%cl(l, coop_healpix_index_BB), map%cl(l, coop_healpix_index_TE)
     map%cl(l, :) = map%cl(l, :)*exp(-l*(l+1.d0)*sigma**2)
     if(il.ne.l) stop "Cls file wrong"
  enddo
  call fp%close()

  cutmap = map

  call imask%open("predx11/predx11_imask_nside"//trim(coop_num2str(nside))//".fits")
  call polmask%open("predx11/predx11_polmask_nside"//trim(coop_num2str(nside))//".fits")

  call map%simulate()
  call map%get_cls()
  cutmap = map
  call cutmap%mask(imask, polmask)
  call cutmap%map2alm()
  call cutmap%get_fullcls(imask, polmask)
  call fp%open("test.txt")
  do l = 2, map%lmax
     write(fp%unit, "(I5, 2E16.7)") l, cutmap%cl(l, :), map%cl(l, :)
  enddo
  call fp%close
     
end program test
