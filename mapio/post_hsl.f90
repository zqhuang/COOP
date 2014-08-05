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
  type(coop_file)::fp
  integer i, nside, id, pix, iminprob
  type(coop_healpix_maps)::map
  COOP_REAL theta, phi, chisq, prob, minprob, l, b
  call map%init(nside = 4, nmaps=1, spin = (/ 0 /))
  map%map = 0.
  minprob = 10.
  iminprob = 0
  do i=0, 95
     call fp%open("predx11/predx111024_T_4id"//trim(coop_num2str(i))//".log", "r")
     read(fp%unit, *) nside, id, theta, phi, chisq, prob
     map%map(i, 1) = log10(prob+1.d-4)
     if(minprob .gt. prob)then
        minprob = prob
        iminprob = i
     endif
     theta = coop_pi - theta
     phi = coop_pi + phi
     call ang2pix_ring(map%nside, theta, phi, pix)
     map%map(pix, 1) = map%map(i, 1)
  enddo
  call map%write("stack_asym.fits")
  call pix2ang_ring(map%nside, iminprob, theta, phi)
  call coop_healpix_ang2lb(theta, phi, l, b)
  write(*,*) "min prob = ", minprob
  write(*,*) "direction l = ", nint(l), " b = ", nint(b)
end program test
