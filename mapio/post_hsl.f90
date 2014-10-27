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
  COOP_REAL theta, phi, chisq(0:95), prob(0:95), minprob, l, b
  call map%init(nside = 4, nmaps=1, spin = (/ 0 /))
  map%map = 0.
  minprob = 10.
  iminprob = 0
  do i=0, 95
     call fp%open("hsloutput/T_on_Tmax_log_4_"//trim(coop_num2str(i))//".txt", "r")
     read(fp%unit, *) nside, id, theta, phi, chisq(i), prob(i)
     map%map(i, 1) = log10(max(prob(i),1.d-4))
     if(minprob .gt. prob(i))then
        minprob = prob(i)
        iminprob = i
     endif
     theta = coop_pi - theta
     phi = coop_pi + phi
     call ang2pix_ring(map%nside, theta, phi, pix)
     call coop_healpix_ang2lb(theta, phi, l, b)
     map%map(pix, 1) = map%map(i, 1)
     print*, i, prob(i), nint(l), nint(b)
  enddo
  call map%write("hsl_T_on_Tmax.fits")
  call pix2ang_ring(map%nside, iminprob, theta, phi)
  call coop_healpix_ang2lb(theta, phi, l, b)
  write(*,*) "min prob = ", minprob
  write(*,*) "direction l = ", nint(l), " b = ", nint(b)
  call coop_asy_histogram(chisq, 10, "chisq_hist.txt")
  call coop_asy_histogram(log(max(prob, 1.d-4)), 10, "prob_hist.txt")
end program test
