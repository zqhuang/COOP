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
  COOP_UNKNOWN_STRING, parameter::map_file = "ffp7/ffp7_smica_cmb_0001_2048_debeam.fits"
!  COOP_UNKNOWN_STRING, parameter::map_file = "ffp7/ffp7_nobpm_smica_cmb_mc_0001_05a_2048_I.fits"
  COOP_UNKNOWN_STRING, parameter::imask_file = ""
  COOP_UNKNOWN_STRING, parameter::polmask_file = ""
  COOP_UNKNOWN_STRING, parameter::output_file = "ffp7_cls.dat"
  type(coop_healpix_maps)::map, imask, polmask
  type(coop_file)::fp
  integer l
  call map%read(map_file)
  if(imask_file .ne. "")then
     call imask%read(imask_file)
     map%map(:, 1) = map%map(:, 1) * imask%map(:, 1)
  endif
  if(map%nmaps .eq. 3)then
     if(polmask_file .ne. "")then
        call polmask%read(polmask_file)
        map%map(:, 2) = map%map(:, 2) * polmask%map(:, 1)
        map%map(:, 3) = map%map(:, 3) * polmask%map(:, 1)
     endif
  endif


  call fp%open(output_file)
  call map%map2alm()
  do l = 0, map%lmax
     write(fp%unit, "(I5, 6E16.7)") l, map%cl(l, :)
  enddo
  call fp%close

end program test
