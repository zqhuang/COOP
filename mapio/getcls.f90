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

  COOP_STRING::imap_file, polmap_file, output_file
  COOP_INT, parameter::lmax = 2000
  COOP_INT, parameter::smooth_delta_ell = 20
  
  COOP_UNKNOWN_STRING, parameter::imask_file = "planck14/dx11_v2_common_int_mask_010a_1024.fits"
  COOP_UNKNOWN_STRING, parameter::polmask_file = "planck14/dx11_v2_common_pol_mask_010a_1024.fits"
  type(coop_healpix_maps)::map, imask, polmask
  type(coop_file)::fp
  integer l, i
  if(iargc().lt.2)then
     write(*,*) "./GetCl imap polmap"
     write(*,*) " or"
     write(*,*) "./GetCl imap polmap outputfile"
     stop
  endif
  imap_file = trim(coop_InputArgs(1))
  polmap_file = trim(coop_InputArgs(2))
  if(iargc() .ge. 3)then
     output_file = trim(coop_InputArgs(3))
  else
     output_file = coop_str_replace(imap_file, ".fits", "_cls.txt")
  endif
  call map%read(imap_file, nmaps_wanted = 3, spin = (/ 0, 2, 2 /) )
  call map%import(polmap_file, index_start = 2, index_end = 3, spin = (/ 2, 2 /))
  call map%map2alm(lmax = lmax)
  call map%get_cls()
  do l=0, map%lmax
     map%cl(l, :) = map%cl(l, :)*(l*(l+1)/coop_2pi*1.e12)
  enddo
  do i=1, 1
     call coop_smooth_data(map%lmax+1, map%cl(0:map%lmax, i), smooth_delta_ell)
  enddo
  call fp%open(output_file, "w")
  do l = 0, map%lmax
     write(fp%unit, "(I5, 6E16.7)") l, map%cl(l,:)
  enddo
  call fp%close()
  write(*,*) "now producing the lmax filtered map"
  call map%alm2map()
  write(*,*) "nmaps = map%nmaps", map%nmaps, size(map%map, 2)
  write(*,*) "npix = ", map%npix, size(map%map, 1)
  write(*,*) "isnan = ", coop_isnan(map%map)
  write(*,*) trim(coop_str_replace(output_file, ".txt", ""))//"_I.fits"
  call map%write(trim(coop_str_replace(output_file, ".txt", ""))//"_I.fits", index_list = (/ 1 /) )
  write(*,*) trim(coop_str_replace(output_file, ".txt", ""))//"_QU.fits"  
  call map%write(trim(coop_str_replace(output_file, ".txt", ""))//"_QU.fits", index_list = (/ 2, 3 /) )
end program test
