program simmaps
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use coop_fitsio_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools

  implicit none
#include "constants.h"
  COOP_REAL::fwhm_arcmin
  COOP_INT::nside
  type(coop_healpix_maps)::map
  COOP_STRING::clfile, fmap
  COOP_INT::hp_l1, hp_l2, i
  type(coop_list_integer)::listinds
  COOP_STRING::strinds
  COOP_INT,dimension(:),allocatable::index_list
  type(coop_cls)::cls
  if(iargc().lt. 2)then
     write(*,*) "Syntax:"
     write(*,*) "./SimuMpas -clfile ... -nside ... -out ... -fwhm_arcmin ... [-highpass_l1 ... -highpass_l2 ... -indlist ... ]"
  endif
  call coop_random_init()
  call coop_get_command_line_argument(key = "clfile", arg = clfile)
  call coop_get_command_line_argument(key = "indlist", arg = strinds, default='')
  call coop_get_command_line_argument(key = "out", arg = fmap)
  call coop_get_command_line_argument(key = "nside", arg = nside)
  call coop_get_command_line_argument(key = "fwhm_arcmin", arg = fwhm_arcmin )
  call coop_get_command_line_argument(key = "highpass_l1", arg = hp_l1, default = 0)
  call coop_get_command_line_argument(key = "highpass_l2", arg = hp_l2, default = hp_l1)
  call cls%load(clfile)
  if(trim(strinds) .ne. '')then
     call coop_string_to_list(strinds, listinds)
     allocate(index_list(listinds%n))
     do i=1, listinds%n
        index_list(i) = listinds%element(i)
     enddo
     call cls%select_maps(index_list)
  endif
  call cls%filter(fwhm_arcmin = fwhm_arcmin, highpass_l1 = hp_l1, highpass_l2 = hp_l2)
  call map%simulate(nside, cls)
  write(*,*) "Simulated map: "//trim(fmap)
  write(*,*) "NSIDE = "//COOP_STR_OF(map%nside)
  write(*,*) "NPIX = "//COOP_STR_OF(map%npix)
  write(*,*) "NMAPS = "//COOP_STR_OF(map%nmaps)
  call map%write(fmap)
  call map%free()
end program simmaps
