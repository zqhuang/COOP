program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  COOP_REAL fwhm_arcmin
  COOP_STRING::qmap, umap, emap, bmap
  type(coop_fits_image_cea)::q,u
  if(iargc().lt.2)then
     write(*,*) "Syntax:"
     write(*,*) "./FQU2EB -qmap QMAP -umap UMAP -emap EMAP -bmap BMAP"
     stop
  endif
  call coop_get_command_line_argument(key = "qmap", arg = qmap)
  call coop_get_command_line_argument(key = "umap", arg = umap)
  call coop_get_command_line_argument(key = "emap", arg = emap)
  call coop_get_command_line_argument(key = "bmap", arg = bmap)
  call q%open(qmap)
  call u%open(umap)
  call coop_fits_image_cea_QU2EB(q, u)
  call q%write(emap)
  call u%write(bmap)
end program test
