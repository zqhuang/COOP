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
  integer i
  real r, g, b
  COOP_LONG_STRING::str
  call fp%open("../include/color_table_planck.h")
  str = "p = Gradient(256"
  do i=1, 256
     read(fp%unit, *) r, g, b
     str = trim(str)//", rgb255("//trim(coop_num2str(r, "(F10.1)"))//", "//trim(coop_num2str(g, "(F10.2)"))//", "//trim(coop_num2str(b, "(F10.1)"))//") "
  enddo
  str = trim(str)//");"
  call fp%close()
  call fp%open("tmp.txt")
  write(fp%unit, "(A)") trim(str)
end program test
