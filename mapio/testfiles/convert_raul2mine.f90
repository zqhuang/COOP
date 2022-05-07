program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  implicit none
#include "constants.h"
  COOP_INT, parameter::n = 100
  COOP_REAL:: f(n, n)
  COOP_STRING::inp, out
  type(coop_file)::fp
  COOP_INT i, j
  if(iargc() .lt. 1) stop "./Test Input Output"
  inp = coop_InputArgs(1)
  out = coop_InputArgs(2)
  if(trim(out) .eq. "" )then
     out = coop_str_replace(inp, ".dat", ".txt")
  endif
  call fp%open(inp, "r")
  do i = 1, n
     read(fp%unit, *) f(:, i)
  enddo
  call fp%close()
  call fp%open(out, "w")
  do i = 11, 90
     write(fp%unit, "(80E16.7)") f(11:90, i)*1.d6
  enddo
  call fp%close()
end program test
