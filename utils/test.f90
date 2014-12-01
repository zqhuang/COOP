program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT, parameter::n = 256, m = 5
  COOP_INT i
  COOP_REAL  r, g, b
  COOP_LONG_STRING::str
  type(coop_file)::fp, fpout
  call fp%open("colors.txt",'r')
  call fpout%open("planckcolor.txt", 'w')
  str = 'p = Gradient(256'
  do i=1, n
     read(fp%unit, *) r, g, b
     if(mod(i-1,m).eq.0)str = trim(str)//",rgbint("//COOP_STR_OF(nint(r))//","//COOP_STR_OF(nint(g))//","//COOP_STR_OF(nint(b))//")"
  enddo
  str = trim(str)//");"
  call fp%close()
  write(fpout%unit, "(A)") trim(str)
  call fpout%close()
end program Test
