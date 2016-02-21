program flatop
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  COOP_INT,parameter::maxn = 10
  type(coop_fits_image_cea)::maps(maxn)
  COOP_STRING, dimension(maxn)::filename
  COOP_REAL::vals(maxn), ans
  COOP_STRING::output, expr
  COOP_INT::i, nmaps, j
  do i=1, maxn
     call coop_get_command_line_argument(key = 'file'//COOP_STR_OF(i), arg = filename(i), default = '' )
     if(trim(filename(i)).eq.'')then
        nmaps = i - 1
        exit
     endif
  enddo
  call coop_get_command_line_argument(key = 'out', arg = output, default = '' )
 if(trim(output) .eq. '')then
     write(*,*) "./FPTOP -file1 ... -file2 ... -out  OUTPUT"
     stop
  endif
  do i=1, nmaps
     call maps(i)%read(filename(i))
     write(*,*) "Reading "//trim(filename(i))
  enddo
  !!write the function here
  maps(1)%image = (maps(1)%image - maps(2)%image)/2.d0
  !!
  write(*,*) "writing to "//trim(output)
  call maps(1)%write(output)

end program flatop


