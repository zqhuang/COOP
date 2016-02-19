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
  call coop_get_command_line_argument(key = 'expr',  arg = expr, default= '' )
 if(trim(expr) .eq. '' .or. output .eq. '')then
     write(*,*) "./FPTOP -file1 ... -file2 ... -out  OUTPUT -expr EXPRESSION"
     stop
  endif
  do i=1, nmaps-1
     call maps(i)%read(filename(i))
  enddo
  do i=0, maps(1)%npix-1
     do j=1, nmaps
        vals(j) = maps(j)%image(i)
     enddo
     call coop_eval_math(expr, ans, vals)
     maps(1)%image(i) = ans
  enddo
  call maps(1)%write(output)

end program flatop


