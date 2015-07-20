program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map
  COOP_STRING::fname, output
  COOP_SHORT_STRING::field, unit
  COOP_INT::i
  if(iargc().lt. 4)then
     write(*,*) "./CCards -inp inputmap -out outputmap -f1 ... -f2 ... -u1 ... -u2 ..."
     stop
  endif
  call coop_get_command_line_argument(key = 'inp', arg = fname)
  call coop_get_command_line_argument(key = 'out', arg = output)  
  call map%read(fname)
  do i= 1, map%nmaps
     call coop_get_command_line_argument(key = 'f'//COOP_STR_OF(i),arg =  field, default = "")
     if(trim(field).ne."")then
        call map%set_field( i, field)
     endif
     call coop_get_command_line_argument(key = 'u'//COOP_STR_OF(i),arg =  unit, default = "")
     if(trim(unit).ne."")then
        call map%set_unit( i, unit)
     endif
     
  enddo
  call map%fields_to_spins()
  call map%write(trim(adjustl(output)))
end program test
