program plot
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_STRING::input, output
  COOP_INT::l
  input = coop_InputArgs(1)
  l = len_trim(input)
  if(l.le.0)then
     write(*,*) "./EBin filename"
  else
     if(l.ge.5)then
        if(input(l-3:l).eq.".cec")then
           output = input(1:l-4)
           call coop_binfile_decrypt(trim(input), trim(output))
        else
           output = trim(input)//".cec"
           call coop_binfile_encrypt(trim(input), trim(output))
        endif
     else
        output = trim(input)//".cec"
        call coop_binfile_encrypt(trim(input), trim(output))
     endif
  endif
  
end program plot
