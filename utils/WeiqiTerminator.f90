program GoTerminator
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  
  type(coop_weiqi_board)::board
  COOP_INT::i
  logical valid
  call board%init( blist = (/ 1, 19, 20, 21,  3 /), wlist=(/ 4, 22,  26, 45, 44, 43, 42, 41, 60, 59, 58  /) , turn = coop_weiqi_white)
  call board%print()
  call board%move(2)
  call board%print()  
  call board%move(3)
  call board%print()  
 ! valid = board%is_alive( coop_weiqi_black)  
 ! print*, valid
end program GoTerminator
