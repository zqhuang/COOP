program test
  use coop_wrapper_utils
  use coop_expint_mod
  implicit none
#include "constants.h"
  character(LEN=*),parameter::postfix = "_h20_s2.txt"
  character(LEN=*),parameter::postfix2 = "_h20_s3.txt"
  COOP_REAL::r, r2
  COOP_INT,parameter::ntab = 10000
  COOP_REAL::tab_lnl(ntab), tab_lnh(ntab), tab_lnh2(ntab)
  COOP_REAL::tab2_lnl(ntab), tab2_lnh(ntab), tab2_lnh2(ntab)
  COOP_REAL::t, tmp1, tmp2
  COOP_INT::i, n
  open(11, file="tab"//postfix)
  read(11, *) r
  do i = 1, ntab
     read(11, *) tmp1, tmp2, tab_lnl(i), tab_lnh(i)
  enddo
  close(11)

  open(11, file="tab"//postfix2)
  read(11, *) r2
  do i = 1, ntab
     read(11, *) tmp1, tmp2, tab2_lnl(i), tab2_lnh(i)
  enddo
  close(11)

  
  call coop_spline(ntab, tab_lnl, tab_lnh, tab_lnh2, -1.5d0, 0.d0)
  call coop_spline(ntab, tab2_lnl, tab2_lnh, tab2_lnh2, -1.5d0, 0.d0)
  
  do i=0, 10000
     t = i*0.005
     write(*, *) t, h_of_t(t) - h2_of_t(t)
  enddo
  
contains

  function h_of_t(t) result(h)
    COOP_REAL::lnl, rl, h, t
    COOP_INT::il
    lnl = -2.d0*abs(t) - log(r)
    if(lnl .lt. tab_lnl(1))then
       h = tab_lnh(1) + 1.5*(lnl - tab_lnl(1))
    else
       call coop_splint(ntab, tab_lnl, tab_lnh, tab_lnh2, lnl, h)
    endif
    h = - exp(h) /sqrt(r)
  end function h_of_t


  function h2_of_t(t) result(h)
    COOP_REAL::lnl, rl, h, t
    COOP_INT::il
    lnl = -2.d0*abs(t) - log(r)
    if(lnl .lt. tab2_lnl(1))then
       h = tab2_lnh(1) + 1.5*(lnl - tab2_lnl(1))
    else
       call coop_splint(ntab, tab2_lnl, tab2_lnh, tab2_lnh2, lnl, h)
    endif
    h = - exp(h) /sqrt(r)
  end function h2_of_t

end program test
