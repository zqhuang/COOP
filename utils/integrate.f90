module coop_integrate_mod
  use coop_wrapper_typedef
  implicit none


#include "constants.h"


  private

  public::coop_integrate, coop_polint

  interface coop_integrate
     module procedure coop_qrombc, coop_qromb_with_arguments
  end interface coop_integrate


contains


  function coop_qrombc(func,a, b, precision) result(integral)
#define QROMB_ARGUMENTS
#include "qromb.h"    
#undef QROMB_ARGUMENTS    
  end function coop_qrombc


  function coop_qromb_with_arguments(func, a, b, args, precision) result(integral)
   type(coop_arguments)::args
#define QROMB_ARGUMENTS  ,args
#include "qromb.h"
#undef QROMB_ARGUMENTS
  end function coop_qromb_with_arguments



  subroutine Coop_polint(xa,ya,n,x,y,dy)
    COOP_INT  :: n,NMAX
    COOP_REAL  :: dy,x,y,xa(n),ya(n)
    parameter (NMAX=10)
    COOP_INT  :: i,m,ns
    COOP_REAL  :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
    ns=1
    dif=ABS(x-xa(1))
    do i=1,n
       dift=ABS(x-xa(i))
       if (dift.LT.dif) then
          ns=i
          dif=dift
       endif
       c(i)=ya(i)
       d(i)=ya(i)
    enddo
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
             write(*,*) 'failure in polint'
             STOP
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       enddo
       if (2*ns.LT.n-m)then
          dy=c(ns+1)
       ELSE
          dy=d(ns)
          ns=ns-1
       endif
       y=y+dy
    enddo
    return
  end subroutine Coop_polint


end module coop_integrate_mod
