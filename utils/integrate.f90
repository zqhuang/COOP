module coop_integrate_mod
  use coop_wrapper_typedef
  implicit none


#include "constants.h"


  private

  public::coop_integrate

  interface coop_integrate
     module procedure coop_qrombc, coop_qromb_with_arguments
  end interface coop_integrate


contains


  function coop_qrombc(func,a, b, precision)
    COOP_INT  :: JMAX,JMAXP,K,KM
    COOP_REAL  :: a,b,func, eps, coop_qrombc
    COOP_REAL, optional::precision
    external func
    parameter (JMAX=25, JMAXP=JMAX+1, K=5, KM=K-1)
    !USES polint,trapzd
    COOP_INT  :: j
    COOP_REAL  :: dss,h(JMAXP),s(JMAXP)
    if(present(precision))then
       eps = precision
    else
       eps = 1.d-6
    endif
    h(1)=1.d0
    do j=1,JMAX
       call D_trapzd(func,a,b,s(j),j)
       if (j.ge.K) then
          if(maxval(s(j-km:j))-minval(s(j-km:j)).le.eps*maxval(abs(s(j-km:j))) .or. maxval(abs(s(km:j))).lt.1.d-40)then
             coop_qrombc=sum(s(j-km:j))/k
             return
          endif
          call d_polint(h(j-KM),s(j-KM),K,0.d0,coop_qrombc,dss)
          if (coop_isnan(coop_qrombc).or.ABS(dss).LE.EPS*ABS(coop_qrombc).or.ABS(dss).LT.1.d-31) return
       endif
       s(j+1)=s(j)
       h(j+1)=0.25D0*h(j)
    enddo
  end function coop_qrombc


  subroutine D_trapzd(func,a,b,s,n)
    COOP_INT  :: n
    COOP_REAL  :: a,b,s,func
    external func
    COOP_INT  :: j,it
    COOP_REAL  :: del,sum,tnm,x
    if (n.eq.1) then
       s=0.5*(b-a)*(func(a)+func(b))
    ELSE
       it=ishft(1,n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5*del
       sum=0.
       do j=1,it
          sum=sum+func(x)
          x=x+del
       enddo
       s=0.5D0*(s+(b-a)*sum/tnm)
    endif
    return
  end subroutine D_trapzd

  subroutine D_polint(xa,ya,n,x,y,dy)
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
  end subroutine D_polint


  function coop_qromb_with_arguments(func, a, b, args, precision) result(integral)
    COOP_INT  :: JMAX,JMAXP,K,KM
    COOP_REAL  :: a,b,integral, eps
    COOP_REAL, optional::precision
    type(coop_arguments) args
    external func
    COOP_REAL func
    parameter (JMAX=25, JMAXP=JMAX+1, K=5, KM=K-1)
    !USES polint,trapzd
    COOP_INT  :: j
    COOP_REAL  :: dss,h(JMAXP),s(JMAXP)
    if(present(precision))then
       eps = precision
    else
       eps = 1.d-6
    endif
    h(1)=1.d0
    do j=1,JMAX
       call trapzd_with_args(args,func,a,b,s(j),j)
       if (j.ge.K) then
          if(maxval(s(j-km:j))-minval(s(j-km:j)).le.eps*maxval(abs(s(j-km:j))) .or. maxval(abs(s(km:j))).lt.1.d-30)then
             integral=sum(s(j-km:j))/k
             return
          endif
          call d_polint(h(j-KM),s(j-KM),K,0.d0,integral,dss)
          if (coop_isnan(integral) .or. ABS(dss) .lt. eps*abs(integral).or.ABS(dss).LT.1.d-31) return
       endif
       s(j+1)=s(j)
       h(j+1)=0.25D0*h(j)
    enddo
  end function coop_qromb_with_arguments


  subroutine trapzd_with_args(args,func,a,b,s,n)
    type(coop_arguments) args
    COOP_INT  :: n
    COOP_REAL  :: a,b,s
    external func
    COOP_REAL func
    COOP_INT  :: j,it
    COOP_REAL  :: del,sum,tnm,x
    if (n.eq.1) then
       s=0.5*(b-a)*(func(a, args)+func(b, args))
    ELSE
       it=ishft(1,n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5*del
       sum=0.
       do j=1,it
          sum=sum+func(x, args)
          x=x+del
       enddo
       s=0.5D0*(s+(b-a)*sum/tnm)
    endif
    return
  end subroutine Trapzd_with_args


  

end module coop_integrate_mod
