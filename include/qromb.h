    external func
    COOP_REAL func
    COOP_INT  :: jmax,jmaxP,k, km
    COOP_REAL  :: a,b,integral, eps
    COOP_REAL, optional::precision
    parameter (jmax=25, jmaxP=jmax+1, k=5, km=k-1)
    !USES polint,trapzd
    COOP_INT  :: j
    COOP_REAL  :: dss,h(jmaxP),s(jmaxP)
    if(present(precision))then
       eps = precision
    else
       eps = 1.d-6
    endif
    h(1)=1.d0
    do j=1,jmax
       call trapzd_local(s(j),j)
       if (j.ge.K) then
          if(maxval(s(j-km:j))-minval(s(j-km:j)).le.eps*maxval(abs(s(j-km:j))) .or. maxval(abs(s(km:j))).lt.1.d-30)then
             integral=sum(s(j-km:j))/k
             return
          endif
          call coop_polint(h(j-KM),s(j-KM),K,0.d0,integral,dss)
          if (coop_isnan(integral) .or. ABS(dss) .lt. eps*abs(integral).or.ABS(dss).LT.1.d-31) return
       endif
       s(j+1)=s(j)
       h(j+1)=0.25D0*h(j)
    enddo

contains

  subroutine trapzd_local(ss, n)
    COOP_INT  :: n
    COOP_REAL  ::ss
    COOP_INT  :: jj, it
    COOP_REAL  :: del,thesum,tnm, x
    if (n.eq.1) then
       ss=0.5*(b-a)*(func(a QROMB_ARGUMENTS)+func(b  QROMB_ARGUMENTS))
    else
       it=ishft(1,n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5*del
       thesum=0.
       do jj=1,it
          thesum=thesum+func(x  QROMB_ARGUMENTS)
          x=x+del
       enddo
       ss=0.5D0*(ss+(b-a)*thesum/tnm)
    endif
    return
  end subroutine trapzd_local

