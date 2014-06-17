module coop_threej_mod
  use coop_wrapper_typedef
  use coop_jl_mod
  use coop_special_function_mod
  implicit none
  private

#include "constants.h"

  public::coop_threej000, coop_threej_fixed_m1, coop_threej_table, coop_threej_generate_table


  COOP_INT ,parameter::dl = kind(1.d0)
  COOP_INT ,parameter::sp = kind(1.)

  type coop_threej_table
     COOP_INT  ntheta
     COOP_REAL  dtheta
     COOP_REAL ,dimension(:),allocatable::mu
     real(sp),dimension(:,:,:),allocatable::s
   contains
     procedure::free => coop_threej_table_free
  end type coop_threej_table


  interface is_triangle
     module procedure is_triangle_real, is_triangle_int
  end interface is_triangle

contains

  function Coop_threej000(l1,l2,l3) result(w3j)
    COOP_REAL  w3j
    COOP_INT  l1,l2,l3, J, g
    J = (l1+l2+l3)
    if(l1+l2.lt.l3 .or. l2+l3.lt.l1 .or. l3+l1.lt.l2 .or. mod(J,2).ne.0)then
       w3j = 0._dl
       return
    endif
    g = J/2
    w3j = dexp(&
         log_gamma(g+1.d0)-log_gamma(g-l1+1.d0) &
         -log_gamma(g-l2+1.d0) - log_gamma(g-l3+1.d0) &
         +( &
         log_gamma(J-2*l1+1.d0) + log_gamma(J-2*l2+1.d0) &
         + log_gamma(J-2*l3+1.d0) - log_gamma(J+2.d0) &
         )/2.d0 &
         )
    if(mod(g,2).ne.0) w3j = -w3j
  end function Coop_threej000

  


!!$for large l' you need to use limber approximation
  function limapp_int(l,func)
    !!return limber approximation: \int_0^\infty func(x) j_l(x) dx
    external func
    COOP_REAL func
    COOP_INT  l
    COOP_REAL  limapp_int
    limapp_int = sqrt(coop_pio2/(l+0.5_dl))*func(l+0.5_dl)
  end function limapp_int


  !!return normalized P_l^m(x) := sqrt(4Pi/(2l+1))Y_l^m(arccos(x),0)
  !!this works for arbitrary l,m
  !!but it is too slow if you want a lot of P_l^m s
  function Normalized_Plm_recur(l,m,x) result(Plm)
    COOP_REAL  x, Plm_prev, Plm, Plm_next, sca
    COOP_INT  l,m, i
    logical need_norm
    if(m.lt.0 .or. m.gt.l .or. dabs(x).gt.1.d0)then
       stop "Invalid argument in Normalized_Plm_recur"
       return
    endif
    if(m.eq.0 .or. l.le.10)then
       Plm = coop_Normalized_Associated_Legendre(l,m,x)
       return
    endif
    if(m.eq.l)then
       Plm = dexp(lngamma_hi_table(m+1)-lngamma_hi_table(1)-log_gamma(2*m+1.d0)/2.d0 + m*(coop_ln2+lnabs(1.d0-x**2)/2.d0))
       if(mod(l,2).ne.0) Plm = -Plm
       return
    endif
    sca = lngamma_hi_table(m+1)-lngamma_hi_table(1)-log_gamma(2*m+1.d0)/2.d0 + m*(coop_ln2+lnabs(1.d0-x**2)/2.d0)
    if(sca .lt. -1.e4_dl)then
       Plm = 0.d0
       return
    endif
    if(sca .gt. -20._dl)then
       if(mod(m,2).eq.0)then
          Plm_prev = dexp(sca)
       else
          Plm_prev = -dexp(sca)
       endif
       Plm =  sqrt(2*m+1.d0)*x*Plm_prev
       do i=m+2,l
          Plm_next = ((2*i-1)*x*Plm-sqrt((i-m-1.d0)*(i+m-1.d0))*Plm_prev)/(sqrt(dble(i-m)*(i+m)))
          Plm_prev = Plm
          Plm = Plm_next
       enddo
    else
       if(mod(m,2).eq.0)then
          Plm_prev = 1.d0
       else
          Plm_prev = -1.d0
       endif
       need_norm = .true.
       Plm =  sqrt(2*m+1.d0)*x*Plm_prev
       do i=m+2,l
          Plm_next = ((2*i-1)*x*Plm-sqrt((i-m-1.d0)*(i+m-1.d0))*Plm_prev)/(sqrt(dble(i-m)*(i+m)))
          Plm_prev = Plm
          Plm = Plm_next
          if(need_norm)then
             if(dabs(Plm).gt.1.e20_dl)then
                if(sca .ge. -30._dl)then
                   call multiply_exp(Plm, sca) 
                   call multiply_exp(Plm_prev, sca) 
                   need_norm = .false.
                else
                   sca = sca + 30._dl
                   call multiply_exp(Plm, -30._dl)
                   call multiply_exp(Plm_prev, -30._dl)
                endif
             endif
          endif
       enddo
       if(need_norm)call  multiply_exp(Plm,sca)
    endif
  end function Normalized_Plm_recur

  !!this only works for small l, m
  subroutine get_all_Normalized_Plms(lmax, mmax, x, Plms)
    COOP_REAL  x, lnsx
    COOP_INT  lmax, mmax
    COOP_REAL  Plms(0:lmax, 0:mmax)
    COOP_INT  l,m
    if(dabs(x).gt.1.d0 .or. m .gt. 100 .or. lmax .gt. 1000) stop "wrong argument in get_all_normalized_Plms"
    Plms = 0.d0
    if(dabs(x).eq.1.d0)then
       Plms(0:lmax:2,0) = 1.d0
       Plms(1:lmax:2,0) = x
    endif
    lnsx = dlog(1.d0-x**2)/2.d0
    !$omp parallel do private(m, l)
    do m=0, min(lmax,mmax)
       if(mod(m,2).eq.0)then
          Plms(m,m)= dexp(lngamma_hi_table(m+1)-lngamma_hi_table(1)-log_gamma(2*m+1)/2.d0+m*(coop_ln2+lnsx))
       else
          Plms(m,m) = - dexp(lngamma_hi_table(m+1)-lngamma_hi_table(1)-log_gamma(2*m+1)/2.d0+m*(coop_ln2+lnsx))
       endif
       if(m .lt. lmax)then
          Plms(m+1,m) = sqrt(2*m+1.d0)*x*Plms(m,m)
          do l=m+2,lmax
             Plms(l, m) = ((2*l-1)*x*Plms(l-1,m)-sqrt((l-m-1.d0)*(l+m-1.d0))*Plms(l-2,m))/(sqrt(dble(l-m)*(l+m)))
          enddo
       endif
    enddo
    !$omp end parallel do
  end subroutine get_all_Normalized_Plms
  


  !!for a given l, obtain all P_l^m(x) for -l<=m<=l;
  !!output: Plms(-l:l)
  subroutine get_all_normalized_Plms_l(l, x, Plms, mmax)
    COOP_REAL  x, A, B, sqrtA, nu, Bmax, xsqrtA, sqrtB, sqrtBmB, Aplus, Aminus
    COOP_INT  l,m, m0, dm0, mstart
    COOP_INT ,optional::mmax
    COOP_INT  mend
    COOP_REAL  Plms(-l:l)
    if(dabs(x).gt.1.d0 .or. l.lt.0) Stop "Invalid arguments in get_all_normalized_Plms_l"
    if(present(mmax))then
       mend = min(mmax, l)
    else
       mend= l
    endif
    Plms(mend+1:l) = 0
    if(l .le. 50)then
       !$omp parallel do
       do m=0,mend
          Plms(m) = Normalized_Plm_recur(l,m,x)
       enddo
       !$omp end parallel do
    else
       if(dabs(x).lt. 0.9d0)then
          Plms(0) = Legendre_P(l,x)
          Plms(1) = Normalized_Plm_recur(l,1,x)
          mstart = 2
       elseif(dabs(x) .lt. 0.99)then
          Plms(0) = Legendre_P(l,x)
          Plms(1) = Normalized_Plm_recur(l,1,x)
          Plms(2) = Normalized_Plm_recur(l,2,x)
          Plms(3) = Normalized_Plm_recur(l,3,x)
          mstart = 4
       elseif(dabs(x) .lt. 0.995)then
          !$omp parallel do
          do m=0, 20
             Plms(m) = Normalized_Plm_recur(l,m,x)
          enddo
          !$omp end parallel do
          mstart = 21 
       else
          !$omp parallel do
          do m=0, 50
             Plms(m) = Normalized_Plm_recur(l,m,x)
          enddo
          !$omp end parallel do
          mstart = 51 
       endif
       if(mstart .le. mend)then
          A = l*(l+1.d0)
          sqrtA = dsqrt(A)
          xsqrtA  = -x*sqrtA
          Aplus = A*(1.d0+x)
          Aminus = A*(1.d0-x)
          Bmax = (1.d0-x**2)*A
          m0 = nint(dsqrt(Bmax+1.d0))
          dm0 = max(ceiling(l*0.032),30)
          Plms(max(m0+dm0+1,mstart):mend) = 0._dl
          !$omp parallel do 
          do m=max(mstart,m0-dm0+1), min(mend, m0+dm0)
             Plms(m) = Normalized_Plm_recur(l,m,x)
          enddo
          !$omp end parallel do
          !$omp parallel do private(B, nu, sqrtB, sqrtBmB)
          do m = mstart, min(m0-dm0, mend)
             B = dble(m*m - 1)
             sqrtB= dsqrt(B)
             sqrtBmB = dsqrt(Bmax-B)
             nu = sqrtB*sqrtBmB
             Plms(m) =  dsin(sqrtA * datan(xsqrtA/sqrtBmB) &
                  + sqrtB/2.d0*(datan((Aplus-B)/nu) - datan((Aminus-B)/nu))+(l+m+1)*coop_pio2)/(1.d0-x**2/(1.d0-B/A))**0.25d0
             if(mod(l+m,2).eq.0)then
                Plms(m) = Plms(m)*dexp((log_gamma(l-m+1)+log_gamma(l+m+1))/2.d0-l*coop_ln2 - log_gamma((l-m)/2+1)-log_gamma((l+m)/2+1))
             else
                Plms(m) = Plms(m)*dexp((log_gamma(l-m+1)+log_gamma(l+m+1))/2.d0-l*coop_ln2 - lngamma_hi_table((l-m)/2+2)-lngamma_hi_table((l+m)/2+2))
             endif
          enddo
          !$omp end parallel do
       endif
    endif
    Plms(-1:-l:-2) = -Plms(1:l:2)
    Plms(-2:-l:-2) = Plms(2:l:2)
  end subroutine get_all_normalized_Plms_l

  !!for given m calculate all P_l^m(x), with l running from m to lmax
  !!outpu Plms(m:lmax)
  subroutine get_all_normalized_Plms_m(lmax, m, x, Plms)
    COOP_INT  lmax, m, l
    COOP_REAL  x, Plms(m:lmax), sca
    logical need_norm
    if(dabs(x).gt.1.d0 .or. m.lt.0)stop "Invalid argument in get_all_normalized_Plms,m"
    if(dabs(x).eq.1.d0)then
       if(m.eq.0)then
          Plms(0:lmax:2) = 1.d0
          Plms(1:lmax:2) = x
       else
          Plms = 0.d0
       endif
       return
    endif
    sca = lngamma_hi_table(m+1)-lngamma_hi_table(1)-log_gamma(2*m+1)/2.d0+m*(coop_ln2+dlog(1.d0-x**2)/2.d0)
    if(sca.lt.-1.e4_dl) then
       Plms = 0.d0
       return
    endif
    if(sca.gt.-20._dl .or. m .eq. lmax)then       
       if(mod(m,2).eq.0)then
          Plms(m) = dexp(sca)
       else
          Plms(m) = -dexp(sca)
       endif
       if(m .lt. lmax)then
          Plms(m+1) = sqrt(2*m+1)*x*Plms(m)
          do l=m+2,lmax
             Plms(l) = ((2*l-1)*x*Plms(l-1)-sqrt(l-m-1)*sqrt(l+m-1)*Plms(l-2))/(sqrt(l-m)*sqrt(l+m))
          enddo
       endif
    else
       if(mod(m,2).eq.0)then
          Plms(m) = 1.d0
       else
          Plms(m) = -1.d0
       endif
       need_norm = .true.
       Plms(m+1) = sqrt(2*m+1)*x*Plms(m)
       do l=m+2,lmax
          Plms(l) = ((2*l-1)*x*Plms(l-1)-sqrt(l-m-1)*sqrt(l+m-1)*Plms(l-2))/(sqrt(l-m)*sqrt(l+m))
          if(need_norm)then
             call multiply_exp(Plms(l-2), sca)
             if(dabs(Plms(l)).gt.1.e20_dl)then
                if(sca.gt.-30._dl)then
                   call multiply_exp(Plms(l-1),sca)
                   call multiply_exp(Plms(l),sca)
                   need_norm = .false.
                else
                   sca = sca+30._dl
                   call multiply_exp(Plms(l-1), -30._dl)
                   call multiply_exp(Plms(l), -30._dl)
                endif
             endif
          endif
       enddo
       if(need_norm)then
          call multiply_exp(Plms(lmax-1), sca)
          call multiply_exp(Plms(lmax), sca)
       endif
    endif
  end subroutine get_all_normalized_Plms_m

  subroutine Accurate_normalized_Plms(lmax, x, Plms, phiisPi)
    COOP_INT  lmax
    COOP_REAL  x
    COOP_REAL  Plms(0:lmax, -lmax:lmax)
    COOP_INT  m
    logical,optional::phiIsPi
    !$omp parallel do
    do m=0,lmax
       Plms(0:m-1, m) = 0._dl
       call get_all_normalized_Plms_m(lmax, m, x, Plms(m:lmax, m))
       if(present(PhiIsPi))then
          if(PhiIsPi .and. mod(m,2).ne.0)then
             Plms(m:lmax,m) = -Plms(m:lmax,m)
          end if
       endif
    enddo
    !$omp end parallel do
    do m=1, lmax, 2
       Plms(0:lmax, -m) = -Plms(0:lmax, m)
    enddo
    do m=2, lmax, 2
       Plms(0:lmax, -m) = Plms(0:lmax, m)
    enddo
  end subroutine Accurate_normalized_Plms

  !!return ThreeJSymbol[{j1,m1},{j2,m},{j3,-m1-m}]
  !!m running from -j2 to j2;
  !!The input: COOP_INT  j1>=0,j2>=0,j3 >=0 and j2<=j3; COOP_INT  m1>=0
  !!           j1+j2>=j3 and j2+j3>=j1 and j3+j1>=j2
  subroutine coop_threej_fixed_m1(j1,j2, j3, m1, threej)
    COOP_INT  j1,j2,j3,m1,m, totjsq
    COOP_REAL  threej(-j2:j2)
    COOP_INT  mmax, mmid
    COOP_REAL  sca
    logical need_norm
    if(j1+j2.lt.j3 .or. j2+j3.lt.j1 .or. j1+j3.lt.j2 .or. m1.gt.j1 .or. m1.lt.0 .or. j2.gt.j3) then
       stop "invalid arguments in coop_threej_fixed_m1"
    endif
    mmax = min(j2, j3-m1)
    threej(mmax+1:j2) = 0._dl
    mmid = min(0, mmax)
    sca = (log_gamma(2*j2+1)+log_gamma(j1+j3-j2+1)+log_gamma(j2+j3-m1+1)+log_gamma(j1+m1+1)-log_gamma(j1+j2-j3+1)-log_gamma(j2+j3-j1+1)-log_gamma(j1-m1+1)-log_gamma(j3+m1-j2+1)-log_gamma(j1+j2+j3+2))/2.d0
    totjsq = j2*(j2+1)+j3*(j3+1)-j1*(j1+1)
    if(mod(j1+m1,2).eq.0)then
       threej(-j2) = 1.d0
    else
       threej(-j2) = -1.d0
    endif
    if(sca .ge. -20._dl .or. -j2+1.ge.mmid)then
       call multiply_exp(threej(-j2),sca)
       need_norm = .false.
    else
       need_norm = .true.
    endif
    if(-j2.lt.mmid)then
       threej(-j2+1) = -coef_D(-j2)/coef_C(-j2+1)*threej(-j2)
       do m=-j2+2, mmid
          threej(m) = (-coef_D(m-1)*threej(m-1)-coef_C(m-1)*threej(m-2))/coef_C(m)
          if(need_norm)then
             call multiply_exp(threej(m-2), sca)
             if(dabs(threej(m)).gt. 1.e30_dl)then
                if(sca .ge. -30._dl)then
                   call multiply_exp(threej(m-1),sca)
                   call multiply_exp(threej(m),sca)
                   need_norm = .false.
                else
                   sca = sca+30._dl
                   call multiply_exp(threej(m-1),-30._dl)
                   call multiply_exp(threej(m),-30._dl)
                endif
             endif
          endif
       enddo
       if(need_norm)then
          call multiply_exp(threej(mmid-1), sca)
          call multiply_exp(threej(mmid),sca)
       endif
    endif
    if(mmid.eq.mmax)  return
    if(mmax.eq.j2)then
       sca = (log_gamma(2*j2+1)+log_gamma(j1+j3-j2+1)+log_gamma(j1-m1+1)+log_gamma(j2+j3+m1+1)-log_gamma(j1+j2-j3+1)-log_gamma(j2+j3-j1+1)-log_gamma(j3-j2-m1+1)-log_gamma(j1+m1+1)-log_gamma(j1+j2+j3+2))/2.d0
    else !!mmax = j3-m1
       sca = (log_gamma(j1+j2-j3+1)+log_gamma(2*j3+1)+log_gamma(j2+j3-m1+1) + log_gamma(j1+m1+1)-log_gamma(j1+j3-j2+1) - log_gamma(j2+j3-j1+1)-log_gamma(j2+j3-j1+1)-log_gamma(j1-m1+1)-log_gamma(j2+m1-j3+1)-log_gamma(j1+j2+j3+2))/2.d0
    endif
    if(mod(j2+j3+m1,2).eq.0)then
       threej(mmax) = 1.d0
    else
       threej(mmax) = -1.d0
    endif
    if(sca.ge.-20._dl .or. mmax-2.le.mmid)then
       call multiply_exp(threej(mmax),sca)
       need_norm = .false.
    else
       need_norm = .true.
    endif
    threej(mmax-1) = -coef_D(mmax)/coef_C(mmax)*threej(mmax)
    do m=mmax-2,mmid+1, -1
       threej(m) = (-coef_C(m+2)*threej(m+2)-coef_D(m+1)*threej(m+1))/coef_C(m+1)
       if(need_norm)then
          call multiply_exp(threej(m+2), sca)
          if(dabs(threej(m)).gt.1.e30_dl)then
             if(sca.ge.-30._dl)then
                call multiply_exp(threej(m+1),sca)
                call multiply_exp(threej(m),sca)
                need_norm = .false.
             else
                sca = sca+30._dl
                call multiply_exp(threej(m+1),-30._dl)
                call multiply_exp(threej(m),-30._dl)
             endif
          endif
       endif
    enddo
    if(need_norm)then
       call multiply_exp(threej(mmid+2),sca)
       call multiply_exp(threej(mmid+1),sca)
    endif
  contains
    
    function coef_C(m)
      COOP_REAL  coef_C
      COOP_INT  m
      coef_C = sqrt((j2-m+1.d0)*(j2+m)*(j3-m-m1+1.d0)*(j3+m+m1))
    end function coef_C


    function coef_D(m)
      COOP_REAL  coef_D
      COOP_INT  m
      coef_D = totjsq - 2*m*(m+m1)
    end function coef_D

  end subroutine coop_threej_fixed_m1


  !!this returns ThreeJSymbol[{j1,m},{j2,-m-m3},{j3,m3}] (m=-j1,-j1+1,...,j1) in an array  threej(-j1:j1)
  !!mmin and mmax are returned such that the 3-j symbol can be nonzero only when mmin<=m<=mmax;
  subroutine get_3js_m3_general(j1,j2, j3, m3, threej, mmin,mmax)
    COOP_INT  j1,j2,j3,m3,mmin,mmax
    COOP_REAL  threej(-j1:j1)
    mmin = max(-j1,-j2-m3)
    mmax = min(j1, j2-m3)
    if(j1.le.j2)then
       call coop_threej_fixed_m1(j3, j1, j2, m3, threej(-j1:j1))
    else
       call coop_threej_fixed_m1(j3, j2, j1, m3, threej(-j2:j2))
       if(mod(j1+j2+j3,2).eq.0)then
          threej(mmin:mmax) = threej(min(j2,j1-m3):-j2:-1)
       else
          threej(mmin:mmax) = -threej(min(j2,j1-m3):-j2:-1)
       endif
       threej(mmax+1:j1) = 0._dl
       threej(-j1:mmin-1)= 0._dl
    endif
  end subroutine get_3js_m3_general


  function coop_threej_approx(l1, l2, l3, m3, mu1, mu2) result(f)
    !!sum_{m1, m2} 4 pi/Sqrt[(2 l1 +1)(2 l2+1)] Y_{l1, m1}(theta1, 0) Y_{l2, m2}(theta, pi) ThreeJSymbol({l1, m1}, {l2, m2}, {l3, m3})
    COOP_REAL  f, mu1, mu2, threej(-l1:l1), Pl1ms(-l1:l1), Pl2ms(-l2:l2)
    COOP_INT  l1, l2, l3, m3, mmin, mmax
    if(is_triangle(l1,l2,l3) .and. m3 .le. l3)then
       call get_3js_m3_general(l1, l2, l3, m3, threej, mmin,mmax)
       call get_all_normalized_Plms_l(l1, mu1, Pl1ms)
       Pl1ms(1:l1:2) = -Pl1ms(1:l1:2)
       Pl1ms(-1:-l1:-2) = -Pl1ms(-1:-l1:-2) !!(-1)^m factor
       call get_all_normalized_Plms_l(l2, mu2, Pl2ms)
       if(abs(-m3-mmin) .gt. l2 .or. abs(-m3-mmax) .gt. l2) stop "Error in sum_YlmYlm3j"
       f= sum(threej(mmin:mmax)*Pl1ms(mmin:mmax)*Pl2ms(-m3-mmin:-m3-mmax:-1))
    else
       f = 0._dl
    endif
  end function coop_threej_approx

  subroutine coop_threej_quick(l1, l2, l3, m3max, Pl1ms, Pl2ms, threej, m1min, m1max, m2min, m2max, f)
    !!sum_{m1, m2} 4 pi/Sqrt[(2 l1 +1)(2 l2+1)] Y_{l1, m1}(theta1, 0) Y_{l2, m2}(theta, pi) ThreeJSymbol({l1, m1}, {l2, m2}, {l3, m3})
    COOP_INT  l1, l2, l3, m3max
    COOP_REAL  threej(-l1:l1, 0:m3max), Pl1ms(-l1:l1),  Pl2ms(-l2:l2)
    real(sp) f(0:m3max)
    COOP_INT  m1min(0:m3max), m1max(0:m3max), m2min(0:m3max), m2max(0:m3max)
    COOP_INT  m3
    do m3 = 0, min(m3max, l3)
       f(m3) = real(sum(threej(m1min(m3):m1max(m3), m3)*Pl1ms(m1min(m3):m1max(m3))*Pl2ms(m2min(m3):m2max(m3))))
    enddo
    if(l3.lt.m3max)f(l3+1:m3max) = 0.
  end subroutine coop_threej_quick


  subroutine coop_threej_generate_table(l1, l2, l3, m3max, st)
    COOP_INT  l1, l2, l3, m3max
    COOP_REAL ,dimension(:,:),allocatable::threej, Pl1ms, Pl2ms
    type(coop_threej_table) st
    COOP_INT  i, im1, im2, m3, m1min(0:m3max),m1max(0:m3max), m2min(0:m3max), m2max(0:m3max)
    st%ntheta = min(max(l1, l2, l3)*4 + 400, 16384) !!make sure this is an even number, since we will use the parity Plm(x) = (-1)^{l+m} Plm(-x)
    st%dtheta = coop_pi/(st%ntheta-1)
    if(allocated(st%s))deallocate(st%s)
    if(allocated(st%mu))deallocate(st%mu)
    allocate(st%s(0:m3max, st%ntheta, st%ntheta), st%mu(st%ntheta))
    call set_uniform(st%ntheta, st%mu, 0.d0, coop_pi)
    st%mu = cos(st%mu)
    allocate(threej(-l1:l1, 0:m3max), Pl1ms(-l1:l1, st%ntheta), Pl2ms(-l2:l2, st%ntheta))

    !$omp parallel do
    do i = 1, st%ntheta/2
       call get_Plms_l1l2(l1, l2, st%mu(i), Pl1ms(-l1:l1, i), Pl2ms(-l2:l2, i))
       Pl1ms(-l1:l1:2, st%ntheta+1-i) = Pl1ms(-l1:l1:2, i)
       Pl1ms(-l1+1:l1-1:2, st%ntheta+1-i) = -Pl1ms(-l1+1:l1-1:2, i)
       Pl2ms(-l2:l2:2, st%ntheta+1-i) = Pl2ms(-l2:l2:2, i)
       Pl2ms(-l2+1:l2-1:2, st%ntheta+1-i) = -Pl2ms(-l2+1:l2-1:2, i)
    enddo
    !$omp end parallel do

    !!prepare the 3j-symbol table
    do m3 = 0, min(m3max, l3)
       call get_3js_m3_general(l1, l2, l3, m3, threej(-l1:l1,m3), m1min(m3), m1max(m3))
       m2min(m3) = m3+m1min(m3)
       m2max(m3) = m3+m1max(m3)
       if(m2min(m3) .lt. -l2 .or. m2max(m3) .gt. l2) stop "error in m2min/m2max"
    enddo


    !$omp parallel do private(im1, im2)
    do im2 = 1, st%ntheta
       do im1 = 1, st%ntheta
          call coop_threej_quick(l1, l2, l3, m3max, Pl1ms(-l1:l1, im1), Pl2ms(-l2:l2, im2), threej(-l1:l1,0:m3max), m1min(0:m3max), m1max(0:m3max), m2min(0:m3max), m2max(0:m3max), st%s(0:m3max, im1, im2))
       enddo
    enddo
    !$omp end parallel do
    deallocate(Pl1ms, Pl2ms, threej)
  end subroutine coop_threej_generate_table

  !!compute (-1)^mP_{l1}^m(mu) (m=-l1, -l1+1, ..., l1) 
  !! save in Pl1ms
  !! and 
  !! P_{l2}^{-m}(mu) (m=-l2, -l2+1, ..., l2)
  !! => save in Pl2ms
  subroutine get_Plms_l1l2(l1, l2, mu, Pl1ms, Pl2ms)
    COOP_INT  l1, l2, lmax, m
    COOP_REAL  Pl1ms(-l1:l1), Pl2ms(-l2:l2)
    COOP_REAL  mu
    COOP_REAL , dimension(:),allocatable::Plms
    !!this will be a bit faster but less accurate
    if(max(l1,l2) .gt. 1000)then
       call get_all_normalized_Plms_l(l1, mu, Pl1ms)
       call get_all_normalized_Plms_l(l2, mu, Pl2ms)
       Pl2ms(-l2:l2) = Pl2ms(l2:-l2:-1)
       Pl1ms(1:l1:2) = -Pl1ms(1:l1:2)
       Pl1ms(-1:-l1:-2) = -Pl1ms(-1:-l1:-2)
       return
    endif
    lmax = max(l1, l2)
    allocate(Plms(0:lmax))
    do m=0, lmax, 2
       call get_all_normalized_Plms_m(lmax, m, mu, Plms(m:lmax))
       if(m.le.l1)then
          Pl1ms(m) = Plms(l1)
          Pl1ms(-m) = Plms(l1)
       endif
       if(m.le.l2)then
          Pl2ms(-m) = Plms(l2)  
          Pl2ms(m) = Plms(l2)
       endif
    enddo
    do m=1, lmax, 2
       call get_all_normalized_Plms_m(lmax, m, mu, Plms(m:lmax))
       if(m.le.l1)then
          Pl1ms(m) = -Plms(l1)
          Pl1ms(-m) = Plms(l1)
       endif
       if(m.le.l2)then
          Pl2ms(-m) = Plms(l2)  
          Pl2ms(m) = -Plms(l2)       
       endif
    enddo
    deallocate(Plms)
  end subroutine get_Plms_l1l2

  subroutine coop_threej_interp(st, m3max, mu1, mu2, ss)
    COOP_INT  m3max
    type(coop_threej_table) st
    COOP_REAL  ss(0:m3max)
    COOP_REAL  mu1, mu2, r1, r2
    COOP_INT  i1, i2
    r1 = acos(mu1)/st%dtheta + 1
    i1 = min(floor(r1), st%ntheta-1)
    r1 = r1 - i1
    r2 = acos(mu2)/st%dtheta + 1
    i2 = min(floor(r2), st%ntheta-1)
    r2 = r2 - i2
    if(i1 .le. 1 .or. i1 .ge. st%ntheta-1 .or. i2 .le. 1 .or. i2 .ge. st%ntheta-1)then
       ss = (st%s(0:m3max, i1, i2)*(1.d0-r1)+st%s(0:m3max, i1+1, i2)*r1)*(1.d0-r2) &
            + (st%s(0:m3max, i1, i2+1)*(1.d0-r1) + st%s(0:m3max, i1+1, i2+1)*r1)*r2
    else
       call bicubic_interp(m3max+1, st%s(0:m3max, i1-1:i1+2, i2-1:i2+2), r1, r2, ss)
    endif
  end subroutine coop_threej_interp


  subroutine coop_threej_table_free(this)
    class(coop_threej_table) this
    if(allocated(this%s))deallocate(this%s)
    if(allocated(this%mu))deallocate(this%mu)
  end subroutine coop_threej_table_free



  subroutine multiply_exp(a, x) !!replace A with A*exp(x), A can be negative
    COOP_REAL  a, x
    if(a .eq. 0.d0)  return
    a = sign(exp(x + log(abs(a))), a)
  end subroutine multiply_exp

  function is_triangle_real(k1, k2, k3)
    COOP_REAL  k1, k2, k3
    logical is_triangle_real
    is_triangle_real = (k1+k2 .gt. k3 .and. k2+k3 .gt. k1 .and. k1+k3 .gt. k2)
  end function is_triangle_real

  function is_triangle_int(l1, l2, l3)
    COOP_INT  l1, l2, l3
    logical is_triangle_int
    is_triangle_int = (l1+l2 .ge. l3 .and. l2+l3.ge.l1 .and. l1+l3.ge.l2)
  end function is_triangle_int



end module Coop_threej_mod
