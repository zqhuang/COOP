module coop_basicutils_mod
  use coop_constants_mod
  use coop_svd_mod
  implicit none
#include "constants.h"

  public

  integer,parameter::sp = kind(1.)
  integer,parameter::dl = kind(1.d0)

  private::sp, dl

  interface coop_swap
     module procedure coop_swap_real, coop_swap_int, coop_swap_real_array, coop_swap_int_array
  end interface coop_swap

  interface coop_isnan
     module procedure coop_isnan_d, coop_isnan_arrd, coop_isnan_arr2d, coop_isnan_s, coop_isnan_arrs, coop_isnan_arr2s
  end interface coop_isnan

  interface coop_smooth_data
     module procedure coop_smooth_data_d, coop_smooth_data_s
  end interface coop_smooth_data


  interface coop_set_uniform
     module procedure coop_set_uniform_d, coop_set_uniform_s
  end interface coop_set_uniform

contains

  Function coop_OuterProd(a, b) result(outerprod)
    COOP_REAL, DIMENSION(:), INTENT(IN) :: a,b
    COOP_REAL, DIMENSION(size(a),size(b)) :: outerprod
    Outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  end function Coop_OuterProd

  function coop_maxloc(arr) result(imaxloc)
    COOP_REAL, DIMENSION(:), INTENT(IN) :: arr
    COOP_INT imaxloc
    COOP_INT,DIMENSION(1)::imax
    imax=maxloc(arr(:))
    imaxloc=imax(1)
  end function coop_maxloc

  function coop_minloc(arr) result(iminloc)
    COOP_REAL, DIMENSION(:), INTENT(IN) :: arr
    COOP_INT iminloc
    COOP_INT,DIMENSION(1)::imin
    imin=minloc(arr(:))
    iminloc=imin(1)
  end function coop_minloc

  subroutine coop_swap_real(x,y)
    COOP_REAL x,y,tmp
    tmp=x
    x=y
    y=tmp
  end subroutine coop_swap_real

  subroutine coop_swap_int(x,y)
    COOP_INT x,y,tmp
    tmp=x
    x=y
    y=tmp
  end subroutine coop_swap_int

  subroutine coop_swap_real_array(x,y)
    COOP_REAL,dimension(:),intent(INOUT):: x,y
    COOP_REAL tmp(size(x))
    tmp=x
    x=y
    y=tmp
  end subroutine coop_swap_real_array

  subroutine coop_swap_int_array(x,y)
    COOP_INT,dimension(:),intent(INOUT):: x,y
    COOP_INT tmp(size(x))
    tmp=x
    x=y
    y=tmp
  end subroutine coop_swap_int_array

  

  function coop_systime_sec(reset) result(sec)
    COOP_REAL sec
    logical,optional::Reset
    COOP_REAL,save::pretime=0.d0 
    COOP_INT nowtime, countrate
    call system_clock(nowtime, countrate)
    if(present(Reset))then
       if(Reset)then
          Pretime = NowTime/dble(CountRate)
          sec = 0.d0
          return
       end if
    end if
    sec = NowTime/dble(CountRate)-PreTime
    return 
  end function coop_systime_sec

  subroutine coop_PrtSysTime(Reset)
    logical,optional::Reset
    COOP_REAL,save::pretime=0.d0
    COOP_INT Nowtime,Countrate
    call system_clock(NOWTIME,COUNTRATE)
    If(Present(Reset))then
       Pretime=Nowtime/Real(Countrate)
       If(Reset)Then
          print*,"===== time label reset to zero. ===="
       endif
    else
       Write(*,'(A18,F14.4,A10)') "===== time label: ",Nowtime/Real(Countrate)-pretime," sec ====="
    endif
  end subroutine Coop_PrtSysTime

  subroutine coop_init_random()
    COOP_INT nowtime, countrate, s(3), si
    call System_clock(NOWTIME,COUNTRATE)
    s(1) = nowtime
    s(2) = countrate
    s(3) = mod(nowtime, 17)+1
    call random_seed(SIZE = si)
    call random_seed(PUT = s(1:si))
  end subroutine coop_init_random


  Subroutine Coop_return_error(name, message, action)
    CHARACTER(LEN=*),optional::name
    Character(Len=*),optional::Message
    Character(Len=*),optional::action
    Character::ans
    if(present(name)) WRITE(*,*) "Error/Warning in "//TRIM(name)
    if(present(Message)) write(*,*) "** "//Trim(Message)//" **"
    if(present(ACTION))THEN
       select case(action)
       case("stop","STOP","Stop", "abort", "Abort", "ABORT")
          stop "program terminated"
       case("return","RETURN","Return", "pass", "Pass", "PASS")
          return
       case("wait","Wait","WAIT", "ask", "ASK", "Ask")
100       write(*,*) "Ignore this error and continue? (Y/N)"
          read(*,*) ans
          if(ans == "y" .or. ans == "Y")then
             return
          else
             goto 100
          endif
       end select
    else
       stop
    endif
  end subroutine Coop_return_error



  function coop_getdim(name, N1,N2,N3,N4,N5,N6) result(getdim)
    COOP_UNKNOWN_STRING name
    COOP_INT,INTENT(IN)::N1,N2
    COOP_INT,OPTIONAL::N3,N4,N5,N6
    COOP_INT getdim
    getdim=N1
    if(N1.NE.N2)THEN
       write(*,"(A,2I8)") "Dimension Error:",N1, N2
       CALL Coop_return_error(name)
       return
    ELSE
       if(present(N3))THEN
          if(N1.NE.N3)THEN
             write(*,"(A,3I8)") "Dimension Error:",N1,N2, N3
             CALL Coop_return_error(name)
             return
          Endif
          if(present(N4))THEN
             if(N1.NE.N4)then
                write(*,"(A,4I8)") "Dimension Error:",N1,N2, N3, N4
                call Coop_return_error(Name)
                return
             endif
             if(Present(N5))then
                if(N5.ne.N1)then
                   write(*,"(A,5I8)") "Dimension Error:",N1,N2, N3, N4, N5
                   call Coop_return_error(Name)
                   return
                endif
                if(present(N6))then
                   if(N1.ne.N6)then
                      write(*,"(A,6I8)") "Dimension Error:",N1,N2, N3, N4, N5, N6
                      call Coop_return_error(Name)
                      return
                   endif
                endif
             endif
          endif
       endif
    endif
  end function coop_getdim

  
  subroutine coop_set_uniform_d(n, x, lower, upper, logscale)
    COOP_INT n, i
    COOP_REAL x(n), lower, upper, rlow, dx
    logical,optional::logscale
    x(1) = lower
    x(n) = upper
    if(present(logscale))then
       if(logscale)then
          dx = (log(upper) - log(lower))/(n-1)
          rlow = log(lower) - dx
          !$omp parallel do
          do i = 2, n-1
             x(i) = exp(rlow + dx*i)
          enddo
          !$omp end parallel do
          return
       endif
    endif
    dx = (upper-lower)/(n-1)
    rlow = lower-dx
    !$omp parallel do
    do i = 2, n-1
       x(i) =rlow + dx*i
    enddo
    !$omp end parallel do

  end subroutine coop_set_uniform_d


    subroutine coop_set_uniform_s(n, x, lower, upper, logscale)
    COOP_INT n, i
    COOP_SINGLE x(n), lower, upper, rlow, dx
    logical,optional::logscale
    x(1) = lower
    x(n) = upper
    if(present(logscale))then
       if(logscale)then
          dx = (log(upper) - log(lower))/(n-1)
          rlow = log(lower) - dx
          !$omp parallel do
          do i = 2, n-1
             x(i) = exp(rlow + dx*i)
          enddo
          !$omp end parallel do
          return
       endif
    endif
    dx = (upper-lower)/(n-1)
    rlow = lower-dx
    !$omp parallel do
    do i = 2, n-1
       x(i) =rlow + dx*i
    enddo
    !$omp end parallel do

  end subroutine coop_set_uniform_s


  subroutine coop_locate(n, x, needle, loc, res)
    COOP_INT n, loc, imin, imax
    COOP_REAL x(n), needle, res
    if(x(1).le. x(n))then
       if(needle .lt. x(1))then
          loc = 0
          return
       endif
       if(needle .gt. x(n))then
          loc = n
          return
       endif
       imin = 1
       imax = n
       do while(imax - imin .gt. 1)
          loc = (imax + imin)/2
          if(x(loc).le. needle)then
             imin = loc
          else
             imax = loc
          endif
       enddo
       loc = imin
       if(needle .gt. x(imin))then
          res = (needle - x(imin))/(x(imax)-x(imin))
       else
          res = 0
       endif
       return
    else
       if(needle .gt. x(1))then
          loc = 0
          return
       endif
       if(needle .lt. x(n))then
          loc = n
          return
       endif
       imin = 1
       imax = n
       do while(imax - imin .gt. 1)
          loc = (imax + imin)/2
          if(x(loc).ge. needle)then
             imin = loc
          else
             imax = loc
          endif
       enddo
       loc = imin
       if(needle .lt. x(imin))then
          res = (needle - x(imin))/(x(imax)-x(imin))
       else
          res = 0
       endif
       return       
    endif
  end subroutine coop_locate


  subroutine coop_spline(n, x, y, y2, ypl, ypr)
    COOP_INT n, i
    COOP_REAL x(n), y(n), y2(n)
    COOP_REAL, optional::ypl,ypr
    COOP_REAL yil, yir, bet, dxr, dxl
    COOP_REAL gam(n-1)
    if(n.le.2)then
       y2 = 0.
       return
    endif
    dxr = x(2) - x(1)
    yir=(y(2)-y(1))/dxr
    if(present(ypl))then
       y2(1)=(yir-ypl)/dxr*3.
       gam(1)= 0.5
    else
       y2(1)=0.
       gam(1)=0.
    endif
    dxr = dxr/6.
    do i=2, n-1
       dxl = dxr
       dxr=x(i+1)-x(i)
       bet=(x(i+1)-x(i-1))/3.-dxl*gam(i-1)
       if(abs(bet) .lt. 1.d-30) stop 'Error in SPLinE.'
       yil=yir
       yir=(y(i+1)-y(i))/dxr
       y2(i)=(yir-yil-dxl*y2(i-1))/bet
       dxr=dxr/6.
       gam(i)=dxr/bet
    enddo
    if(present(ypr))then
       bet=(x(n)-x(n-1))/3.-dxr*gam(n-1)
       if(abs(bet) .lt. 1.d-30) stop 'Error in SPLinE.'
       y2(n)=(ypr-yir-dxr*y2(n-1))/bet
    else
       y2(n)=0.
    endif
    do i=n-1, 1 , -1
       y2(i)=y2(i)-gam(i)*y2(i+1)
    enddo
  end subroutine coop_spline


  subroutine coop_splint(n, x, y, y2, xs, ys)
    COOP_INT n, l, r
    COOP_REAL x(n), y(n), y2(n)
    COOP_REAL xs, ys, a, b
    call coop_locate(n, x, xs, l, b)
    if(l .lt. 1)then
       ys=y(1)
       return
    endif
    if( l .ge. n)then
       ys=y(n)
       return
    endif
    r = l + 1
    a =  1.d0 - b
    ys=y(l)*a+y(r)*b+  &
         (y2(l)*(a*a-1.)*a+y2(r)*(b*b-1.)*b)/6.*(x(r)-x(l))**2
  end subroutine coop_splint

  subroutine coop_splint_derv(n, x, y, y2, xs, ys)
    COOP_REAL,DIMENSION(:),INTENT(IN)::x, y, y2
    COOP_REAL xs, ys
    COOP_INT n, j, l, r
    COOP_REAL::a, b
    call coop_locate(n, x, xs, l, b)
    If(l.lt.1)then
       ys = (y(2)-y(1))/(x(2)-x(1)) 
    endif
    if(l.ge. n)Then
       ys = (y(n)-y(n-1))/(x(n)-x(n-1))
       return
    endif
    r = l + 1
    a = 1.d0 - b
    a=(x(r)-xs)/(X(r)-X(r-1))
    b=1.D0-a
    ys=(y(r)-y(l))/(x(r)-x(l)) + &
         (-y2(l)*(0.5d0*a*a-1.d0/6.d0) + y2(r)*(0.5d0*b*b-1.d0/6.d0))*(x(r)-x(r-1))
  end subroutine Coop_splint_derv



  subroutine coop_spline_uniform(n, y, y2)
    COOP_INT n, i
    COOP_REAL  y(n), y2(n)
    COOP_REAL yil, yir, bet, ypl, ypr
    COOP_REAL gam(n-1)
    if(n.le.2)then
       y2 = 0.
       return
    endif
    yir=(y(2)-y(1))
    ypl=(2.d0*y(2)-1.5d0*y(1)-0.5d0*y(3))
    ypr = -(2.d0*y(n-1)-1.5d0*y(n)-0.5d0*y(n-2))
    y2(1)=(yir-ypl)*3.
    gam(1)= 0.5
    do i=2, n-1
       bet=2.d0/3.d0-gam(i-1)/6.d0
       if(abs(bet) .lt. 1.d-30) stop 'Error in spline.'
       yil=yir
       yir=(y(i+1)-y(i))
       y2(i)=(yir-yil-y2(i-1)/6.d0)/bet
       gam(i)=(1.d0/6.d0)/bet
    enddo
    bet=1.d0/3.-gam(n-1)/6.d0
    if(abs(bet) .lt. 1.d-30) stop 'Error in spline.'
    y2(n)=(ypr-yir-y2(n-1)/6.d0)/bet
    do i=n-1, 1 , -1
       y2(i)=y2(i)-gam(i)*y2(i+1)
    enddo
    y2 =  y2/6.
  end subroutine coop_spline_uniform


  subroutine coop_naturalspline_uniform(n, y, y2)
    COOP_INT n, i
    COOP_REAL  y(n), y2(n)
    COOP_REAL yil, yir, bet
    COOP_REAL gam(n-1)
    if(n.le.2)then
       y2 = 0.
       return
    endif
    yir=(y(2)-y(1))
    y2(1)=0.
    gam(1)=0.
    do i=2, n-1
       bet = 2.d0/3.-gam(i-1)/6.d0
       if(abs(bet) .lt. 1.d-30) stop 'Error in SPLinE.'
       yil=yir
       yir=(y(i+1)-y(i))
       y2(i)=(yir-yil-y2(i-1)/6.d0)/bet
       gam(i)=(1.d0/6.d0)/bet
    enddo
    y2(n)=0.
    do i=n-1, 1 , -1
       y2(i)=y2(i)-gam(i)*y2(i+1)
    enddo
    y2 = y2/6.d0
  end subroutine coop_naturalspline_uniform


  subroutine coop_splint_uniform(n, xmin, dx, y, y2, xs, ys)
    COOP_INT n, l, r
    COOP_REAL xmin, dx, y(n), y2(n), xs, ys, a, b
    b = (xs - xmin)/dx + 1.
    l = floor(b)
    if(l.lt.1)then
       ys = y(1)
       return
    endif
    if(l.ge.n)then
       ys = y(n)
       return
    endif
    b=b-l
    r = l+1
    a=1.-b
    ys=y(l)*a+y(r)*b+  &
         (y2(l)*(a*a-1.)*a+y2(r)*(b*b-1.)*b)
  end subroutine coop_splint_uniform




  subroutine coop_splint_derv_uniform(n, xmin, dx, y, y2, xs, ys)
    COOP_INT n, l, r
    COOP_REAL xmin, dx, y(n), y2(n), xs, ys, a, b
    b = (xs - xmin)/dx + 1.d0
    l = floor(b)
    if(l.lt.1)then
       ys = (2.d0*y(2)-1.5d0*y(1)+0.5d0*y(3))/dx
       return
    endif
    if(l.ge.n)then
       ys = (1.5d0*y(n)-2.d0*y(n-1)+0.5d0*y(n-2))/dx
       return
    endif
    b=b-l
    r = l+1
    a=1.-b
    ys=(y(r)-y(l))/dx + &
         (-y2(l)*(3.d0*a*a-1.d0) + y2(r)*(3.d0*b*b-1.d0))/dx
  end subroutine coop_splint_derv_uniform



  subroutine coop_naturalsplint_derv_uniform(n, xmin, dx, y, y2, xs, ys)
    COOP_INT n, l, r
    COOP_REAL xmin, dx, y(n), y2(n), xs, ys, a, b
    b = (xs - xmin)/dx + 1.d0
    l = floor(b)
    if(l.lt.1)then
       ys = (y(2)-y(1))/dx
       return
    endif
    if(l.ge.n)then
       ys = (y(n)-y(n-1))/dx
       return
    endif
    b=b-l
    r = l+1
    a=1.-b
    ys=(y(r)-y(l))/dx + &
         (-y2(l)*(3.d0*a*a-1.d0) + y2(r)*(3.d0*b*b-1.d0))/dx
  end subroutine coop_naturalsplint_derv_uniform
  

  subroutine coop_cheb_eval_all(n, x, y)
    COOP_INT,intent(IN):: n
    COOP_REAL,intent(IN):: x
    COOP_REAL,intent(OUT)::y(0:n)
    COOP_REAL twox
    COOP_INT i
    y(0) = 1.d0
    if(n.eq.0) return
    y(1) = x
    twox = 2.*x
    do i = 2, n
       y(i) = twox * y(i-1) - y(i-2)
    enddo
  end subroutine coop_cheb_eval_all


  subroutine coop_get_cheb_value(n, c, x, y)
    COOP_INT n
    COOP_REAL c(n), x, y
    COOP_REAL twox, y1, y2, y3
    COOP_INT i
    twox = x+x
    y3 = 0.d0
    y2 = 0.d0
    do i = n, 2, -1
       y1 = c(i) + twox * y2 - y3
       y3 = y2
       y2 = y1
    enddo
    y = c(1) + x* y2 - y3    
  end subroutine coop_get_cheb_value


  subroutine coop_fit_template(n, m, y, tpls, c)
    COOP_INT,intent(IN):: n, m
    COOP_REAL,intent(IN)::y(n)
    COOP_REAL,intent(INOUT)::tpls(n, m)
    COOP_REAL,intent(OUT):: c(m)
    call coop_svd_least_square_one(n, m, tpls, y, c)
  end subroutine coop_fit_template


  subroutine coop_fit_templates(n, m, nd, y, tpls, c)
    COOP_INT,intent(IN):: n, m, nd
    COOP_REAL,intent(IN)::y(n, nd)
    COOP_REAL,intent(INOUT)::tpls(n, m)
    COOP_REAL,intent(OUT):: c(m, nd)
    call coop_svd_least_square_all(n, m, nd, tpls, y, c)
  end subroutine coop_fit_templates


  subroutine coop_chebfit(n, x,y, m, a,b,c)
    COOP_INT,intent(in):: n, m
    COOP_REAL,intent(in)::a, b
    COOP_REAL,intent(in)::x(n), y(n)
    COOP_REAL,intent(out)::c(m)
    COOP_REAL fx(n, m ),t(n)
    COOP_INT i
    if(m.gt.n)then
       call coop_return_error("coop_chebfit", "Not enough data.", "stop")
    endif
    t=2.d0*(x-a)/(b-a)-1.d0
    do i=1,n
       call coop_cheb_eval_all(m-1,t(i), Fx(i,1:m))
    enddo
    call coop_fit_template(n, m, y, fx, c)
  end subroutine coop_chebfit

  subroutine coop_chebfit_uniform(n, y, m, c)
    COOP_INT,intent(in):: n, m
    COOP_REAL,intent(in)::y(n)
    COOP_REAL,intent(out)::c(m)
    COOP_REAL fx(n, m ),dt
    COOP_INT i
    if(m.gt.n)then
       call coop_return_error("coop_chebfit_uniform", "Not enough data", "stop")
    endif
    dt = 2.d0/(n-1)
    do i=1,n
       call coop_cheb_eval_all(m-1, -1.+ (i-1) * dt, Fx(i,1:m))
    enddo
    call coop_fit_template(n, m, y, fx, c)
  end subroutine coop_chebfit_uniform


  !!suppose f(x) is expanded with chebyshev polynomial c (with boundaries a, b mapped to -1, 1)
  !!this subroutine expand f'(x) with chebyshev polynomials (same boundaries)
  subroutine coop_chebfit_derv(n, a, b, c, cder)
    COOP_INT n, j
    COOP_REAL c(n), cder(n), a, b
    cder(n) = 0.d0
    cder(n-1) = 2*(n-1)*c(n)
    do j = n-2,1,-1
       cder(j)=cder(j+2)+2*j*c(j+1)
    enddo
    cder(1) = cder(1)/2.d0
    cder = cder*(2.d0/(b-a))
    return
  end subroutine coop_chebfit_derv

  subroutine coop_chebeval(n, a, b, c, x, y)
    COOP_INT n
    COOP_REAL a, b, c(n), x, y, t
    t = 2.*(x-a)/(b-a)-1.
    call  coop_get_cheb_value(n, c, t, y)
  end subroutine coop_chebeval

  function coop_isnan_d(x)  result(isnan)
    real(dl) x
    logical isnan
    isnan = .not. (x.gt.0.d0 .or. x.le.0.d0)
  end function coop_isnan_d

  function coop_isnan_arrd(x) result(isnan)
    real(dl),dimension(:):: x
    logical isnan
    isnan = .not. (all(x.gt.0.d0 .or. x.le.0.d0))
  end function coop_isnan_arrd

  function coop_isnan_arr2d(x) result(isnan)
    real(dl),dimension(:,:):: x
    logical isnan
    isnan = .not. (all(x.gt.0.d0 .or. x.le.0.d0))
  end function coop_isnan_arr2d


  function coop_isnan_s(x)  result(isnan)
    real(sp) x
    logical isnan
    isnan = .not. (x.gt.0.d0 .or. x.le.0.d0)
  end function coop_isnan_s

  function coop_isnan_arrs(x) result(isnan)
    real(sp),dimension(:):: x
    logical isnan
    isnan = .not. (all(x.gt.0.d0 .or. x.le.0.d0))
  end function coop_isnan_arrs

  function coop_isnan_arr2s(x) result(isnan)
    real(sp),dimension(:,:):: x
    logical isnan
    isnan = .not. (all(x.gt.0.d0 .or. x.le.0.d0))
  end function coop_isnan_arr2s


  subroutine coop_vector_cross_product(v1, v2, v3)
    COOP_REAL,dimension(3)::v1,v2,v3
    v3 = (/ v1(2)*v2(3) - v1(3)*v2(2),  v1(3)*v2(1)-v1(1)*v2(3), v1(1)*v2(2)-v1(2)*v2(1) /)
  end subroutine coop_vector_cross_product

  !!this is a simple version assuming that n is not too big
  function coop_polyvalue(n, p, x) result(y)
    COOP_INT n, j
    COOP_REAL p(n)
    COOP_REAL x, y
    y = p(n)
    do j= n-1,1,-1
       y= y*x + p(j)
    enddo    
  end function coop_polyvalue

  function coop_InputArgs(i) result(args)
    COOP_STRING args
    COOP_INT, intent(in) :: i
    if (iargc() < i) then
       args = ''
    else
       call getarg(i, args)
    end if
  end function Coop_InputArgs


  subroutine coop_smooth_data_d(n, y, sigma)
    COOP_INT::n
    COOP_REAL::y(n)
    COOP_INT::sigma
    COOP_REAL::w(-3*sigma:3*sigma), ycopy(1-3*sigma:n+3*sigma)
    COOP_INT i, m
    w(0) = 1.d0
    m = 3*sigma
    !$omp parallel do
    do i = 1, m
       w(i) = exp(-(dble(i)/sigma)**2/2.d0)
       w(-i) = w(i)
    enddo
    !$omp end parallel do    
    w = w/sum(w)
    ycopy(1:n) = y
    ycopy(1-m:0) = y(1)    
    ycopy(n+1:n+m) = y(n)
    !$omp parallel do
    do i=1, n
       y(i) = sum(w*ycopy(i-m:i+m))
    enddo
    !$omp end parallel do
  end subroutine coop_smooth_data_d

  subroutine coop_smooth_data_s(n, y, sigma)
    COOP_INT::n
    COOP_SINGLE::y(n)
    COOP_INT::sigma
    COOP_SINGLE::w(-3*sigma:3*sigma), ycopy(1-3*sigma:n+3*sigma)
    COOP_INT i, m
    w(0) = 1.
    m = 3*sigma
    !$omp parallel do
    do i = 1, m
       w(i) = exp(-(dble(i)/sigma)**2/2.d0)
       w(-i) = w(i)
    enddo
    !$omp end parallel do
    w = w/sum(w)
    ycopy(1:n) = y
    ycopy(1-m:0) = y(1)    
    ycopy(n+1:n+m) = y(n)
    !$omp parallel do
    do i=1, n
       y(i) = sum(w*ycopy(i-m:i+m))
    enddo
    !$omp end parallel do
  end subroutine coop_smooth_data_s
  


end module coop_basicutils_mod
