module coop_interpolation_mod
  use coop_wrapper_typedef
  use coop_matrix_mod
  implicit none

#include "constants.h"
  private

  public::coop_bilinear_interp, coop_bicubic_interp, coop_linear_interp, coop_linear_least_square_fit, coop_spline_fill, coop_ratint, coop_ratval, coop_pade, coop_ratlsq, coop_smooth_fit

  Interface coop_ratlsq
     module procedure coop_ratlsq_f, coop_ratlsq_arr
  end Interface coop_ratlsq

  Interface coop_linear_least_square_fit
     module procedure coop_linear_least_square_fit_s, coop_linear_least_square_fit_d
  end Interface coop_linear_least_square_fit

  Interface coop_linear_interp
     module procedure linear_interp_s, linear_interp_d
  end Interface coop_linear_interp

  Interface coop_bilinear_interp
     module procedure bilinear_interp_s, bilinear_interp_v, bilinear_interp_s_sp, bilinear_interp_v_sp
  end Interface coop_bilinear_interp


  Interface coop_bicubic_interp
     module procedure bicubic_interp_s, bicubic_interp_v, bicubic_interp_s_sp, bicubic_interp_v_sp, bicubic_interp_v_sp2sp
  end Interface coop_bicubic_interp

  Interface coop_spline_fill
     module procedure coop_spline_fill_d, coop_spline_fill_s, coop_spline_fill_vd, coop_spline_fill_vs, coop_spline_fill_md, coop_spline_fill_ms
  end Interface coop_spline_fill


  type coop_smooth_fit
     logical::initialized = .false.
     logical::xlog = .false.
     logical::ylog = .false.
     COOP_REAL::dev=0.d0
     COOP_INT::n_up = 0
     COOP_INT::n_down = 0
     COOP_REAL,dimension(:),allocatable::c
     COOP_REAL::center = 0.d0
     COOP_REAL::norm = 1.d0
   contains
     procedure::free => coop_smooth_fit_free
     procedure::fit => coop_smooth_fit_fit
     procedure::eval => coop_smooth_fit_eval     
  end type coop_smooth_fit

contains

  subroutine coop_linear_least_square_fit_d(n, x, y, k, b, r)
    integer n
    COOP_REAL x(n), y(n), k, b, xbar, ybar, sxx, syy, sxy
    COOP_REAL, optional::r
    xbar  = sum(x)/n
    ybar = sum(y)/n
    sxx = sum((x-xbar)**2)
    syy = sum((y-ybar)**2)
    sxy = sum((x-xbar)*(y-ybar))
    if(sxx .le. 0.d0)then
       k = coop_infinity
    else
       k = sxy/sxx
    endif
    b = ybar - k*xbar
    if(present(r))then
       if(sxx .lt. 1.d-30 .or. syy .lt. 1.d-30)then
          r = 1.d0
       else
          r = sxy/sqrt(sxx*syy)
       endif
    endif
  end subroutine coop_linear_least_square_fit_d


  subroutine coop_linear_least_square_fit_s(n, x, y, k, b, r)
    integer n
    COOP_SINGLE x(n), y(n), k, b, xbar, ybar, sxx, syy, sxy
    COOP_SINGLE, optional::r
    xbar  = sum(x)/n
    ybar = sum(y)/n
    sxx = sum((x-xbar)**2)
    syy = sum((y-ybar)**2)
    sxy = sum((x-xbar)*(y-ybar))
    if(sxx .le. 0.)then
       k = coop_infinity
    else
       k = sxy/sxx
    endif
    b = ybar - k*xbar
    if(present(r))then
       if(sxx .le. 0. .or. syy .le. 0.)then
          r = 1.
       else
          r = sxy**2/sxx/syy
       endif
    endif
  end subroutine coop_linear_least_square_fit_s

  subroutine linear_interp_s(n, x, y, xs, ys)
    COOP_SINGLE xs, ys
    integer n, j ,l , r
    COOP_SINGLE a
    COOP_SINGLE x(n), y(n)
    if(x(n).gt.x(1))then
       if(xs .lt. x(1))then
          ys = y(1)
          return
       elseif(xs .gt. x(n))then
          ys = y(n)
          return
       end if
       l = 1
       r = n
       do while(r - l .gt. 1)
          j = (r + l) /2
          if(x(j) .gt. xs)then
             r = j
          else
             l = j
          end if
       end do
    else
       if(xs .ge. x(1))then
          ys = y(1)
          return
       elseif(xs .le. x(n))then
          ys = y(n)
          return
       end if
       l = 1
       r = n
       do while(r- l .gt.1)
          j=(r+l)/2
          if(x(j).lt.xs)then
             r= j
          else
             l = j
          endif
       enddo
    endif
    if(x(r).ne.x(r-1))then
       a=(xs-x(r-1))/(x(r)-x(r-1))
       ys=y(r)*a + y(r-1)*(1.d0-a)
    else
       ys = y(r)
    endif
  end subroutine linear_interp_s

  subroutine linear_interp_d(n, x, y, xs, ys)
    COOP_REAL xs, ys
    integer n, j ,l , r
    COOP_REAL a
    COOP_REAL x(n), y(n)
    if(x(n).gt.x(1))then
       if(xs .lt. x(1))then
          ys = y(1)
          return
       elseif(xs .gt. x(n))then
          ys = y(n)
          return
       end if
       l = 1
       r = n
       do while(r - l .gt. 1)
          j = (r + l) /2
          if(x(j) .gt. xs)then
             r = j
          else
             l = j
          end if
       end do
    else
       if(xs .ge. x(1))then
          ys = y(1)
          return
       elseif(xs .le. x(n))then
          ys = y(n)
          return
       end if
       l = 1
       r = n
       do while(r- l .gt.1)
          j=(r+l)/2
          if(x(j).lt.xs)then
             r= j
          else
             l = j
          endif
       enddo
    endif
    if(x(r).ne.x(r-1))then
       a=(xs-x(r-1))/(x(r)-x(r-1))
       ys=y(r)*a + y(r-1)*(1.d0-a)
    else
       ys = y(r)
    endif
  end subroutine linear_interp_d


  subroutine bilinear_interp_s(f, rx, ry, z)
    COOP_REAL ,intent(IN)::f(2,2)
    COOP_REAL rx, ry
    COOP_REAL z
    z = (f(1,1)*(1.d0-rx)+f(2,1)*rx)*(1.d0-ry) + (f(1,2)*(1.d0-rx)+f(2,2)*rx)*ry
  end subroutine bilinear_interp_s


  subroutine bilinear_interp_v(n, f, rx, ry, z)
    integer n
    COOP_REAL ,intent(IN)::f(n, 2,2)
    COOP_REAL rx, ry
    COOP_REAL z(n)
    z = (f(:,1,1)*(1.d0-rx)+f(:,2,1)*rx)*(1.d0-ry) + (f(:,1,2)*(1.d0-rx)+f(:,2,2)*rx)*ry
  end subroutine bilinear_interp_v


  subroutine bicubic_interp_s(f, rx, ry, z)
    COOP_REAL,intent(IN)::f(4,4)
    COOP_REAL rx, ry
    COOP_REAL :: x(16)
    COOP_REAL, INTENT(OUT) :: z
    COOP_REAL, DIMENSION(16,16),parameter:: wt = reshape( (/ &
         0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0, &
         0,-2,0,0,0,0,0,0,0,2,0,0,0,0,0,0, &
         0,4,0,0,0,-10,0,0,0,8,0,0,0,-2,0,0, &
         0,-2,0,0,0,6,0,0,0,-6,0,0,0,2,0,0, &
         0,0,0,0,-2,0,2,0,0,0,0,0,0,0,0,0, &
         1,0,-1,0,0,0,0,0,-1,0,1,0,0,0,0,0, &
         -2,0,2,0,5,0,-5,0,-4,0,4,0,1,0,-1,0, &
         1,0,-1,0,-3,0,3,0,3,0,-3,0,-1,0,1,0, &
         0,0,0,0,4,-10,8,-2,0,0,0,0,0,0,0,0, &
         -2,5,-4,1,0,0,0,0,2,-5,4,-1,0,0,0,0, &
         4,-10,8,-2,-10,25,-20,5,8,-20,16,-4,-2,5,-4,1, &
         -2,5,-4,1,6,-15,12,-3,-6,15,-12,3,2,-5,4,-1, &
         0,0,0,0,-2,6,-6,2,0,0,0,0,0,0,0,0, &
         1,-3,3,-1,0,0,0,0,-1,3,-3,1,0,0,0,0, &
         -2,6,-6,2,5,-15,15,-5,-4,12,-12,4,1,-3,3,-1, &
         1,-3,3,-1,-3,9,-9,3,3,-9,9,-3,-1,3,-3,1 &
         /), (/ 16,16 /)) / 4.d0
    x = reshape(matmul(reshape(f, (/1, 16/) ), wt), (/ 16 /))
    z=((x(16 )*ry+x(15 ))*ry+x(14 ))*ry+x(13 )
    z=rx*z+((x(12 )*ry+x(11 ))*ry+x(10 ))*ry+x(9 )
    z=rx*z+((x(8 )*ry+x(7 ))*ry+x(6 ))*ry+x(5 )
    z=rx*z+((x(4 )*ry+x(3 ))*ry+x(2 ))*ry+x(1 )
  End subroutine bicubic_interp_s

  Subroutine bicubic_interp_v(n, f, rx, ry, z)
    integer n
    COOP_REAL,intent(IN)::f(n,4,4)
    COOP_REAL rx, ry
    COOP_REAL :: x(n, 16)
    COOP_REAL, INTENT(OUT) :: z(n)
    COOP_REAL, DIMENSION(16,16),parameter:: wt = reshape( (/ &
         0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0, &
         0,-2,0,0,0,0,0,0,0,2,0,0,0,0,0,0, &
         0,4,0,0,0,-10,0,0,0,8,0,0,0,-2,0,0, &
         0,-2,0,0,0,6,0,0,0,-6,0,0,0,2,0,0, &
         0,0,0,0,-2,0,2,0,0,0,0,0,0,0,0,0, &
         1,0,-1,0,0,0,0,0,-1,0,1,0,0,0,0,0, &
         -2,0,2,0,5,0,-5,0,-4,0,4,0,1,0,-1,0, &
         1,0,-1,0,-3,0,3,0,3,0,-3,0,-1,0,1,0, &
         0,0,0,0,4,-10,8,-2,0,0,0,0,0,0,0,0, &
         -2,5,-4,1,0,0,0,0,2,-5,4,-1,0,0,0,0, &
         4,-10,8,-2,-10,25,-20,5,8,-20,16,-4,-2,5,-4,1, &
         -2,5,-4,1,6,-15,12,-3,-6,15,-12,3,2,-5,4,-1, &
         0,0,0,0,-2,6,-6,2,0,0,0,0,0,0,0,0, &
         1,-3,3,-1,0,0,0,0,-1,3,-3,1,0,0,0,0, &
         -2,6,-6,2,5,-15,15,-5,-4,12,-12,4,1,-3,3,-1, &
         1,-3,3,-1,-3,9,-9,3,3,-9,9,-3,-1,3,-3,1 &
         /), (/ 16,16 /)) / 4.d0
    x = matmul(reshape(f, (/ n , 16 /)), wt) 
    z=((x(:,16)*ry+x(:,15))*ry+x(:,14))*ry+x(:,13)
    z=rx*z+((x(:,12)*ry+x(:,11))*ry+x(:,10))*ry+x(:,9)
    z=rx*z+((x(:,8)*ry+x(:,7))*ry+x(:,6))*ry+x(:,5)
    z=rx*z+((x(:,4)*ry+x(:,3))*ry+x(:,2))*ry+x(:,1)
  End subroutine bicubic_interp_v

  subroutine bilinear_interp_s_sp(f, rx, ry, z)
    COOP_SINGLE ,intent(IN)::f(2,2)
    COOP_REAL rx, ry
    COOP_REAL z
    z = (f(1,1)*(1.d0-rx)+f(2,1)*rx)*(1.d0-ry) + (f(1,2)*(1.d0-rx)+f(2,2)*rx)*ry
  end subroutine bilinear_interp_s_sp


  subroutine bilinear_interp_v_sp(n, f, rx, ry, z)
    integer n
    COOP_SINGLE ,intent(IN)::f(n, 2,2)
    COOP_REAL rx, ry
    COOP_REAL z(n)
    z = (f(:,1,1)*(1.d0-rx)+f(:,2,1)*rx)*(1.d0-ry) + (f(:,1,2)*(1.d0-rx)+f(:,2,2)*rx)*ry
  end subroutine bilinear_interp_v_sp


  subroutine bicubic_interp_s_sp(f, rx, ry, z)
    COOP_SINGLE,intent(IN)::f(4,4)
    COOP_REAL rx, ry
    COOP_REAL :: x(16)
    COOP_REAL, INTENT(OUT) :: z
    COOP_REAL, DIMENSION(16,16),parameter:: wt = reshape( (/ &
         0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0, &
         0,-2,0,0,0,0,0,0,0,2,0,0,0,0,0,0, &
         0,4,0,0,0,-10,0,0,0,8,0,0,0,-2,0,0, &
         0,-2,0,0,0,6,0,0,0,-6,0,0,0,2,0,0, &
         0,0,0,0,-2,0,2,0,0,0,0,0,0,0,0,0, &
         1,0,-1,0,0,0,0,0,-1,0,1,0,0,0,0,0, &
         -2,0,2,0,5,0,-5,0,-4,0,4,0,1,0,-1,0, &
         1,0,-1,0,-3,0,3,0,3,0,-3,0,-1,0,1,0, &
         0,0,0,0,4,-10,8,-2,0,0,0,0,0,0,0,0, &
         -2,5,-4,1,0,0,0,0,2,-5,4,-1,0,0,0,0, &
         4,-10,8,-2,-10,25,-20,5,8,-20,16,-4,-2,5,-4,1, &
         -2,5,-4,1,6,-15,12,-3,-6,15,-12,3,2,-5,4,-1, &
         0,0,0,0,-2,6,-6,2,0,0,0,0,0,0,0,0, &
         1,-3,3,-1,0,0,0,0,-1,3,-3,1,0,0,0,0, &
         -2,6,-6,2,5,-15,15,-5,-4,12,-12,4,1,-3,3,-1, &
         1,-3,3,-1,-3,9,-9,3,3,-9,9,-3,-1,3,-3,1 &
         /), (/ 16,16 /)) / 4.d0
    x = reshape(matmul(reshape(f, (/1, 16/) ), wt), (/ 16 /))
    z=((x(16 )*ry+x(15 ))*ry+x(14 ))*ry+x(13 )
    z=rx*z+((x(12 )*ry+x(11 ))*ry+x(10 ))*ry+x(9 )
    z=rx*z+((x(8 )*ry+x(7 ))*ry+x(6 ))*ry+x(5 )
    z=rx*z+((x(4 )*ry+x(3 ))*ry+x(2 ))*ry+x(1 )
  End subroutine bicubic_interp_s_sp

  Subroutine bicubic_interp_v_sp(n, f, rx, ry, z)
    integer n
    COOP_SINGLE,intent(IN)::f(n,4,4)
    COOP_REAL rx, ry
    COOP_REAL :: x(n, 16)
    COOP_REAL, INTENT(OUT) :: z(n)
    COOP_REAL, DIMENSION(16,16),parameter:: wt = reshape( (/ &
         0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0, &
         0,-2,0,0,0,0,0,0,0,2,0,0,0,0,0,0, &
         0,4,0,0,0,-10,0,0,0,8,0,0,0,-2,0,0, &
         0,-2,0,0,0,6,0,0,0,-6,0,0,0,2,0,0, &
         0,0,0,0,-2,0,2,0,0,0,0,0,0,0,0,0, &
         1,0,-1,0,0,0,0,0,-1,0,1,0,0,0,0,0, &
         -2,0,2,0,5,0,-5,0,-4,0,4,0,1,0,-1,0, &
         1,0,-1,0,-3,0,3,0,3,0,-3,0,-1,0,1,0, &
         0,0,0,0,4,-10,8,-2,0,0,0,0,0,0,0,0, &
         -2,5,-4,1,0,0,0,0,2,-5,4,-1,0,0,0,0, &
         4,-10,8,-2,-10,25,-20,5,8,-20,16,-4,-2,5,-4,1, &
         -2,5,-4,1,6,-15,12,-3,-6,15,-12,3,2,-5,4,-1, &
         0,0,0,0,-2,6,-6,2,0,0,0,0,0,0,0,0, &
         1,-3,3,-1,0,0,0,0,-1,3,-3,1,0,0,0,0, &
         -2,6,-6,2,5,-15,15,-5,-4,12,-12,4,1,-3,3,-1, &
         1,-3,3,-1,-3,9,-9,3,3,-9,9,-3,-1,3,-3,1 &
         /), (/ 16,16 /)) / 4.d0
    x = matmul(reshape(f, (/ n , 16 /)), wt) 
    z=((x(:,16)*ry+x(:,15))*ry+x(:,14))*ry+x(:,13)
    z=rx*z+((x(:,12)*ry+x(:,11))*ry+x(:,10))*ry+x(:,9)
    z=rx*z+((x(:,8)*ry+x(:,7))*ry+x(:,6))*ry+x(:,5)
    z=rx*z+((x(:,4)*ry+x(:,3))*ry+x(:,2))*ry+x(:,1)
  End subroutine bicubic_interp_v_sp


  Subroutine bicubic_interp_v_sp2sp(n, f, rx, ry, z)
    integer n
    COOP_SINGLE,intent(IN)::f(n,4,4)
    COOP_REAL rx, ry
    COOP_REAL :: x(n, 16)
    COOP_SINGLE, INTENT(OUT) :: z(n)
    COOP_REAL, DIMENSION(16,16),parameter:: wt = reshape( (/ &
         0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0, &
         0,-2,0,0,0,0,0,0,0,2,0,0,0,0,0,0, &
         0,4,0,0,0,-10,0,0,0,8,0,0,0,-2,0,0, &
         0,-2,0,0,0,6,0,0,0,-6,0,0,0,2,0,0, &
         0,0,0,0,-2,0,2,0,0,0,0,0,0,0,0,0, &
         1,0,-1,0,0,0,0,0,-1,0,1,0,0,0,0,0, &
         -2,0,2,0,5,0,-5,0,-4,0,4,0,1,0,-1,0, &
         1,0,-1,0,-3,0,3,0,3,0,-3,0,-1,0,1,0, &
         0,0,0,0,4,-10,8,-2,0,0,0,0,0,0,0,0, &
         -2,5,-4,1,0,0,0,0,2,-5,4,-1,0,0,0,0, &
         4,-10,8,-2,-10,25,-20,5,8,-20,16,-4,-2,5,-4,1, &
         -2,5,-4,1,6,-15,12,-3,-6,15,-12,3,2,-5,4,-1, &
         0,0,0,0,-2,6,-6,2,0,0,0,0,0,0,0,0, &
         1,-3,3,-1,0,0,0,0,-1,3,-3,1,0,0,0,0, &
         -2,6,-6,2,5,-15,15,-5,-4,12,-12,4,1,-3,3,-1, &
         1,-3,3,-1,-3,9,-9,3,3,-9,9,-3,-1,3,-3,1 &
         /), (/ 16,16 /)) / 4.d0
    x = matmul(reshape(f, (/ n , 16 /)), wt) 
    z=real(((x(:,16)*ry+x(:,15))*ry+x(:,14))*ry+x(:,13))
    z=real(rx*z+((x(:,12)*ry+x(:,11))*ry+x(:,10))*ry+x(:,9))
    z=real(rx*z+((x(:,8)*ry+x(:,7))*ry+x(:,6))*ry+x(:,5))
    z=real(rx*z+((x(:,4)*ry+x(:,3))*ry+x(:,2))*ry+x(:,1))
  End subroutine bicubic_interp_v_sp2sp



  subroutine coop_spline_fill_d(n, x, y, computed, logx, logy)
    COOP_INT n
    COOP_REAL x(n), y(n)
    logical computed(n), logx, logy
    COOP_INT nc, i, j
    COOP_REAL,dimension(:),allocatable::xc, yc, yc2
    nc = count(computed)
    if(nc.eq.n)return
    allocate(xc(nc), yc(nc), yc2(nc))
    j = 1
    do i = 1, n
       if(computed(i))then
          xc(j) = x(i)
          yc(j) = y(i)
          j = j + 1
       endif
    enddo
    if(logx) xc = log(xc)
    if(logy) yc = log(yc)
    call coop_spline(nc, xc, yc, yc2)
    if(logx)then
       if(logy)then
          !$omp parallel do
          do i = 1, n
             if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, log(x(i)), y(i))
             y(i) = exp(y(i))
          enddo
          !$omp end parallel do
       else
          !$omp parallel do
          do i = 1, n
             if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, log(x(i)), y(i))
          enddo
          !$omp end parallel do
       endif
    else
       if(logy)then
          !$omp parallel do
          do i = 1, n
             if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, x(i), y(i))
             y(i) = exp(y(i))
          enddo
          !$omp end parallel do
       else
          !$omp parallel do
          do i = 1, n
             if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, x(i), y(i))
          enddo
          !$omp end parallel do
       endif
    endif
    deallocate(xc, yc, yc2)
  end subroutine coop_spline_fill_d



  subroutine coop_spline_fill_s(n, x, y, computed, logx, logy)
    COOP_INT n
    COOP_SINGLE x(n), y(n)
    logical computed(n), logx, logy
    COOP_INT nc, i, j
    COOP_REAL,dimension(:),allocatable::xc, yc, yc2
    COOP_REAL ytmp
    nc = count(computed)
    if(nc.eq.n)return
    allocate(xc(nc), yc(nc), yc2(nc))
    j = 1
    do i = 1, n
       if(computed(i))then
          xc(j) = COOP_REAL_OF(x(i))
          yc(j) = y(i)
          j = j + 1
       endif
    enddo
    if(logx) xc = log(xc)
    if(logy) yc = log(yc)
    call coop_spline(nc, xc, yc, yc2)
    if(logx)then
       if(logy)then
          !$omp parallel do private(ytmp)
          do i = 1, n
             if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, log(COOP_REAL_OF(x(i))), ytmp)
             y(i) = exp(ytmp)
          enddo
          !$omp end parallel do
       else
          !$omp parallel do private(ytmp)
          do i = 1, n
             if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, log(COOP_REAL_OF(x(i))), ytmp)
             y(i) = ytmp
          enddo
          !$omp end parallel do 
       endif
    else
       if(logy)then
          !$omp parallel do private(ytmp)
          do i = 1, n
             if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, COOP_REAL_OF(x(i)), ytmp)
             y(i) = exp(ytmp)
          enddo
          !$omp end parallel do
       else
          !$omp parallel do private(ytmp)
          do i = 1, n
             if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, COOP_REAL_OF(x(i)), ytmp)
             y(i) = ytmp
          enddo
          !$omp end parallel do
       endif
    endif
    deallocate(xc, yc, yc2)
  end subroutine coop_spline_fill_s



  subroutine coop_spline_fill_vd(m, n, x, y, computed, logx, logy)
    COOP_INT m, n
    COOP_REAL x(n), y(m,n)
    logical computed(n), logx, logy
    COOP_INT nc, i, j, im
    COOP_REAL,dimension(:),allocatable::xc, yc, yc2
    nc = count(computed)
    if(nc.eq.n)return
    allocate(xc(nc), yc(nc), yc2(nc))
    j = 1
    do i = 1, n
       if(computed(i))then
          xc(j) = x(i)
          j = j + 1
       endif
    enddo
    if(logx) xc = log(xc)

    do im = 1, m
       j = 1
       do i = 1, n
          if(computed(i))then
             yc(j) = y(im, i)
             j = j + 1
          endif
       enddo
       if(logy) yc = log(yc)
       call coop_spline(nc, xc, yc, yc2)
       if(logx)then
          if(logy)then
             !$omp parallel do
             do i = 1, n
                if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, log(x(i)), y(im, i))
                y(im, i) = exp(y(im, i))
             enddo
             !$omp end parallel do
          else
             !$omp parallel do
             do i = 1, n
                if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, log(x(i)), y(im, i))
             enddo
             !$omp end parallel do
          endif
       else
          if(logy)then
             !$omp parallel do
             do i = 1, n
                if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, x(i), y(im, i))
                y(im, i) = exp(y(im, i))
             enddo
             !$omp end parallel do
          else
             !$omp parallel do
             do i = 1, n
                if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, x(i), y(im, i))
             enddo
             !$omp end parallel do
          endif
       endif
    enddo
    deallocate(xc, yc, yc2)
  end subroutine coop_spline_fill_vd



  subroutine coop_spline_fill_vs(m, n, x, y, computed, logx, logy)
    COOP_INT m,n
    COOP_SINGLE x(n), y(m,n)
    logical computed(n), logx, logy
    COOP_INT nc, i, j
    COOP_REAL,dimension(:),allocatable::xc, yc, yc2
    COOP_REAL ytmp
    COOP_INT im
    nc = count(computed)
    if(nc.eq.n)return
    allocate(xc(nc), yc(nc), yc2(nc))
    j = 1
    do i = 1, n
       if(computed(i))then
          xc(j) = x(i)
          j = j + 1
       endif
    enddo
    if(logx) xc = log(xc)

    do im = 1, m
       j = 1
       do i = 1, n
          if(computed(i))then
             yc(j) = y(im, i)
             j = j + 1
          endif
       enddo
       if(logy) yc = log(yc)
       call coop_spline(nc, xc, yc, yc2)
       if(logx)then
          if(logy)then
             !$omp parallel do private(ytmp)
             do i = 1, n
                if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, log(COOP_REAL_OF(x(i))), ytmp)
                y(im, i) = exp(ytmp)
             enddo
             !$omp end parallel do
          else
             !$omp parallel do private(ytmp)
             do i = 1, n
                if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, log(COOP_REAL_OF(x(i))), ytmp)
                y(im, i) = ytmp
             enddo
             !$omp end parallel do 
          endif
       else
          if(logy)then
             !$omp parallel do private(ytmp)
             do i = 1, n
                if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, COOP_REAL_OF(x(i)), ytmp)
                y(im, i) = exp(ytmp)
             enddo
             !$omp end parallel do
          else
             !$omp parallel do private(ytmp)
             do i = 1, n
                if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, COOP_REAL_OF(x(i)), ytmp)
                y(im, i) = ytmp
             enddo
             !$omp end parallel do
          endif
       endif
    enddo
    deallocate(xc, yc, yc2)
  end subroutine coop_spline_fill_vs




  subroutine coop_spline_fill_md(m1, m2, n, x, y, computed, logx, logy)
    COOP_INT m1, m2, n
    COOP_REAL x(n), y(m1, m2,n)
    logical computed(n), logx, logy
    COOP_INT nc, i, j, im1, im2
    COOP_REAL,dimension(:),allocatable::xc, yc, yc2
    nc = count(computed)
    if(nc.eq.n)return
    allocate(xc(nc), yc(nc), yc2(nc))
    j = 1
    do i = 1, n
       if(computed(i))then
          xc(j) = x(i)
          j = j + 1
       endif
    enddo
    if(logx) xc = log(xc)

    do im1 = 1, m1; do im2=1,m2
       j = 1
       do i = 1, n
          if(computed(i))then
             yc(j) = y(im1, im2, i)
             j = j + 1
          endif
       enddo
       if(logy) yc = log(yc)
       call coop_spline(nc, xc, yc, yc2)
       if(logx)then
          if(logy)then
             !$omp parallel do
             do i = 1, n
                if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, log(x(i)), y(im1, im2, i))
                y(im1, im2, i) = exp(y(im1, im2, i))
             enddo
             !$omp end parallel do
          else
             !$omp parallel do
             do i = 1, n
                if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, log(x(i)), y(im1, im2, i))
             enddo
             !$omp end parallel do
          endif
       else
          if(logy)then
             !$omp parallel do
             do i = 1, n
                if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, x(i), y(im1, im2, i))
                y(im1, im2, i) = exp(y(im1, im2, i))
             enddo
             !$omp end parallel do
          else
             !$omp parallel do
             do i = 1, n
                if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, x(i), y(im1, im2, i))
             enddo
             !$omp end parallel do
          endif
       endif
    enddo; enddo
    deallocate(xc, yc, yc2)
  end subroutine coop_spline_fill_md



  subroutine coop_spline_fill_ms(m1, m2, n, x, y, computed, logx, logy)
    COOP_INT m1, m2, n
    COOP_SINGLE x(n), y(m1, m2,n)
    logical computed(n), logx, logy
    COOP_INT nc, i, j
    COOP_REAL,dimension(:),allocatable::xc, yc, yc2
    COOP_REAL ytmp
    COOP_INT im1, im2
    nc = count(computed)
    if(nc.eq.n)return
    allocate(xc(nc), yc(nc), yc2(nc))
    j = 1
    do i = 1, n
       if(computed(i))then
          xc(j) = x(i)
          j = j + 1
       endif
    enddo
    if(logx) xc = log(xc)

    do im1 = 1, m1; do im2=1, m2
       j = 1
       do i = 1, n
          if(computed(i))then
             yc(j) = y(im1, im2, i)
             j = j + 1
          endif
       enddo
       if(logy) yc = log(yc)
       call coop_spline(nc, xc, yc, yc2)
       if(logx)then
          if(logy)then
             !$omp parallel do private(ytmp)
             do i = 1, n
                if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, log(COOP_REAL_OF(x(i))), ytmp)
                y(im1, im2, i) = exp(ytmp)
             enddo
             !$omp end parallel do
          else
             !$omp parallel do private(ytmp)
             do i = 1, n
                if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, log(COOP_REAL_OF(x(i))), ytmp)
                y(im1, im2, i) = ytmp
             enddo
             !$omp end parallel do 
          endif
       else
          if(logy)then
             !$omp parallel do private(ytmp)
             do i = 1, n
                if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, COOP_REAL_OF(x(i)), ytmp)
                y(im1, im2, i) = exp(ytmp)
             enddo
             !$omp end parallel do
          else
             !$omp parallel do private(ytmp)
             do i = 1, n
                if(.not. computed(i)) call coop_splint(nc, xc, yc, yc2, COOP_REAL_OF(x(i)), ytmp)
                y(im1, im2, i) = ytmp
             enddo
             !$omp end parallel do
          endif
       endif
    enddo; enddo
    deallocate(xc, yc, yc2)
  end subroutine coop_spline_fill_ms

  !!input n, xa(n), ya(n), x
  !!return the interpolated/extrapolated value y and error estimation dy
  subroutine coop_ratint(n, xa,ya,x,y,dy)
    integer n
    COOP_REAL:: xa(n), ya(n), x, y, dy
    INTEGER :: m, ns
    COOP_REAL :: c(n),d(n),dd(n),h(n),t(n)
    COOP_REAL, parameter :: eps=1.d-25
    h=xa-x
    ns=coop_minloc(abs(h))
    y=ya(ns)
    if (x == xa(ns)) then
       dy=0.0
       RETURN
    end if
    c=ya
    d=ya+eps
    ns=ns-1
    do m=1,n-1
       t(1:n-m)=(xa(1:n-m)-x)*d(1:n-m)/h(1+m:n)
       dd(1:n-m)=t(1:n-m)-c(2:n-m+1)
       if (any(dd(1:n-m) == 0.0)) &
            stop 'failure in ratint'
       dd(1:n-m)=(c(2:n-m+1)-d(1:n-m))/dd(1:n-m)
       d(1:n-m)=c(2:n-m+1)*dd(1:n-m)
       c(1:n-m)=t(1:n-m)*dd(1:n-m)
       if (2*ns < n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       end if
       y=y+dy
    end do
  end subroutine coop_ratint

  

  function coop_ratval(x,cof,mm,kk) result(ratval)
    COOP_INT kk,mm, j
    COOP_REAL ratval,x,cof(mm+kk+1)
    !!From NR: Given mm, kk, and cof(1:mm+kk+1), evaluate and return the rational function (cof(1)+ cof(2)x + · · · + cof(mm+1)x^mm)/(1 + cof(mm+2)x + · · · + cof(mm+kk+1)xkk). INTEGER j
    COOP_REAL:: sumd,sumn
    sumn=cof(mm+1)
    do j=mm,1,-1
       sumn=sumn*x+cof(j)
    enddo
    sumd=0.d0
    do j=mm+kk+1,mm+2,-1
       sumd=(sumd+cof(j))*x
    enddo
    ratval=sumn/(1.d0+sumd)
    return
  end function coop_ratval

  subroutine coop_pade(cof,resid)
    COOP_REAL, DIMENSION(:), INTENT(INOUT) :: cof
    COOP_REAL, INTENT(OUT) :: resid
    COOP_INT :: k,n
    COOP_INT, DIMENSION((size(cof)-1)/2) :: indx
    COOP_REAL, PARAMETER ::big=1.d30
    COOP_REAL :: d,rr,rrold
    COOP_REAL, DIMENSION((size(cof)-1)/2) :: x,y,z
    COOP_REAL, DIMENSION((size(cof)-1)/2,(size(cof)-1)/2) :: q,qlu
    n=(size(cof)-1)/2
    x=cof(n+2:2*n+1)
    y=x
    do k=1,n
       q(:,k)=cof(n+2-k:2*n+1-k)
    end do
    qlu=q
    call coop_matrix_ludcmp(qlu,indx,d, n)
    call coop_matrix_lubksb(qlu,indx,x, n)
    rr=BIG
    do
       rrold=rr
       z=x
       call coop_mprove(q,qlu,indx,y,x, n)
       rr=sum((z-x)**2)
       if (rr >= rrold) exit
    end do
    resid=sqrt(rrold)
    do k=1,n
       y(k)=cof(k+1)-dot_product(z(1:k),cof(k:1:-1))
    end do
    cof(2:n+1)=y
    cof(n+2:2*n+1)=-z
  end subroutine coop_pade


  subroutine coop_mprove(a,alud,indx,b,x, n)
    COOP_INT::n
    COOP_REAL:: a(n,n),alud(n,n)
    COOP_INT :: indx(n)
    COOP_REAL::b(n), x(n), r(n)
    r=matmul(a,x)-b
    call coop_matrix_lubksb(alud,indx,r, n)
    x=x-r
  end subroutine coop_mprove

  subroutine coop_ratlsq_f(func, a, b, mm, kk, cof, dev)
    external func
    COOP_REAL::func
    COOP_REAL, INTENT(IN) :: a,b
    COOP_INT, INTENT(IN) :: mm,kk
    COOP_REAL :: cof(mm+kk+1)
    COOP_REAL :: dev
    COOP_INT::ncof, npt, npth
    COOP_INT::i
    COOP_INT, PARAMETER :: NPFAC=8
    COOP_REAL, DIMENSION((mm+kk+1)*NPFAC) :: xs, fs
    COOP_REAL::theta
    ncof=mm+kk+1
    npt=NPFAC*ncof
    theta=coop_pio2/(npt-1)    
    npth=npt/2
    xs(1:npth-1)=a+(b-a)*sin(theta*arth(0,1,npth-1))**2
    xs(npth:npt)=b-(b-a)*sin(theta*arth(npt-npth,-1,npt-npth+1))**2
    !$omp parallel do
    do i=1, npt
       fs(i)=func(xs(i))
    enddo
    !$omp end parallel do
    call coop_ratlsq_arr(npt, xs, fs, mm, kk, cof, dev)
  contains
    FUNCTION arth(first,increment,n)
      COOP_INT, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8      
      COOP_INT, INTENT(IN) :: first,increment
      COOP_INT, INTENT(IN) :: n
      COOP_REAL, DIMENSION(n) :: arth
      COOP_INT :: k,k2
      COOP_REAL :: temp
      if (n > 0) arth(1)=first
      if (n <= NPAR_ARTH) then
         do k=2,n
            arth(k)=arth(k-1)+increment
         end do
      else
         do k=2,NPAR2_ARTH
            arth(k)=arth(k-1)+increment
         end do
         temp=increment*NPAR2_ARTH
         k=NPAR2_ARTH
         do
            if (k >= n) exit
            k2=k+k
            arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
            temp=temp+temp
            k=k2
         end do
      end if
    END FUNCTION arth
  end subroutine coop_ratlsq_f



  subroutine coop_ratlsq_arr(npt, xs, fs, mm,kk,cof,dev)
    COOP_INT::npt
    COOP_INT, INTENT(IN) :: mm,kk
    COOP_REAL,intent(IN) :: xs(npt), fs(npt)
    COOP_REAL,intent(OUT) :: cof(mm+kk+1)
    COOP_REAL::coff(mm+kk+1)
    COOP_REAL, INTENT(OUT) :: dev
    COOP_INT, PARAMETER :: MAXIT=5
    COOP_REAL, PARAMETER :: BIG=1.d30
    COOP_INT :: it,ncof,i
    COOP_REAL :: devmax, e, dev_tiny
    COOP_REAL, DIMENSION(npt) :: bb,ee,wt
    COOP_REAL, DIMENSION(npt,mm+kk+1) :: u,temp
    COOP_REAL, parameter::eps = 1.d-10
    COOP_REAL::xsmax, fac
    xsmax = 1.d2**(-1.d0/max(kk, mm))*maxval(abs(xs))
    ncof=mm+kk+1
    dev=BIG
    wt=1.0
    ee=1.0
    e=0.0
    dev_tiny = maxval(abs(fs))*eps
    do it=1,MAXIT
       bb=wt*(fs+sign(e,ee))
       temp=geop(spread(1.d0,1,npt),xs/xsmax,ncof)
       u(:,1:mm+1)=temp(:,1:mm+1)*spread(wt,2,mm+1)
       u(:,mm+2:ncof)=-temp(:,2:ncof-mm)*spread(bb,2,ncof-mm-1)
       call coop_svd_least_square_one(npt, ncof, u, bb, coff)
       do i = 1, npt
          ee(i)=coop_ratval(xs(i)/xsmax,coff,mm,kk)-fs(i)
       enddo
       wt=abs(ee)
       devmax=maxval(wt)
       e=sum(wt)/npt
       if (devmax <= dev) then
          cof=coff
          dev=devmax
       end if
       if(dev .lt. dev_tiny)exit
    end do
    fac = 1.d0/xsmax
    do i=2, mm+1
       cof(i) = cof(i)*fac
       fac = fac/xsmax
    enddo
    fac = 1.d0/xsmax
    do i = mm+2, mm+kk+1
       cof(i) = cof(i)*fac
       fac = fac/xsmax
    enddo
!10  format (' ratlsq iteration=',i2,' max error=',1p,e10.3)
  contains
    

    FUNCTION geop(first,factor,n)
      COOP_INT, PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2      
      COOP_REAL, DIMENSION(:), INTENT(IN) :: first,factor
      COOP_INT, INTENT(IN) :: n
      COOP_REAL, DIMENSION(size(first),n) :: geop
      COOP_INT :: k,k2
      COOP_REAL, DIMENSION(size(first)) :: temp
      if (n > 0) geop(:,1)=first(:)
      if (n <= NPAR_GEOP) then
         do k=2,n
            geop(:,k)=geop(:,k-1)*factor(:)
         end do
      else
         do k=2,NPAR2_GEOP
            geop(:,k)=geop(:,k-1)*factor(:)
         end do
         temp=factor**NPAR2_GEOP
         k=NPAR2_GEOP
         do
            if (k >= n) exit
            k2=k+k
            geop(:,k+1:min(k2,n))=geop(:,1:min(k,n-k))*&
                 spread(temp,2,size(geop(:,1:min(k,n-k)),2))
            temp=temp*temp
            k=k2
         end do
      end if
    END FUNCTION geop
  end subroutine coop_ratlsq_arr

  subroutine coop_smooth_fit_free(this)
    class(coop_smooth_fit)::this
    this%n_up = -1
    this%n_down = -1
    if(allocated(this%c))deallocate(this%c)
    this%initialized = .false.
    this%center = 0.d0
    this%norm = 1.d0
    this%xlog = .false.
    this%ylog = .false.
    this%dev = 0.d0
  end subroutine coop_smooth_fit_free

  subroutine coop_smooth_fit_fit(this, n, x, y, dof, xlog, ylog)
    class(coop_smooth_fit)::this
    COOP_INT::dof    
    COOP_INT::n, i, nd,nu
    COOP_REAL,intent(IN)::x(n), y(n)
    COOP_REAL::ctry(dof)
    COOP_REAL::xs(n), ys(n)    
    logical,optional::xlog, ylog
    COOP_REAL::dev
    if(dof .lt. 1) stop "smooth_fit: dof must be > 0"
    if(dof .ge. n) stop "smooth_fit: dof must be < n"
    call this%free()
    allocate(this%c(dof))
    
    if(present(xlog))then
       this%xlog = xlog
    endif
    if(present(ylog))then
       this%ylog = ylog
    endif
    if(this%xlog)then
       if(any(x.le.0.d0))stop "smooth_fit_fit: xlog option conflits with x<0"
       xs = log(x)
    else
       xs = x
    endif
    if(this%ylog)then
       if(any(y.le.0.d0))stop "smooth_fit_fit: ylog option conflits with y<0"
       ys = log(y)
    else
       ys = y
    endif
    if(dof .eq. 1)then
       this%n_down = 0
       this%n_up = 0
       this%c = sum(ys)/n
       return
    endif
    this%center = sum(xs)/n
    xs = xs - this%center
    this%norm = 1.d2**(2.d0/dof)/maxval(abs(xs))
    xs = xs * this%norm
    this%dev = 1.d99
    do i = 0, dof-1
       nd = i
       nu = dof-1-nd
       call coop_ratlsq(n, xs, ys, nu, nd, ctry, dev)
       if(abs(dev).lt. abs(this%dev))then
          this%dev = dev
          this%n_up = nu
          this%n_down = nd
          this%c = ctry
       endif
    enddo
    this%initialized = .true.
  end subroutine coop_smooth_fit_fit

  function coop_smooth_fit_eval(this, x) result(y)
    class(coop_smooth_fit)::this    
    COOP_REAL::x, y, xs
    if(this%xlog)then
       xs = (log(x) - this%center)*this%norm
    else
       xs = (x - this%center)*this%norm
    endif
    y =  coop_ratval(xs, this%c, this%n_up, this%n_down)
    if(this%ylog)y=exp(y)
  end function coop_smooth_fit_eval

  

end module coop_interpolation_mod



