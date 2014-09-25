module coop_interpolation_mod
  use coop_wrapper_typedef
  use coop_matrix_mod
  implicit none

#include "constants.h"
  private


  public::coop_bilinear_interp, coop_bicubic_interp, coop_linear_interp, coop_linear_least_square_fit, coop_spline_fill

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
       if(sxx .eq. 0.d0 .or. syy .eq. 0.d0)then
          r = 1.d0
       else
          r = sxy**2/sxx/syy
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


end module coop_interpolation_mod



