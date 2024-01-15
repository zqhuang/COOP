!!this module does low-spin hankel transform
!!for a given 2D picture f(r, \phi) with certain rotational symmetry (low m), you want to extract the Hankel transform
!! F_{lm} = \int_0^\infty r dr f_m(r),
!! where f_m(r) = \frac{1}{2\pi} \int_0^{\2pi} f(r,\phi) e^{-im\phi} d\phi is the spin-m mode

module coop_hankel_mod
  use coop_wrapper_typedef
  use coop_special_function_mod
  use coop_list_mod
  implicit none
#include "constants.h"

  private

#define BESSEL_J coop_bessj  
!! change to besjn for compilers that support besjn
#include "jmaxmin.h"
#include "jzeros.h"  

  public::coop_hankel_kernel, coop_2D_radial_decompose
  
  type coop_hankel_kernel
     COOP_INT::m
     COOP_INT::n = 0
     COOP_INT::index_cut1, index_cut2
     COOP_REAL,dimension(:),allocatable::kernel
     COOP_REAL,dimension(:),allocatable::x     
   contains
     procedure::free => coop_hankel_kernel_free
     procedure::init => coop_hankel_kernel_init
     procedure::integrate => coop_hankel_kernel_integrate
  end type coop_hankel_kernel

contains

  subroutine coop_hankel_kernel_init(this, m)
    COOP_INT, parameter::nbins = 101
    class(coop_hankel_kernel)::this
    COOP_INT::m, i, j
    type(coop_list_double)::xlist, klist
    call xlist%init()
    call klist%init()
    if(m.lt. 0 .or. m .gt. coop_jmaxmin_mmax) &
         call coop_return_error("COOP only support hankel transformation for 0<=m<="//COOP_STR_OF(coop_jmaxmin_mmax))
    if(m.eq.0)then
       call add_points(0.d0, coop_jzeros(1, m), nbins*5)
    endif
    do i=1, coop_jzeros_size - 3
       call add_points(coop_jzeros(i, m), coop_jzeros(i+1, m), nbins + nbins/i)
    enddo
    this%index_cut1 = xlist%n
    call add_points(coop_jzeros(coop_jzeros_size-2, m), coop_jzeros(coop_jzeros_size-1, m), nbins)    
    this%index_cut2 = xlist%n
    call add_points(coop_jzeros(coop_jzeros_size-1, m), coop_jzeros(coop_jzeros_size, m), nbins)
    this%n = xlist%n
    allocate(this%x(this%n), this%kernel(this%n))
    do i=1, this%n
       call xlist%get_element(i, this%x(i))
       call klist%get_element(i, this%kernel(i))       
    enddo
    call xlist%free()
    call klist%free()

  contains

    
    subroutine add_points(a, b, n)
      COOP_REAL::a, b, c, w
      COOP_INT,optional::n
      COOP_INT::i
      COOP_REAL::dx, x
      if(present(n))then
         dx = (b-a)/n
         x = a
         do i=1, n
            call get_weight(x, x+dx, c, w)
            call xlist%push(c)
            call klist%push(w)
            x = x + dx
         enddo
      else
         call get_weight(a, b, c, w)
         call xlist%push(c)
         call klist%push(w)
      endif
    end subroutine add_points

    subroutine  get_weight(a, b, c, weight)
      COOP_REAL::a, b, c, weight, x, dx, wx, val
      COOP_INT::i
      COOP_INT, parameter::n = 101
      x = a
      dx = (b-a)/n
      weight  = (sqrt(a)*BESSEL_J(m, a) + sqrt(b)*BESSEL_J(m,b))/2.d0
      wx = sqrt(b)*BESSEL_J(m,b)*(b-a)/2.d0
      do i = 1, n-1
         x = x+ dx
         val = sqrt(x)*BESSEL_J(m, x)
         weight  = weight + val
         wx = wx + val*(x-a)
      end do
      c = a + wx/weight
      weight = weight*dx*sqrt(c)
    end subroutine get_weight
    
  end subroutine coop_hankel_kernel_init
  
  subroutine coop_hankel_kernel_free(this)
    class(coop_hankel_kernel)::this
    COOP_DEALLOC(this%kernel)
    COOP_DEALLOC(this%x)    
    this%m = 0
    this%n = 0
  end subroutine coop_hankel_kernel_free
  

  function coop_hankel_kernel_integrate(this, f, l) result(s)
    class(coop_hankel_kernel)::this
    external f    
    COOP_REAL::f
    !!f gives f_m(r)
    COOP_REAL::ds1, ds2, s, l
    COOP_INT::i
    s = 0.d0
    ds1 = 0.d0
    ds2 = 0.d0
    !$omp parallel do reduction(+:s)
    do i=1, this%index_cut1
       s = s + this%kernel(i)*f(this%x(i)/l) 
    enddo
    !$omp end parallel do
    
    !$omp parallel do reduction(+:s)
    do i=this%index_cut1 + 1, this%index_cut2
       ds1 = ds1 + this%kernel(i)*f(this%x(i)/l) 
    enddo
    !$omp end parallel do

    !$omp parallel do reduction(+:s)
    do i=this%index_cut2 + 1, this%n
       ds2 = ds2 + this%kernel(i)*f(this%x(i)/l) 
    enddo
    !$omp end parallel do
    if(abs(ds1 - ds2) .gt. abs(ds1)*0.01)then
       s = s + ds1**2/(ds1-ds2)
    else
       s = s + ds1 + ds2 !!failure to do extrapolation
    endif
    s = s/l**2
  end function coop_hankel_kernel_integrate

  !!decompose f(r, phi) = \sum_{m=0}^{\infty} C_m(r)\cos{m\phi} + S_m(r)\sin{m\phi}
  subroutine coop_2D_radial_decompose(n, f, mmax, Cr, Sr, Ck, Sk)
    COOP_INT::n, mmax, m
    COOP_REAL::f(-n:n, -n:n)
    COOP_REAL::Cr(1:n, 0:mmax), Sr(1:n, 0:mmax)
    COOP_REAL,optional::Ck(1:n, 0:mmax), Sk(1:n, 0:mmax)
    COOP_INT::i, ix, iy, nsteps, j
    COOP_REAL::dtheta, r, theta, rx, ry, fv, l
    type(coop_hankel_kernel)::hankel
    nsteps = n*20
    dtheta = coop_2pi/nsteps
    Cr = 0.d0
    Sr = 0.d0
    do i=1, n
       r = dble(i)
       do j=0, nsteps-1
          theta = dtheta*j
          rx = r*cos(theta)
          ry = r*sin(theta)
          ix = min(floor(rx), n-1)
          iy = min(floor(ry), n-1)
          rx = rx - ix
          ry = ry - iy
          fv = (1.d0-rx)*(f(ix, iy)*(1.d0-ry) + f(ix, iy+1)*ry) &
               + rx * ( f(ix+1, iy)*(1.d0-ry) + f(ix+1, iy+1)*ry)
          Cr(i, 0) = Cr(i, 0) + fv
          do m = 1, mmax
             Cr(i, m) = Cr(i, m) + fv * cos(m*theta)
             Sr(i, m) = Sr(i, m) + fv * sin(m*theta)             
          enddo
       enddo       
    enddo
    Cr(1:n,0) = Cr(1:n, 0)/nsteps 
    Cr(1:n,1:mmax) =  Cr(1:n,1:mmax) * (2.d0/nsteps)
    Sr(1:n,1:mmax) =  Sr(1:n,1:mmax) * (2.d0/nsteps)
    if(present(Ck) .and. present(Sk))then
       Sk(:,0) = 0.d0
       do m=0, mmax
          call hankel%init(m)
          do i=1, n
             l = dble(i)/n
             Ck(i,m) = hankel%integrate(crtemp, l)
             if(m.ne.0)Sk(i, m)= hankel%integrate(srtemp, l)
          enddo
          call hankel%free()
       enddo
    endif
  contains
    function crtemp(r)
      COOP_REAL::r, crtemp
      COOP_INT::ir
      COOP_REAL::s
      ir = floor(r)
      if(ir .lt. n)then
         s = r - ir
         crtemp = cr(ir, m)*(1.d0-s) + cr(ir+1, m)*s - cr(n, m)
      else
         crtemp = 0.d0
      endif
    end function crtemp

    function srtemp(r)
      COOP_REAL::r, srtemp
      COOP_INT::ir
      COOP_REAL::s
      ir = floor(r)      
      if(ir .lt. n)then
         s = r - ir
         srtemp = sr(ir, m)*(1.d0-s) + sr(ir+1, m)*s - sr(n, m)
      else
         srtemp = 0.d0
      endif
    end function srtemp
    
  end subroutine coop_2D_radial_decompose


  
end module coop_hankel_mod


