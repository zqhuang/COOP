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
  
#include "jmaxmin.h"
#include "jzeros.h"  

  public::coop_hankel_kernel
  
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
      weight  = (sqrt(a)*besjn(m, a) + sqrt(b)*besjn(m,b))/2.d0
      wx = sqrt(b)*besjn(m,b)*(b-a)/2.d0
      do i = 1, n-1
         x = x+ dx
         val = sqrt(x)*besjn(m, x)
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
  

end module coop_hankel_mod


