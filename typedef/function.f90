module coop_type_function
  use coop_constants
  use coop_basicutils
  use coop_type_arguments
  implicit none
#include "constants.h"

  private

  public:: coop_function

  type coop_function
     COOP_INT method
     logical xlog, ylog
     COOP_INT n
     COOP_REAL xmin, xmax, dx, fleft, fright
     COOP_REAL, dimension(:), allocatable::f, f2
   contains
     procedure::init => coop_function_initialize
     procedure::eval => coop_function_evaluate
     procedure::free => coop_function_free
  end type coop_function


  interface coop_function
     procedure coop_function_constructor
  end interface coop_function

contains


  function coop_function_constructor(f, xmin, xmax, xlog, ylog, arg) result(cf)
    external f
    COOP_REAL f, xmin, xmax, dx, lnxmin, lnxmax
    logical,optional::xlog
    logical,optional::ylog
    type(coop_arguments),optional::arg
    COOP_REAL_ARRAY::y
    COOP_INT i
    type(coop_function) :: cf
    if(present(xlog))then
       cf%xlog = xlog
    else
       cf%xlog = .false.
    endif
    if(present(ylog))then
       cf%ylog = ylog
    else
       cf%ylog = .false.
    endif
    if(cf%xlog)then
       lnxmin = log(xmin)
       lnxmax = log(xmax)
       dx = (lnxmax - lnxmin)/(coop_default_array_size - 1)
       if(present(arg))then
          y(1) = f(xmin, arg)
          y(coop_default_array_size) = f(xmax, arg)
          !$omp parallel do
          do i=2, coop_default_array_size-1
             y(i) = f(exp(lnxmin + dx*(i-1)), arg)
          enddo
          !$omp end parallel do
       else
          y(1) = f(xmin)
          y(coop_default_array_size) = f(xmax)
          !$omp parallel do
          do i=2, coop_default_array_size-1
             y(i) = f(exp(lnxmin + dx*(i-1)))
          enddo
          !$omp end parallel do
       endif
    else
       dx = (xmax - xmin)/(coop_default_array_size - 1)
       if(present(arg))then
          y(1) = f(xmin, arg)
          y(coop_default_array_size) = f(xmax, arg)
          !$omp parallel do
          do i=2, coop_default_array_size-1
             y(i) = f(xmin + dx*(i-1), arg)
          enddo
          !$omp end parallel do
       else
          y(1) = f(xmin)
          y(coop_default_array_size) = f(xmax)
          !$omp parallel do
          do i=2, coop_default_array_size-1
             y(i) = f(xmin + dx*(i-1))
          enddo
          !$omp end parallel do
       endif
    endif
    call cf%init(coop_default_array_size, xmin, xmax, y, method = COOP_INTERPOLATE_SPLINE, xlog = xlog, ylog = ylog)
  end function coop_function_constructor

  subroutine coop_function_free(this)
    class(coop_function):: this
    if(allocated(this%f))deallocate(this%f)
    if(allocated(this%f2))deallocate(this%f2)
  end subroutine coop_function_free

  subroutine coop_function_initialize(this, n, xmin, xmax, f, method, fleft, fright, chebyshev_order, xlog, ylog)
    class(coop_function):: this
    logical, optional::xlog, ylog
    COOP_INT,intent(IN):: n
    COOP_REAL,intent(IN):: xmin, xmax, f(n)
    COOP_REAL, optional::fleft, fright
    COOP_INT, optional::method
    COOP_INT, optional::chebyshev_order
    call this%free()
    if(present(xlog))then
       this%xlog = xlog
    else
       this%xlog = .false.
    endif
    if(present(ylog))then
       this%ylog = ylog
    else
       this%ylog = .false.
    endif
    if(n .le. 1 .or. xmin .eq. xmax)then
       write(*,*) "coop function cannot be initialized for xmin = xmax"
       stop 
    endif
    if(present(method))then
       this%method = method
    else
       this%method = COOP_INTERPOLATE_LINEAR
    endif
    if(this%method .eq. COOP_INTERPOLATE_CHEBYSHEV)then
       if(present(chebyshev_order))then
          this%n = chebyshev_order + 1
       else
          this%n = coop_default_chebyshev_fit_order + 1
       endif
    else
       this%n = n
    endif
    allocate(this%f(this%n), this%f2(this%n))
    if(this%method .ne. COOP_INTERPOLATE_CHEBYSHEV)then
       if(this%ylog)then
          this%f = log(f)
       else
          this%f = f
       endif
    endif
    if(this%xlog)then
       this%xmin = log(xmin)
       this%xmax = log(xmax)
    else
       this%xmin = xmin
       this%xmax = xmax
    endif
    this%dx = (this%xmax - this%xmin)/(this%n-1)
    if(present(fleft))then
       this%fleft = fleft
    else
       this%fleft = f(1)
    endif
    if(present(fright))then
       this%fright = fright
    else
       this%fright = f(n)
    endif
    select case(this%method)
    case(COOP_INTERPOLATE_LINEAR)
       this%f2 = 0.
    case(COOP_INTERPOLATE_QUDRATIC)
       this%f2(2:n-1) = this%f(3:n) + this%f(1:n-2) - 2.*this%f(2:n-1)
       this%f2(1) = this%f2(2)
       this%f2(n) = this%f2(n-1)
       this%f2 = this%f2/6.
    case(COOP_INTERPOLATE_SPLINE)
       call coop_spline_uniform(this%n, this%f, this%f2)
    case(COOP_INTERPOLATE_CHEBYSHEV)
       call coop_chebfit_uniform(n, f, this%n, this%f)
    end select
  end subroutine coop_function_initialize

  function coop_function_evaluate(this, x) result(f)
    class(coop_function):: this
    COOP_REAL x, f, a, b
    COOP_INT l, r
    if(this%xlog)then
       b = (log(x) - this%xmin)/this%dx + 1.
    else
       b = (x - this%xmin)/this%dx + 1.
    endif
    l = floor(b)
    if(l .lt. 1)then
       f = this%fleft
       return
    endif
    if(l.ge. this%n)then
       f = this%fright
       return
    endif
    select case(this%method)
    case(COOP_INTERPOLATE_LINEAR, COOP_INTERPOLATE_QUDRATIC, COOP_INTERPOLATE_SPLINE)
       b = b - l
       r = l + 1
       a = 1. - b
       f = this%f(l) * a + this%f(r) * b + this%f2(l) * (a**2-1.)*a + this%f2(r)*(b**2-1.)*b
    case(COOP_INTERPOLATE_CHEBYSHEV)
       if(this%xlog)then
          call coop_chebeval(this%n, this%xmin, this%xmax, this%f, log(x), f)
       else
          call coop_chebeval(this%n, this%xmin, this%xmax, this%f, x, f)
       endif
    end select
    if(this%ylog)then
       f = exp(f)
    endif
  end function coop_function_evaluate


end module coop_type_function
