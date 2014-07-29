module coop_function_mod
  use coop_constants_mod
  use coop_basicutils_mod
  use coop_arguments_mod
  use coop_sort_mod
  implicit none
#include "constants.h"

  private

  public:: coop_function

  type coop_function
     COOP_INT method
     logical xlog, ylog, check_boundary
     COOP_INT n
     COOP_REAL xmin, xmax, dx, fleft, fright, slopeleft, sloperight
     COOP_REAL, dimension(:), allocatable::f, f2, f1
   contains
     procedure::set_boundary => coop_function_set_boundary
     procedure::init => coop_function_initialize
     procedure::init_NonUniform => coop_function_initialize_NonUniform
     procedure::eval => coop_function_evaluate  
     procedure::eval_bare => coop_function_evaluate_bare !!without log scaling
     procedure::integrate => coop_function_integrate
     procedure::derivative => coop_function_derivative
     procedure::derivative2 => coop_function_derivative2
     procedure::free => coop_function_free
  end type coop_function


  interface coop_function
     module procedure coop_function_constructor
  end interface coop_function

contains

  subroutine coop_function_set_boundary(this, fleft, fright, slopeleft, sloperight)
    class(coop_function)::this
    COOP_REAL,optional::fleft, fright, slopeleft, sloperight
    if(present(fleft))this%fleft = fleft
    if(present(fright))this%fright = fright
    if(present(slopeleft))this%slopeleft = slopeleft
    if(present(sloperight))this%sloperight = sloperight
  end subroutine coop_function_set_boundary


  function coop_function_constructor(f, xmin, xmax, xlog, ylog, args, method, check_boundary) result(cf)
    external f
    COOP_REAL f, xmin, xmax, dx, lnxmin, lnxmax
    logical,optional::xlog
    logical,optional::ylog
    type(coop_arguments),optional::args
    COOP_INT, optional::method
    logical,optional::check_boundary
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
       if(xmin .le. 0.d0 .or. xmax .le. 0.d0) stop "when xmin or xmax is negative you cannot set xlog = .true."
       lnxmin = log(xmin)
       lnxmax = log(xmax)
       dx = (lnxmax - lnxmin)/(coop_default_array_size - 1)
       if(present(args))then
          y(1) = f(xmin, args)
          y(coop_default_array_size) = f(xmax, args)
          !$omp parallel do
          do i=2, coop_default_array_size-1
             y(i) = f(exp(lnxmin + dx*(i-1)), args)
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
       if(present(args))then
          y(1) = f(xmin, args)
          y(coop_default_array_size) = f(xmax, args)
          !$omp parallel do
          do i=2, coop_default_array_size-1
             y(i) = f(xmin + dx*(i-1), args)
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
    if(coop_isnan(y))then
       write(*,*) "Cannot construct the function type: found f = NAN within the specified range."
       stop
    endif
    if(present(check_boundary))then
       if(present(method))then
          call cf%init(n=coop_default_array_size, xmin=xmin, xmax=xmax, f=y, method = method, xlog = cf%xlog, ylog = cf%ylog, check_boundary = check_boundary)
       else
          call cf%init(n=coop_default_array_size, xmin=xmin, xmax=xmax, f=y, method = COOP_INTERPOLATE_SPLINE, xlog = cf%xlog, ylog = cf%ylog, check_boundary = check_boundary)
       endif
    else
       if(present(method))then
          call cf%init(n=coop_default_array_size, xmin=xmin, xmax=xmax, f=y, method = method, xlog = cf%xlog, ylog = cf%ylog)
       else
          call cf%init(n=coop_default_array_size, xmin=xmin, xmax=xmax, f=y, method = COOP_INTERPOLATE_SPLINE, xlog = cf%xlog, ylog = cf%ylog)
       endif
    endif
  end function coop_function_constructor

  subroutine coop_function_free(this)
    class(coop_function):: this
    if(allocated(this%f))deallocate(this%f)
    if(allocated(this%f1))deallocate(this%f1)
    if(allocated(this%f2))deallocate(this%f2)
  end subroutine coop_function_free


  subroutine coop_function_initialize_NonUniform(this, x, f, xlog, ylog)
    class(coop_function):: this
    logical, optional::xlog, ylog
    COOP_REAL,dimension(:),intent(IN):: x, f
    COOP_INT i
    COOP_INT m, n, loc
    COOP_REAL_ARRAY::xx, ff
    COOP_REAL :: xc(size(x)), fc(size(f)), fc2(size(f)), res
    logical do_spline
    if(coop_isnan(f))then
       write(*,*) "Cannot construct the function type: found f = NAN within the specified range."
       stop
    endif
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
    m  = coop_getdim("coop_function%init_NonUniform", size(x), size(f))
    this%method = COOP_INTERPOLATE_SPLINE
    this%n = coop_default_array_size
    allocate(this%f(this%n), this%f2(this%n))
    if(this%xlog)then
       if(any(x.le.0.d0))stop "Error: cannot set xlog = .true. for x<0"
       xc = log(x)
    else
       xc  = x
    endif
    if(this%ylog)then
       fc = log(f)
    else
       fc = f  
    endif
    call coop_quicksortacc(xc, fc)
    this%xmin = xc(1)
    this%xmax = xc(m)
    this%dx = (this%xmax - this%xmin)/(this%n-1)
    if(this%dx  .eq. 0.) stop "All x's are equal. Cannot generate function for a single point."
    this%fleft = fc(1)    
    this%fright = fc(m)
    this%f(1) = fc(1)
    this%f(this%n) = fc(m)
    this%slopeleft= 0.d0
    this%sloperight = 0.d0
    if(m .lt. 16384)then
       if( minval(xc(2:m) - xc(1:m-1)) .gt. 1.d-2 * maxval(xc(2:m)-xc(1:m-1)))then
          do_spline = .true.
       else
          do_spline = .false.
       endif
    else
       do_spline = .false.
    endif
    if(do_spline)then
       call coop_spline(m, xc, fc, fc2)
       !$omp parallel do
       do i=2, coop_default_array_size-1
          call coop_splint(m, xc, fc, fc2, this%xmin + (i-1)*this%dx, this%f(i))
       enddo
       !$omp end parallel do
    else
       !$omp parallel do private(loc, res)
       do i = 2, coop_default_array_size - 1
          call coop_locate(m, xc, this%xmin + (i-1)*this%dx, loc, res)
          this%f(i) = fc(loc)*(1.d0-res) + fc(loc+1)*res
       enddo
       !$omp end parallel do
    endif
    call coop_spline_uniform(this%n, this%f, this%f2)
    this%check_boundary = .true.
  end subroutine coop_function_initialize_NonUniform

  subroutine coop_function_initialize(this, n, xmin, xmax, f, method, fleft, fright, slopeleft, sloperight, chebyshev_order, xlog, ylog, check_boundary)
    class(coop_function):: this
    logical, optional::xlog, ylog
    COOP_INT,intent(IN):: n
    COOP_REAL,intent(IN):: xmin, xmax
    COOP_REAL,intent(IN)::f(n)
    COOP_REAL, optional::fleft, fright, slopeleft, sloperight
    COOP_INT, optional::method
    COOP_INT, optional::chebyshev_order
    logical,optional::check_boundary
    if(coop_isnan(f))then
       write(*,*) "Cannot construct the function type: found f = NAN within the specified range."
       stop
    endif
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
    else
       allocate(this%f1(this%n))
    endif
    if(this%xlog)then
       if(xmin .le. 0.d0 .or. xmax .le. 0.d0) stop "Error: cannot set xlog = .true. for xmin<0 or xmax<0"
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
    if(this%ylog)then
       this%fleft = log(this%fleft)
       this%fright = log(this%fright)
    endif
    if(present(slopeleft))then
       this%slopeleft = slopeleft
    else
       this%slopeleft= 0.d0
    endif
    if(present(sloperight))then
       this%sloperight = sloperight
    else
       this%sloperight = 0.d0
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
       if(this%ylog)then
          call coop_chebfit_uniform(n, log(f), this%n, this%f)
       else
          call coop_chebfit_uniform(n, f, this%n, this%f)
       endif
       call coop_chebfit_derv(this%n, this%xmin, this%xmax, this%f, this%f1)
       call coop_chebfit_derv(this%n, this%xmin, this%xmax, this%f1, this%f2)
    end select
    if(present(check_boundary))then
       this%check_boundary = check_boundary
    else
       this%check_boundary = .true.
    endif
  end subroutine coop_function_initialize

  function coop_function_evaluate(this, x) result(f)
    class(coop_function):: this
    COOP_REAL x, f
    if(this%xlog)then
       f = coop_function_evaluate_bare(this, log(x))
    else
       f = coop_function_evaluate_bare(this, x)
    endif
    if(this%ylog)then
       f = exp(f)
    endif
  end function coop_function_evaluate

  function coop_function_evaluate_bare(this, x) result(f)
    class(coop_function):: this
    COOP_REAL x, f, a, b, xdiff
    COOP_INT l, r
    xdiff = x - this%xmin
    b = xdiff/this%dx + 1.d0
    l = floor(b)
    if(l .lt. 1)then
       if(this%check_boundary)then
          if(b.gt. 0.999999)then
             f = this%fleft
             return
          endif
          write(*,*) this%xlog, this%ylog, this%n, b
          write(*,*) x, ":", this%xmin, " -- ", this%xmax, " ---", this%dx
          write(*,*) "coop_function cannot be evaluated out of its boundary"
          stop
       endif
       f = this%fleft + this%slopeleft*xdiff
    elseif(l.ge. this%n)then
       if(this%check_boundary)then
          if(b.le.this%n+1.e-6)then
             f = this%fright
             return
          endif
          write(*,*) this%xlog, this%ylog, this%n, b
          write(*,*) x, ":", this%xmin, " -- ", this%xmax, " ---", this%dx
          write(*,*) "coop_function cannot be evaluated out of its boundary"
          stop
       endif
       xdiff = x - this%xmax
       f = this%fright + this%sloperight*xdiff
    else
       select case(this%method)
       case(COOP_INTERPOLATE_LINEAR)
          b = b - l
          f = this%f(l) * (1.d0-b) + this%f(l+1) * b 
       case(COOP_INTERPOLATE_QUDRATIC, COOP_INTERPOLATE_SPLINE)
          b = b - l
          r = l + 1
          a = 1. - b
          f = this%f(l) * a + this%f(r) * b + this%f2(l) * (a**2-1.)*a + this%f2(r)*(b**2-1.)*b
       case(COOP_INTERPOLATE_CHEBYSHEV)
          call coop_chebeval(this%n, this%xmin, this%xmax, this%f, x, f)
       case default
          stop "UNKNOWN interpolation method"
       end select
    endif
  end function coop_function_evaluate_bare

  
  !!simple integration 
  !!to be optimized
  function coop_function_integrate(this, a, b) result(integral)
    class(coop_function)::this
    integer,parameter::nsteps = 8192
    COOP_REAL a, b, integral, dx, y1, y2, ym, x1, x2, xm, lna, lnb, dxby2
    integer i
    if(this%xlog)then
       lna = log(a)
       lnb = log(b)
       dx = (lnb - lna)/nsteps
       dxby2 = dx/2.d0
       x1 = lna
       y1 = this%eval(exp(x1))*exp(x1)
       integral = y1
       do i=1, nsteps-1
          xm = x1 + dxby2
          x2 = x1 + dx
          ym = this%eval(exp(xm))*exp(xm)
          y2 = this%eval(exp(x2))*exp(x2)
          integral = integral + (2.d0*y2+ym*4.d0)
          x1 = x2
          y2 = y1
       enddo
       xm = x1 + dxby2
       x2 = x1 + dx
       ym = this%eval(exp(xm))*exp(xm)
       y2 = this%eval(exp(x2))*exp(x2)
       integral = (integral + (y2+ym*4.d0))*(dx/6.d0)
    else
       dx = (b - a)/nsteps
       dxby2 = dx/2.d0
       x1 = a
       y1 = this%eval(x1)
       integral = y1
       do i=1, nsteps-1
          xm = x1 + dxby2
          x2 = x1 + dx
          ym = this%eval(xm)
          y2 = this%eval(x2)
          integral = integral + (2.d0*y2+ym*4.d0)
          x1 = x2
          y2 = y1
       enddo
       xm = x1 + dxby2
       x2 = x1 + dx
       ym = this%eval(xm)
       y2 = this%eval(x2)
       integral = (integral + (y2+ym*4.d0))*(dx/6.d0)
    endif
  end function coop_function_integrate


!!to be optimized
  function coop_function_derivative(this, x) result(fp)
    class(coop_function)::this
    COOP_REAL x, fp, dx, xbare
    if(this%xlog)then
       xbare  = log(x)
    else
       xbare = x
    endif
    select case(this%method)
    case(COOP_INTERPOLATE_CHEBYSHEV)
       call coop_chebeval(this%n, this%xmin, this%xmax, this%f1, xbare, fp)
    case default
       dx = this%dx/4.d0
       fp = (this%eval_bare(xbare+dx) - this%eval_bare(xbare - dx))/(2.d0*dx)
    end select
    if(this%xlog) fp = fp/x
    if(this%ylog) fp = fp * this%eval(x)
  end function coop_function_derivative


  function coop_function_derivative2(this, x) result(fpp)
    class(coop_function)::this
    COOP_REAL x, fpp, f, dx, fp, fplus, fminus, xbare
    if(this%xlog)then
       xbare = log(x)
    else
       xbare = x
    endif
    select case(this%method)
    case(COOP_INTERPOLATE_CHEBYSHEV)
       if(this%ylog)call coop_chebeval(this%n, this%xmin, this%xmax, this%f, xbare, f)
       if(this%ylog .or. this%xlog) call coop_chebeval(this%n, this%xmin, this%xmax, this%f1, xbare, fp)
       call coop_chebeval(this%n, this%xmin, this%xmax, this%f2, xbare, fpp)
    case default
       dx = this%dx*0.9d0
       fplus = this%eval_bare(xbare+dx)
       fminus = this%eval_bare(xbare-dx)
       f = this%eval_bare(xbare)
       if(this%xlog .or. this%ylog) fp = (fplus - fminus)/(2.d0*dx)
       fpp = (fplus + fminus - 2.d0*f)/(dx**2)
    end select
    if(this%xlog)then
       fpp = (fpp  - fp)/x**2
       fp = fp/x
    endif
    if(this%ylog)then
       fpp = (fpp + fp**2) * exp(f)
    endif
  end function coop_function_derivative2


end module coop_function_mod
