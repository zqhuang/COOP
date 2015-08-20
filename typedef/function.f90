module coop_function_mod
  use coop_constants_mod
  use coop_basicutils_mod
  use coop_arguments_mod
  use coop_sort_mod
  implicit none
#include "constants.h"

  private

  public:: coop_function, coop_function_multeval, coop_function_multeval_bare, coop_function_constructor

  type coop_function
     COOP_SHORT_STRING::name="NoName"
     logical::initialized = .false.
     COOP_INT::method = COOP_INTERPOLATE_QUADRATIC  
     logical::xlog =.false.
     logical::ylog = .false.
     logical::check_boundary = .true.
     COOP_INT::n = 0
     COOP_REAL xmin, xmax, dx, fleft, fright, slopeleft, sloperight
     COOP_REAL, dimension(:), allocatable::f, f2, f1
   contains
     procedure::set_boundary => coop_function_set_boundary
     procedure::init => coop_function_init
     procedure::init_polynomial => coop_function_init_polynomial
     procedure::init_NonUniform => coop_function_init_NonUniform
     procedure::eval => coop_function_evaluate  
     procedure::eval_bare => coop_function_evaluate_bare !!without log scaling
     procedure::derivative_bare => coop_function_derivative_bare !!without log scaling
     procedure::derivative2_bare => coop_function_derivative2_bare !!without logc scaling
     procedure::integrate => coop_function_integrate
     procedure::derivative => coop_function_derivative
     procedure::derivative2 => coop_function_derivative2
     procedure::maxloc => coop_function_maxloc
     procedure::minloc => coop_function_minloc
     procedure::free => coop_function_free
     procedure::monotonic_solution => coop_function_monotonic_solution
  end type coop_function

!!$  interface coop_function
!!$     module procedure coop_function_constructor
!!$  end interface coop_function  

contains


  function coop_function_constructor(f, xmin, xmax, xlog, ylog, args, method, check_boundary,name) result(cf)
    external f
    COOP_UNKNOWN_STRING,optional::name
    COOP_REAL f, xmin, xmax, dx, lnxmin, lnxmax
    logical,optional::xlog
    logical,optional::ylog
    type(coop_arguments),optional::args
    COOP_INT, optional::method
    logical,optional::check_boundary
    COOP_REAL_ARRAY::y
    COOP_INT i
    type(coop_function) :: cf
    if(present(name))then
       cf%name = trim(adjustl(name))
    endif
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
       write(*,*) "Cannot construct the function  "//trim(cf%name)//": found f = NAN within the specified range."
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
    cf%initialized = .true.
  end function coop_function_constructor
  

  subroutine coop_function_set_boundary(this, fleft, fright, slopeleft, sloperight)
    class(coop_function)::this
    COOP_REAL,optional::fleft, fright, slopeleft, sloperight
    if(present(fleft))this%fleft = fleft
    if(present(fright))this%fright = fright
    if(present(slopeleft))this%slopeleft = slopeleft
    if(present(sloperight))this%sloperight = sloperight
    this%check_boundary = .false.
  end subroutine coop_function_set_boundary



  subroutine coop_function_free(this)
    class(coop_function):: this
    if(allocated(this%f))deallocate(this%f)
    if(allocated(this%f1))deallocate(this%f1)
    if(allocated(this%f2))deallocate(this%f2)
    this%initialized = .false.
  end subroutine coop_function_free


  subroutine coop_function_init_polynomial(this, p, xlog, ylog,  name)
    class(coop_function)::this
    logical,optional::xlog, ylog
    COOP_REAL,dimension(:)::p
    COOP_UNKNOWN_STRING,optional::name    
    COOP_INT::i
    call this%free()
    this%method = COOP_INTERPOLATE_POLYNOMIAL
    if(present(name))this%name = trim(adjustl(name))
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
    this%check_boundary = .false.
    this%n = size(p)
    allocate(this%f(this%n), this%f1(this%n), this%f2(this%n))
    this%f = p
    this%f1(this%n) = 0.d0
    this%f2(this%n) = 0.d0
    if(this%n .gt.1) &
         this%f2(this%n-1) = 0.d0        
    do i = 1, this%n-1
       this%f1(i) = this%f(i+1)*i
    enddo
    do i=1, this%n - 2
       this%f2(i) = this%f1(i+1)*i
    enddo
    this%xmax = 1.d99**(1.d0/max(this%n-1,1))
    this%xmin = - this%xmax
    this%initialized = .true.        
  end subroutine coop_function_init_polynomial

  subroutine coop_function_init_NonUniform(this, x, f, xlog, ylog, check_boundary, smooth, name)
    class(coop_function):: this
    logical, optional::xlog, ylog
    COOP_REAL,dimension(:),intent(IN):: x, f
    COOP_INT i
    COOP_INT m, n, loc
    COOP_REAL_ARRAY :: xx, ff
    COOP_REAL :: xc(size(x)), fc(size(f)), fc2(size(f)), res
    logical do_spline
    logical,optional::check_boundary
    logical,optional::smooth
    COOP_UNKNOWN_STRING,optional::name
    if(present(name))then
       this%name = trim(adjustl(name))
    endif    
    if(coop_isnan(f))then
       write(*,*) "Cannot construct the function "//trim(this%name)//": found f = NAN within the specified range."
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
    this%method = COOP_INTERPOLATE_QUADRATIC
    m  = coop_getdim("coop_function%init_NonUniform", size(x), size(f))
    this%n = min(coop_default_array_size, max(m, 128))
    allocate(this%f(this%n), this%f2(this%n))
    if(this%xlog)then
       if(any(x.le.0.d0))stop "Error: cannot set xlog = .true. for x<0"
       xc = log(x)
    else
       xc  = x
    endif
    if(this%ylog)then
       if(any(f .le. 0.d0)) stop "Error: cannot set ylog = .true. for f(x) < 0"
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
    if(m .lt. 8192)then
       if( minval(xc(2:m) - xc(1:m-1)) .gt. 2.d-2 * maxval(xc(2:m)-xc(1:m-1)))then
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
       do i=2, this%n - 1
          call coop_splint(m, xc, fc, fc2, this%xmin + (i-1)*this%dx, this%f(i))
       enddo
       !$omp end parallel do
    else
       !$omp parallel do private(loc, res)
       do i = 2, this%n - 1
          call coop_locate(m, xc, this%xmin + (i-1)*this%dx, loc, res)
          this%f(i) = fc(loc)*(1.d0-res) + fc(loc+1)*res
       enddo
       !$omp end parallel do
    endif
    this%f2(2:this%n-1) = this%f(3:this%n) + this%f(1:this%n-2) - 2.*this%f(2:this%n-1)
    this%f2(1) = this%f2(2)
    this%f2(this%n) = this%f2(this%n-1)
    this%f2 = this%f2/6.d0
    if(present(smooth))then
       if(smooth)then
          if(this%n .ge. 200)then  !!check f2 is smooth
             call coop_smooth_data(this%n, this%f2, min(this%n/200, 50))
          endif
       endif
    endif
    if(present(check_boundary))then
       this%check_boundary = check_boundary
    endif
    this%initialized = .true.
  end subroutine coop_function_init_NonUniform

  subroutine coop_function_init(this, n, xmin, xmax, f, method, fleft, fright, slopeleft, sloperight, chebyshev_order, xlog, ylog, check_boundary, smooth, name)
    class(coop_function):: this
    logical, optional::xlog, ylog
    COOP_INT,intent(IN):: n
    COOP_REAL,intent(IN):: xmin, xmax
    COOP_REAL,intent(IN)::f(n)
    COOP_REAL, optional::fleft, fright, slopeleft, sloperight
    COOP_INT, optional::method
    COOP_INT, optional::chebyshev_order
    logical,optional::check_boundary
    logical,optional::smooth
    COOP_UNKNOWN_STRING, optional::name
    COOP_INT i, count_tiny, count_small
    COOP_REAL::fmean, ftiny, curv, flarge, fsmall
    if(present(name))then
       this%name = trim(adjustl(name))
    endif    
    if(coop_isnan(f))then
       write(*,*) "Cannot construct the function "//trim(this%name)//": found f = NAN within the specified range."
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
    if(this%xlog .and. (xmin .le. 0.d0 .or. xmax .le. 0.d0)) call coop_return_error("function:init", "cannot do xlog for x<=0", "stop")
    if(this%ylog)then
       if(any(f.le.0.d0)) call coop_return_error("function:init", "cannot do ylog for y<=0", "stop")
    endif
    if(n .le. 1 .or. xmin .eq. xmax)then
       write(*,*) "coop function cannot be initialized for xmin = xmax"
       stop 
    endif
    if(n.eq.2)then
       this%method = COOP_INTERPOLATE_LINEAR
    elseif(present(method))then
       this%method = method
    else
       !!automatically choose the best method
       allocate(this%f(n))
       if(this%ylog)then
          this%f = log(f)
       else
          this%f = f
       endif
       fmean = sum(this%f)/n
       flarge = sqrt(sum((this%f-fmean)**2)/n)/n
       ftiny = flarge*1.d-7
       fsmall = flarge*1.d-4

       count_small = 0
       count_tiny = 0
       do i=1, n - 3
          curv = abs(this%f(i) - 2.d0*this%f(i+1) + this%f(i+2))+ abs( this%f(i+1) -2.d0 * this%f(i+2) + this%f(i+3))
          if(curv .gt. flarge)then
             this%method = COOP_INTERPOLATE_LINEAR
             goto 100
          endif
          if(curv .lt. fsmall)then
             count_small = count_small + 1
             if(curv .lt.ftiny)then
                count_tiny = count_tiny + 1
                if(count_tiny .gt. 5)then
                   this%method = COOP_INTERPOLATE_QUADRATIC
                   goto 100
                endif
             endif
             if(count_small .gt. 20)then
                this%method = COOP_INTERPOLATE_QUADRATIC
                goto 100
             endif
          endif
       enddo
       this%method = COOP_INTERPOLATE_SPLINE
    endif
100 continue
    if(this%method .eq. COOP_INTERPOLATE_CHEBYSHEV)then
       if(present(chebyshev_order))then
          this%n = chebyshev_order + 1
       else
          this%n = coop_default_chebyshev_fit_order + 1
       endif
    else
       this%n = n
    endif
    if(.not. allocated(this%f))then
       allocate(this%f(this%n))
       if(this%method .ne. COOP_INTERPOLATE_CHEBYSHEV)then
          if(this%ylog)then
             this%f = log(f)
          else
             this%f = f
          endif
       endif
    endif
    allocate(this%f2(this%n))
    if(this%method .eq. COOP_INTERPOLATE_CHEBYSHEV) allocate(this%f1(this%n))

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
    case(COOP_INTERPOLATE_QUADRATIC)
       this%f2(2:n-1) = this%f(3:n) + this%f(1:n-2) - 2.*this%f(2:n-1)
       this%f2(1) = this%f2(2)
       this%f2(n) = this%f2(n-1)
       this%f2 = this%f2/6.
       if(present(smooth))then
          if(smooth)then
             if(this%n .ge. 200)then  !!check f2 is smooth
                call coop_smooth_data(this%n, this%f2, min(this%n/200, 50))
             endif
          endif
       endif
    case(COOP_INTERPOLATE_SPLINE)
       call coop_spline_uniform(this%n, this%f, this%f2)
       if(present(smooth))then
          if(smooth)then
             if(this%n .ge. 200)then  !!check f2 is smooth
                call coop_smooth_data(this%n, this%f2, min(this%n/200, 50))
             endif
          endif
       endif
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
    this%initialized = .true.
  end subroutine coop_function_init





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
    if(this%method .eq. COOP_INTERPOLATE_POLYNOMIAL)then
       f = coop_polyvalue(this%n, this%f, x)
       return
    endif
    xdiff = x - this%xmin
    b = xdiff/this%dx + 1.d0
    l = floor(b)
    if(l .lt. 1)then
       if(this%check_boundary)then
          if(b.gt. 0.9999d0)then
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
          if(b .le. dble(this%n)+1.d-4)then
             f = this%fright
             return
          endif
          write(*,*) "function "//trim(this%name)//": boundary overflow (from the left)"
          stop
       endif
       xdiff = x - this%xmax
       f = this%fright + this%sloperight*xdiff
    else
       select case(this%method)
       case(COOP_INTERPOLATE_LINEAR)
          b = b - l
          f = this%f(l) * (1.d0-b) + this%f(l+1) * b 
       case(COOP_INTERPOLATE_QUADRATIC, COOP_INTERPOLATE_SPLINE)
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


  function coop_function_multeval(nf, funcs, x) result(f)
    COOP_INT nf
    type(coop_function)::funcs(nf)
    COOP_REAL x, f(nf)
    if(funcs(1)%xlog)then
       f = coop_function_multeval_bare(nf, funcs, log(x))
    else
       f = coop_function_multeval_bare(nf, funcs, x)
    endif
    if(funcs(1)%ylog)then
       f = exp(f)
    endif
  end function coop_function_multeval


  function coop_function_multeval_bare(nf, funcs, x) result(f)
    COOP_INT nf
    type(coop_function)::funcs(nf)
    COOP_REAL x, f(nf), a, b, xdiff
    COOP_INT l, r, i
    xdiff = x - funcs(1)%xmin
    b = xdiff/funcs(1)%dx + 1.d0
    l = floor(b)
    if(l .lt. 1)then
       if(funcs(1)%check_boundary)then
          if(b.gt. 0.9999d0)then
             do i=1, nf
                f(i) = funcs(i)%fleft
             enddo
             return
          endif
          write(*,*) funcs(1)%xlog, funcs(1)%ylog, funcs(1)%n, b
          write(*,*) x, ":", funcs(1)%xmin, " -- ", funcs(1)%xmax, " ---", funcs(1)%dx
          write(*,*) "coop_function cannot be evaluated out of its boundary"
          stop
       endif
       do i=1, nf
          f(i) = funcs(i)%fleft + funcs(i)%slopeleft*xdiff
       enddo
       return
    elseif(l.ge. funcs(1)%n)then
       if(funcs(1)%check_boundary)then
          if(b .le. dble(funcs(1)%n)+1.d-4)then
             do i=1, nf
                f(i) = funcs(i)%fright
             enddo
             return
          endif
          write(*,*) funcs(1)%xlog, funcs(1)%ylog, funcs(1)%n, b
          write(*,*) x, ":", funcs(1)%xmin, " -- ", funcs(1)%xmax, " ---", funcs(1)%dx
          write(*,*) "coop_function cannot be evaluated out of its boundary"
          stop
       endif
       xdiff = x - funcs(1)%xmax
       do i = 1, nf
          f(i) = funcs(i)%fright + funcs(1)%sloperight*xdiff
       enddo
    else
       select case(funcs(1)%method)
       case(COOP_INTERPOLATE_LINEAR)
          b = b - l
          do i=1, nf
             f(i) = funcs(i)%f(l) * (1.d0-b) + funcs(i)%f(l+1) * b 
          enddo
       case(COOP_INTERPOLATE_QUADRATIC, COOP_INTERPOLATE_SPLINE)
          b = b - l
          r = l + 1
          a = 1. - b
          do i=1, nf
             f(i) = funcs(i)%f(l) * a + funcs(i)%f(r) * b + funcs(i)%f2(l) * (a**2-1.)*a + funcs(i)%f2(r)*(b**2-1.)*b
          enddo
       case(COOP_INTERPOLATE_CHEBYSHEV)
          do i=1, nf
             call coop_chebeval(funcs(i)%n, funcs(i)%xmin, funcs(i)%xmax, funcs(i)%f, x, f(i))
          enddo
       case default
          stop "UNKNOWN interpolation method"
       end select
    endif
  end function coop_function_multeval_bare

  function coop_function_derivative_bare(this, x) result(f)
    class(coop_function):: this
    COOP_REAL x, f, a, b, xdiff
    COOP_INT l, r
    if(this%method .eq. COOP_INTERPOLATE_POLYNOMIAL)then
       f = coop_polyvalue(this%n-1, this%f1, x)
       return
    endif
    xdiff = x - this%xmin
    b = xdiff/this%dx + 1.d0
    l = floor(b)
    if(l .lt. 1)then
       if(this%check_boundary)then
          if(b .gt. 0.9999d0)then
             b = 1.d0
             l = 1
          else
             write(*,*) this%xlog, this%ylog, this%n, b
             write(*,*) x, ":", this%xmin, " -- ", this%xmax, " ---", this%dx
             write(*,*) "coop_function cannot be evaluated out of its boundary"
             stop
          endif
       else
          f = this%slopeleft
          return
       endif
    elseif(l.ge. this%n)then
       if(this%check_boundary)then
          if(b .le. dble(this%n)+1.d-4)then
             b = this%n
             l = this%n-1
          else
             write(*,*) this%xlog, this%ylog, this%n, b
             write(*,*) x, ":", this%xmin, " -- ", this%xmax, " ---", this%dx
             write(*,*) "coop_function cannot be evaluated out of its boundary"
             stop
          endif
       else
          f = this%sloperight
          return
       endif
    endif
    
    select case(this%method)
    case(COOP_INTERPOLATE_LINEAR)
       if(l.eq.1.or.l.eq.this%n-1)then
          f = (this%f(l+1)-this%f(l))/this%dx
       else
          f = (this%eval_bare(x+0.5*this%dx)-this%eval_bare(x-0.5*this%dx))/this%dx
       endif
    case(COOP_INTERPOLATE_QUADRATIC, COOP_INTERPOLATE_SPLINE)
       b = b - l
       r = l + 1
       a = 1. - b
       f =  ((this%f(r) - this%f(l)) + ( this%f2(r)*(3.d0*b**2-1.) - this%f2(l) * (3.d0*a**2-1.) ))/this%dx
    case(COOP_INTERPOLATE_CHEBYSHEV)
       call coop_chebeval(this%n, this%xmin, this%xmax, this%f1, x, f)
    case default
       call coop_return_error("coop_functiion_derivative_bare", "UNKNOWN interpolation method", "stop")
    end select

  end function coop_function_derivative_bare


  function coop_function_derivative2_bare(this, x) result(f)
    class(coop_function):: this
    COOP_REAL x, f, a, b, xdiff
    COOP_INT l, r
    if(this%method .eq. COOP_INTERPOLATE_POLYNOMIAL)then
       f = coop_polyvalue(this%n-2, this%f2, x)
       return
    endif
    xdiff = x - this%xmin
    b = xdiff/this%dx + 1.d0
    l = floor(b)
    if(l .lt. 1)then
       if(this%check_boundary)then
          if(b .gt. 0.9999d0)then
             b = 1.d0
             l = 1
          else
             write(*,*) this%xlog, this%ylog, this%n, b
             write(*,*) x, ":", this%xmin, " -- ", this%xmax, " ---", this%dx
             write(*,*) "coop_function cannot be evaluated out of its boundary"
             stop
          endif
       else
          f = 0.d0
          return
       endif
    elseif(l.ge. this%n)then
       if(this%check_boundary)then
          if(b .le. dble(this%n)+1.d-4)then
             b = this%n
             l = this%n-1
          else
             write(*,*) this%xlog, this%ylog, this%n, b
             write(*,*) x, ":", this%xmin, " -- ", this%xmax, " ---", this%dx
             write(*,*) "coop_function cannot be evaluated out of its boundary"
             stop
          endif
       else
          f = 0.d0
          return
       endif
    endif
    
    select case(this%method)
    case(COOP_INTERPOLATE_LINEAR)
       if(l .eq. 1)then
          f = (this%f(3)+this%f(1) - 2.d0*this%f(2))/this%dx**2
          return
       elseif(l.eq.this%n-1)then
          f = (this%f(this%n)+this%f(this%n-2) - 2.d0*this%f(this%n-1))/this%dx**2
          return
       else
          f = (this%eval_bare(x+this%dx)+this%eval_bare(x-this%dx)-2.d0*this%eval_bare(x))/this%dx**2
       endif

    case(COOP_INTERPOLATE_QUADRATIC, COOP_INTERPOLATE_SPLINE)
       b = b - l
       r = l + 1
       a = 1. - b
       f =  ( this%f2(r)*b + this%f2(l) * a )*6.d0/this%dx**2
    case(COOP_INTERPOLATE_CHEBYSHEV)
       call coop_chebeval(this%n, this%xmin, this%xmax, this%f2, x, f)
    case default
       call coop_return_error("coop_functiion_derivative2_bare", "UNKNOWN interpolation method", "stop")
    end select

  end function coop_function_derivative2_bare


  function coop_function_maxloc(this) result(xm)
    class(coop_function)::this
    COOP_REAL xm
    integer iloc, n, iloop
    COOP_REAL xmin, xmax, x, fm, dx, f
    select case(this%method)
    case(COOP_INTERPOLATE_CHEBYSHEV)
       xmin = this%xmin
       xmax = this%xmax
       n = 512
    case default
       iloc = coop_maxloc(this%f)
       xmin = this%xmin + this%dx*max(iloc-2, 0)
       xmax = this%xmin + this%dx*min(iloc, this%n-1)
       n = 16
    end select
    do iloop = 1, 3
       dx = (xmax-xmin)/n
       xm = xmin+dx 
       fm = this%eval_bare(xm)
       do iloc = 2, n-1
          x = xmin + iloc*dx
          f= this%eval_bare(x)
          if(f .gt. fm)then
             fm = f
             xm = x
          endif
       enddo
       xmin = xm - dx
       xmax = xm + dx
       n = 8
    enddo
    if(this%xlog) xm = exp(xm)
  end function coop_function_maxloc


  function coop_function_minloc(this) result(xm)
    class(coop_function)::this
    COOP_REAL xm
    integer iloc, n, iloop
    COOP_REAL xmin, xmax, x, fm, dx, f
    select case(this%method)
    case(COOP_INTERPOLATE_CHEBYSHEV)
       xmin = this%xmin
       xmax = this%xmax
       n = 512
    case default
       iloc = coop_minloc(this%f)
       xmin = this%xmin + this%dx*max(iloc-2, 0)
       xmax = this%xmin + this%dx*min(iloc, this%n-1)
       n = 16
    end select
    do iloop = 1, 3
       dx = (xmax-xmin)/n
       xm = xmin+dx 
       fm = this%eval_bare(xm)
       do iloc = 2, n-1
          x = xmin + iloc*dx
          f= this%eval_bare(x)
          if(f .lt. fm)then
             fm = f
             xm = x
          endif

       enddo
       xmin = xm - dx
       xmax = xm + dx
       n = 8
    enddo
    if(this%xlog) xm = exp(xm)
  end function coop_function_minloc
  
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
    COOP_REAL x, fp
    if(this%xlog)then
       fp  = this%derivative_bare(log(x))/x
    else
       fp = this%derivative_bare(x)
    endif
    if(this%ylog) fp = fp * this%eval(x)
  end function coop_function_derivative


  function coop_function_derivative2(this, x) result(fpp)
    class(coop_function)::this
    COOP_REAL x, fpp, xbare, fp, f
    if(this%xlog)then
       xbare = log(x)
    else
       xbare = x
    endif
    fpp = this%derivative2_bare(xbare)
    if(this%ylog)f = this%eval_bare(xbare)
    if(this%ylog .or. this%xlog) fp = this%derivative_bare(xbare)
    if(this%xlog)then
       fpp = (fpp  - fp)/x**2
       fp = fp/x
    endif
    if(this%ylog)then
       fpp = (fpp + fp**2) * exp(f)
    endif
  end function coop_function_derivative2

  subroutine coop_function_monotonic_solution(this, f, x)
    class(coop_function)::this
    COOP_REAL,intent(in)::f
    COOP_REAL, intent(out)::x
    COOP_REAL xmin, xmax, fmin, fmax, xmid, fs, dx
    COOP_INT imin, imax, imid
    if(this%ylog)then
       if(f.le.0.d0)return
       fs = log(f)
    else
       fs = f
    endif
    select case(this%method)
    case(COOP_INTERPOLATE_CHEBYSHEV, COOP_INTERPOLATE_POLYNOMIAL)
       xmin = this%xmin
       xmax = this%xmax
       fmin = this%eval_bare(xmin)
       fmax = this%eval_bare(xmax)
       if((fmin .lt. fs .and. fmax .lt. fs) .or. (fmin .gt. fs .and. fmax .gt. fs))then
          call coop_return_error("coop_function_monotonic_solution", "solution does not exist or the function is not monotonic", "stop")
       endif
       dx = abs(xmax - xmin)*1.d-12
       if(fmax .ge. fmin)then
          do while(abs(xmax - xmin) .gt. dx)
             xmid = (xmax+xmin)/2.d0
             if(this%eval_bare(xmid) .gt. fs)then
                xmax = xmid
             else
                xmin = xmid
             endif
          enddo
       else
          do while(abs(xmax - xmin) .gt. dx)
             xmid = (xmax+xmin)/2.d0
             if(this%eval_bare(xmid) .lt. fs)then
                xmax = xmid
             else
                xmin = xmid
             endif
          enddo
       endif
       x = (xmin + xmax)/2.d0
    case(COOP_INTERPOLATE_LINEAR, COOP_INTERPOLATE_QUADRATIC, COOP_INTERPOLATE_SPLINE)
       imin = 1
       imax = this%n
       if((this%f(imin) .lt. fs .and. this%f(imax) .lt. fs) .or. (this%f(imin) .gt. fs .and. this%f(imax) .gt. fs))then
          call coop_return_error("coop_function_monotonic_solution", "solution does not exist or the function is not monotonic", "stop")
       endif
       if(this%f(imin) .lt. this%f(imax))then
          do while(imax - imin .gt. 1)
             imid = (imax+imin)/2
             if(this%f(imid) .gt. fs)then
                imax = imid
             else
                imin = imid
             endif
          enddo
          if(this%method .eq. COOP_INTERPOLATE_LINEAR)then
             dx = (this%f(imax)- this%f(imin))
             if(dx .ne. 0.d0)then
                x = this%xmin + (imin-1 + (fs - this%f(imin))/dx)*this%dx
             else
                x = this%xmin + (imin-0.5d0)*this%dx
             endif
          else
             xmax = this%xmin + (imax-1)*this%dx
             xmin = this%xmin + (imin-1)*this%dx
             dx = abs(this%dx) * 1.d-8
             do while(abs(xmax - xmin) .gt. dx)
                xmid = (xmax+xmin)/2.d0
                if(this%eval_bare(xmid) .gt. fs)then
                   xmax = xmid
                else
                   xmin = xmid
                endif
             enddo
             x = (xmin + xmax)/2.d0
          endif
       else
          do while(imax - imin .gt. 1)
             imid = (imax+imin)/2
             if(this%f(imid) .gt. fs)then
                imin = imid
             else
                imax = imid
             endif
          enddo
          if(this%method .eq. COOP_INTERPOLATE_LINEAR)then
             dx = (this%f(imax)- this%f(imin))
             if(dx .ne. 0.d0)then
                x = this%xmin + (imin-1 + (fs - this%f(imin))/dx)*this%dx
             else
                x = this%xmin + (imin-0.5d0)*this%dx
             endif
          else
             xmax = this%xmin + (imax-1)*this%dx
             xmin = this%xmin + (imin-1)*this%dx
             dx = abs(this%dx) * 1.d-8
             do while(abs(xmax - xmin) .gt. dx)
                xmid = (xmax+xmin)/2.d0
                if(this%eval_bare(xmid) .lt. fs)then
                   xmax = xmid
                else
                   xmin = xmid
                endif
             enddo
             x = (xmin + xmax)/2.d0
          endif
       endif
    case default
       call coop_return_error("coop_function_find_solution", "Unknown interpolation method", "stop")
    end select
    if(this%xlog) x = exp(x)
  end subroutine coop_function_monotonic_solution


end module coop_function_mod


