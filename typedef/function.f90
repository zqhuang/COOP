module coop_function_mod
  use coop_constants_mod
  use coop_basicutils_mod
  use coop_arguments_mod
  use coop_sort_mod
  implicit none
#include "constants.h"

  private

  public:: coop_function, coop_function_multeval, coop_function_multeval_bare, coop_function_constructor, coop_function_polynomial, coop_2dfunction

  type coop_function
     COOP_SHORT_STRING::name="NoName"
     logical::initialized = .false.
     COOP_INT::method = COOP_INTERPOLATE_LINEAR
     logical::xlog =.false.
     logical::ylog = .false.
     logical::check_boundary = .true.
     COOP_INT::n = 0
     COOP_INT::n_down = 0
     COOP_REAL::scale = 1.d0   !!f -> scale * f + shift
     COOP_REAL::shift = 0.d0
     COOP_REAL xmin, xmax, dx, fleft, fright, slopeleft, sloperight
     COOP_REAL, dimension(:), allocatable::f, f1, f2, f3
   contains
     procedure::set_boundary => coop_function_set_boundary
     procedure::init => coop_function_init
     procedure::mult_const => coop_function_mult_const
     procedure::add_const => coop_function_add_const
     procedure::init_polynomial => coop_function_init_polynomial
     procedure::init_rational => coop_function_init_rational
     procedure::init_powerlaw => coop_function_init_powerlaw
     procedure::init_NonUniform => coop_function_init_NonUniform
     procedure::init_zigzag => coop_function_init_zigzag
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


  type coop_2dfunction
     COOP_SHORT_STRING::name = "NoName2D"
     logical::initialized = .false.
     COOP_INT::nx = 0
     COOP_INT::ny = 0
     COOP_REAL::xmin=0.d0
     COOP_REAL::xmax=1.d0
     COOP_REAL::ymin=0.d0
     COOP_REAL::ymax=1.d0
     COOP_REAL::dx = 0.d0
     COOP_REAL::dy = 0.d0
     logical::xlog = .false.
     logical::ylog = .false.
     logical::zlog = .false.
     COOP_REAL,dimension(:,:),allocatable::f, fxx, fyy
   contains
     procedure::init => coop_2dfunction_init
     procedure::init_symmetric => coop_2dfunction_init_symmetric
     procedure::free => coop_2dfunction_free
     procedure::eval_bare => coop_2dfunction_evaluate_bare
     procedure::eval => coop_2dfunction_evaluate
  end type coop_2dfunction

contains
  
  function coop_2dfunction_evaluate_bare(this, xbare, ybare) result(zbare)
    class(coop_2dfunction)::this
    COOP_REAL::xbare, ybare, zbare
    COOP_INT::ix, iy
    COOP_REAL::rx, ry
    rx = (xbare - this%xmin)/this%dx+1.d0
    ry = (ybare - this%ymin)/this%dy+1.d0
    ix = min(max(floor(rx), 1), this%nx-1)
    rx = max(min(rx - ix, 1.d0), 0.d0)
    iy = min(max(floor(ry), 1), this%ny-1)
    ry = max(min(ry - iy, 1.d0), 0.d0)
    zbare =(this%f(ix, iy)*(1.d0-ry) + this%f(ix, iy+1)*ry)*(1.d0-rx) &
         + (this%f(ix+1, iy)*(1.d0-ry) + this%f(ix+1, iy+1)*ry)*rx &
         + ((this%fxx(ix, iy)*(1.d0-ry) + this%fxx(ix, iy+1)*ry)*(1.d0-rx) &
         + (this%fxx(ix+1, iy)*(1.d0-ry) + this%fxx(ix+1, iy+1)*ry)*rx)*rx*(rx-1.d0) &
         + ((this%fyy(ix, iy)*(1.d0-ry) + this%fyy(ix, iy+1)*ry)*(1.d0-rx) &
         + (this%fyy(ix+1, iy)*(1.d0-ry) + this%fyy(ix+1, iy+1)*ry)*rx)*ry*(ry-1.d0)
  end function coop_2dfunction_evaluate_bare


  function coop_2dfunction_evaluate(this, x, y) result(z)
    class(coop_2dfunction)::this
    COOP_REAL::x, y, z
    if(this%xlog)then
       if(this%ylog)then
          z = this%eval_bare(log(x), log(y))
       else
          z = this%eval_bare(log(x), y)
       endif
    else
       if(this%ylog)then
          z = this%eval_bare(x, log(y))
       else
          z = this%eval_bare(x, y)
       endif
    endif
    if(this%zlog) z = exp(z)
  end function coop_2dfunction_evaluate


  subroutine coop_2dfunction_free(this)
    class(coop_2dfunction)::this
    COOP_DEALLOC(this%f)
    COOP_DEALLOC(this%fxx)
    COOP_DEALLOC(this%fyy)
    this%xlog = .false.
    this%ylog = .false.
    this%zlog = .false.
    this%nx  = 0
    this%ny = 0
    this%initialized = .false.
    this%name = "NoName2D"
  end subroutine coop_2dfunction_free

  subroutine coop_2dfunction_init(this, f, xmin, xmax, ymin, ymax, nx, ny, xlog, ylog, zlog, name)
    !!save f(x, y) into coop_2dfunction object this
    class(coop_2dfunction)::this
    external f
    COOP_REAL::f
    COOP_INT::nx, ny
    COOP_REAL::xmin, xmax, ymin, ymax
    logical,optional::xlog,ylog,zlog
    COOP_UNKNOWN_STRING,optional::name
    COOP_INT::ix, iy
    COOP_REAL::dx, dy
    call this%free()
    if(present(name))this%name = trim(adjustl(name))
    if(nx .lt. 2 .or. ny .lt. 2) stop "2d function init assumes nx >=2 and ny >=2,"
    this%nx = nx
    this%ny = ny
    allocate(this%f(nx, ny), this%fxx(nx, ny), this%fyy(nx, ny))
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
    if(present(zlog))then
       this%zlog = zlog
    else
       this%zlog = .false.
    endif
    if(this%xlog)then
       if(xmin .le. 0.d0 .or. xmax .le. 0.d0) stop "2dfunction_init: xlog = T but xmin<=0 or xmax<=0"
       this%xmin = log(xmin)
       this%xmax = log(xmax)
    else
       this%xmin = xmin
       this%xmax = xmax
    endif
    if(this%ylog)then
       if(ymin .le. 0.d0 .or. ymax .le. 0.d0) stop "2dfunction_init: ylog = T but ymin<=0 or ymax<=0"
       this%ymin = log(ymin)
       this%ymax = log(ymax)
    else
       this%ymin = ymin
       this%ymax = ymax
    endif
    this%dx = (this%xmax - this%xmin)/(this%nx-1)
    this%xmin = this%xmin + this%dx/1.d10
    this%xmax = this%xmax - this%dx/1.d10
    this%dx = (this%xmax - this%xmin)/(this%nx-1)

    this%dy = (this%ymax - this%ymin)/(this%ny-1)
    this%ymin = this%ymin + this%dy/1.d10
    this%ymax = this%ymax - this%dy/1.d10
    this%dy = (this%ymax - this%ymin)/(this%ny-1)
    if(abs(this%dx) .le. 0.d0) stop "2dfunction_init: xmax == xmin?"
    if(abs(this%dy) .le. 0.d0) stop "2dfunction_init: ymax == ymin?"

    if(this%xlog)then
       if(this%ylog)then
          !$omp parallel do private(ix, iy)
          do iy = 1, this%ny
             do ix = 1, this%nx
                this%f(ix, iy) = f(exp(this%xmin + this%dx*(ix-1)), exp(this%ymin + this%dy*(iy-1)))
             enddo
          enddo
          !$omp end parallel do
       else
          !$omp parallel do private(ix, iy)
          do iy = 1, this%ny
             do ix = 1, this%nx
                this%f(ix, iy) = f(exp(this%xmin + this%dx*(ix-1)), this%ymin + this%dy*(iy-1))
             enddo
          enddo
          !$omp end parallel do
       endif
    else
       if(this%ylog)then
          !$omp parallel do private(ix, iy)
          do iy = 1, this%ny
             do ix = 1, this%nx
                this%f(ix, iy) = f(this%xmin + this%dx*(ix-1), exp(this%ymin + this%dy*(iy-1)))
             enddo
          enddo
          !$omp end parallel do
       else
          !$omp parallel do private(ix, iy)
          do iy = 1, this%ny
             do ix = 1, this%nx
                this%f(ix, iy) = f(this%xmin + this%dx*(ix-1), this%ymin + this%dy*(iy-1))
             enddo
          enddo
          !$omp end parallel do
       endif
    endif
    if(this%zlog)then
       if(any(this%f .le. 0.d0)) stop "2dfunction_init: negative values"
       this%f = log(this%f)
    endif
    !$omp parallel do private(ix, iy)
    do iy = 2, this%ny-1
       do ix = 1, this%nx
          this%fyy(ix, iy) = (this%f(ix, iy+1) + this%f(ix, iy-1) - 2.d0*this%f(ix, iy))
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(ix, iy)
    do iy = 1, this%ny
       do ix = 2, this%nx-1
          this%fxx(ix, iy) = (this%f(ix+1, iy) + this%f(ix-1, iy) - 2.d0*this%f(ix, iy))
       enddo
    enddo
    !$omp end parallel do


    if(this%nx .gt. 2)then
       this%fxx(1, :) = this%fxx(2, :)
       this%fxx(this%nx, :) = this%fxx(this%nx-1, :)
    else
       this%fxx(1, :) = 0.d0
       this%fxx(this%nx, :) = 0.d0
    endif
    if(this%ny .gt. 2)then
       this%fyy(:, 1) = this%fyy(:, 2)
       this%fyy(:, this%ny) = this%fyy(:, this%ny-1)
    else
       this%fyy(:,1) = 0.d0
       this%fyy(:,this%ny) = 0.d0
    endif
    this%fxx = this%fxx/2.d0
    this%fyy = this%fyy/2.d0
    this%initialized = .true.
    return
  end subroutine coop_2dfunction_init


  subroutine coop_2dfunction_init_symmetric(this, f, xmin, xmax, nx, xlog, zlog, name)
    !!save f(x, y) = f(y, x) into coop_2dfunction object this
    class(coop_2dfunction)::this
    external f
    COOP_REAL::f
    COOP_INT::nx
    COOP_REAL::xmin, xmax
    logical,optional::xlog, zlog
    COOP_UNKNOWN_STRING,optional::name
    COOP_INT::ix, iy
    COOP_REAL::dx, dy
    call this%free()
    if(present(name))this%name = trim(adjustl(name))
    if(nx .lt. 2 ) stop "2d function init assumes nx >=2"
    this%nx = nx
    this%ny = nx
    allocate(this%f(nx, nx), this%fxx(nx, nx), this%fyy(nx, nx))
    if(present(xlog))then
       this%xlog = xlog
    else
       this%xlog = .false.
    endif
    this%ylog = this%xlog
    if(present(zlog))then
       this%zlog = zlog
    else
       this%zlog = .false.
    endif
    if(this%xlog)then
       if(xmin .le. 0.d0 .or. xmax .le. 0.d0) stop "2dfunction_init_symmetric: xlog = T but xmin<=0 or xmax<=0"
       this%xmin = log(xmin)
       this%xmax = log(xmax)
    else
       this%xmin = xmin
       this%xmax = xmax
    endif

    this%dx = (this%xmax - this%xmin)/(this%nx-1)
    this%xmin = this%xmin + this%dx/1.d10
    this%xmax = this%xmax - this%dx/1.d10
    this%dx = (this%xmax - this%xmin)/(this%nx-1)
    if(abs(this%dx) .le. 0.d0) stop "2dfunction_init_symmetric: xmax == xmin?"

    this%dy = this%dx
    this%ymin = this%xmin
    this%ymax = this%xmax



    if(this%xlog)then
       !$omp parallel do private(ix, iy)
       do iy = 1, this%ny
          do ix = 1, iy
             this%f(ix, iy) = f(exp(this%xmin + this%dx*(ix-1)), exp(this%ymin + this%dy*(iy-1)))
          enddo
       enddo
       !$omp end parallel do
    else
       !$omp parallel do private(ix, iy)
       do iy = 1, this%ny
          do ix = 1, iy
             this%f(ix, iy) = f(this%xmin + this%dx*(ix-1), this%ymin + this%dy*(iy-1))
          enddo
       enddo
       !$omp end parallel do
    endif
    !$omp parallel do private(ix, iy)
    do iy = 1, this%ny
       do ix = 1, iy-1
          this%f(iy, ix) = this%f(ix, iy)
       enddo
    enddo
    !$omp end parallel do

    if(this%zlog)then
       if(any(this%f .le. 0.d0)) stop "2dfunction_init: negative values"
       this%f = log(this%f)
    endif
    !$omp parallel do private(ix, iy)
    do iy = 2, this%ny-1
       do ix = 1, this%nx
          this%fyy(ix, iy) = (this%f(ix, iy+1) + this%f(ix, iy-1) - 2.d0*this%f(ix, iy))
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(ix, iy)
    do iy = 1, this%ny
       do ix = 2, this%nx-1
          this%fxx(ix, iy) = (this%f(ix+1, iy) + this%f(ix-1, iy) - 2.d0*this%f(ix, iy))
       enddo
    enddo
    !$omp end parallel do
    this%fxx(1, :) = this%fyy(:, 1)
    this%fxx(this%nx, :) = this%fyy(:, this%ny)
    this%fyy(:,1) = this%fxx(1,:)
    this%fyy(:,this%ny) = this%fxx(this%nx, :)

    this%fxx(1, 1) = (this%fxx(2, 1)+this%fxx(1, 2))/2.d0
    this%fxx(1, this%nx) = (this%fxx(2, this%nx) + this%fxx(1, this%nx-1))/2.d0
    this%fxx(this%nx, 1) = (this%fxx(this%nx, 2)+this%fxx(this%nx-1, 1))/2.d0
    this%fxx(this%nx, this%nx) = (this%fxx(this%nx-1, this%nx) + this%fxx(this%nx, this%nx-1))/2.d0


    this%fyy(1, 1) = (this%fyy(2, 1)+this%fyy(1, 2))/2.d0
    this%fyy(1, this%nx) = (this%fyy(2, this%nx) + this%fyy(1, this%nx-1))/2.d0
    this%fyy(this%nx, 1) = (this%fyy(this%nx, 2)+this%fyy(this%nx-1, 1))/2.d0
    this%fyy(this%nx, this%nx) = (this%fyy(this%nx-1, this%nx) + this%fyy(this%nx, this%nx-1))/2.d0

    this%fxx = this%fxx/2.d0
    this%fyy = this%fyy/2.d0
    this%initialized = .true.
    return
  end subroutine coop_2dfunction_init_symmetric

  subroutine coop_function_mult_const(this, c)
    class(coop_function)::this
    COOP_REAL::c
    this%scale = this%scale * c
    this%shift = this%shift * c
  end subroutine coop_function_mult_const

  subroutine coop_function_add_const(this, c)
    class(coop_function)::this
    COOP_REAL::c
    this%shift = this%shift + c
  end subroutine coop_function_add_const



  function coop_function_constructor(f, xmin, xmax, xlog, ylog, args, method, check_boundary,name) result(cf)
    external f
    COOP_UNKNOWN_STRING,optional::name
    COOP_REAL f, xmin, xmax, dx, lnxmin, lnxmax
    logical,optional::xlog
    logical,optional::ylog
    logical::doxlog, doylog
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
       doxlog = xlog
    else
       doylog = .false.
    endif
    if(present(ylog))then
       doxlog = ylog
    else
       doylog = .false.
    endif
    if(doxlog)then
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
          call cf%init(n=coop_default_array_size, xmin=xmin, xmax=xmax, f=y, method = method, xlog = doxlog, ylog = doylog, check_boundary = check_boundary)
       else
          call cf%init(n=coop_default_array_size, xmin=xmin, xmax=xmax, f=y, method = COOP_INTERPOLATE_SPLINE, xlog = doxlog, ylog = doylog, check_boundary = check_boundary)
       endif
    else
       if(present(method))then
          call cf%init(n=coop_default_array_size, xmin=xmin, xmax=xmax, f=y, method = method, xlog = doxlog, ylog = doylog)
       else
          call cf%init(n=coop_default_array_size, xmin=xmin, xmax=xmax, f=y, method = COOP_INTERPOLATE_SPLINE, xlog = doxlog, ylog = doylog)
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
    COOP_DEALLOC(this%f)
    COOP_DEALLOC(this%f1)
    COOP_DEALLOC(this%f2)
    COOP_DEALLOC(this%f3)
    this%method = COOP_INTERPOLATE_LINEAR
    this%name = "NoName"
    this%check_boundary = .true.
    this%n = 0
    this%n_down  = 0
    this%scale = 1.d0
    this%shift = 0.d0
    this%xlog = .false.
    this%ylog = .false.
    this%initialized = .false.
  end subroutine coop_function_free


  !!initialize a function
  !! y = c_1 x^(alpha_1) + c_2 x^(alpha_2) + ...
  subroutine coop_function_init_powerlaw(this, c, alpha, xlog, ylog, name)
    class(coop_function)::this
    logical,optional::xlog, ylog
    COOP_REAL,dimension(:)::alpha, c
    COOP_UNKNOWN_STRING,optional::name    
    COOP_INT::i
    call this%free()
    this%method = COOP_INTERPOLATE_POWERLAW
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
    this%n = coop_getdim("init_powerlaw", size(c), size(alpha))
    allocate(this%f(this%n), this%f1(this%n), this%f2(this%n))
    this%f = c
    this%f1 = alpha
    this%f2 = c*alpha
    this%xmax = 1.d99
    this%xmin = -1.d99
    this%initialized = .true.        
  end subroutine coop_function_init_powerlaw


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

!!this is generalized rational function
!!f = sum(c_up * x**alpha_up)/sum(c_down * x**alpha_down)
  subroutine coop_function_init_rational(this, c_up, alpha_up, c_down, alpha_down,  xlog, ylog,  name)
    class(coop_function)::this
    logical,optional::xlog, ylog
    COOP_REAL,dimension(:)::c_up, alpha_up, c_down, alpha_down
    COOP_UNKNOWN_STRING,optional::name    
    COOP_INT::i
    call this%free()
    this%method = COOP_INTERPOLATE_RATIONAL
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
    this%n = coop_GetDim("init_rational", size(c_up), size(alpha_up))
    this%n_down = coop_getDim("init_rational", size(c_down), size(alpha_down))
    allocate(this%f(this%n), this%f1(this%n), this%f2(this%n_down), this%f3(this%n_down))
    this%f = c_up
    this%f1 = alpha_up
    this%f2 = c_down
    this%f3 = alpha_down
    this%xmax = 1.d99
    this%xmin = - 1.d99
    this%initialized = .true.        
  end subroutine coop_function_init_rational


  subroutine coop_function_init_zigzag(this, x, f, xlog, ylog, name)
    class(coop_function):: this    
    COOP_REAL,dimension(:),intent(IN):: x, f
    COOP_UNKNOWN_STRING,optional::name
    logical, optional::xlog, ylog
    call this%free()
    if(present(name))then
       this%name = trim(adjustl(name))
    endif    
    if(coop_isnan(f))then
       write(*,*) "Cannot construct the function "//trim(this%name)//": found f = NAN within the specified range."
       stop
    endif
    this%n  = coop_getdim("coop_function%init_ZigZag", size(x), size(f))
    if(this%n.eq.1)then
       call this%init_polynomial(f)
       return
    endif
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
    this%method = COOP_INTERPOLATE_ZIGZAG
    allocate(this%f(this%n), this%f2(this%n))
    if(this%xlog)then
       if(any(x.le.0.d0))stop "Error: cannot set xlog = .true. for x<0"
       this%f2 = log(x)
    else
       this%f2  = x
    endif
    if(this%ylog)then
       if(any(f .le. 0.d0)) stop "Error: cannot set ylog = .true. for f(x) < 0"
       this%f = log(f)
    else
       this%f = f  
    endif
    this%xmin = minval(this%f2)
    this%xmax = maxval(this%f2)
    this%initialized = .true.
  end subroutine coop_function_init_zigzag

  
  subroutine coop_function_init_NonUniform(this, x, f, xlog, ylog, check_boundary, name)
    class(coop_function):: this
    logical, optional::xlog, ylog
    COOP_REAL,dimension(:),intent(IN):: x, f
    COOP_INT i
    logical,optional::check_boundary
    COOP_UNKNOWN_STRING,optional::name
    call this%free()
    if(present(name))then
       this%name = trim(adjustl(name))
    endif    
    if(coop_isnan(f))then
       write(*,*) "Cannot construct the function "//trim(this%name)//": found f = NAN within the specified range."
       stop
    endif
    this%n = coop_getdim("coop_function_init_NonUniform", size(x), size(f))
    if(this%n .eq. 1)then
       call this%init_polynomial( f )
       return
    endif
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
    this%method = COOP_INTERPOLATE_NONUNIFORM
    allocate(this%f(this%n), this%f2(this%n), this%f1(this%n))
    if(this%xlog)then
       if(any(x.le.0.d0))stop "Error: cannot set xlog = .true. for x<0"
       this%f1 = log(x)
    else
       this%f1 = x
    endif
    if(this%ylog)then
       if(any(f .le. 0.d0)) stop "Error: cannot set ylog = .true. for f(x) < 0"
       this%f = log(f)
    else
       this%f = f  
    endif
    call coop_quicksortacc(this%f1, this%f)
    this%xmin = this%f1(1)
    this%xmax = this%f1(this%n)
    this%dx = (this%xmax - this%xmin)/(this%n-1)
    if(this%dx  .eq. 0.) stop "All x's are equal. Cannot generate function for a single point."
    this%fleft = this%f(1) 
    this%fright = this%f(this%n)
    if(this%n .gt. 3)then
       if(this%f1(2) .gt. this%f1(1) .and. this%f1(3) .gt. this%f1(2) )then
          this%slopeleft = ( (this%f(2) - this%f(1))/(this%f1(2) - this%f1(1)) *(this%f1(3) - this%f1(1)) - (this%f(3) - this%f(1))/(this%f1(3) - this%f1(1))*(this%f1(2) - this%f1(1)) )/(this%f1(3) - this%f1(2))
       elseif(this%f1(2) .gt. this%f1(1))then
          this%slopeleft =  (this%f(2) - this%f(1))/(this%f1(2) - this%f1(1))
       else
          this%slopeleft = 0.d0
       endif
       if( this%f1(this%n) .gt. this%f1(this%n-1) .and. this%f1(this%n-1) .gt. this%f1(this%n-2))then
          this%sloperight = ( (this%f(this%n-1) - this%f(this%n))/(this%f1(this%n-1) - this%f1(this%n)) *(this%f1(this%n-2) - this%f1(this%n)) - (this%f(this%n-2) - this%f(this%n))/(this%f1(this%n-2) - this%f1(this%n))*(this%f1(this%n-1) - this%f1(this%n)) )/(this%f1(this%n-2) - this%f1(this%n-1))
       elseif( this%f1(this%n) .gt. this%f1(this%n-1))then
          this%sloperight = (this%f(this%n-1) - this%f(this%n))/(this%f1(this%n-1) - this%f1(this%n))
       else
          this%sloperight = 0.d0
       endif
       call coop_spline(this%n, this%f1, this%f, this%f2, this%slopeleft, this%sloperight)
    else
       this%f2 = 0.d0
    endif
100 if(present(check_boundary))then
       this%check_boundary = check_boundary
    else
       this%check_boundary = .false.
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
    call this%free()
    if(present(name))then
       this%name = trim(adjustl(name))
    endif    
    if(coop_isnan(f))then
       write(*,*) "Cannot construct the function "//trim(this%name)//": found f = NAN within the specified range."
       stop
    endif
    if(n.eq.1)then
       call this%init_polynomial (f)
       return
    endif
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
    f = f*this%scale + this%shift
  end function coop_function_evaluate

  function coop_function_evaluate_bare(this, x) result(f)
    class(coop_function):: this
    COOP_REAL x, f, a, b, xdiff
    COOP_INT l, r
    select case(this%method)
    case(COOP_INTERPOLATE_POWERLAW)
       f = sum(this%f * x**this%f1, this%f .ne. 0.d0)
       return
    case(COOP_INTERPOLATE_RATIONAL)
       f = sum(this%f * x**this%f1, this%f .ne. 0.d0)/sum(this%f2 * x**this%f3, this%f2 .ne. 0.d0)
       return
    case(COOP_INTERPOLATE_POLYNOMIAL)
       f = coop_polyvalue(this%n, this%f, x)
       return
    case(COOP_INTERPOLATE_ZIGZAG)
       l = coop_left_index(this%n, this%f2, x)
       if(l.le.0)then
          f = this%f(1)
          return
       endif
       if(l.ge.this%n)then
          f = this%f(this%n)
          return
       endif
       f = (this%f(l+1) * (x - this%f2(l)) + this%f(l) * (this%f2(l+1) - x))/(this%f2(l+1) - this%f2(l))
       return
    case(COOP_INTERPOLATE_NONUNIFORM)
       call coop_splint(this%n, this%f1, this%f, this%f2, x, f)
       return
    end select
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
    if(funcs(1)%ylog) &
       f = exp(f)
    f = f*funcs%scale + funcs%shift
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
          do i=1, nf
             f(i) = funcs(i)%eval(x)
          enddo
       end select
    endif
  end function coop_function_multeval_bare

  function coop_function_derivative_bare(this, x) result(f)
    class(coop_function):: this
    COOP_REAL x, f, a, b, xdiff
    COOP_INT l, r
    select case(this%method)
    case(COOP_INTERPOLATE_POWERLAW)
       f = sum(this%f2*x**(this%f1-1.d0), this%f2 .ne. 0.d0)
       return
    case(COOP_INTERPOLATE_RATIONAL)
       a = sum(this%f*x**this%f1, this%f .ne. 0.d0)
       b = sum(this%f2*x**this%f3, this%f2 .ne. 0.d0)
       f = (sum(this%f*this%f1*x**(this%f1-1.d0), this%f .ne. 0.d0 .and. this%f1 .ne. 0.d0)*b - sum(this%f2*this%f3*x**(this%f3-1.d0), this%f2 .ne. 0.d0 .and. this%f3 .ne. 0.d0)*a)/b**2
       return
    case(COOP_INTERPOLATE_POLYNOMIAL)
       f = coop_polyvalue(this%n-1, this%f1, x)
       return
    case(COOP_INTERPOLATE_ZIGZAG)
       if(x .lt. this%xmin .or. x .gt. this%xmax .or. this%n .le. 1)then
          f = 0.d0
          return
       endif
       l = min(max(coop_left_index(this%n, this%f2, x), 1), this%n-1)
       f = (this%f(l+1) - this%f(l)) / (this%f2(l+1) - this%f2(l))
       return
    case(COOP_INTERPOLATE_NONUNIFORM)
       call coop_splint_derv(this%n, this%f1, this%f, this%f2, x, f)
       return
    end select
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
    COOP_REAL x, f, a, b, xdiff, ap, bp
    COOP_INT l, r
    select case(this%method)
    case(COOP_INTERPOLATE_POWERLAW)
       f = sum(this%f2*(this%f1-1.d0)*x**(this%f1-2.d0), this%f2 .ne. 0.d0 .and. this%f1.ne. 1.d0)
       return
    case(COOP_INTERPOLATE_RATIONAL)
       a = sum(this%f*x**this%f1, this%f .ne. 0.d0)
       b = sum(this%f2*x**this%f3, this%f2 .ne. 0.d0)
       ap = sum(this%f*this%f1*x**(this%f1-1.d0), this%f .ne. 0.d0 .and. this%f1 .ne. 0.d0)
       bp = sum(this%f2*this%f3*x**(this%f3-1.d0), this%f2 .ne. 0.d0 .and. this%f3 .ne. 0.d0)
       f =  (sum(this%f*this%f1*(this%f1-1.d0)*x**(this%f1-2.d0), this%f .ne. 0.d0 .and. this%f1 .ne. 0.d0 .and. this%f1 .ne. 1.d0) - 2.d0*ap*bp/b + a/b*(2.d0*bp**2/b-sum(this%f2*this%f3*(this%f3-1.d0)*x**(this%f3-2.d0), this%f2.ne.0.d0 .and. this%f3 .ne. 0.d0 .and. this%f3.ne. 1.d0)))/b
       return
    case(COOP_INTERPOLATE_POLYNOMIAL)
       f = coop_polyvalue(this%n-2, this%f2, x)
       return
    case(COOP_INTERPOLATE_ZIGZAG)
       f = 0.d0
       return
    case(COOP_INTERPOLATE_NONUNIFORM)
       call coop_splint_derv2(this%n, this%f1, this%f, this%f2, x, f)
       return
    end select
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
    case(COOP_INTERPOLATE_CHEBYSHEV, COOP_INTERPOLATE_POWERLAW, COOP_INTERPOLATE_POLYNOMIAL, COOP_INTERPOLATE_RATIONAL, COOP_INTERPOLATE_NONUNIFORM)
       xmin = this%xmin
       xmax = this%xmax
       n = 512
    case default
       if(this%scale .ge. 0.d0)then
          iloc = coop_maxloc(this%f)
       else
          iloc = coop_minloc(this%f)
       endif
       select case(this%method)
       case(COOP_INTERPOLATE_ZIGZAG)
          xm = this%f1(iloc)
          goto 100
       case(COOP_INTERPOLATE_LINEAR)
          xm = this%xmin + this%dx * (iloc-1)
          goto 100
       end select
       xmin = this%xmin + this%dx*max(iloc-2, 0)
       xmax = this%xmin + this%dx*min(iloc, this%n-1)
       n = 16
    end select
    if(this%scale .ge. 0.d0)then
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
    else
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
    endif
100 if(this%xlog) xm = exp(xm)
  end function coop_function_maxloc


  function coop_function_minloc(this) result(xm)
    class(coop_function)::this
    COOP_REAL xm
    integer iloc, n, iloop
    COOP_REAL xmin, xmax, x, fm, dx, f
    select case(this%method)
    case(COOP_INTERPOLATE_CHEBYSHEV, COOP_INTERPOLATE_POWERLAW, COOP_INTERPOLATE_POLYNOMIAL, COOP_INTERPOLATE_RATIONAL, COOP_INTERPOLATE_NONUNIFORM)
       xmin = this%xmin
       xmax = this%xmax
       n = 512
    case default
       if(this%scale .ge. 0.d0)then
          iloc = coop_minloc(this%f)
       else
          iloc = coop_maxloc(this%f)
       endif
       select case(this%method)
       case(COOP_INTERPOLATE_ZIGZAG)
          xm = this%f1(iloc)
          goto 100
       case(COOP_INTERPOLATE_LINEAR)
          xm = this%xmin + this%dx * (iloc-1)
          goto 100
       end select
       xmin = this%xmin + this%dx*max(iloc-2, 0)
       xmax = this%xmin + this%dx*min(iloc, this%n-1)
       n = 16
    end select
    if(this%scale .ge. 0.d0)then
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
    else
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
    endif
100 if(this%xlog) xm = exp(xm)
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
       if(this%ylog)then
          fp = this%derivative_bare(log(x)) * exp(this%eval_bare(log(x))-log(x))*this%scale
       else
          fp  = this%derivative_bare(log(x))/x * this%scale
       endif
    else
       if(this%ylog)then
          fp = this%derivative_bare(x) * exp(this%eval_bare(x)) * this%scale
       else
          fp = this%derivative_bare(x) * this%scale
       endif
    endif
  end function coop_function_derivative


  function coop_function_derivative2(this, x) result(fpp)
    class(coop_function)::this
    COOP_REAL x, fpp, xbare, fp
    if(this%xlog)then
       xbare = log(x)
       fpp = this%derivative2_bare(xbare)
       fp = this%derivative_bare(xbare)
       if(this%ylog)then
          fpp = (fpp  - fp)
          fpp = (fpp + fp**2) * exp(this%eval_bare(xbare) - 2.d0*xbare)
       else
          fpp = (fpp  - fp)/x**2
       endif
    else
       fpp = this%derivative2_bare(x)
       if(this%ylog)then
          fpp = (fpp + this%derivative_bare(x)**2) * exp(this%eval_bare(x))
       endif
    endif
    fpp = fpp*this%scale
  end function coop_function_derivative2

  subroutine coop_function_monotonic_solution(this, f, x)
    class(coop_function)::this
    COOP_REAL,intent(in)::f
    COOP_REAL, intent(out)::x
    COOP_REAL xmin, xmax, fmin, fmax, xmid, fs, dx
    COOP_INT imin, imax, imid
    fs = (f - this%shift)/this%scale
    if(this%ylog)then  
       if(fs .lt. 0.d0) return
       fs = log(fs)
    endif
    select case(this%method)
    case(COOP_INTERPOLATE_CHEBYSHEV, COOP_INTERPOLATE_POLYNOMIAL, COOP_INTERPOLATE_ZIGZAG, COOP_INTERPOLATE_NONUNIFORM, COOP_INTERPOLATE_RATIONAL, COOP_INTERPOLATE_POWERLAW)
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

  function coop_function_polynomial( p ) result(f)
    COOP_REAL,dimension(:)::p
    type(coop_function)::f
    call f%init_polynomial(p)
  end function coop_function_polynomial


end module coop_function_mod


