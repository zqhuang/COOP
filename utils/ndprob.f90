module coop_nd_prob_mod
  use coop_list_mod
  use coop_matrix_mod
  use coop_file_mod
  use coop_wrapper_typedef
  implicit none
#include "constants.h"  


  
  type coop_nd_prob
     COOP_STRING::name
     COOP_REAL::pow = 2.d0  !!weight = 1/r^pow
     COOP_INT:: range = 128 !!search # samples in one dimension
     COOP_INT::dim = 0
     COOP_INT::n = 0
     COOP_INT::capacity
     COOP_REAL::total_mult = 0.d0
     logical::need_update = .true.     
     COOP_REAL,dimension(:,:),allocatable::x, u, x2u
     COOP_REAL,dimension(:),allocatable::f, r2, mult, mean, umin, umax
     !!f saves - log(prob), r2 saves (x C^{-1} x )
     COOP_INT,dimension(:,:),allocatable::ind
   contains
     procedure::report_err => coop_nd_prob_report_err
     procedure::load => coop_nd_prob_load
     procedure::load_chains => coop_nd_prob_load_chains     
     procedure::free => coop_nd_prob_free
     procedure::alloc => coop_nd_prob_alloc
     procedure::eval => coop_nd_prob_eval
     procedure::add => coop_nd_prob_add
     procedure::update => coop_nd_prob_update
  end type coop_nd_prob


contains


  subroutine coop_nd_prob_update(this)
    class(coop_nd_prob)::this
    COOP_INT i1, i2, dim, i
    COOP_REAL, dimension(:),allocatable::f
    COOP_REAL:: worst_lnlike
    this%total_mult = sum(this%mult(1:this%n))
    if(this%total_mult .le. 0.d0) stop "nd_prob: zero total multiplicity"
    do dim = 1, this%dim
       this%mean(dim) = sum(this%x(dim, 1:this%n)*this%mult)/this%total_mult
    enddo
    do i1 = 1, this%dim
       do i2 = 1, i1
          this%x2u(i1, i2) = sum((this%x(i1, 1:this%n) - this%mean(i1))*(this%x(i2, 1:this%n) - this%mean(i2))*this%mult(1:this%n))/this%total_mult
          this%x2u(i2, i1) = this%x2u(i1, i2)
       enddo
    enddo
    call coop_matsym_cholesky(this%dim, this%x2u)
    call coop_matsym_cholesky_invert(this%dim, this%x2u)
    do i = 1, this%n
       this%u(:, i) = matmul(this%x2u, this%x(:,i)-this%mean)
    enddo
    !!now sort and make index list
    allocate(f(this%n))
    do dim = 1, this%dim
       f= this%u(dim, 1:this%n)
       call coop_quicksort_index(f, this%ind(1:this%n, dim))
       this%umin(dim) = this%u(dim, this%ind(1, dim))
       this%umax(dim) = this%u(dim, this%ind(this%n, dim))
    enddo
    deallocate(f)
    !$omp parallel do
    do i=1, this%n
       this%r2(i) = sum(this%u(:, i)**2)
       if(this%f(i) .lt. coop_logZero) &
            this%f(i) = this%f(i) - this%r2(i)/2.d0
    enddo
    !$omp end parallel do
    worst_lnlike = maxval(this%f(1:this%n), mask = (this%f(1:this%n) .lt. coop_logZero))
    where (this%f(1:this%n) .ge. coop_logZero)
       this%f(1:this%n) = worst_lnlike + 10.d0
       !!no event detection: assuming probability e^{-10} times smaller than gaussian interpolation.
    end where
  end subroutine coop_nd_prob_update

  subroutine coop_nd_prob_alloc(this, dim, capacity)
    class(coop_nd_prob)::this
    COOP_INT::capacity, dim
    if(capacity .le. 0 .or. dim .le.0) stop "cannot allocate nd_prob with 0 data or 0 dimension"
    call this%free()    
    this%dim = dim
    this%n = 0
    this%total_mult = 0.d0
    this%capacity = capacity
    this%need_update = .true.
    allocate(this%x(dim, capacity), this%u(dim, capacity), this%x2u(dim, dim), this%f(capacity), this%mult(capacity), this%mean(dim), this%ind(capacity, dim), this%umin(dim), this%umax(dim), this%r2(capacity))    
  end subroutine coop_nd_prob_alloc

  subroutine coop_nd_prob_add(this, line)
    class(coop_nd_prob)::this    
    COOP_SINGLE::line(:)
    if(this%n.eq. this%capacity) return
    this%n = this%n + 1
    this%mult(this%n) = line(1)
    this%f(this%n) = line(2)    
    this%x(:, this%n) = line(3:this%dim+2)
    this%need_update = .true.
  end subroutine coop_nd_prob_add
  
  subroutine coop_nd_prob_free(this)
    class(coop_nd_prob)::this
    if(allocated(this%x))deallocate(this%x)
    if(allocated(this%u))deallocate(this%u)    
    if(allocated(this%x2u))deallocate(this%x2u)    
    if(allocated(this%f))deallocate(this%f)
    if(allocated(this%r2))deallocate(this%r2)    
    if(allocated(this%mult))deallocate(this%mult)    
    if(allocated(this%mean))deallocate(this%mean)            
    if(allocated(this%ind))deallocate(this%ind)
    if(allocated(this%umin))deallocate(this%umin)
    if(allocated(this%umax))deallocate(this%umax)
    this%n = 0
    this%dim = 0
    this%total_mult = 0.d0
    this%need_update = .true.    
  end subroutine coop_nd_prob_free

  subroutine coop_nd_prob_load_chains(this, prefix, nvars)
    COOP_UNKNOWN_STRING::prefix
    class(coop_nd_prob)::this
    type(coop_file)::fp
    COOP_STRING::filename
    COOP_INT::nvars, n, nchains, i
    type(coop_list_realarr)::rl
    if(.not. coop_file_exists(trim(prefix)//"_1.txt"))then
       write(*,*) "cannot find chain "//trim(prefix)
       stop
    endif
    n = 0
    i = 1
    filename = trim(prefix)//"_"//COOP_STR_OF(i)//".ndf"
    do while(coop_file_exists(filename))
       call coop_file_load_realarr(filename, rl, 2+nvars)
       i = i + 1
       filename = trim(prefix)//"_"//COOP_STR_OF(i)//".ndf"       
    enddo
    call this%alloc(dim = nvars, capacity = rl%n)
    do i=1, rl%n
       call this%add(rl%element(i))
    enddo
    call rl%free()
    call this%update()
  end subroutine coop_nd_prob_load_chains

  subroutine coop_nd_prob_load(this, filename, form, nvars, name)
    COOP_UNKNOWN_STRING::filename
    COOP_UNKNOWN_STRING::form
    class(coop_nd_prob)::this
    type(coop_file)::fp
    COOP_INT, optional::nvars
    COOP_UNKNOWN_STRING,optional::name
    COOP_REAL,dimension(:),allocatable::f
    COOP_INT i, dim, i1, i2, n
    COOP_REAL::worst_lnlike, best_lnlike
    if(.not. coop_file_exists(filename))then
       write(*,*) "nd_prob_load: "//trim(filename)//" does not exit"
       stop
    endif
    call this%free()
    if(present(name))then
       this%name = trim(name)
    else
       this%name = "NoName"
    endif
    select case(trim(form))
    case("MCMC", "mcmc", "CHAIN", "chain", "Chain")
       dim = coop_file_numcolumns(filename) - 2
    case("table", "TABLE", "Table", "NormalizedTable")
       dim = coop_file_numcolumns(filename) - 1
    case default
       write(*,*) "nd_prob_load: "//trim(form)
       stop "unknown format"
    end select
    if(present(nvars))then
       if(nvars .gt. dim) stop "gd_function_load: dimension overflow"       
       dim = nvars
    endif
    n = coop_file_numlines(filename)
    if(dim .le. 0 .or. n .le. 0 .or. dim .gt. 64)then
       write(*,*) "nd_prob_load: dim = "//COOP_STR_OF(dim)
       write(*,*) "nd_prob_load: n = "//COOP_STR_OF(n)
       write(*,*) "cannot initialize (1<=dim<=64; n>=1 not satisfied)"
       stop
    endif
    call this%alloc(dim = dim, capacity = n)
    this%n = n  
    if(this%dim .eq. 1)then
       this%range = 2
    else
       this%range = ceiling(dble(this%n)**(1.d0/this%dim))*2  !!search roughly two columns
    endif
    
    call fp%open(filename, "r")
    select case(trim(form))
    case("mcmc", "MCMC", "CHAIN", "chain")
       do i=1, this%n
          read(fp%unit, *, ERR=100, END=100) this%mult(i), this%f(i), this%x(:, i)
       enddo
    case("NormalizedTable", "Table", "TABLE", "table")
       do i=1, this%n
          read(fp%unit, *, ERR=100, END=100) this%x(:, i), this%f(i)
          if(this%f(i) .lt. 0.d0)then
             stop "nd_prob_load: with table format the probability must be semi-positive definite"
          endif
          if(this%f(i) .gt. 0.d0)then
             this%f(i) = -log(this%f(i))
          else
             this%f(i) = coop_logZero
          endif
       enddo
       this%mult(1:this%n) = 1.d0
       if(trim(form).ne."NormalizedTable")then
          best_lnlike = minval(this%f(1:this%n), mask = this%f(1:this%n) .lt. coop_logZero) 
          where (this%f(1:this%n) .lt. coop_logZero)
             this%f(1:this%n) = this%f(1:this%n) - best_lnlike
          end where
       endif
    case default
       stop "nd_prob_load: unknown format"
    end select
    call fp%close()    
    call this%update()
    return
100 write(*,*) "nd_prob_load: "//trim(filename)//" is broken at line"//COOP_STR_OF(i)
    stop
    
    
  end subroutine coop_nd_prob_load

  function coop_nd_prob_eval(this, x) result(y)
    class(coop_nd_prob)::this
    COOP_REAL,dimension(:),intent(IN)::x
    COOP_REAL::y
    COOP_INT,dimension(:),allocatable::imin
    COOP_INT::lower(size(x)), upper(size(x))
    COOP_REAL::u(size(x))
    COOP_INT::dim, ib, it, im, i, nmin
    COOP_REAL::r2,  far(size(x))
    COOP_REAL,dimension(:),allocatable::dis2
    if(size(x).ne.this%dim)then
       call this%report_err(x)
    endif
    u = matmul(this%x2u, x - this%mean)
    r2 = sum(u**2)
    select case(this%dim)
    case(1)
       nmin = 2
    case(2)
       nmin = 4
    case default
       nmin = 8
    end select
    allocate(imin(nmin))
    if(this%n .lt. 1024)then  !!do brute force search, O(dim * n) operations
       allocate(dis2(this%n))
       !$omp parallel do
       do i=1, this%n
          dis2(i) = max(sum( (u - this%u(:, i))**2 ), 1.d-30)
       enddo
       !$omp end parallel do
       call coop_find_minlocs(dis2, imin)
       y = sum(this%f(imin)/dis2(imin))/sum(1.d0/dis2(imin)) + r2/2.d0
       goto 100
    endif
    
    !!for large n, do sub-sample search, O(dim log_2 n) operations
    do dim = 1, this%dim
       if(u(dim) .le. this%u(dim, this%ind(1, dim)))then
          lower(dim) = 1
          upper(dim) = 1+this%range
       elseif(u(dim).ge. this%u(dim, this%ind(this%n, dim)))then
          lower(dim) = this%n - this%range
          upper(dim) = this%n
       else
          ib = 1
          it = this%n
          do while(it - ib .gt. 1)
             im = (ib + it)/2
             if(this%u(dim, this%ind(im, dim)).gt. u(dim))then
                it = im
             else
                ib = im
             endif
          enddo
          upper(dim) = min(it + this%range, this%n)
          lower(dim) = max(1, ib - this%range)
       endif
       far(dim) = min( abs(this%u(dim, this%ind(lower(dim), dim))-u(dim)), &
            abs(this%u(dim, this%ind(upper(dim), dim))-u(dim)) )
    enddo
    dim = coop_maxloc(far)  !!choose the dimension with sparse samples
    allocate(dis2(lower(dim):upper(dim)))
    !$omp parallel do
    do i = lower(dim), upper(dim)
       dis2(i) = max(sum( (u - this%u(:, this%ind(i, dim)))**2), 1.d-30)
    enddo
    !$omp end parallel do
    call coop_find_minlocs(dis2, imin)
    imin = imin + lower(dim)-1
    y = sum(this%f(this%ind(imin, dim))/dis2(imin))/sum(1.d0/dis2(imin)) + r2/2.d0
100 deallocate(dis2, imin)
    return
    
  end function coop_nd_prob_eval

  subroutine coop_nd_prob_report_err(this, x)
    class(coop_nd_prob)::this
    COOP_REAL,dimension(:)::x
    write(*,*) "coop_nd_prob: "//trim(this%name)
    write(*,*) "dimension: ", this%dim
    write(*,*) "saved data points: ", this%n
    write(*,*) "failed to evaluate at x = ", x
    stop
  end subroutine coop_nd_prob_report_err

  
end module coop_nd_prob_mod
    

