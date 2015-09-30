module coop_random_mod
  use coop_wrapper_typedef
  use coop_MPI_mod
  implicit none
#include "constants.h"

  type coop_random_cycl
     COOP_INT::n = 0
     COOP_INT::loc = 0     
     COOP_INT, dimension(:), allocatable::ind
   contains
     procedure::init => coop_random_cycl_init
     procedure::free => coop_random_cycl_free
     procedure::next => coop_random_cycl_next
  end type coop_random_cycl

contains

  subroutine coop_random_cycl_free(this)
    class(coop_random_cycl)::this
    if(allocated(this%ind))deallocate(this%ind)
    this%n = 0
    this%loc = 1
  end subroutine coop_random_cycl_free

  subroutine coop_random_cycl_init(this, n)
    class(coop_random_cycl)::this
    COOP_INT, optional::n
    COOP_INT::i, j, tmp
    if(present(n))then
       if(this%n .ne. n)then
          if(allocated(this%ind))deallocate(this%ind)
          this%n = n
          allocate(this%ind(n))
       endif
    endif
    this%ind = (/ (i, i=1, this%n) /)
    do i=1, this%n-1
       tmp = this%ind(i)
       j = coop_random_index(this%n - i) + i
       this%ind(i) = this%ind(j)
       this%ind(j) = tmp
    enddo
    this%loc = 1    
  end subroutine coop_random_cycl_init

  function coop_random_cycl_next(this) result(next)
    class(coop_random_cycl)::this
    COOP_INT next
    next = this%ind(this%loc)
    this%loc = this%loc + 1
    if(this%loc .gt. this%n)call this%init()  !!initialize new random indices
  end function coop_random_cycl_next

  subroutine Coop_random_Init(i)
    COOP_INT,optional::i
    COOP_INT n, k
    COOP_INT,dimension(:),allocatable::seed
    call random_seed(size=n)
    allocate(seed(n))
    if(present(i))then
       seed(1) = i
    else
       seed(1)=coop_MPI_Rank()+ mod(floor(coop_systime_sec()), 99999)
    endif
    do k=2, n
       seed(k) = mod(seed(k-1)**2 + seed(k-1)*98 + 1, 999983)
    enddo
    call random_seed(put=seed)
    deallocate(seed)
  end subroutine Coop_random_Init


  function coop_random_vector(n) result(vec)
    COOP_INT n, i
    COOP_REAL vec(n)
    call coop_random_get_Gaussian_vector(n, vec)
    vec = vec/sqrt(sum(vec**2))
  end function coop_random_vector

  subroutine Coop_random_rotation(r, n)
    COOP_INT, intent(in) :: n
    COOP_REAL r(n, n), vec(N), norm
    COOP_INT i,j
    do j = 1,  n 
       do
          call coop_random_get_gaussian_vector(n, vec)
          do i = 1, j-1
             vec = vec - sum(vec*R(:,i))*R(:,i)
          end do
          norm = sum(vec**2)
          if (norm .gt. 1.d-3) exit
       end do
       R(:, j) = vec / sqrt(norm)
    end do

  end subroutine Coop_random_rotation

  subroutine coop_random_get_gaussian_vector(n, s)
    COOP_INT n
    COOP_REAL s(n), v(2), r
    COOP_INT i
    do i = 1, n-1, 2
       call random_number(v)
       v = 2.d0*v - 1.d0
       r=v(1)**2 + v(2)**2
       do while (R >= 1.d0)
          call random_number(v)
          v = 2.d0*v - 1.d0
          r=v(1)**2 + v(2)**2
       end do
       s(i:i+1) =  v * dsqrt(-2.d0*dlog(r)/r)
    enddo
    if(mod(n,2).ne.0)then
       call random_number(v)
       v = 2.d0*v - 1.d0
       r=v(1)**2 + v(2)**2
       do while (R >= 1.d0)
          call random_number(v)
          v = 2.d0*v - 1.d0
          r=v(1)**2 + v(2)**2
       end do
       s(n) = v(1) * dsqrt(-2.d0*dlog(r)/r)
    endif
  end subroutine coop_random_get_gaussian_vector

  function coop_random_gaussian_vector(n) result(s)
    COOP_INT n
    COOP_REAL s(n)
    call coop_random_get_gaussian_vector(n, s)
  end function coop_random_gaussian_vector

  !!thread-safe
  function coop_random_complex_Gaussian(real_imaginary) result(c)
    COOP_COMPLEX c    
    COOP_REAL r, v(2)
    logical,optional::real_imaginary
    call random_number(v)
    v = 2.d0*v - 1.d0
    r=v(1)**2 + v(2)**2
    do while (R >= 1.d0)
       call random_number(v)
       v = 2.d0*v - 1.d0
       r=v(1)**2 + v(2)**2
    end do
    if(present(real_imaginary))then
       if(real_imaginary)then
          c = cmplx(v(1), 0.d0)*dsqrt(-2.d0*dlog(r)/r)
       else
          c = cmplx(0.d0, v(2))*dsqrt(-2.d0*dlog(r)/r)
       endif
    else
       c = cmplx(v(1), v(2))*dsqrt(-dlog(r)/r)
    endif
  end function coop_random_complex_Gaussian

  function coop_random_complex_Gaussian_vector(n) result(c)
    COOP_INT::n
    COOP_COMPLEX c(n)
    COOP_INT::i
    do i=1, n
       c(i) = coop_random_complex_Gaussian()
    enddo
  end function coop_random_complex_Gaussian_vector
  

  function coop_random_Gaussian() result(s)
    COOP_REAL r, v(2), s
    call random_number(v)
    v = 2.d0*v - 1.d0
    r=v(1)**2 + v(2)**2
    do while (R >= 1.d0)
       call random_number(v)
       v = 2.d0*v - 1.d0
       r=v(1)**2 + v(2)**2
    end do
    s = v(1) * dsqrt(-2.d0*dlog(r)/r)
  end function coop_random_Gaussian


  subroutine coop_random_get_Gaussian_pair(g1, g2)
    COOP_REAL r, v(2), g1, g2
    call random_number(v)
    v = 2.d0*v - 1.d0
    r=v(1)**2 + v(2)**2
    do while (R >= 1.d0)
       call random_number(v)
       v = 2.d0*v - 1.d0
       r=v(1)**2 + v(2)**2
    end do
    g2 = dsqrt(-2.d0*dlog(r)/r)
    g1 = g2*v(1)
    g2 = g2*v(2)
  end subroutine coop_random_get_Gaussian_pair



  function Coop_random_exp()
    !
    !     Coop_random-number generator for the exponential distribution
    !     Algorithm EA from J. H. Ahrens and U. Dieter,
    !     Communications of the ACM, 31 (1988) 1330--1337.
    !     Coded by K. G. Hamilton, December 1996, with corrections.
    !
    COOP_REAL coop_random_exp
    COOP_REAL u, up, g, y
    COOP_REAL,parameter ::   alog2= 0.6931471805599453
    COOP_REAL,parameter ::      a = 5.7133631526454228
    COOP_REAL,parameter ::      b = 3.4142135623730950
    COOP_REAL,parameter ::     c = -1.6734053240284925
    COOP_REAL,parameter ::      p = 0.9802581434685472
    COOP_REAL,parameter ::     aa = 5.6005707569738080
    COOP_REAL,parameter ::     bb = 3.3468106480569850
    COOP_REAL,parameter ::     hh = 0.0026106723602095
    COOP_REAL,parameter ::     dd = 0.0857864376269050

    u = coop_random_unit()
    do while (u.le.0)                 ! Comment out this block 
       u = coop_random_unit()                    ! if your RNG can never
    enddo                             ! return exact zero
    g = c
    u = u+u
    do while (u.lt.1.0)
       g = g + alog2
       u = u+u
    enddo
    u = u-1.0
    if (u.le.p) then
       Coop_random_exp = g + aa/(bb-u)
       return
    endif
    do
       u = coop_random_unit()
       y = a/(b-u)
       up = coop_random_unit()
       if ((up*hh+dd)*(b-u)**2 .le. exp(-(y+c))) then
          Coop_random_exp = g+y
          return
       endif
    enddo

  end function Coop_random_exp

  function coop_random_unit() result(r)
    COOP_REAL r
    call random_number(r)
  end function coop_random_unit

  subroutine coop_random_pick_index(n, a, b)
    COOP_INT n, a, b
    if(b.gt.a)then
       n = a + floor((b-a+0.999999999d0)*coop_random_unit())
    else
       n = b + floor((a-b+0.999999999d0)*coop_random_unit())
    endif
  end subroutine coop_random_pick_index

  function coop_random_index(n)
    COOP_INT n
    COOP_INT coop_random_index
    coop_random_index = min(floor(coop_random_unit()*n) + 1, n)
  end function coop_random_index

  function coop_random_index_paral(n) result(ind)
    COOP_INT n, ind
    COOP_REAL x
    call random_number(x)
    ind = min(floor(x*n) + 1, n)
  end function coop_random_index_paral

  function coop_random_uniform(a, b)
    COOP_REAL a, b, coop_random_uniform
    coop_random_uniform = a  + coop_random_unit()*(b-a)
  end function coop_random_uniform

  function coop_rand01() result(i01)
    COOP_INT i01
    COOP_REAL r
    call random_number(r)
    if(r .gt. 0.5d0)then
       i01 = 1
    else
       i01 = 0
    end if
  end function coop_rand01



  subroutine coop_fit_gaussian(x, nbins, xbar, sigma, A)
    COOP_INT nbins    
    COOP_REAL::x(:)
    COOP_REAL::xcopy(size(x))
    COOP_REAL c(nbins), xb(nbins), err(nbins), dx, newxbar, newsigma, newA, newchisq, chisq, xbar, sigma, A
    COOP_INT::i, n, nstep, ntail, j, next, nstep1, nstep2, nstep3, fail
    n = size(x)
    if(n.lt. 48 .or. nbins.lt. 6) stop "fit_gaussian assumes at least 6 bins and 48 samples"
    xcopy = x
    call coop_quicksort(xcopy)
    nstep1 = n/nbins/4
    if(nstep1 .lt. 2)stop "too many bins for fit_gaussian: max number of bins is n/8"        
    nstep2 = nstep1*2
    nstep3 = nstep1*3
    if(nbins.gt.6)then
       nstep = (n- (nstep1+nstep2+nstep3)*2)/(nbins-6)
    endif
    ntail = n - nstep*(nbins-6) - (nstep1+nstep2+nstep3)*2
    nstep1 = nstep1 + ntail/6
    nstep2 = nstep2 + ntail/6
    nstep3 = (n - (nstep1+nstep2)*2 - nstep*(nbins-6) ) / 2
    j = nstep1
    xb(1) = sum(xcopy(1:j))/nstep1
    dx = ((xcopy(j+1)+xcopy(j))/2.d0 - xcopy(1) - (xcopy(2)-xcopy(1))/2.d0)
    c(1) = nstep1/dx
    err(1) =1.d0/dx
    do i=2, nbins-1
       next = j+1
       if(i.eq.2 .or. i.eq. nbins-1)then
          j = j + nstep2
       elseif(i.eq.3.or.i.eq.nbins-2)then
          j = j + nstep3
       else
          j = j + nstep
       endif
       xb(i) = sum(xcopy(next:j))/(j-next+1)
       dx = ((xcopy(j+1)+xcopy(j))/2.d0-(xcopy(next)+xcopy(next-1))/2.d0)
       c(i) = (j-next+1)/dx
       err(i) = 1.d0/dx
    enddo
    next = j+1
    xb(nbins) = sum(xcopy(next:n))/(n-j)
    dx = ((xcopy(n)+(xcopy(n)-xcopy(n-1))/2.d0)-(xcopy(j)+xcopy(next))/2.d0)
    c(nbins) = (n-j)/dx
    err(nbins) = 1.d0/dx
    xbar = sum(x)/n
    sigma = sqrt(sum((x-xbar)**2)/n)
    call get_chisq(xbar, sigma, chisq)
    dx = sigma/20.d0
    do i = 1, 30
       call random_walk(xbar, sigma, dx, newxbar, newsigma)
       call get_chisq(newxbar, newsigma, newchisq)
       if(newchisq .lt. chisq)then
          xbar = newxbar
          sigma = newsigma
          chisq = newchisq
       endif
    enddo
    dx = sigma/20.d0
    fail = 0
    do i = 1, 50000
       call random_walk(xbar, sigma, dx, newxbar, newsigma)
       call get_chisq(newxbar, newsigma, newchisq)
       if(newchisq .lt. chisq)then
          xbar = newxbar
          sigma = newsigma
          chisq = newchisq
          fail = 0
       else
          fail = fail + 1
       endif
       if(fail .gt. 10)then
          dx = dx * 0.8d0
          fail  = 0
          if(dx/sigma .lt. 1.d-5)exit
       endif
    enddo
    call get_A(xbar, sigma, A)
    A = A*sigma*sqrt(coop_2pi)
  contains

    subroutine random_walk(xbar, sigma, dx, newxbar, newsigma)
      COOP_REAL::r, xbar, sigma, dx, newxbar, newsigma
      call random_number(r)
      newxbar = xbar + dx*(r-0.5d0)
      call random_number(r)      
      newsigma = exp(log(sigma)+dx/sigma*(r-0.5d0))
    end subroutine random_walk

    subroutine get_A(xbar, sigma, A)
      COOP_REAL::xbar, sigma, A
      A = sum( c*exp(-((xb-xbar)/sigma)**2/2.d0)/err) / sum(exp(-((xb-xbar)/sigma)**2)/err)
    end subroutine get_A
    
    subroutine get_chisq(xbar, sigma, chisq)
      COOP_REAL::xbar, sigma, A, chisq
      call get_A(xbar, sigma, A)
      chisq  = sum((A*exp(-((xb - xbar)/sigma)**2/2.d0) - c)**2/err)
    end subroutine get_chisq
    
  end subroutine coop_fit_gaussian  
  
end module Coop_random_mod






