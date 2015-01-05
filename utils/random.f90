module coop_random_mod
  use coop_wrapper_typedef
  use coop_MPI_mod
  implicit none
#include "constants.h"

contains


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
    do j = 1, N
       do
          call coop_random_get_gaussian_vector(n, vec)
          do i = 1, j-1
             vec = vec - sum(vec*R(:,i))*R(:,i)
          end do
          norm = sum(vec**2)
          if (norm > 1e-3) exit
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
    COOP_REAL r, v(2)
    COOP_COMPLEX c
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

end module Coop_random_mod






