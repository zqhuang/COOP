module coop_cholesky_mod
  use coop_wrapper_typedef
  use coop_sort_mod
  use coop_sortrev_mod
  use coop_MPI_mod
  implicit none

#include "constants.h"

  private


  public::coop_solve_constrained, coop_cholesky, coop_cholesky_inv, coop_cholesky_sq, coop_cholesky_solve, coop_cholesky_solve_transpose, coop_cholesky_solve_mult, coop_sympos_inverse

contains

  !! for constrained problems
  !!suppose you have a covariance matrix C for known vector x_I and unknown vector x_J, you want to construct a mapping matrix Fmean (n_J, n_I) such that mean <x_J> = Fmean x_I and a matrix Ffluc (n_J, n_J) such that \delta x_J = Ffluc (Gaussian random vector y_J)
  !!after calling this subroutine C is destroyed
  subroutine coop_solve_constrained(m, n_known, n_unknown, dim_fmean, dim_ffluc, C, Fmean, Ffluc)
    COOP_INT,intent(IN)::m, n_known, n_unknown, dim_fmean, dim_ffluc
    COOP_INT::n, info
    COOP_REAL,intent(INOUT)::C(m, n_known+n_unknown)
    COOP_REAL,intent(OUT)::Fmean(dim_fmean, n_known), FFluc(dim_ffluc, n_unknown)
    COOP_REAL::cI(n_known, n_known)
    COOP_INT::i, j
    !!do cholesky for  n_known x n_known
    n = n_unknown + n_known
    do i=1, n_known
       do j=1, i
          cI(i, j) = c(i, j)
          cI(j, i) = c(i, j)
       enddo
    enddo
    call coop_cholesky(m, n_known, C, info)
    if(info .ne. 0) stop "solve_constrained: matrix not positive definite"
    Fmean(1:n_unknown, 1:n_known) = C(n_known+1:n, 1:n_known)
    call coop_cholesky_solve_transpose(m = m, n = n_known, mb = dim_fmean, neqs = n_unknown, a = C,  b = Fmean)
    Ffluc(1:n_unknown, 1:n_unknown) = C(n_known+1:n, n_known+1:n) - matmul(Fmean(1:n_unknown, 1:n_known), transpose(C(n_known+1:n, 1:n_known)))
    call coop_cholesky(dim_ffluc, n_unknown, FFluc, info)
    if(info .ne. 0) stop "solve_constrained: matrix not positive definite"    
  end subroutine coop_solve_constrained

  subroutine coop_fillzero(uplo, m, n, a)
    character uplo
    COOP_INT m, n, i, j
    COOP_REAL a(m, n)
    select case(uplo)
    case('U', 'u')
       do i=2, n
          a(1:i-1, i) = 0.d0
       enddo
    case('L', 'l')
       do i=1, n-1
          a(i+1:n, i) = 0.d0
       enddo
    case default
       stop 'fillzero: unknown uplo option'
    end select
  end subroutine coop_fillzero
  
  subroutine coop_cholesky(m, n, a, info)
    COOP_INT m, n, i, j, info
    COOP_REAL a(m, n), mineig
    info = 0
#ifdef HAS_LAPACK
    call dpotrf('L', n, a, m, info)
#else    
    if(a(1,1).le.0.d0)then
       info = 1
       return
    endif
    a(1,1)= sqrt(a(1,1))
    a(2:n,1) = a(2:n,1)/a(1,1)
    do  i = 2, n
       a(i, i)= a(i,i) - sum(a(i,1:i-1)**2)
       if(a(i,i).gt. 0.d0)then
          a(i,i) = sqrt(a(i,i))
          do j = i+1, n 
             a(j, i)=(a(j, i)-sum(a(j,1:i-1)*a(i,1:i-1)))/a(i,i)
          enddo
       else
          info = i
          return
       endif
       a(1:i-1, i) = 0.d0       
    enddo
#endif    
  end Subroutine coop_cholesky



    subroutine coop_cholesky_solve(m, n, a, b)
    COOP_INT m, n, mb, neqs, i, ieq
    COOP_REAL  a(m, n), b(n)  !m>=n; mb>=n
    b(1)=b(1)/a(1,1)
    do i = 2, n
       b(i)=(b(i)-sum(a(i,1:i-1)*b(1:i-1)))/a(i,i)
    enddo
    b(n)=b(n)/a(n, n)
    do i=n-1,1,-1
       b(i)=(b(i)-sum(b(i+1:n)*a(i+1:n,i)))/a(i,i)
    enddo
  end subroutine coop_cholesky_solve

  subroutine coop_cholesky_solve_transpose(m, n, mb, neqs, a, b)
    COOP_INT m, n, mb, neqs, i, ieq
    COOP_REAL  a(m, n), b(mb, n)  !!m >= n; mb >= neqs
    b(1:neqs,1)=b(1:neqs,1)/a(1,1)
    do i = 2, n
       do ieq = 1, neqs
          b(ieq, i)=(b(ieq, i)-sum(a(i,1:i-1)*b(ieq, 1:i-1)))/a(i,i)
       enddo
    enddo
    b(1:neqs, n)=b(1:neqs, n)/a(n, n)
    do i=n-1,1,-1
       do ieq = 1, neqs
          b(ieq, i)=(b(ieq, i)-sum(b(ieq, i+1:n)*a(i+1:n,i)))/a(i,i)
       enddo
    enddo
  end subroutine coop_cholesky_solve_transpose


  subroutine coop_cholesky_solve_mult(m, n, mb, neqs, a, b)
    COOP_INT m, n, mb, neqs, i, ieq
    COOP_REAL  a(m, n), b(mb, neqs)  !m>=n; mb>=n
    b(1, 1:neqs)=b(1, 1:neqs)/a(1,1)
    do i = 2, n
       do ieq = 1, neqs
          b(i, ieq)=(b(i, ieq)-sum(a(i,1:i-1)*b(1:i-1, ieq)))/a(i,i)
       enddo
    enddo
    b(n, 1:neqs)=b(n, 1:neqs)/a(n, n)
    do i=n-1,1,-1
       do ieq = 1, neqs
          b(i, ieq)=(b(i, ieq)-sum(b(i+1:n, ieq)*a(i+1:n,i)))/a(i,i)
       enddo
    enddo
  end subroutine coop_cholesky_solve_mult
  
  subroutine coop_cholesky_inv(m, n, a)
    COOP_INT,intent(IN)::m, n
    COOP_REAL::a(m, n)
    COOP_INT i,j
    do  i=1,n
       a(i,i)=1.d0/a(i, i)
       do  j=i+1,n
          a(j,i)=-dot_product(a(j,i:j-1),a(i:j-1,i))/a(j,j)
       enddo
    enddo
  end subroutine coop_cholesky_inv

  subroutine coop_cholesky_sq(m, n, a)
    COOP_INT,intent(IN)::m, n
    COOP_REAL::a(m, n)
    COOP_INT::i, j
    do i=1,n
       do j=i+1,n
          a(i,j)=dot_product(a(j:n,i),a(j:n,j))
       enddo
    enddo
    do i=1,n-1
       a(i,i)=dot_product(a(i:n,i),a(i:n,i))
       a(i+1:n,i)=a(i,i+1:n)
    enddo
    a(n,n)=a(n,n)**2
  end subroutine coop_cholesky_sq

  !!a(1:n, 1:n) is a positive definite symmetric matrix, replace a(1:n, 1:n) with a^{-1}(1:n, 1:n)
  !!Input: lower left triangle of a is required
  !! m is the first physical dimension of a (i.e., a is actually stored as a(m, n); m>=n
  subroutine coop_sympos_inverse(m, n, a)
    COOP_INT,intent(IN)::m, n
    COOP_REAL::a(m, n) !m>=n
    COOP_INT info
    select case(n)
    case(1)
       a=1.d0/a
    case(2)
       a=ReShape( (/ a(2,2), -a(2,1), -a(1,2), a(1,1) /) / (a(1,1)*a(2,2)-a(1,2)*a(2,1)), (/ 2, 2 /) )
    case(3)
       a = Reshape( (/ a(2,2)*a(3,3) - a(2,3)*a(3,2), a(2,3)*a(3,1)-a(2,1)*a(3,3), a(2,1)*a(3,2)- a(2,2)*a(3,1),  &
            a(3,2)*a(1,3)-a(1,2)*a(3,3), a(3,3)*a(1,1)-a(3,1)*a(1,3), a(1,2)*a(3,1)-a(1,1)*a(3,2), &
            a(1,2)*a(2,3)- a(2,2)*a(1,3), a(2,1)*a(1,3)-a(1,1)*a(2,3),  a(1,1)*a(2,2)-a(1,2)*a(2,1) /), (/ 3, 3 /) ) &
            / ( a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2)) + a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) + a(1,3)*(a(2,1)*a(3,2)- a(2,2)*a(3,1)) )
    case default
       call coop_cholesky(m, n, a, info)
       if(info.ne.0) call Coop_return_error("coop_sympos_inverse", "the matrix is not positive definite", "stop")
#ifdef HAS_LAPACK
       call dpotri("L", n, a, m, info)
       if(info.ne.0) call Coop_return_error("coop_sympos_inverse", "the matrix is not positive definite", "stop")       
#else       
       call coop_cholesky_inv(m, n, a)
       call coop_cholesky_sq(m, n, a)
#endif       
    end select
  end subroutine coop_sympos_inverse

  
end module coop_cholesky_mod
