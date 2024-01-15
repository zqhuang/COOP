module coop_matsym_mod
  use coop_matrix_mod
  use coop_wrapper_typedef
  use coop_sort_mod
  use coop_sortrev_mod
#include "constants.h"
  private

  type coop_symmetric_matrix
     COOP_INT::n = 0
     logical::has_UR = .false.
     logical::diag_normalized = .false.
     logical::posdef = .false.
     COOP_REAL,dimension(:,:),allocatable::c
   contains

  end type coop_symmetric_matrix
  
  
  
  
contains

  !! for constrained problems
  !!suppose you have a covariance matrix C for known vector x_I and unknown vector x_J, you want to construct a mapping matrix Fmean (n_J, n_I) such that mean <x_J> = Fmean x_I and a matrix Ffluc (n_J, n_J) such that \delta x_J = Ffluc (Gaussian random vector y_J)
  !!after calling this subroutine C is destroyed
  subroutine coop_solve_constrained(m, n_known, n_unknown, C, Fmean, Ffluc, epsilon)
    COOP_INT,intent(IN)::m, n_known, n_unknown
    COOP_INT::n
    COOP_REAL,intent(INOUT)::C(m, n_known+n_unknown)
    COOP_REAL,intent(OUT)::Fmean(n_unknown, n_known), FFluc(n_unknown, n_unknown)
    COOP_REAL::cI(n_known, n_known)
    COOP_INT::i, j
    COOP_REAL::epsilon
    !!do cholesky for  n_known x n_known
    n = n_unknown + n_known
    do i=1, n_known
       do j=1, i
          cI(i, j) = c(i, j)
          cI(j, i) = c(i, j)
       enddo
    enddo
    call coop_cholesky(m, n_known, C, epsilon)
    Fmean = C(n_known+1:n, 1:n_known)
    call coop_cholesky_solve(m, n_known, n_unknown, C,  Fmean)
    Ffluc = C(n_known+1:n, n_known+1:n) - matmul(Fmean, transpose(C(n_known+1:n, 1:n_known)))
    call coop_cholesky(n_unknown, n_unknown, FFluc, epsilon)
  end subroutine coop_solve_constrained
    
  
  subroutine coop_cholesky(m, n, a, epsilon)
    COOP_INT m, n, i, j
    COOP_REAL a(m, n), mineig
    COOP_REAL,optional::epsilon
    if(present(epsilon))then
       mineig = a(1,1)
       do i=2,n
          if(a(i,i).gt.mineig)mineig = a(i,i)
       enddo
       mineig = epsilon*mineig       
    else
       mineig = 0.d0
    endif
    if(a(1,1).gt. mineig)then
       a(1,1)= sqrt(a(1,1))
       a(2:n,1) = a(2:n,1)/a(1,1)
    else
       a(1:n,1) = 0.d0
    endif
    do  i = 2, n
       a(i, i)= a(i,i) - sum(a(i,1:i-1)**2)
       if(a(i,i).gt. mineig)then
          a(i,i) = sqrt(a(i,i))
          do j = i+1, n 
             a(j, i)=(a(j, i)-sum(a(j,1:i-1)*a(i,1:i-1)))/a(i,i)
          enddo
       else
          a(i:n,i) = 0.d0
       endif
    enddo
  end Subroutine coop_cholesky



  subroutine coop_cholesky_solve(m, n, neqs, a, b)
    COOP_INT m, n, mb, neqs, i, ieq
    COOP_REAL  a(m, n), b(neqs, n)
    COOP_REAL, dimension(:, :),allocatable::bcopy
    if(a(1,1).ne.0.d0)then
       b(1:neqs,1)=b(1:neqs,1)/a(1,1)
    else
       b(1:neqs,1) = 0.d0
    endif
    do i = 2, n
       if(a(i, i) .ne. 0.d0)then
          do ieq = 1, neqs
             b(ieq, i)=(b(ieq, i)-sum(a(i,1:i-1)*b(ieq, 1:i-1)))/a(i,i)
          enddo
       else
          b(1:neqs, i) = 0.d0
       endif
    enddo
    if(a(n,n).ne.0.d0)then
       b(1:neqs, n)=b(1:neqs, n)/a(n, n)
    else
       b(1:neqs, n) = 0.d0
    endif
    do i=n-1,1,-1
       if(a(i,i).ne.0.d0)then
          do ieq = 1, neqs
             b(ieq, i)=(b(ieq, i)-sum(b(ieq, i+1:n)*a(i+1:n,i)))/a(i,i)
          enddo
       endif
    enddo
  end subroutine coop_cholesky_solve
  


  subroutine coop_cholesky_inv(m, n, a)
    COOP_INT,intent(IN)::m, n
    COOP_REAL::a(m, n)
    COOP_INT i,j
    do  i=1,n
       if(a(i,i) .ne. 0.d0)then
          a(i,i)=1.d0/a(i, i)
          do  j=i+1,n
             if(a(j, j).ne.0.d0)then
                a(j,i)=-dot_product(a(j,i:j-1),a(i:j-1,i))/a(j,j)
             else
                a(j, i) = 0.d0
             endif
          enddo
       endif          
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
  !!epsilon (optional) is the toleratable upper bound of  diagonal element
  subroutine coop_sympos_inverse(m, n, a, epsilon)
    COOP_INT,intent(IN)::m, n
    COOP_REAL::a(m, n)
    COOP_REAL, optional::epsilon
    if(present(epsilon))then
       call coop_cholesky(m, n, a, epsilon)
    else
       call coop_cholesky(m, n, a)
    endif
    call coop_cholesky_inv(m, n, a)
    call coop_cholesky_sq(m, n, a)
  end subroutine coop_sympos_inverse
  

  

end module coop_matsym_mod
