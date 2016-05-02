module coop_cholesky_mod
  use coop_wrapper_typedef
  use coop_sort_mod
  use coop_sortrev_mod
  use coop_MPI_mod
  implicit none

#include "constants.h"

  private


  public:: coop_cholesky, coop_cholesky_inv, coop_cholesky_sq, coop_cholesky_solve, coop_cholesky_solve_transpose, coop_cholesky_solve_mult, coop_sympos_inverse, coop_sympos_clean, coop_fillzero, coop_fill_symmetrize

contains

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


  subroutine coop_fill_symmetrize(uplo, m, n, a)
    character uplo
    COOP_INT m, n, i, j
    COOP_REAL a(m, n)
    select case(uplo)
    case('U', 'u')
       do i=2, n
          a(1:i-1, i) = a(i, 1:i-1)
       enddo
    case('L', 'l')
       do i=1, n-1
          a(i+1:n, i) = a(i, i+1:n)
       enddo
    case default
       stop 'fillzero: unknown uplo option'
    end select
  end subroutine coop_fill_symmetrize

 
 

  subroutine coop_sympos_clean(m, n, a, epsilon)
    COOP_INT::m, n, i, j
    COOP_REAL::a(m, n), onemeps, epsilon, tmp, incr
    COOP_REAL::sigma(n)
    onemeps = 1.d0-epsilon
    do i=1, n
       sigma(i) = sqrt(a(i,i))
       a(i, 1:i-1) =  a(i, 1:i-1) /sigma(i)
       a(i+1:n, i) =  a(i+1:n, i) /sigma(i)
       a(i, i) = 1.d0
    enddo
    do j=1, n-1
       do i=j+1, n
          if(abs(a(i, j)) .lt. epsilon)then  !!add a positive definite matrix to eliminate a(i, j)
             a(i, i) = a(i, i) + abs(a(i, j))
             a(j, j) = a(j, j) + abs(a(i, j))
             a(i, j)=0.d0
          elseif(abs(a(i,j)) .gt. onemeps)then
             tmp = sign(onemeps, a(i, j))
             a(i, i) = a(i, i) + abs(a(i, j) - tmp)
             a(j, j) = a(j, j) + abs(a(i, j) - tmp)             
             a(i, j) = tmp
          endif
       enddo
    enddo
    do i=1, n
       a(i, 1:i) =  a(i, 1:i) * sigma(i)
       a(i:n, i) =  a(i:n, i) * sigma(i)
    enddo
  end subroutine coop_sympos_clean

  subroutine coop_cholesky(m, n, a, info)
    COOP_INT m, n, i, j, info
    COOP_REAL a(m, n)
    info = 0
    if(a(1,1).le.0.d0)then
       info = 1
       return
    endif
    if(n.eq.1)then
       a(1,1) = sqrt(a(1,1))
       return
    endif
#ifdef HAS_LAPACK
    call dpotrf('L', n, a, m, info)
    call coop_fillzero('U',m, n, a)
#else    
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
    COOP_INT info, i
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
       if(info.ne.0)then
          if(m.lt. 10 .and. n.lt. 10)then
             do i=1, n
                write(*,"("//COOP_STR_OF(n)//"G14.5)") a(i, 1:n)
             enddo
          endif
          call Coop_return_error("coop_sympos_inverse", "the matrix is not positive definite", "stop")
       endif
#ifdef HAS_LAPACK
       call dpotri("L", n, a, m, info)
       if(info.ne.0) call Coop_return_error("coop_sympos_inverse", "the matrix is not positive definite", "stop")
       call coop_fill_symmetrize('U', m, n, a)
#else       
       call coop_cholesky_inv(m, n, a)
       call coop_cholesky_sq(m, n, a)
#endif       
    end select
  end subroutine coop_sympos_inverse

  
end module coop_cholesky_mod
