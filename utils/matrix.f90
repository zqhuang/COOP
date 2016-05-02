!Matrix utility routines. Uses BLAS/LAPACK. Mostly wrapper routines.

module coop_matrix_mod
  use coop_wrapper_typedef
  use coop_sort_mod
  use coop_sortrev_mod
  use coop_MPI_mod
  use coop_cholesky_mod
  implicit none

#include "constants.h"

  private

  public::coop_write_matrix, coop_print_matrix, coop_read_matrix, coop_set_identity_matrix, coop_identity_matrix, coop_diagonal_matrix, coop_matrix_add_diagonal, coop_matrix_sum_columns,  coop_matrix_solve_small,  Coop_matrix_Solve, Coop_matrix_Inverse, coop_matsym_xCx, coop_matrix_det_small,coop_matsym_mat2vec, coop_matsym_vec2mat, coop_matsym_Diagonalize,  Coop_matsym_Sqrt,  coop_matsym_sqrt_small,  coop_matsym_power_small,  Coop_matsym_power,  coop_matsym_function, Coop_matsym_LnDet, Coop_matsym_Solve,  coop_matrix_sorted_svd,  coop_covmat, coop_matrix_dominant_eigen_value, coop_solve_constrained, coop_matsym_diag, coop_matrix_ludcmp, coop_matrix_lubksb, coop_matsym_get_indices

  type coop_covmat   !!assume only lower triangle is saved in C; !!invC contains full matrix; L is lower triangle Cholesky, with zero filled.
     COOP_INT::n = 0
     COOP_REAL::mult = 0.d0
     COOP_REAL,dimension(:,:),allocatable::C, invC, L
     COOP_REAL,dimension(:),allocatable::mean, sigma
   contains
     procedure::free => coop_covmat_free
     procedure::alloc => coop_covmat_alloc
     procedure::diagonal => coop_covmat_diagonal
     procedure::write => coop_covmat_write
     procedure::read => coop_covmat_read
     procedure::import => coop_covmat_import
     procedure::export => coop_covmat_export
     procedure::normalize => coop_covmat_normalize     
     procedure::chi2 => coop_covmat_chi2
     procedure::invert => coop_covmat_invert
     procedure::Cholesky => coop_covmat_Cholesky
     procedure::MPI_sync => coop_covmat_MPI_sync
  end type coop_covmat



  interface coop_write_matrix
     module procedure coop_write_matrix_s, coop_write_matrix_d
  end interface coop_write_matrix

  
  interface coop_print_matrix
     module procedure coop_print_matrix_s, coop_print_matrix_d
  end interface coop_print_matrix
  
  interface coop_read_matrix
     module procedure coop_read_matrix_s, coop_read_matrix_d
  end interface coop_read_matrix

contains

!!****section 1: lapack-independent routines for small matrices **********!!
!!========================================================================!!
!!**slower and probably also less accurate than lapack lib ***************!!
!!****** ok for small matrices with rank ~ a few tens ********************!!
!!******* should not use for huge matrices *******************************!!
!!========================================================================!!

  subroutine coop_set_identity_matrix(n, a)
    COOP_INT i, n
    COOP_REAL  a(n,n)
    a = 0
    !$omp parallel do
    do i=1,n
       a(i,i) = 1.d0
    enddo
    !$omp end parallel do
  end subroutine coop_set_identity_matrix

  function coop_identity_matrix(n) result(a)
    COOP_INT i, n
    COOP_REAL  a(n,n)
    a = 0
    !$omp parallel do
    do i=1,n
       a(i,i) = 1.d0
    enddo
    !$omp end parallel do
  end function coop_identity_matrix

  function coop_diagonal_matrix(diag) result(a)
    COOP_REAL ,dimension(:),intent(IN)::diag
    COOP_REAL  a(size(diag),size(diag))
    COOP_INT i
    !$omp parallel do
    do i=1,size(diag)
       a(i,i)=diag(i)
    enddo
    !$omp end parallel do
  end function coop_diagonal_matrix

  subroutine coop_matrix_add_diagonal(n, mat, diag)
    COOP_INT n, i
    COOP_REAL  mat(n,n), diag(n)
    !$omp parallel do
    do i=1,n
       mat(i,i) = mat(i,i)+diag(i)
    enddo
    !$omp end parallel do
  end subroutine coop_matrix_add_diagonal

  Subroutine Coop_write_matrix_d(funit, mat,nx,ny)
    COOP_INT funit
    COOP_INT,optional::nx,ny
    COOP_REAL  mat(:,:)
    COOP_INT i
    COOP_SHORT_STRING Form
    if(present(nx).and.present(ny))then
       Form="("//Trim(Coop_num2str(ny))//"E16.7)"
       do i=1,nx
          write(funit,trim(form))mat(i,1:ny)
       enddo
    else
       Form="("//Trim(Coop_num2str(size(mat,2)))//"E13.4)"
       do i=1,size(mat,1)
          write(funit,trim(form))mat(i,:)
       enddo
    end if
  End Subroutine Coop_write_matrix_d

  Subroutine coop_print_matrix_d(mat, nx, ny)
    COOP_INT,optional:: nx,ny
    COOP_REAL  mat(:,:)
    COOP_INT i
    character(Len=128)Form
    if(present(nx).and.present(ny))then
       Form="("//Trim(Coop_num2str(ny))//"E13.4)"
       do i=1,nx
          write(*,trim(form))mat(i,1:ny)
       enddo
    else
       Form="("//Trim(Coop_num2str(size(mat,2)))//"E13.4)"
       do i=1,size(mat,1)
          write(*,trim(form))mat(i,:)
       enddo
    endif
  End Subroutine Coop_print_matrix_d

  Subroutine Coop_read_matrix_d(funit, mat, nx, ny, success)
    logical,optional::success
    COOP_INT funit,nx,ny
    COOP_REAL  mat(:, :)
    COOP_INT i
    COOP_LONG_STRING line
    if(ny.gt. 256)then
       call coop_feedback("warning: reading a huge matrix: will not check comment lines")
       do i=1,nx
          read(funit,*) mat(i,1:ny)
       enddo
    else
       i = 1
       do
          read(funit, '(A)', ERR=100, end=100) line
          line =  adjustl(trim(line))
          if(line(1:1) .eq. "#" .or. trim(line) .eq. "")cycle
          read(line, *, ERR=100) mat(i, 1:ny)
          i = i+1
          if(i.gt. nx) exit
       enddo
    endif
    if(present(success)) success = .true.
    return
100 if(present(success))then
       success = .false.
       return
    else
       stop "Coop_read_matrix: matrix data error"
    endif
  End Subroutine Coop_read_matrix_d


  Subroutine Coop_write_matrix_s(funit, mat,nx,ny)
    COOP_INT funit
    COOP_INT,optional::nx,ny
    COOP_SINGLE mat(:,:)
    COOP_INT i
    COOP_SHORT_STRING form
    if(present(nx).and.present(ny))then
       Form="("//Trim(Coop_num2str(ny))//"E16.7)"
       do i=1,nx
          write(funit,trim(form))mat(i,1:ny)
       enddo
    else
       Form="("//Trim(Coop_num2str(size(mat,2)))//"E13.4)"
       do i=1,size(mat,1)
          write(funit,trim(form))mat(i,:)
       enddo
    end if
  End Subroutine Coop_write_matrix_s

  Subroutine coop_print_matrix_s(mat, nx, ny)
    COOP_INT,optional:: nx,ny
    COOP_SINGLE mat(:,:)
    COOP_INT i
    COOP_SHORT_STRING form
    if(present(nx).and.present(ny))then
       Form="("//Trim(Coop_num2str(ny))//"E13.4)"
       do i=1,nx
          write(*,trim(form))mat(i,1:ny)
       enddo
    else
       Form="("//Trim(Coop_num2str(size(mat,2)))//"E13.4)"
       do i=1,size(mat,1)
          write(*,trim(form))mat(i,:)
       enddo
    endif
  End Subroutine Coop_print_matrix_s

  Subroutine Coop_read_matrix_s(funit,mat, nx,ny, success)
    logical,optional::success
    COOP_INT funit,nx,ny
    COOP_SINGLE mat(:,:)
    COOP_INT i
    COOP_LONG_STRING line
    if(ny.gt. 256)then
       call coop_feedback( "warning: reading a huge matrix: will not check comment lines")
       do i=1,nx
          read(funit,*) mat(i,1:ny)
       enddo
    else
       i = 1
       do
          read(funit, '(A)', ERR=100, end=100) line
          line =  adjustl(trim(line))
          if(line(1:1) .eq. "#" .or. trim(line) .eq. "")cycle
          read(line, *, ERR=100) mat(i, 1:ny)
          i = i+1
          if(i.gt. nx) exit
       enddo
    endif
    if(present(success)) success = .true.
    return
100 if(present(success))then
       success = .false.
       return
    else
       stop "Coop_read_matrix: matrix data error"
    endif
  End Subroutine Coop_read_matrix_s


  function coop_matrix_sum_columns(a, m, n) result(s)
    COOP_INT m, n, i
    COOP_REAL  a(m, n), s(n)
    do i = 1, n
       s(i) = sum(a(1:m,i))
    enddo
  end function coop_matrix_sum_columns

  !!return Tr(A)
  function coop_matrix_trace(A) result(tr)
    COOP_REAL ,dimension(:,:),intent(in)::A
    COOP_INT i
    COOP_REAL  tr
    tr = A(1,1)
    do i=2, size(A, 1)
       tr = tr + A(i,i)
    enddo
  end function coop_matrix_trace

  function coop_matrix_product_trace(A, B) result(tr)
    COOP_REAL, dimension(:, :)::A, B
    COOP_REAL::tr
    COOP_INT::i
    if(size(A, 1) .ne. size(B, 2) .or. size(A, 2) .ne. size(B, 1)) stop "wrong input in matrix_product_trace"
    tr = dot_product(A(1, :), B(:, 1))
    do i = 2, size(A, 1)
       tr = tr + dot_product(A(i, :), B(:, i))
    enddo
  end function coop_matrix_product_trace


  !!solve Ax=b and save x in b; A is not destroyed (while coop_matrix_solve does destroy a)
  subroutine coop_matrix_solve_small(n, a, b)
    COOP_INT n
    COOP_REAL  a(n,n), b(n), acopy(n,n)
    COOP_INT indx(n), info
#ifndef HAS_LAPACK
    COOP_REAL  d
#endif
    select case(n)
    case(1)
       b = b/a(1,1)
    case(2)
       b = (/ A(2,2)*b(1)-A(1,2)*b(2), A(1,1)*b(2)-A(2,1)*b(1) /) &
            / ( A(1,1)*A(2,2)-A(1,2)*A(2,1) )
    case(3)
       b = (/ (b(1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))+a(1,2)*(a(2,3)*b(3)-b(2)*a(3,3))+a(1,3)*(b(2)*a(3,2) - a(2,2)*b(3))) , &
            (a(1,1)*(b(2)*a(3,3)-a(2,3)*b(3))+b(1)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+a(1,3)*(a(2,1)*b(3) - b(2)*a(3,1))) , &
            (a(1,1)*(a(2,2)*b(3)-b(2)*a(3,2))+a(1,2)*(b(2)*a(3,1)-a(2,1)*b(3))+b(1)*(a(2,1)*a(3,2) - a(2,2)*a(3,1))) /) &
            / (a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))+a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+a(1,3)*(a(2,1)*a(3,2) - a(2,2)*a(3,1)))
    case default
       acopy = a
#ifdef HAS_LAPACK
       call dgesv(n, 1, acopy, n, indx, b, n,info)
#else
       call coop_matrix_LUdcmp(acopy, indx, d, n)
       call coop_matrix_LUbksb(acopy, indx, b, n)
#endif
    end select
  end subroutine coop_matrix_solve_small

  !!solve A x = B
  !!INPUT: A(n,n), B(n,m)
  !!OUTPUT: solution returned in B
  subroutine Coop_matrix_Solve(a,b)
    COOP_INT n,m
    COOP_REAL ,dimension(:,:)::a,b
    COOP_INT indx(size(a,1))
    COOP_INT i
#ifndef HAS_LAPACK
    COOP_REAL  d
#endif
    n = Coop_getdim("Coop_matrix_Solve", size(a,1), size(a,2), size(b,1))
    m = size(B,2)
#ifdef HAS_LAPACK
    call dgesv(n, m, a, n, indx, b, n, i)
#else
    call coop_matrix_LUdcmp(a,indx,d,n)
    Do i = 1,m
       call coop_matrix_LUbksb(a, indx,b(:,i),n)
    Enddo
#endif
  End subroutine Coop_matrix_Solve

  !!return Inverse of A
  Subroutine Coop_matrix_Inverse(A)
    COOP_REAL ,dimension(:,:)::A
    COOP_INT i
#ifdef HAS_LAPACK
    COOP_INT n
    COOP_INT indx(size(A,1))
    COOP_REAL  work(4*size(A,1))
    n = Coop_getdim("Coop_matrix_Inverse", size(a,1), size(a,2))
    call dgetrf(n, n, a, n, indx, i)
    if(i.ne.0)then
       call coop_return_error( "Coop_matrix_Inverse", "singular matrix cannot be inverted", "stop")
    endif
    call dgetri(n, a, n, indx, work, 4*n, i)
#else
    COOP_REAL  acopy(size(A,1), size(A,2))
    Acopy=a
    a=0.d0
    do i=1,size(a)
       a(i,i)=1.d0
    enddo
    call Coop_matrix_Solve(acopy,a)
#endif
  end Subroutine Coop_matrix_Inverse


  function coop_matsym_xCx(C, x)  result(xCx)
    COOP_REAL ,dimension(:,:)::C
    COOP_REAL ,dimension(:)::x
    COOP_REAL  xCx
    COOP_INT n, i, j
    n =  Coop_getdim("coop_matrix_xCx", size(x), size(C,1), size(C,2))
    xCx = C(1,1)*x(1)**2/2.d0
    do i=2,n
       xCx = xCx + C(i,i)*x(i)**2/2.d0
       do j=1,i-1
          xCx = xCx + C(j, i)*x(i)*x(j)
       enddo
    enddo
    xCx = xCx * 2.d0
  end function coop_matsym_xCx

  subroutine coop_matrix_det_small(n, A, det)
    COOP_INT n
    COOP_REAL  A(n, n), det
    select case(n)
    case(1)
       det = A(1,1)
    case(2)
       det = A(1,1)*A(2,2)-A(1,2)*A(2,1)
    case(3)
       det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) + A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3)) + A(1,3)*(A(2,1)*A(3,2)- A(2,2)*A(3,1))
    case(4)
       det =  A(1,1)*(A(2,2)*(A(3,3)*A(4,4) - A(3,4)*A(4,3)) + A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4)) + A(2,4)*(A(3,2)*A(4,3)- A(3,3)*A(4,2))) &
            + A(1,2)*(A(2,3)*(A(3,4)*A(4,1) - A(3,1)*A(4,4)) + A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)) + A(2,1)*(A(3,3)*A(4,4)- A(3,4)*A(4,3))) &
            + A(1,3)*(A(2,4)*(A(3,1)*A(4,2) - A(3,2)*A(4,1)) + A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2)) + A(2,2)*(A(3,4)*A(4,1)- A(3,1)*A(4,4))) &
            + A(1,4)*(A(2,1)*(A(3,2)*A(4,3) - A(3,3)*A(4,2)) + A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)) + A(2,3)*(A(3,1)*A(4,2)- A(3,2)*A(4,1)))
    case default
       stop "coop_matrix_det_small only works for n<=4"
    end select
  end subroutine coop_matrix_det_small

  !!====================== for positive definite symetric matrix ================
  subroutine coop_matsym_mat2vec(n,mat,vec)
    COOP_INT n
    COOP_REAL ,intent(IN):: mat(n,n)
    COOP_REAL ,intent(OUT)::vec(n*(n+1)/2)
    COOP_INT i, j, k
    k = 1
    do i = 0, n-1
       do j=1, n-i
          vec(k) = mat(j, j+i)
          k = k + 1
       enddo
    enddo
  end subroutine coop_matsym_mat2vec

  subroutine coop_matsym_vec2mat(n, vec, mat)
    COOP_INT n
    COOP_REAL ,intent(IN)::vec(n*(n+1)/2)
    COOP_REAL ,intent(OUT):: mat(n,n)
    COOP_INT i, j, k
    k = 1
    do i = 0, n-1
       do j=1, n-i
          mat(j, j+i) = vec(k) 
          k = k + 1
       enddo
    enddo
    do i = 1, n
       do j=1, i-1
          mat(i, j) = mat(j, i)
       enddo
    enddo
  end subroutine coop_matsym_vec2mat

  !!find R such that R^T H R = diag(e), and R^T R = 1
  !!EIGEN VALUES ARE STORED IN E WITH E(1)<=E(2)...(from ground state..)
  !!EIGEN VECTORS ARE STORED IN COLUMNS OF H(:,1), H(:,2),...
  Subroutine Coop_matsym_Diagonalize(m, n, H, e)
    COOP_INT::m, n
    COOP_REAL::h(m, n)
    COOP_REAL::e(n)
#ifdef HAS_LAPACK
    call coop_matrix_diagonalize(m,n,h,e)
#else
    COOP_REAL  psi(n, n)
    COOP_INT Indx(n)
    COOP_INT i
    call Coop_matsym_diag(m, n, h, psi, 1.d-8)
    do i = 1,N
       e(i) = h(i,i)
    enddo
    call Coop_quicksortacc(E,Indx)
    do i=1,n
       h(:,I)=psi(:,indx(I))
    enddo
#endif

  End Subroutine Coop_matsym_Diagonalize


  subroutine coop_matsym_sqrt_small(n, a)
    COOP_INT n
    COOP_REAL  a(n,n)
    COOP_REAL  s, t
    select case(n)
    case(1)
       a = sqrt(a)
    case(2)
       s = sqrt(a(1,1)*a(2,2)-a(1,2)*a(2,1))
       t = sqrt(a(1,1)+a(2,2)+2*s)
       a = reshape( (/ a(1,1)+s, a(2,1), a(1,2), a(2,2) + s /) / t, (/ 2, 2 /))
    case default
       call coop_matsym_sqrt(a)
    end select
  end subroutine coop_matsym_sqrt_small

  subroutine Coop_matsym_Sqrt(A, mineig)
    !!get the square root of a positive definite symmetric matrix A
    COOP_REAL ,dimension(:,:)::A
    COOP_REAL,optional::mineig
    if(present(mineig))then
       call coop_matsym_power(A, 0.5d0, mineig)
    else
       call coop_matsym_power(A, 0.5d0)
    endif
  End subroutine Coop_matsym_Sqrt

  subroutine coop_matsym_power_small(n, A, alpha)
    COOP_INT n
    COOP_REAL  A(n, n)
    COOP_REAL  alpha, l1, l2, s, t, u(2,2), trans(2,2), norm
    select case(n)
    case(1)
       A = A**alpha
    case(2)
       if(A(1,2).eq.0.d0)then
          A(1,1) = A(1,1)**alpha
          A(2,2) = A(2,2)**alpha
          return
       endif
       s = A(1,1) + A(2,2)
       t = sqrt((A(1,1)-A(2,2))**2+4.d0*A(1,2)**2)
       l1 = (s+t)/2.d0
       l2 = (s-t)/2.d0
       norm = sqrt(A(1,2)**2+(l1 - A(1,1))**2)
       u(:,1) = (/ A(1,2)/norm, (l1-A(1,1))/norm /)
       norm = sqrt(A(1,2)**2+(l2 - A(1,1))**2)
       u(:,2) = (/ A(1,2)/norm, (l2-A(1,1))/norm /)
       trans(1,:) = u(:,1)*l1**alpha
       trans(2,:) = u(:,2)*l2**alpha
       a = matmul(u, trans)
    case default
       call coop_matsym_power(a, alpha)
    end select
  end subroutine coop_matsym_power_small


  subroutine Coop_matsym_power(A,alpha, mineig)
    !!get A^alpha for a positive definite symmetric matrix A
    COOP_REAL ,dimension(:,:)::A
    COOP_REAL alpha
    COOP_REAL  E(size(A,1)), trans(size(A,1), size(A,1))
    COOP_INT n, i
    COOP_REAL, optional::mineig
    COOP_REAL::me
    n= Coop_getdim("coop_matsym_sqrt", size(A,1),size(A,2))
    if(present(mineig))then
       me = mineig
    else
       me = 1.d-10
    endif
#ifdef HAS_LAPACK
    call coop_matsym_diagonalize(n, n, a,e)
    me = maxval(abs(e))*me
    where (e.lt. me)
       e = 0.d0
    elsewhere
       e = e**alpha       
    end where

    do i=1,n
       trans(i,:) = A(:,i)*e(i)
    enddo
    A = matmul(A, trans)
#else
    call coop_svd_decompose(n,n,a,e,trans)
    me = maxval(abs(e))*me    
    where (e.lt. me)
       e = 0.d0
    elsewhere
       e = e**alpha       
    endwhere
    do i=1,n
       trans(:,i) = trans(:,i)*e(i)
    enddo
    a = matmul(a, transpose(trans))
#endif
  End subroutine Coop_matsym_power

  subroutine coop_matsym_function(n, A, f)
    COOP_INT n,i
    COOP_REAL  E(n), trans(n, n)
    COOP_REAL  A(n, n)
    external f
    COOP_REAL  f
    if(n.eq.1)then
       A(1,1) = f(A(1,1))
       return
    endif
#ifdef HAS_LAPACK
    call coop_matsym_diagonalize(n, n, a, e)
    do i=1,n
       trans(i,:) = A(:,i)*f(e(i))
    enddo
    A = matmul(A, trans)
#else
    call coop_svd_decompose(n,n,a,e,trans)
    do i=1,n
       trans(:,i) = trans(:,i)*f(e(i))
    enddo
    A = matmul(A, transpose(trans))
#endif    
  end subroutine coop_matsym_function
 


  function Coop_matsym_LnDet(n,A)  result(LnSymMatDet)
    COOP_INT n,i,j
    COOP_REAL ,intent(in)::A(n,n)
    COOP_REAL  LnSymMatDet
    COOP_REAL  Acopy(n,n)
    select case(n)
    case(1)
       LnSymMatDet = log(A(1,1))
       return
    case(2)
       LnSymMatDet = log(A(1,1)*A(2,2)-A(2,1)**2)
       return
    case default
       Acopy(1,1)=SQRT(A(1,1))
       LnSymMatDet = log(Acopy(1,1))
       Acopy(2:N,1)=A(2:N,1)/Acopy(1,1)
       do I=2,N-1
          Acopy(I,I)=SQRT(A(I,I)-sum(Acopy(I,1:I-1)**2)) !!L(I,I) doNE
          do J=I+1,N  !! NOW CALCULATE L(J,I), J>=I
             Acopy(J,I)=(A(J,I)-sum(Acopy(J,1:I-1)*Acopy(I,1:I-1)))/Acopy(I,I)
          enddo
          LnSymMatDet = LnSymMatDet+log(Acopy(i,i))
       enddo
       LnSymMatDet = LnSymMatDet*2+log(A(n,n)-sum(Acopy(n,1:n-1)**2))
    end select
  end function Coop_matsym_LnDet

  !!this solves Ax=b, where A is postive symmetric matrix.
  !! the result is saved in b, A is NOT destroyed (while coop_matsym_solve does destroy A)
  subroutine coop_matsym_solve_small(n, A, b)
    COOP_INT n, info
    COOP_REAL  a(n,n), b(n), acopy(n,n)
    select case(n)
    case(1)
       b = b/A(1,1)
    case(2)
       b = (/ A(2,2)*b(1)-A(1,2)*b(2), A(1,1)*b(2)-A(2,1)*b(1) /) &
            / ( A(1,1)*A(2,2)-A(1,2)*A(2,1) )
    case(3)
       b = (/ (b(1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))+a(1,2)*(a(2,3)*b(3)-b(2)*a(3,3))+a(1,3)*(b(2)*a(3,2) - a(2,2)*b(3))) , &
            (a(1,1)*(b(2)*a(3,3)-a(2,3)*b(3))+b(1)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+a(1,3)*(a(2,1)*b(3) - b(2)*a(3,1))) , &
            (a(1,1)*(a(2,2)*b(3)-b(2)*a(3,2))+a(1,2)*(b(2)*a(3,1)-a(2,1)*b(3))+b(1)*(a(2,1)*a(3,2) - a(2,2)*a(3,1))) /) &
            / (a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))+a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+a(1,3)*(a(2,1)*a(3,2) - a(2,2)*a(3,1)))
    case default
       acopy = a
#ifdef HAS_LAPACK
       call dposv('L', n, 1, acopy, n, b, n, info)
       if(info.ne.0) call Coop_return_error("coop_matsym_solve_small", "the matrix is not positive definite", "stop")
#else
       call Coop_cholesky(n, n, acopy, info)
       if(info.ne.0) call Coop_return_error("coop_matsym_solve_small", "the matrix is not positive definite", "stop")       
       call Coop_cholesky_solve(n, n, acopy, b)
#endif
    end select
  end subroutine coop_matsym_solve_small

  Subroutine Coop_matsym_Solve(A,b)
    COOP_REAL ,dimension(:,:),INTENT(INOUT)::A
    COOP_REAL ,dimension(:,:),INTENT(INOUT)::b
    COOP_INT n, m, info, i
    n=COOP_GETDIM("Coop_matsym_Solve",SIZE(A,1),SIZE(A,2),SIZE(b,1))
    m = size(b, 2)
#ifdef HAS_LAPACK
    call dposv('L', n, m, a, n, b, n, info)
    if(info.ne.0) call Coop_return_error("coop_matsym_solve", "the matrix is not positive definite", "stop")           
#else
    call Coop_cholesky(n, n, a, info)
    if(info.ne.0) call Coop_return_error("coop_matsym_solve", "the matrix is not positive definite", "stop")               
    do i=1,m
       call Coop_cholesky_solve_mult(n, n, n, m, a, b)
    enddo
#endif
  End Subroutine Coop_matsym_Solve

!!******section 2: lapack library routines *********************************!!
!!==========================================================================!!
!!**************compiled only if lapack library is available****************!!
!!==========================================================================!!
#ifdef HAS_LAPACK

  subroutine coop_matrix_least_square_all(m, n, nrhs, A, b, x)
    !!find x that minimize ||Ax-B||
    COOP_INT,intent(IN)::m,n, nrhs
    COOP_REAL ,intent(INOUT)::A(m, n)
    COOP_REAL ,intent(INOUT)::b(m, nrhs)
    COOP_REAL ,intent(OUT)::x(n, nrhs)
    COOP_REAL ,dimension(:),allocatable::work
    COOP_INT lwork, info
    lwork = min(m,n) + max(1,m,n,nrhs)
    allocate(work(lwork))
    if(m.le.n)then
       x(1:m, 1:nrhs) = b
       call dgels('N',m, n, nrhs, A, m, x, n, work, lwork, info)
    else
       call dgels('N',m, n, nrhs, A, m, b, m, work, lwork, info)
       x = b(1:n, 1:nrhs)
    endif
    deallocate(work)
    if(info .ne. 0) then
       call coop_return_error("coop_matrix_least_square_all", "Error info = "//trim(coop_num2str(info)), "stop")
    endif
  end subroutine coop_matrix_least_square_all

  subroutine coop_matrix_least_square_one(m, n, A, b, x)
    !!find x that minimize || A x - B ||
    COOP_INT,intent(IN) :: m, n
    COOP_REAL ,intent(IN) :: a(m, n)
    COOP_REAL ,intent(INOUT) :: b(m)
    COOP_REAL ,intent(OUT) :: x(n)
    COOP_REAL ,dimension(:),allocatable :: work
    COOP_INT lwork, info
    COOP_REAL  acpy(m, n)
    acpy = a
    lwork = min(m,n) + max(1,m,n)
    allocate(work(lwork))
    if(m.le.n)then
       x(1:m) = b
       call dgels('N', m, n, 1, acpy, m, x, n, work, lwork, info)
    else
       call dgels('N', m, n, 1, acpy, m, b, m, work, lwork, info)
       x = b(1:n)
    endif
    deallocate(work)
    if(info .ne. 0) then
       call coop_return_error("coop_matrix_least_square_one", "Error info = "//trim(coop_num2str(info)), "stop")
    endif
  end subroutine coop_matrix_least_square_one


  subroutine coop_matrix_least_square_solve(A, x, b)
    !!find x that minimize ||Ax-b||
    !!on exit:
    !!A(m, n) is destroyed
    COOP_REAL ,dimension(:,:)::A
    COOP_REAL ,dimension(:)::b
    COOP_REAL ,dimension(:),intent(out)::x
    COOP_REAL ,dimension(:),allocatable::work, bcopy
    COOP_INT m,n, lwork, info
    m = Coop_getdim("matrix_least_square_solve", size(A, 1), size(b))
    n = Coop_getdim("matrix_least_square_solve", size(A, 2), size(x))
    if(m .lt. n) call Coop_return_error("matrix_least_square_solve", "A is not full rank", "stop")
    lwork = max(1, min(m,n)*2)
    allocate(work(lwork), bcopy(m))
    bcopy = b
    call dgels('N',m, n, 1, A, m, bcopy, m, work, lwork, info)
    if(info .ne. 0) call Coop_return_error("matrix_least_square_solve","solution cannot be found","return")
    x = bcopy(1:n)
    deallocate(work, bcopy)
  end subroutine coop_matrix_least_square_solve


  !Does m = U diag U^T, returning U in real symmetric matrix M
  subroutine coop_Matrix_Diagonalize(m, n, c, diag)
    COOP_INT::m, n
    COOP_REAL , intent(inout):: c(m, n )
    COOP_REAL , intent(out) :: diag(n)
    COOP_INT ierr, tmpsize
    COOP_REAL , allocatable, dimension(:) :: tmp
    tmpsize =  9*n+1
    allocate(tmp(tmpsize));
    call DSYEV('V','L',n,c,m,diag,tmp,tmpsize,ierr) !evalues and vectors of symmetric matrix
    if (ierr .ne. 0) call Coop_return_error("matrix_diagonalize", "cannot diagonalize", "return")
    deallocate(tmp)
  end subroutine Coop_Matrix_Diagonalize

#endif

!!auxilliary routines
  Subroutine coop_matrix_ludcmp(A,INDX,D,n)
    COOP_INT n
    COOP_REAL  a(n,n)
    COOP_INT  indx(n)
    COOP_REAL ,INTENT(OUT) :: d
    COOP_REAL  vv(n)
    COOP_REAL ,parameter::TINY=1.E-20
    COOP_INT  j,imax
    D=1.D0
    vv=maxval(abs(a),dim=2)
    if (any(vv .EQ. 0.0)) STOP 'singular matrix in coop_matrix_ludcmp'
    vv=1.d0/vv
    do j=1,n
       imax=(j-1)+coop_maxloc(vv(j:n)*Dabs(a(j:n,j)))
       if (j .NE. imax) then
          call coop_swap(a(imax,1:N),a(j,1:N))
          d=-d
          vv(imax)=vv(j)
       endif
       indx(j)=imax
       if (a(j,j) .EQ. 0.0) a(j,j)=TINY
       a(j+1:n,j)=a(j+1:n,j)/a(j,j)
       a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-coop_outerprod(a(j+1:n,j),a(j,j+1:n))
    enddo
  end Subroutine coop_matrix_ludcmp

  Subroutine coop_matrix_lubksb(A,indx,b,n)
    COOP_INT n
    COOP_REAL  A(n,n)
    COOP_INT indx(n)
    COOP_REAL  b(n)
    COOP_INT:: i,ii,ll
    COOP_REAL  :: summ
    ii=0
    do i=1,n
       ll=indx(i)
       summ=b(ll)
       b(ll)=b(i)
       if (ii .ne. 0) then
          summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
       else if (summ .ne. 0.0) then
          ii=i
       end if
       b(i)=summ
    end do
    do i=n,1,-1
       b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
    end do
  end Subroutine coop_matrix_lubksb

  !!********************************************************
  !!matrix diagonalization; diag(eig1; eig2; ..; eign ) = R^T A R
  !!return R;
  !!return A=diag(eig1, eig2, ..., eign)
  Subroutine Coop_matsym_diag(m, n, a, R, precision)
    COOP_INT m,n, i    
    COOP_INT,parameter::MAXLOOP=100000
    COOP_REAL::a(m, n), R(n, n), l1, l2
    COOP_REAL, optional::precision
    COOP_REAL::prec
    Logical Tag
    if(n.eq.1)then
       R=1.d0
       return
    endif
    if(n.eq.2)then
       if(A(1,2).eq.0.d0)then
          R(1,1) = 1.d0
          R(2,2) = 1.d0
          R(1,2) = 0.d0
          R(2,1) = 0.d0
          return
       endif
       prec = max(sqrt((A(1,1)-A(2,2))**2/4.d0+A(1,2)**2), 1.d-99)
       l1 = (A(1,1)+A(2,2))/2.d0 + prec
       l2 = (A(1,1)+A(2,2))/2.d0 - prec
       prec = max(sqrt((A(1,1)- l1)**2 + A(1, 2)**2), 1.d-99)
       R(1, 1) = A(1, 2)/prec
       R(2, 1) = (l1 - A(1,1))/prec
       prec = max(sqrt((A(1,1)-l2)**2 + A(1, 2)**2), 1.d-99)
       R(1, 2) = A(1,2)/prec
       R(2,2) = (l2 - A(1,1))/prec
       A(1,1) = l1
       A(2,2) = l2
       A(1,2)=0.d0
       A(2,1) = 0.d0
       return
    endif
    if(present(precision))then
       prec = precision*sum(abs(A(1:n, 1:n)))/N/N
    else
       prec = 1.d-6*sum(abs(A(1:n, 1:n)))/N/N       
    endif
    R=0.d0
    do i=1,n
       R(i,i)=1.d0
    enddo
    i = 0
    Tag=.true.
    do while(i .lt.MAXLOOP .and. tag)
       call Coop_matsym_diag_step(m, A, R, N, prec, tag)
       I=I+1
    enddo
  end Subroutine Coop_matsym_diag

  Subroutine Coop_matsym_rotang(m, n, A, m1, n1, tant,cost,sint)
    COOP_INT m, n, m1, n1    !!m1 > n1
    COOP_REAL ::A(m, n)
    COOP_REAL ,parameter::EPS=1.D-6
    COOP_REAL ,parameter::LARGEPAR=1.d0/EPS

    COOP_REAL  tant,cost,sint,S,S2
    s = (A(n1,n1)-A(m1,m1))*0.5d0/A(m1,n1)
    s2 = s**2
    if(S2.lt.EPS)then
       tant=sign(1.d0+S2*(0.5d0-S2/8.d0),S)-S
       S2=(tant*tant-1.d0)/2.d0
       cost=(1.d0+S2*(-0.5d0+S2*(0.375D0-0.3125D0*S2)))/coop_sqrt2
    elseif(S2.gt.LARGEPAR)then
       tant=0.5d0/S*(1.d0-0.25D0/S2)
       S2=tant*tant
       cost=1.d0+S2*(-0.5d0+S2*(0.375D0-0.3125D0*S2))
    else
       tant=sign(dsqrt(S2+1.d0),S)-S
       cost=1.d0/dsqrt(1.d0 + tant**2)
    endif
    sint  = tant * cost
  end Subroutine Coop_matsym_rotang

  Subroutine Coop_matsym_findmax(m, A,N,m1,n1,precision,tag)
    COOP_INT m, m1, n1, n, i, j
    COOP_REAL  A(m, n)
    COOP_REAL  precision,temp,Y
    logical tag
    temp = precision
    do I=2,N
       do J=1,I-1
          Y=abs(A(I,J))
          if(Y.gt.temp)then
             temp=Y
             m1=I
             n1=J
          endif
       enddo
    enddo
    tag = (temp .gt. precision)
  end  Subroutine Coop_matsym_findmax

  Subroutine Coop_matsym_rot(R,N,m1,n1,cost,sint)
    COOP_INT n, m1, n1, i
    COOP_REAL  cost,sint
    COOP_REAL  R(N,N),Rm1(N)
    do i = 1, n
       Rm1(i)=R(i,m1)*cost-R(i,n1)*sint
       R(i,n1)=R(i,n1)*cost+R(i,m1)*sint
    enddo
    R(:,m1)=Rm1
  end Subroutine Coop_matsym_rot

  Subroutine Coop_matsym_diag_step(m, A, R, N, precision, tag)
    COOP_INT m, n, n1, m1  !!n1<m1
    COOP_REAL  A(m, n),R(n, n), precision, cost,sint,tant,temp
    logical tag
    call Coop_matsym_findmax(m,A,N,m1,n1,precision,tag)
    if(tag)then
       call Coop_matsym_rotang(m, n, A,m1,n1,tant,cost,sint)
       call Coop_matsym_rot(R,N,m1,n1,cost,sint)
       A(1:N1-1,M1)=cost*A(M1,1:N1-1)-sint*A(N1,1:N1-1)
       A(N1+1:M1-1,M1)=cost*A(M1,N1+1:M1-1)-sint*A(N1,N1+1:M1-1)
       A(M1+1:N,M1)=cost*A(M1,M1+1:N)-sint*A(N1,M1+1:N)
       A(1:N1-1,N1)=cost*A(N1,1:N1-1)+sint*A(M1,1:N1-1)
       A(N1+1:M1-1,N1)=cost*A(N1,N1+1:M1-1)+sint*A(M1,N1+1:M1-1)
       A(M1+1:N,N1)=cost*A(N1,M1+1:N)+sint*A(M1,M1+1:N)
       A(M1,1:N1-1)=A(1:N1-1,M1)
       A(M1,N1+1:M1-1)=A(N1+1:M1-1,M1)
       A(M1,M1+1:N)=A(M1+1:N,M1)
       A(N1,1:N1-1)=A(1:N1-1,N1)
       A(N1,N1+1:M1-1)=A(N1+1:M1-1,N1)
       A(N1,M1+1:N)=A(M1+1:N,N1)
       TEMP=tant*A(M1,N1)
       A(M1,M1)=A(M1,M1)-TEMP
       A(N1,N1)=A(N1,N1)+TEMP
       A(M1,N1)=0.d0
       A(N1,M1)=0.d0       
    endif
  end Subroutine Coop_matsym_diag_step

  subroutine coop_matsym_get_indices(n, ind, i, j)
    COOP_INT n, ind, i, j
    COOP_REAL  x
    COOP_INT imj, ix
    i = (2*n+1)**2-8*ind
    x = sqrt(dble(i)+1.d-8)
    imj = ceiling((2*n-1-x)/2.d0)
    j = ind - (2*n+1-imj)*imj/2
    i = j + imj
  end subroutine coop_matsym_get_indices

  subroutine coop_matrix_sorted_svd(n, a, w, v)
    COOP_INT n, i
    COOP_INT ind(n)
    COOP_REAL a(n, n), w(n), v(n, n), tmp(n, n)
    call coop_svd_decompose(n, n, a, w, v)
    call coop_quicksortrevacc(w, ind)
    tmp = a
    !$omp parallel do
    do i=1, n
       a(:, i) = tmp(:, ind(i))
    enddo
    !$omp end parallel do
    tmp = v
    !$omp parallel do
    do i=1, n
       v(:, i) = tmp(:, ind(i))
    enddo
    !$omp end parallel do
  end subroutine coop_matrix_sorted_svd



!!! covmat
  subroutine coop_covmat_free(this)
    class(coop_covmat)::this
    this%n = 0
    this%mult = 0.d0
    if(allocated(this%mean))deallocate(this%mean)
    if(allocated(this%sigma))deallocate(this%sigma)
    if(allocated(this%C))deallocate(this%C)
    if(allocated(this%invC))deallocate(this%invC)
    if(allocated(this%L))deallocate(this%L)            
  end subroutine coop_covmat_free

  subroutine coop_covmat_alloc(this, n)
    class(coop_covmat)::this
    COOP_INT::n
    if(this%n .eq. n)goto 100
    call this%free()
    this%n = n
    allocate(this%mean(n), this%sigma(n), this%C(n, n), this%invC(n, n), this%L(n, n))
100 this%mean = 0.d0
    this%sigma = 1.d0
    this%c = 0.d0
    this%L = 0.d0
    this%invc = 0.d0
    this%mult = 0.d0
  end subroutine coop_covmat_alloc
  
  subroutine coop_covmat_export(this, filename)
    class(coop_covmat)::this
    COOP_UNKNOWN_STRING::filename
    COOP_INT i
    open(coop_tmp_file_unit, file = trim(adjustl(filename)), form='unformatted')
    write(coop_tmp_file_unit) this%n
    write(coop_tmp_file_unit) this%mean, this%sigma
    do i=1, this%n
       write(coop_tmp_file_unit) this%C(i+1:this%n, i)
    enddo
    close(coop_tmp_file_unit)
  end subroutine coop_covmat_export

  subroutine coop_covmat_import(this, filename)
    class(coop_covmat)::this
    COOP_UNKNOWN_STRING::filename    
    COOP_INT::n, i
    open(coop_tmp_file_unit, file = trim(adjustl(filename)), form="unformatted")
    read(coop_tmp_file_unit, ERR=100, END=100) n
    call this%alloc(n)
    read(coop_tmp_file_unit, ERR=100, END=100)this%mean, this%sigma
    do i=1, n
       read(coop_tmp_file_unit, ERR = 100, END=100) this%C(i+1:n, i)
       this%C(i, i) = 1.d0
    enddo
    close(coop_tmp_file_unit)
    do i=2, this%n
       this%c(1:i-1,i) = this%c(i, 1:i-1)
    enddo    
    call this%cholesky()
    call this%invert()
    return
100 write(*,*) "file: "//trim(adjustl(filename))
    stop "covmat_read: failed"
  end subroutine coop_covmat_import


  subroutine coop_covmat_write(this, unit)
    class(coop_covmat)::this
    COOP_INT::unit
    COOP_INT i
    write(unit,"(I10)") this%n
    write(unit,"("//COOP_STR_OF(this%n)//"E16.7)") this%mean
    write(unit, "("//COOP_STR_OF(this%n)//"E16.7)")this%sigma
    do i=1, this%n
       write(unit, "("//COOP_STR_OF(this%n-i)//"E16.7)") this%C(i+1:this%n, i)
    enddo
  end subroutine coop_covmat_write

  subroutine coop_covmat_read(this, unit)
    class(coop_covmat)::this
    COOP_INT::unit
    COOP_INT::n, i
    read(unit,*, ERR=100, END=100) n
    call this%alloc(n)
    read(unit,*, ERR=100, END=100)this%mean, this%sigma
    do i=1, n
       read(unit,*, ERR = 100, END=100) this%C(i+1:n, i)
       this%c(i,i) = 1.d0
    enddo
    do i=2, this%n
       this%c(1:i-1,i) = this%c(i, 1:i-1)
    enddo    
    call this%cholesky()
    call this%invert()    
    return
100 stop "covmat_read: failed"
  end subroutine coop_covmat_read
  

  subroutine coop_covmat_invert(this)
    class(coop_covmat)::this
    COOP_INT i,j    
    this%invC = this%L
    !!invert L
    call coop_cholesky_inv(this%n, this%n, this%invC)
    call coop_cholesky_sq(this%n, this%n, this%invC)
  end subroutine coop_covmat_invert

  subroutine coop_covmat_cholesky(this)
    class(coop_covmat)::this
    COOP_INT::i, j
    COOP_REAL::s
    this%L(1, 1) = 1.d0
    this%L(2:this%n, 1) = this%C(2:this%n, 1)
    do  i =2 , this%n
       s = 1.d0 - sum(this%L(i,1:i-1)**2)
       if(s .lt. -1.d-14)then
          stop "covmat_cholesky: matrix is not positive definite"
       endif
       this%L(i, i)=sqrt(max(s, 1.d-14)) 
       do j = i+1, this%n  
          this%L(j, i)=(this%C(j, i)-sum(this%L(j,1:i-1)*this%L(i,1:i-1)))/this%L(i,i)
       enddo
       this%L(1:i-1,i)=0.d0
    enddo
  end subroutine coop_covmat_cholesky

  function coop_covmat_chi2(this, x) result(chisq)
    class(coop_covmat)::this
    COOP_REAL::x(this%n), mu(this%n), chisq
    COOP_INT::i, j
    mu = (x - this%mean)/this%sigma
    chisq = dot_product(mu, matmul(this%invC, mu))/2.d0
  end function coop_covmat_chi2

  subroutine coop_covmat_normalize(this)
    class(coop_covmat)::this
    COOP_INT::i
    COOP_REAL::rescale(this%n)
    !$omp parallel do
    do i=1, this%n
       rescale(i) = sqrt(this%c(i, i))
    end do
    !$omp end parallel do    
    do i=1, this%n
       if(rescale(i).ne.1.d0)then
          this%c(i+1:this%n, i) = this%c(i+1:this%n, i)/rescale(i)
          this%c(i, 1:i-1) = this%c(i, 1:i-1)/rescale(i)
          this%c(i,i) = 1.d0
          this%sigma(i) = this%sigma(i)*rescale(i)
       endif
    enddo
    !!symmetrize
    do i=2, this%n
       this%c(1:i-1,i) = this%c(i, 1:i-1)
    enddo
    call this%cholesky()
    call this%invert()
  end subroutine coop_covmat_normalize

  subroutine coop_covmat_diagonal(this, sigma)
    class(coop_covmat)::this
    COOP_REAL::sigma(this%n)
    this%sigma = sigma
    call coop_set_identity_matrix(this%n, this%c)
    call coop_set_identity_matrix(this%n, this%invc)
    call coop_set_identity_matrix(this%n, this%L)    
  end subroutine coop_covmat_diagonal


  subroutine coop_covmat_MPI_sync(this, converge_R, weight_B)
    class(coop_covmat)::this
    COOP_REAL, optional::converge_R, weight_B !!convergence test
    COOP_REAL, dimension(:),allocatable::info, covinfo
    COOP_REAL,dimension(:,:),allocatable::cov, meanscov
    COOP_REAL:: R, bfac
    COOP_INT:: i, j, ii, ms, num_proc
    if(present(converge_R))converge_R = 0.d0    
    num_proc = coop_MPI_NumProc()
    if(num_proc .eq.1)goto 100
    ms = this%n*(this%n+1)/2
    allocate(info(0:this%n), covinfo(ms*2))
    info(0) = this%mult
    info(1:this%n) = this%mean*this%mult
    call coop_MPI_sum(info(0:this%n))
    if(info(0) .le. 0.d0)then
       return !!does not change anything for mult = 0
    endif
    info(1:this%n) = info(1:this%n)/info(0)
    ii = 0
    do i=1, this%n
       do j = 1, i
          ii = ii + 1
          covinfo(ii)  = this%c(i, j)*this%sigma(i)*this%sigma(j) *this%mult !!mean of covariance 
          covinfo(ms+ii) = (this%mean(i)-info(i))*(this%mean(j)-info(j))*this%mult  !!covariance of mean
       enddo
    enddo
    call coop_MPI_Sum(covinfo)
    !!normalize
    covinfo = covinfo/info(0)
    !!correction
    covinfo(ms+1:2*ms) = covinfo(ms+1:2*ms)*(dble(num_proc)/(num_proc-1.d0))
    
    this%mult = info(0)
    this%sigma = 1.d0
    this%mean = info(1:this%n)
    ii = 0
    allocate(cov(this%n, this%n), meanscov(this%n, this%n))
    do i=1, this%n
       do j = 1, i
          ii = ii + 1          
          cov(i, j) = covinfo(ii)
          meanscov(i, j) = covinfo(ii+ms)
          cov(j, i) = cov(i, j)
          meanscov(j, i) = meanscov(i, j)
       enddo
    enddo
    if(present(converge_R))then
       converge_R = coop_GelmanRubin_R(this%n, cov, meanscov)
    endif
    if(present(weight_B))then
       bfac = weight_B
    else
       if(present(converge_R))then
          bfac = converge_R/(converge_R + 3.d2)
       else
          bfac = 1.d0/info(0) !!default weight
       endif
    endif
    
    this%c = cov  !+ meanscov * bfac
    deallocate(info, covinfo, cov, meanscov)

100 call this%normalize()    
  end subroutine coop_covmat_MPI_sync


  function coop_matrix_dominant_eigen_value(n, a, accuracy) result(lambda)
    COOP_INT, parameter::max_loops = 150
    COOP_INT::n, iloop
    COOP_REAL::a(n, n)
    COOP_REAL,optional::accuracy
    COOP_REAL::eps
    COOP_REAL::x(n), y(n), last_lambda, lambda, sigma
    lambda = 0.d0
    sigma = maxval(abs(a))    
    if(sigma .eq. 0.d0)then
       return
    endif
    x = 1.d0/sqrt(dble(n))
    do iloop = 1, max(5, n/10)  !!do 5 iterations
       y = matmul(a, x)
       sigma = sqrt(dot_product(y, y))
       if(sigma .le. 1.d-30)then
          write(*,*) "Warning: dominant eigen too small"
          return
       endif
       x =  y/sigma
    enddo
    y = matmul(a, x)
    last_lambda = dot_product(x, y)
    if(abs(last_lambda) .lt. 1.d-30)then
       write(*,*) "Warning: dominant eigen too small"
       return       
    endif
    sigma = sqrt(dot_product(y, y))
    if(sigma .le. 1.d-30)then
       write(*,*) "Warning: dominant eigen too small"
       return
    endif
    x =  y/sigma
    iloop = 0
    if(present(accuracy))then
       eps = accuracy/2.d0
    else
       eps = 1.d-4
    endif
    do 
       y = matmul(a, x)
       lambda = dot_product(x, y)
       if(abs(lambda/last_lambda - 1.d0) .le. eps)return
       sigma = sqrt(dot_product(y, y))
       if(sigma .le. 1.d-30)then
          write(*,*) "Warning: dominant eigen too small"
          return
       endif
       x =  y/sigma
       iloop = iloop+1
       if(iloop .gt. max_loops)then
          write(*,*) "Warning: dominant eigen value does not converge"
          call coop_print_matrix(a, n, n)
          write(*,*) "last iteration: ", last_lambda, lambda
          return
       endif
       last_lambda = lambda
       if(abs(last_lambda) .lt. 1.d-30)then
          write(*,*) "Warning: dominant eigen too small"
          return       
       endif
    enddo
  end function coop_matrix_dominant_eigen_value



  function coop_GelmanRubin_R(n, cov, meanscov) result(R)
    COOP_INT, intent(in) :: n
    COOP_REAL,intent(in) :: cov(n,n), meanscov(n,n)
    COOP_REAL:: evals(n), R
    COOP_INT::i, info
    COOP_REAL:: rot(n, n), rotmeans(n, n)
    COOP_REAL::sc
    rot = cov
    rotmeans = meanscov
    do i=1,n
       sc = sqrt(cov(i,i))
       rot(i,:) = rot(i,:) / sc
       rot(:,i) = rot(:,i) / sc
       rotmeans(i,:) = rotmeans(i,:) /sc
       rotmeans(:,i) = rotmeans(:,i) /sc
    end do
    call coop_cholesky(n, n, rot, info)
    if(info .ne. 0) call coop_return_error("GelmanRubin_R", "matrix is not positive definite", "stop")
    call coop_cholesky_inv(n, n, rot)
    rotmeans =  matmul(matmul(rot, rotmeans), transpose(rot))
    R = coop_matrix_dominant_eigen_value(n, rotmeans, 1.d-3)
  end function coop_GelmanRubin_R



  !! for constrained problems
  !!suppose you have a covariance matrix C for known vector x_I and unknown vector x_J, you want to construct a mapping matrix Fmean (n_J, n_I) such that mean <x_J> = Fmean x_I and a matrix Ffluc (n_J, n_J) such that \delta x_J = Ffluc (Gaussian random vector y_J)
  !!after calling this subroutine C is destroyed
  subroutine coop_solve_constrained(m, n_known, n_unknown, dim_fmean, dim_ffluc, C, Fmean, Ffluc, epsilon)
    COOP_INT,intent(IN)::m, n_known, n_unknown, dim_fmean, dim_ffluc
    COOP_INT::n
    COOP_REAL,intent(INOUT)::C(m, n_known+n_unknown)
    COOP_REAL,intent(OUT)::Fmean(dim_fmean, n_known), FFluc(dim_ffluc, n_unknown)
    COOP_REAL::eig(n_known)
    COOP_REAL::epsilon, mineig
    COOP_INT::i, j
    !!do cholesky for  n_known x n_known
    n = n_unknown + n_known
    call coop_matsym_diagonalize(m, n_known, C, eig)
    mineig = maxval(abs(eig))*epsilon
    where (eig.lt.mineig)
       eig = 0.d0
    elsewhere
       eig = 1.d0/eig
    end where
    Fmean(1:n_unknown, 1:n_known) = matmul(C(n_known+1:n, 1:n_known), C(1:n_known, 1:n_known))
    do i = 1, n_known
       Fmean(1:n_unknown, i) = Fmean(1:n_unknown, i)*eig(i)
    enddo
    Fmean(1:n_unknown, 1:n_known) = matmul( Fmean(1:n_unknown, 1:n_known), transpose(C(1:n_known, 1:n_known)))

    Ffluc(1:n_unknown, 1:n_unknown) = C(n_known+1:n, n_known+1:n) - matmul(Fmean(1:n_unknown, 1:n_known), transpose(C(n_known+1:n, 1:n_known)))
    call coop_matsym_power(Ffluc, 0.5d0, mineig = epsilon)
  end subroutine coop_solve_constrained
  

end module Coop_matrix_mod


