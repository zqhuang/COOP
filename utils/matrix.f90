!Matrix utility routines. Uses BLAS/LAPACK. Mostly wrapper routines.

module coop_matrix_mod
  use coop_wrapper_typedef
  use coop_sort_mod
  implicit none

#include "constants.h"

  private

  public::coop_write_matrix, coop_print_matrix, coop_read_matrix, coop_set_identity_matrix, coop_identity_matrix, coop_diagonal_matrix, coop_matrix_add_diagonal, coop_matrix_sum_columns,  coop_matrix_solve_small,  Coop_matrix_Solve, Coop_matrix_Inverse, coop_matsym_xCx, coop_matrix_det_small,coop_matsym_mat2vec, coop_matsym_vec2mat, coop_matsym_inverse_small, Coop_matsym_Inverse, Coop_matsym_Diagonalize, coop_matsymdiag_small,  Coop_matsym_Sqrt,  coop_matsym_sqrt_small,  coop_matsym_power_small,  Coop_matsym_power,  coop_matsym_function, Coop_matsym_LnDet, Coop_matsym_Solve, coop_matsym_index

  COOP_INT,parameter::dl = kind(1.d0)
  COOP_INT,parameter::sp = kind(1.)

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
    do i=1,n
       a(i,i) = 1.d0
    enddo
  end subroutine coop_set_identity_matrix

  function coop_identity_matrix(n) result(a)
    COOP_INT i, n
    COOP_REAL  a(n,n)
    a = 0
    do i=1,n
       a(i,i) = 1.d0
    enddo
  end function coop_identity_matrix

  function coop_diagonal_matrix(diag) result(a)
    COOP_REAL ,dimension(:),intent(IN)::diag
    COOP_REAL  a(size(diag),size(diag))
    COOP_INT i
    do i=1,size(diag)
       a(i,i)=diag(i)
    enddo
  end function coop_diagonal_matrix

  subroutine coop_matrix_add_diagonal(n, mat, diag)
    COOP_INT n, i
    COOP_REAL  mat(n,n), diag(n)
    do i=1,n
       mat(i,i) = mat(i,i)+diag(i)
    enddo
  end subroutine coop_matrix_add_diagonal

  Subroutine Coop_write_matrix_d(funit, mat,nx,ny)
    COOP_INT funit
    COOP_INT,optional::nx,ny
    real(dl)  mat(:,:)
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
    real(dl)  mat(:,:)
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

  Subroutine Coop_read_matrix_d(funit,nx,ny,mat, success)
    logical,optional::success
    COOP_INT funit,nx,ny
    real(dl)  mat(nx,ny)
    COOP_INT i
    COOP_LONG_STRING line
    if(ny.gt. 1000)then
       write(*,*) "warning: reading a huge matrix: will not check comment lines"
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
    real(sp) mat(:,:)
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
    real(sp) mat(:,:)
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

  Subroutine Coop_read_matrix_s(funit,nx,ny,mat, success)
    logical,optional::success
    COOP_INT funit,nx,ny
    real(sp) mat(nx,ny)
    COOP_INT i
    COOP_LONG_STRING line
    if(ny.gt. 1000)then
       write(*,*) "warning: reading a huge matrix: will not check comment lines"
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
    do i=2, size(A)
       tr = tr + A(i,i)
    enddo
  end function coop_matrix_trace


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
       write(*,*) "Coop_matrix_Inverse: singular matrix cannot be inverted"
       stop
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

  subroutine coop_matsym_inverse_small(n,A)
    COOP_INT n
    COOP_REAL  A(n,n)
#ifndef HAS_LAPACK
    COOP_INT i, j
    COOP_REAL  sigma(n), p(n)
#endif
    select case(n)
    case(1)
       A=1.d0/A
    case(2)
       A=ReShape( (/ A(2,2), -A(2,1), -A(1,2), A(1,1) /) / (A(1,1)*A(2,2)-A(1,2)*A(2,1)), (/ 2, 2 /) )
    case(3)
       A = Reshape( (/ A(2,2)*A(3,3) - A(2,3)*A(3,2), A(2,3)*A(3,1)-A(2,1)*A(3,3), A(2,1)*A(3,2)- A(2,2)*A(3,1),  &
            A(3,2)*A(1,3)-A(1,2)*A(3,3), A(3,3)*A(1,1)-A(3,1)*A(1,3), A(1,2)*A(3,1)-A(1,1)*A(3,2), &
            A(1,2)*A(2,3)- A(2,2)*A(1,3), A(2,1)*A(1,3)-A(1,1)*A(2,3),  A(1,1)*A(2,2)-A(1,2)*A(2,1) /), (/ 3, 3 /) ) &
            / ( A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) + A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3)) + A(1,3)*(A(2,1)*A(3,2)- A(2,2)*A(3,1)) )
    case default
#ifdef HAS_LAPACK
       call coop_matrix_sympos_inverse(a)
#else
       do i=1, n
          sigma(i) = sqrt(A(i,i))
       enddo
       do i=1,n
          do j=1, n
             A(i,j) = A(i,j)/sigma(i)/sigma(j)
          enddo
       enddo
       call coop_matsym_choldc(n,a,p)
       call coop_matsym_cholinv(n,a,p)
       do i=1,n
          do j=i+1,n
             A(i,j)=dot_product(A(j:n,i),A(j:n,j))
          enddo
       enddo
       do i=1,n-1
          A(i,i)=dot_product(A(i:n,i),A(i:n,i))
          A(i+1:n,i)=A(i,i+1:n)
       enddo
       A(n,n)=A(n,n)**2
       do i=1,n
          do j=1, n
             A(i,j) = A(i,j)/sigma(i)/sigma(j)
          enddo
       enddo
#endif
    end select

  end subroutine coop_matsym_inverse_small


  subroutine Coop_matsym_Inverse(A)
    !! Invert a positive definite symetric matrix
    COOP_INT n
    COOP_REAL ,dimension(:,:),intent(inout):: A
#ifndef HAS_LAPACK
    COOP_REAL  p(size(A,1)), sigma(size(A,1))
    COOP_INT i,j
#endif
    n = Coop_getdim("Coop_matsym_Inverse", size(A,1), size(A,2))
    select case(n)
    case(1)
       A=1.d0/A
    case(2)
       A=ReShape( (/ A(2,2),-A(2,1), -A(1,2),A(1,1) /) / (A(1,1)*A(2,2)-A(1,2)*A(2,1)), (/ 2, 2 /) )
    case default
#ifdef HAS_LAPACK
       call coop_matrix_sympos_inverse(a)
#else
       do i=1, n
          sigma(i) = sqrt(A(i,i))
       enddo
       do i=1,n
          do j=1, n
             A(i,j) = A(i,j)/sigma(i)/sigma(j)
          enddo
       enddo
       call coop_matsym_choldc(n,a,p)
       call coop_matsym_cholinv(n,a,p)
       do i=1,n
          do j=i+1,n
             A(i,j)=dot_product(A(j:n,i),A(j:n,j))
          enddo
       enddo
       do i=1,n-1
          A(i,i)=dot_product(A(i:n,i),A(i:n,i))
          A(i+1:n,i)=A(i,i+1:n)
       enddo
       A(n,n)=A(n,n)**2
       do i=1,n
          do j=1, n
             A(i,j) = A(i,j)/sigma(i)/sigma(j)
          enddo
       enddo
#endif
    end select

  end subroutine Coop_matsym_Inverse


  !!find R such that R^T H R = diag(e), and R^T R = 1
  !!EIGEN VALUES ARE STORED IN E WITH E(1)<=E(2)...(from ground state..)
  !!EIGEN VECTORS ARE STORED IN COLUMNS OF H(:,1), H(:,2),...
  Subroutine Coop_matsym_Diagonalize(H, e, sort)
    COOP_REAL ,dimension(:,:)::h
    COOP_REAL ,dimension(:)::e
    logical,optional:: sort
    COOP_REAL  psi(size(e),size(e))
    COOP_INT Indx(size(E))
    COOP_INT i, n
    n = Coop_getdim("Coop_matsym_diagonalize",SIZE(H,1),SIZE(H,2),SIZE(E))
#ifdef HAS_LAPACK
    call coop_matrix_diagonalize(h, e)
    if(present(sort))then
       if(sort)then
          psi = H
          call coop_quicksortacc(e,indx)
          do i=1,n
             h(:,I)=psi(:,indx(I))
          enddo
       endif
    endif
#else
    call Coop_matsymDIAG(h,psi,1.d-8)
    do i = 1,N
       e(i) = h(i,i)
    enddo
    call Coop_quicksortacc(E,Indx)
    do i=1,n
       h(:,I)=psi(:,indx(I))
    enddo
#endif

  End Subroutine Coop_matsym_Diagonalize

  subroutine coop_matsymdiag_small(n, a, psi)
    COOP_INT n
    COOP_REAL  a(n,n), psi(n,n)
    COOP_REAL  s, d, theta
    select case(n)
    case(1)
       psi = 1.d0
       return
    case(2)
       if(a(1,2).eq.0.d0)then
          psi = reshape( (/ 1.d0, 0.d0, 0.d0, 1.d0 /), (/ 2, 2 /) )
          return
       endif
       d = sign(sqrt((a(1,1)-a(2,2))**2/4.d0+a(1,2)*a(2,1)), a(1,1)-a(2,2))
       s = (a(1,1) + a(2,2))/2.d0
       theta = atan2((a(1,1)-a(2,2))/2.d0 - d, a(1,2))
       a(1,1) = (s+d)
       a(2,1) = 0.d0
       a(1,2) = 0.d0
       a(2,2) = (s-d)
       psi(:,1) = (/ cos(theta), - sin(theta) /)
       psi(:,2) = (/ -psi(2,1), psi(1,1) /)
    case default
       call coop_matsymdiag(a, psi, 1.d-10)
    end select
  end subroutine coop_matsymdiag_small


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

  subroutine Coop_matsym_Sqrt(A)
    !!get the square root of a positive definite symmetric matrix A
    COOP_REAL ,dimension(:,:)::A
    COOP_REAL  E(size(A,1)), trans(size(A,1), size(A,1))
    COOP_INT n, i
    n= Coop_getdim("coop_matsym_sqrt", size(A,1),size(A,2))
#ifdef HAS_LAPACK
    call coop_matsym_diagonalize(a,e)
    if(any(E.lt. -1.d-30))then
       call Coop_return_error("Coop_matsym_Sqrt", "the matrix is not positive definite", "stop")
    endif
    e = sqrt(max(e,0.d0))
    do i=1,n
       trans(i,:) = A(:,i)*e(i)
    enddo
    A = matmul(A, trans)
#else
    call coop_svd_decompose(n,n,a,e,trans)
    do i=1,n
       trans(:,i) = trans(:,i)*sqrt(e(i))
    enddo
    A = matmul(A, transpose(trans))
#endif
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


  subroutine Coop_matsym_power(A,alpha)
    !!get A^alpha for a positive definite symmetric matrix A
    COOP_REAL ,dimension(:,:)::A
    COOP_REAL alpha
    COOP_REAL  E(size(A,1)), trans(size(A,1), size(A,1))
    COOP_INT n, i
    n= Coop_getdim("coop_matsym_sqrt", size(A,1),size(A,2))
#ifdef HAS_LAPACK
    call coop_matsym_diagonalize(a,e)
    if(any(E.lt. -1.d-30))then
       call Coop_return_error("Coop_matsym_Sqrt", "the matrix is not positive definite", "stop")
    endif
    e = (max(e,0.d0))**alpha
    do i=1,n
       trans(i,:) = A(:,i)*e(i)
    enddo
    A = matmul(A, trans)
#else
    call coop_svd_decompose(n,n,a,e,trans)
    if(any(E.lt. -1.d-30))then
       call Coop_return_error("Coop_matsym_Sqrt", "the matrix is not positive definite", "stop")
    endif
    e = (max(e,0.d0))**alpha
    do i=1,n
       trans(:,i) = trans(:,i)*e(i)
    enddo
    A = matmul(A, transpose(trans))
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
    call coop_matsym_diagonalize(a, e)
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
#else
       call Coop_matrix_Cholesky(n,acopy)
       call Coop_matrix_cholesky_solve(n,acopy,b)
#endif
    end select
  end subroutine coop_matsym_solve_small

  Subroutine Coop_matsym_Solve(A,b)
    COOP_REAL ,dimension(:,:),INTENT(INOUT)::A
    COOP_REAL ,dimension(:,:),INTENT(INOUT)::b
    COOP_INT n, m, i
    n=COOP_GETDIM("Coop_matsym_Solve",SIZE(A,1),SIZE(A,2),SIZE(b,1))
    m = size(b, 2)
#ifdef HAS_LAPACK
    call dposv('L', n, m, a, n, b, n, i)
#else
    call Coop_matrix_Cholesky(n,a)
    do i=1,m
       call Coop_matrix_cholesky_solve(n,a,b(:,i))
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
       write(*,*) "Warning: linear least square problem failed"
       write(*,*) "Error info = ", info
       x = 0
       return
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
       write(*,*) "Warning: linear least square problem failed"
       write(*,*) "Error info = ", info
       x = 0
       return
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

  !!calculate the inverse of a symmetric positive definite A
  subroutine coop_matrix_sympos_inverse(A)
    COOP_REAL ,dimension(:,:) :: A
    COOP_INT N, info, i, j
    n = coop_getdim("matrix_sympose_inverse", size(A,1), size(A,2))
    call DPOTRF("L", n, A, n, info)
    if(info.ne.0) call Coop_return_error("matrix_sympose_inverse", "the matrix is not positive definite", "stop")
    call DPOTRI("L", n, A, n, info)
    if(info.ne.0) call Coop_return_error("matrix_sympose_inverse", "cannot find the inverse of the matrix", "stop")
    do i=2,n
       do j=1,i-1
          A(j,i)= A(i,j)
       enddo
    enddo
  end subroutine coop_matrix_sympos_inverse



  !Does m = U diag U^T, returning U in real symmetric matrix M
  subroutine coop_Matrix_Diagonalize(M, diag)
    COOP_REAL , intent(inout):: m(:,:)
    COOP_REAL , intent(out) :: diag(:)
    COOP_INT n
    COOP_INT ierr, tmpsize
    COOP_REAL , allocatable, dimension(:) :: tmp
    n = coop_getdim("matrix_diagonalize", size(M, 1), size(M,2), size(diag))
    tmpsize =  3*n -1
    allocate(tmp(tmpsize));
    call DSYEV('V','U',n,m,n,diag,tmp,tmpsize,ierr) !evalues and vectors of symmetric matrix
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

  Subroutine coop_matsym_choldc(n,a,p) !!a=LL^T,
    COOP_INT n
    COOP_REAL  a(n,n)
    COOP_REAL  p(n)
    COOP_INT i
    COOP_REAL  summ
    do i=1,n
       summ=a(i,i)-dot_product(a(i,1:i-1),a(i,1:i-1))
       if (summ <= 0.0)then
          write(*,*) 'coop_matsym_choldc failed;', i,  summ
          write(*,*) 'try using LAPACK library'
          stop 
       endif
       p(i)=sqrt(summ)
       a(i+1:n,i)=(a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
    end do
  end Subroutine coop_matsym_choldc

  Subroutine coop_matsym_cholsl(n,a,p,b,x) !!solve LL^T x = b
    COOP_INT n
    COOP_REAL  a(n,n)
    COOP_REAL  p(n),b(n)
    COOP_REAL  x(n)
    COOP_INT i
    do i=1,n
       x(i)=(b(i)-dot_product(a(i,1:i-1),x(1:i-1)))/p(i) !!solve L y = b
    end do
    do i=n,1,-1
       x(i)=(x(i)-dot_product(a(i+1:n,i),x(i+1:n)))/p(i) !! solve L^T x = y
    end do
  end Subroutine coop_matsym_cholsl

  subroutine coop_matsym_cholinv(n,a,p)  !!invert L
    COOP_INT n
    COOP_REAL  a(n,n),p(n)
    COOP_INT i,j
    do  i=1,n
       a(i,i)=1./p(i)
       do  j=i+1,n
          a(j,i)=-dot_product(a(j,i:j-1),a(i:j-1,i))/p(j)
       enddo
    enddo
  end subroutine coop_matsym_cholinv

  !!********************************************************
  !!matrix diagonalization; diag(eig1; eig2; ..; eign ) = R^T A R
  !!return R;
  !!return A=diag(eig1, eig2, ..., eign)
  Subroutine Coop_matsymDIAG(A,R,PRECISION)
    COOP_INT,parameter::MAXLOOP=200000
    COOP_REAL ,dimension(:,:),INTENT(INOUT)::A
    COOP_REAL ,dimension(:,:),INTENT(INOUT)::R
    COOP_REAL  PRECISION,PREC
    COOP_INT N,I
    Logical Tag
    N=COOP_GETDIM("Coop_matsymDIAG",SIZE(A,1),SIZE(A,2),SIZE(R,1),SIZE(R,2))
    PREC=PRECISION*sum(abs(A))/N/N
    R=0.d0
    do i=1,n
       R(i,i)=1.d0
    enddo
    I=0
    Tag=.true.

    do while(I.lt.MAXLOOP.and.TAG)
       call Coop_matsymDIAGSTEP_D(A,R,N,PREC,tag)
       I=I+1
    enddo
  end Subroutine Coop_matsymDIAG

  Subroutine Coop_matsymROTANGLE_D(A,M1,N1,tant,cost,sint)
    COOP_REAL ,dimension(:,:),INTENT(IN)::A
    COOP_REAL ,parameter::EPS=1.D-6
    COOP_REAL ,parameter::LARGEPAR=1.d0/EPS
    COOP_INT M1,N1
    COOP_REAL  tant,cost,sint,S,S2
    s = (A(N1,N1)-A(M1,M1))*0.5d0/A(M1,N1)
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
  end Subroutine Coop_matsymROTANGLE_D

  Subroutine Coop_matsymFINDMAX_D(A,N,M1,N1,PRECISION,TAG)
    COOP_INT M1,N1,I,N,J
    COOP_REAL  A(N,N)
    COOP_REAL  PRECISION,TEMP,Y
    logical TAG
    temp = precision
    do I=2,N
       do J=1,I-1
          Y=abs(A(I,J))
          if(Y.gt.TEMP)then
             TEMP=Y
             M1=I
             N1=J
          endif
       enddo
    enddo
    tag = (temp .gt. precision)
  end  Subroutine Coop_matsymFINDMAX_D

  Subroutine Coop_matsymROT_D(R,N,M1,N1,cost,sint)
    COOP_INT N,M1,N1,I
    COOP_REAL  cost,sint
    COOP_REAL  R(N,N),RM1(N)
    do I=1,N
       RM1(i)=R(i,M1)*cost-R(i,N1)*sint
       R(i,N1)=R(i,N1)*cost+R(i,M1)*sint
    enddo
    R(:,M1)=RM1
  end Subroutine Coop_matsymROT_D

  Subroutine Coop_matsymDIAGSTEP_D(A,R,N,PRECISION,TAG)
    COOP_INT N,N1,M1  !!N1<M1
    COOP_REAL  A(N,N),R(N,N),PRECISION,cost,sint,tant,TEMP
    logical TAG
    call Coop_matsymFINDMAX_D(A,N,M1,N1,PRECISION,TAG)
    if(TAG)then
       call Coop_matsymROTANGLE_D(A,M1,N1,tant,cost,sint)
       call Coop_matsymROT_D(R,N,M1,N1,cost,sint)
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
  end Subroutine Coop_matsymDIAGSTEP_D

  !!%%%%%%%%%%%%%% CHOLESKY DECOMPOSITION %%%%%%%%%%%%%%%%%%%%%%%
  !! decompose positive definit real symetric matrix A=LL^T,
  !! where elements of L satisfy L(i,j)=0 if j>i
  !! L is saved in A, A is destroyed
  Subroutine Coop_matrix_cholesky(n, a)
    COOP_INT n, i, j
    COOP_REAL  a(n,n)
    a(1,1)=sqrt(a(1,1))
    a(2:n,1)=A(2:n,1)/A(1,1)
    do I=2,N
       A(I,I)=sqrt(A(I,I)-sum(A(I,1:I-1)**2)) !!L(I,I) doNE
       do J=I+1,N  !! NOW CALCULATE L(J,I), J>=I
          A(J,I)=(A(J,I)-sum(A(J,1:I-1)*A(I,1:I-1)))/A(I,I)
       enddo
       A(1:I-1,I)=0.d0
    enddo
  end Subroutine Coop_matrix_cholesky

  Subroutine Coop_matrix_cholesky_solve(n,a,b)
    COOP_INT N,i
    COOP_REAL  A(N,N),b(N)
    b(1)=b(1)/A(1,1)
    do I=2,N
       b(i)=(b(i)-sum(A(i,1:i-1)*b(1:i-1)))/A(i,i)
    enddo
    b(N)=b(N)/A(N,N)
    do i=n-1,1,-1
       b(i)=(b(i)-sum(b(i+1:n)*A(i+1:n,i)))/A(i,i)
    enddo
  end Subroutine Coop_matrix_cholesky_solve

  function coop_matsym_index(n, i, j) result(ind)
    COOP_INT n, i, j
    COOP_INT ind
    ind = (2*n+1 - abs(j-i))*abs(j-i)/2 + min(i, j)
  end function coop_matsym_index

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

end module Coop_matrix_mod
