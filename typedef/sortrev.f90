module coop_sortrev_mod
  use coop_basicutils_mod
  implicit none
#include "constants.h"

  private

  COOP_INT ,parameter::sp = kind(1.)
  COOP_INT ,parameter::dl = kind(1.d0)

  public::coop_quicksortrev, coop_quicksortrevAcc, coop_quicksortrev_index

  !!coop_quicksortrev(array x) will sortrev x such that x(1)<=x(2)<...<=x(n)
  interface Coop_quicksortrev
     Module Procedure coop_quicksortrev_int,coop_quicksortrev_d, coop_quicksortrev_s,  coop_quicksortrev2d_int,coop_quicksortrev2d_d, coop_quicksortrev2d_s
  End interface



  interface Coop_quicksortrevAcc
     Module Procedure coop_quicksortrevacc_d, coop_quicksortrevAcc_s, Coop_quicksortrevacc_dd, coop_quicksortrevacc_ss
  End interface

  interface coop_quicksortrev_index
     module procedure coop_quicksortrev_index_d, coop_quicksortrev_index_s
  end interface coop_quicksortrev_index


contains

  !!%%%%%%%%%%%%%%% Coop_quicksortrev %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !!non-recursive QUiCK Sortrev r-->r(1)<r(2)<...
  subroutine Coop_quicksortrev_d(r)
    real(dl),dimension(:),intent(inout)::r
    COOP_INT  n
    n=size(r)
    call quicksortrev_double(r, n)
  end subroutine Coop_quicksortrev_d


  subroutine Coop_quicksortrev2d_d(r)
    real(dl),dimension(:,:),intent(inout)::r
    COOP_INT  n
    n=size(r)
    call quicksortrev_double(r, n)
  end subroutine Coop_quicksortrev2d_d



  subroutine Coop_quicksortrev_s(r)
    real(sp),dimension(:),intent(inout)::r
    COOP_INT  n
    n = size(r)
    call quicksortrev_float(r, n)
  end subroutine Coop_quicksortrev_s

  subroutine Coop_quicksortrev2d_s(r)
    real(sp),dimension(:,:),intent(inout)::r
    COOP_INT  n
    call quicksortrev_float(r, n)
  end subroutine Coop_quicksortrev2d_s


  subroutine coop_quicksortrev_int(r)
    COOP_INT ,dimension(:),intent(inout)::r
    COOP_INT  n
    n=size(r)
    call quicksortrev_int(r, n)
  end subroutine coop_quicksortrev_int

  subroutine coop_quicksortrev2d_int(r)
    COOP_INT ,dimension(:,:),intent(inout)::r
    COOP_INT  n
    n=size(r)
    call quicksortrev_int(r, n)
  end subroutine coop_quicksortrev2d_int



  subroutine Coop_quicksortrevacc_d(r, indices)
    real(dl),dimension(:),intent(inOUT)::r
    COOP_INT,dimension(:),intent(OUT)::indices
    COOP_INT n, i
    n = Coop_getdim("Coop_quicksortrevAcc",SiZE(r),Size(indices))
    indices = (/ (i, i=1, n) /)
    call quicksortrev_double_with_indices(r, indices, n)
  end subroutine Coop_quicksortrevacc_d

  subroutine Coop_quicksortrevacc_s(r, indices)
    real(sp),dimension(:),intent(inOUT)::r
    COOP_INT,dimension(:),intent(OUT)::indices
    COOP_INT n, i
    n = Coop_getdim("Coop_quicksortrevAcc",SiZE(r),Size(indices))
    indices = (/ (i, i=1, n) /)
    call quicksortrev_float_with_indices(r, indices, n)
  end subroutine Coop_quicksortrevacc_s

  subroutine Coop_quicksortrevacc_dd(r, indices)
    real(dl),dimension(:),intent(inOUT)::r
    real(dl),dimension(:),intent(OUT)::indices
    COOP_INT  n
    n =  Coop_getdim("Coop_quicksortrevAcc",SiZE(r),Size(indices))
    call quicksortrev_double_with_double(r, indices, n)
  end subroutine Coop_quicksortrevacc_dd


  subroutine Coop_quicksortrevacc_ss(r, indices)
    real(sp),dimension(:),intent(inOUT)::r
    real(sp),dimension(:),intent(OUT)::indices
    COOP_INT  n
    n =  Coop_getdim("Coop_quicksortrevAcc",SiZE(r),Size(indices))
    call quicksortrev_float_with_float(r, indices, n)
  end subroutine Coop_quicksortrevacc_ss


  Subroutine Coop_quicksortrev_index_d(r, ind)
    real(dl),dimension(:),intent(in)::r
    COOP_INT,dimension(:),intent(out)::ind
    COOP_INT n, i
    real(dl),dimension(:),allocatable::x
    n = Coop_getdim("coop_quicksortrev_index", size(r), size(ind))
    allocate(x(n))
    x = r
    ind = (/ (i , i=1,n) /)
    call coop_quicksortrevacc(x, ind)
    deallocate(x)
  End Subroutine Coop_quicksortrev_index_d

  Subroutine Coop_quicksortrev_index_s(r, ind)
    real(sp),dimension(:),intent(in)::r
    COOP_INT,dimension(:),intent(out)::ind
    COOP_INT n, i
    real(sp),dimension(:),allocatable::x
    n = Coop_getdim("coop_quicksortrev_index", size(r), size(ind))
    allocate(x(n))
    x = r
    ind = (/ (i , i=1,n) /)
    call coop_quicksortrevacc(x, ind)
    deallocate(x)
  End Subroutine Coop_quicksortrev_index_s

end module coop_sortrev_mod
