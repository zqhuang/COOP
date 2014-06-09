module coop_sort_mod
  use coop_wrapper_typedef
  implicit none
#include "constants.h"

  private

  COOP_INT ,parameter::sp = kind(1.)
  COOP_INT ,parameter::dl = kind(1.d0)

  public::coop_quicksort, coop_quicksortAcc, coop_quicksort_index

  !!coop_quicksort(array x) will sort x such that x(1)<=x(2)<...<=x(n)
  interface Coop_quicksort
     Module Procedure coop_quicksort_int,coop_quicksort_d, coop_quicksort_s,  coop_quicksort2d_int,coop_quicksort2d_d, coop_quicksort2d_s
  End interface



  interface Coop_quicksortAcc
     Module Procedure coop_quicksortacc_d, coop_quicksortAcc_s, Coop_quicksortacc_dd, coop_quicksortacc_ss
  End interface

  interface coop_quicksort_index
     module procedure coop_quicksort_index_d, coop_quicksort_index_s
  end interface coop_quicksort_index

contains

  function left_location_of(n, arr, x) result(l)
    COOP_INT,intent(in):: n
    COOP_REAL,intent(in):: arr(n), x
    COOP_INT::l
    COOP_INT r, m
    l = 1
    r = n
    if(arr(1).gt.arr(n))then
       do while(r - l .gt. 1)
          m = (l+r)/2
          if(arr(m) .gt. x)then
             l  = m
          else
             r = m
          endif
       enddo
    else
       do while(r-l.gt.1)
          m = (l+r)/2
          if(arr(m) .le. x)then
             l = m
          else
             r = m
          endif
       enddo
    endif
  end function left_location_of

  function right_location_of(n, arr, x) result(r)
    COOP_INT,intent(in):: n
    COOP_REAL,intent(in):: arr(n), x
    COOP_INT::r
    COOP_INT l, m
    l = 1
    r = n
    if(arr(1).gt.arr(n))then
       do while(r - l .gt. 1)
          m = (l+r)/2
          if(arr(m) .gt. x)then
             l  = m
          else
             r = m
          endif
       enddo
    else
       do while(r-l.gt.1)
          m = (l+r)/2
          if(arr(m) .le. x)then
             l = m
          else
             r = m
          endif
       enddo
    endif
  end function right_location_of


  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  !!%%%%%%%%%%%%%%% Coop_quicksort %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !!non-recursive QUiCK SOrT r-->r(1)<r(2)<...
  subroutine Coop_quicksort_d(r)
    real(dl),dimension(:),intent(inout)::r
    COOP_INT  n
    n=size(r)
    call quicksort_double(r, n)
  end subroutine Coop_quicksort_d


  subroutine Coop_quicksort2d_d(r)
    real(dl),dimension(:,:),intent(inout)::r
    COOP_INT  n
    n=size(r)
    call quicksort_double(r, n)
  end subroutine Coop_quicksort2d_d



  subroutine Coop_quicksort_s(r)
    real(sp),dimension(:),intent(inout)::r
    COOP_INT  n
    n = size(r)
    call quicksort_float(r, n)
  end subroutine Coop_quicksort_s

  subroutine Coop_quicksort2d_s(r)
    real(sp),dimension(:,:),intent(inout)::r
    COOP_INT  n
    n=size(r)
    call quicksort_float(r, n)
  end subroutine Coop_quicksort2d_s


  subroutine coop_quicksort_int(r)
    COOP_INT ,dimension(:),intent(inout)::r
    COOP_INT  n
    n=size(r)
    call quicksort_int(r, n)
  end subroutine coop_quicksort_int

  subroutine coop_quicksort2d_int(r)
    COOP_INT ,dimension(:,:),intent(inout)::r
    COOP_INT  n
    n=size(r)
    call quicksort_int(r, n)
  end subroutine coop_quicksort2d_int



  subroutine Coop_quicksortacc_d(r, indices)
    real(dl),dimension(:),intent(inOUT)::r
    COOP_INT,dimension(:),intent(OUT)::indices
    COOP_INT n, i
    n = Coop_getdim("Coop_quicksortAcc",SiZE(r),Size(indices))
    indices = (/ (i, i=1, n) /)
    call quicksort_double_with_indices(r, indices, n)
  end subroutine Coop_quicksortacc_d

  subroutine Coop_quicksortacc_s(r, indices)
    real(sp),dimension(:),intent(inOUT)::r
    COOP_INT,dimension(:),intent(OUT)::indices
    COOP_INT n, i
    n = Coop_getdim("Coop_quicksortAcc",SiZE(r),Size(indices))
    indices = (/ (i, i=1, n) /)
    call quicksort_float_with_indices(r, indices, n)
  end subroutine Coop_quicksortacc_s

  subroutine Coop_quicksortacc_dd(r, indices)
    real(dl),dimension(:),intent(inOUT)::r
    real(dl),dimension(:),intent(OUT)::indices
    integer n
    n =  Coop_getdim("Coop_quicksortAcc",SiZE(r),Size(indices))
    call quicksort_double_with_double(r, indices, n)
  end subroutine Coop_quicksortacc_dd


  subroutine Coop_quicksortacc_ss(r, indices)
    real(sp),dimension(:),intent(inOUT)::r
    real(sp),dimension(:),intent(OUT)::indices
    integer n
    n =  Coop_getdim("Coop_quicksortAcc",SiZE(r),Size(indices))
    call quicksort_float_with_float(r, indices, n)
  end subroutine Coop_quicksortacc_ss


  Subroutine Coop_quicksort_index_d(r, ind)
    real(dl),dimension(:),intent(in)::r
    COOP_INT,dimension(:),intent(out)::ind
    COOP_INT n, i
    real(dl),dimension(:),allocatable::x
    n = Coop_getdim("coop_quicksort_index", size(r), size(ind))
    allocate(x(n), source=r)
    ind = (/ (i , i=1,n) /)
    call coop_quicksortacc(x, ind)
    deallocate(x)
  End Subroutine Coop_quicksort_index_d

  Subroutine Coop_quicksort_index_s(r, ind)
    real(sp),dimension(:),intent(in)::r
    COOP_INT,dimension(:),intent(out)::ind
    COOP_INT n, i
    real(sp),dimension(:),allocatable::x
    n = Coop_getdim("coop_quicksort_index", size(r), size(ind))
    allocate(x(n), source=r)
    ind = (/ (i , i=1,n) /)
    call coop_quicksortacc(x, ind)
    deallocate(x)
  End Subroutine Coop_quicksort_index_s



end module coop_sort_mod
