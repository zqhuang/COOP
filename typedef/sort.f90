module coop_sort_mod
  use coop_basicutils_mod
  implicit none
#include "constants.h"

  private

  public::coop_quicksort, coop_quicksortAcc, coop_quicksort_index, coop_bin_data, coop_array_get_threshold, coop_get_bounds

  !!coop_quicksort(array x) will sort x such that x(1)<=x(2)<...<=x(n)
  interface Coop_quicksort
     Module Procedure coop_quicksort_int,coop_quicksort_d, coop_quicksort_s,  coop_quicksort2d_int,coop_quicksort2d_d, coop_quicksort2d_s
  End interface

    interface coop_get_bounds
     module procedure coop_get_bounds_d, coop_get_bounds_s
  end interface coop_get_bounds


  interface Coop_quicksortAcc
     Module Procedure coop_quicksortacc_d, coop_quicksortAcc_s, Coop_quicksortacc_dd, coop_quicksortacc_ss
  End interface

  interface coop_quicksort_index
     module procedure coop_quicksort_index_d, coop_quicksort_index_s
  end interface coop_quicksort_index

  interface coop_bin_data
     module procedure coop_bin_data_d, coop_bin_data_s, coop_bin_data2d_d, coop_bin_data2d_s
  end interface coop_bin_data

  interface coop_array_get_threshold
     module procedure coop_array_get_threshold_d, coop_array_get_threshold_s, coop_array_get_threshold2d_d, coop_array_get_threshold2d_s
  end interface coop_array_get_threshold

contains

  function coop_left_location_of(n, arr, x) result(l)
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
  end function coop_left_location_of

  function coop_right_location_of(n, arr, x) result(r)
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
  end function coop_right_location_of


  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  !!%%%%%%%%%%%%%%% Coop_quicksort %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !!non-recursive QUiCK SOrT r-->r(1)<r(2)<...
  subroutine Coop_quicksort_d(r)
    COOP_REAL ,dimension(:),intent(inout)::r
    COOP_INT  n
    n=size(r)
    call quicksort_double(r, n)
  end subroutine Coop_quicksort_d


  subroutine Coop_quicksort2d_d(r)
    COOP_REAL ,dimension(:,:),intent(inout)::r
    COOP_INT  n
    n=size(r)
    call quicksort_double(r, n)
  end subroutine Coop_quicksort2d_d



  subroutine Coop_quicksort_s(r)
    COOP_SINGLE,dimension(:),intent(inout)::r
    COOP_INT  n
    n = size(r)
    call quicksort_float(r, n)
  end subroutine Coop_quicksort_s

  subroutine Coop_quicksort2d_s(r)
    COOP_SINGLE,dimension(:,:),intent(inout)::r
    COOP_INT  n
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
    COOP_REAL ,dimension(:),intent(inOUT)::r
    COOP_INT,dimension(:),intent(OUT)::indices
    COOP_INT n, i
    n = Coop_getdim("Coop_quicksortAcc",SiZE(r),Size(indices))
    indices = (/ (i, i=1, n) /)
    call quicksort_double_with_indices(r, indices, n)
  end subroutine Coop_quicksortacc_d

  subroutine Coop_quicksortacc_s(r, indices)
    COOP_SINGLE,dimension(:),intent(inOUT)::r
    COOP_INT,dimension(:),intent(OUT)::indices
    COOP_INT n, i
    n = Coop_getdim("Coop_quicksortAcc",SiZE(r),Size(indices))
    indices = (/ (i, i=1, n) /)
    call quicksort_float_with_indices(r, indices, n)
  end subroutine Coop_quicksortacc_s

  subroutine Coop_quicksortacc_dd(r, indices)
    COOP_REAL ,dimension(:),intent(inOUT)::r
    COOP_REAL ,dimension(:),intent(OUT)::indices
    COOP_INT  n
    n =  Coop_getdim("Coop_quicksortAcc",SiZE(r),Size(indices))
    call quicksort_double_with_double(r, indices, n)
  end subroutine Coop_quicksortacc_dd


  subroutine Coop_quicksortacc_ss(r, indices)
    COOP_SINGLE,dimension(:),intent(inOUT)::r
    COOP_SINGLE,dimension(:),intent(OUT)::indices
    COOP_INT  n
    n =  Coop_getdim("Coop_quicksortAcc",SiZE(r),Size(indices))
    call quicksort_float_with_float(r, indices, n)
  end subroutine Coop_quicksortacc_ss


  Subroutine Coop_quicksort_index_d(r, ind)
    COOP_REAL ,dimension(:),intent(in)::r
    COOP_INT,dimension(:),intent(out)::ind
    COOP_INT n, i
    COOP_REAL ,dimension(:),allocatable::x
    n = Coop_getdim("coop_quicksort_index", size(r), size(ind))
    allocate(x(n), source=r)
    ind = (/ (i , i=1,n) /)
    call coop_quicksortacc(x, ind)
    deallocate(x)
  End Subroutine Coop_quicksort_index_d

  Subroutine Coop_quicksort_index_s(r, ind)
    COOP_SINGLE,dimension(:),intent(in)::r
    COOP_INT,dimension(:),intent(out)::ind
    COOP_INT n, i
    COOP_SINGLE,dimension(:),allocatable::x
    n = Coop_getdim("coop_quicksort_index", size(r), size(ind))
    allocate(x(n), source=r)
    ind = (/ (i , i=1,n) /)
    call coop_quicksortacc(x, ind)
    deallocate(x)
  End Subroutine Coop_quicksort_index_s

  subroutine coop_bin_data_d(x, nbins, center, density)
    COOP_REAL  x(:), center(:), density(:)
    COOP_INT  nbins, n, nb
    n = size(x)
    nb = min(size(center), size(density))
    call get_binned_data_double(x, n, nb, nbins, center, density)
  end subroutine coop_bin_data_d

  subroutine coop_bin_data_s(x, nbins, center, density)
    COOP_SINGLE x(:), center(:), density(:)
    COOP_INT  nbins, n, nb
    n = size(x)
    nb = min(size(center), size(density))
    call get_binned_data_float(x, n, nb, nbins, center, density)
  end subroutine coop_bin_data_s
  
  subroutine coop_bin_data2d_d(x, nbins, center, density)
    COOP_REAL  x(:,:), center(:), density(:)
    COOP_INT  nbins, n, nb
    n = size(x)
    nb = min(size(center), size(density))
    call get_binned_data_double(x, n, nb, nbins, center, density)
  end subroutine coop_bin_data2d_d

  subroutine coop_bin_data2d_s(x, nbins, center, density)
    COOP_SINGLE x(:,:), center(:), density(:)
    COOP_INT  nbins, n, nb
    n = size(x)
    nb = min(size(center), size(density))
    call get_binned_data_float(x, n, nb, nbins, center, density)
  end subroutine coop_bin_data2d_s


  subroutine coop_array_get_threshold_d(x, perc, threshold)
    COOP_REAL  x(:)
    COOP_REAL  perc,threshold
    COOP_INT n
    n = size(x)
    call array_get_threshold_double(x, n, perc, threshold)
  end subroutine coop_array_get_threshold_d


  subroutine coop_array_get_threshold_s(x, perc, threshold)
    COOP_SINGLE x(:)
    COOP_SINGLE perc,threshold
    COOP_INT n
    n = size(x)
    call array_get_threshold_float(x, n, perc, threshold)
  end subroutine coop_array_get_threshold_s

  subroutine coop_array_get_threshold2d_d(x, perc, threshold)
    COOP_REAL  x(:,:)
    COOP_REAL  perc,threshold
    COOP_INT n
    n = size(x)
    call array_get_threshold_double(x, n, perc, threshold)
  end subroutine coop_array_get_threshold2d_d


  subroutine coop_array_get_threshold2d_s(x, perc, threshold)
    COOP_SINGLE x(:,:)
    COOP_SINGLE perc,threshold
    COOP_INT n
    n = size(x)
    call array_get_threshold_float(x, n, perc, threshold)
  end subroutine coop_array_get_threshold2d_s

  subroutine coop_get_bounds_d(x, cls, bounds)
    COOP_REAL x(:), cls(:), bounds(:)
    COOP_REAL,dimension(:),allocatable::xcopy
    COOP_INT n, i, j, m
    COOP_REAL ir
    n = size(x)
    m = coop_getdim("get_bounds", size(cls), size(bounds))
    allocate(xcopy(n))
    xcopy = x
    call coop_quicksort(xcopy)
    do j=1, m
       ir = cls(j)*(n-1.d0) + 1.d0
       i = max(1, min(floor(ir), n-1))
       ir = ir - i
       bounds(j) = xcopy(i) * (1.d0 - ir) + xcopy(i+1)*ir
    enddo
    deallocate(xcopy)
  end subroutine coop_get_bounds_d


subroutine coop_get_bounds_s(x, cls, bounds)
    COOP_SINGLE x(:), cls(:), bounds(:)
    COOP_SINGLE,dimension(:),allocatable::xcopy
    COOP_INT n, i, j, m
    COOP_REAL ir
    n = size(x)
    m = coop_getdim("get_bounds", size(cls), size(bounds))
    allocate(xcopy(n))
    xcopy = x
    call coop_quicksort(xcopy)
    do j=1, m
       ir = cls(j)*(n-1.d0) + 1.d0
       i = max(1, min(floor(ir), n-1))
       ir = ir - i
       bounds(j) = xcopy(i) * (1.d0 - ir) + xcopy(i+1)*ir
    enddo
    deallocate(xcopy)
  end subroutine coop_get_bounds_s
  


end module coop_sort_mod
