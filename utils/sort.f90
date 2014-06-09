module coop_sort_mod
  use coop_wrapper_typedef
  implicit none
#include "constants.h"

  private


  public::coop_quicksort, coop_quicksortAcc, coop_quicksort_index

  !!coop_quicksort(array x) will sort x such that x(1)<=x(2)<...<=x(n)
  interface Coop_quicksort
     Module Procedure coop_quicksort_int,coop_quicksort_d, coop_quicksort_s
  End interface



  interface Coop_quicksortAcc
     Module Procedure coop_quicksort_intAcc,coop_quicksort_dAcc,Coop_quicksort_d_dAcc
  End interface

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
    COOP_REAL,dimension(:),intent(inout)::r
    COOP_INT,dimension(:,:),allocatable::stack  !!STACK
    COOP_REAL X
    COOP_INT n,l1,l2,lpoint,i,j
    n=size(r)
    allocate(stack(2,n))
    l1=1
    l2=N
    lpoint=0
    do
       do while(l1.lt.l2)
          i=l1
          j=l2
          X=r(l1)
          do while(i.lt.j)
             do while(i.lt.j.and.X.lE.r(j))
                j=j-1
             enddo
             if(i.lt.j)then
                r(i)=r(j)
                i=i+1
                do while(i.lt.j.and.r(i).lt.X)
                   i=i+1
                enddo
                if(i.lt.j)then
                   r(j)=r(i)
                   j=j-1
                endif
             endif
          enddo
          r(i)=X
          lpoint=lpoint+1
          stack(1,lpoint)=i+1
          stack(2,lpoint)=l2
          l2=i-1
       enddo
       if(lpoint.eq.0)EXiT
       l1=stack(1,lpoint)
       l2=stack(2,lpoint)
       lpoint=lpoint-1
    enddo
    deallocate(stack)
  end subroutine Coop_quicksort_d


  subroutine Coop_quicksort_s(r)
    real,dimension(:),intent(inout)::r
    COOP_INT,dimension(:,:),allocatable::stack  !!STACK
    real X
    COOP_INT n,l1,l2,lpoint,i,j
    n=size(r)
    allocate(stack(2,n))
    l1=1
    l2=N
    lpoint=0
    do
       do while(l1.lt.l2)
          i=l1
          j=l2
          X=r(l1)
          do while(i.lt.j)
             do while(i.lt.j.and.X.lE.r(j))
                j=j-1
             enddo
             if(i.lt.j)then
                r(i)=r(j)
                i=i+1
                do while(i.lt.j.and.r(i).lt.X)
                   i=i+1
                enddo
                if(i.lt.j)then
                   r(j)=r(i)
                   j=j-1
                endif
             endif
          enddo
          r(i)=X
          lpoint=lpoint+1
          stack(1,lpoint)=i+1
          stack(2,lpoint)=l2
          l2=i-1
       enddo
       if(lpoint.eq.0)EXiT
       l1=stack(1,lpoint)
       l2=stack(2,lpoint)
       lpoint=lpoint-1
    enddo
    deallocate(stack)
  end subroutine Coop_quicksort_s


  subroutine coop_quicksort_int(r)
    COOP_INT,dimension(:),intent(inOUT)::r
    COOP_INT,dimension(:,:),allocatable::stack  !!STACK
    COOP_INT X
    COOP_INT N,l1,l2,lpoint,i,j
    n=size(r)
    allocate(stack(2,n))
    l1=1
    l2=N
    lpoint=0
    do
       do while(l1.lt.l2)
          i=l1
          j=l2
          X=r(l1)
          do while(i.lt.j)
             do while(i.lt.j.and.X.lE.r(j))
                j=j-1
             enddo
             if(i.lt.j)then
                r(i)=r(j)
                i=i+1
                do while(i.lt.j.and.r(i).lt.X)
                   i=i+1
                enddo
                if(i.lt.j)then
                   r(j)=r(i)
                   j=j-1
                endif
             endif
          enddo
          r(i)=X
          lpoint=lpoint+1
          stack(1,lpoint)=i+1
          stack(2,lpoint)=l2
          l2=i-1
       enddo
       if(lpoint.eq.0)EXiT
       l1=stack(1,lpoint)
       l2=stack(2,lpoint)
       lpoint=lpoint-1
    enddo
    deallocate(stack)
  end subroutine coop_quicksort_int

  subroutine Coop_quicksort_dACC(r,sortindeX,ini)
    COOP_REAL,dimension(:),intent(inOUT)::r
    logical,optional::ini
    COOP_INT,dimension(:),intent(OUT)::sortindex
    COOP_INT,dimension(:,:),allocatable::stack  !!STACK
    COOP_REAL X
    COOP_INT N,l1,l2,lpoint,i,j
    COOP_INT iX
    n = Coop_getdim("Coop_quicksortAcc",SiZE(r),Size(Sortindex))
    allocate(stack(2,N))
    if(.not.PrESENT(ini))then
       do i=1,N
          sortindex(i)=i
       enddo
    ElSEif(ini)then
       do i=1,N
          sortindex(i)=i
       enddo
    endif
    l1=1
    l2=N
    lpoint=0
    do
       do while(l1.lt.l2)
          i=l1
          j=l2
          X=r(l1)
          iX=sortindex(l1)
          do while(i.lt.j)
             do while(i.lt.j.and.X.lE.r(j))
                j=j-1
             enddo
             if(i.lt.j)then
                r(i)=r(j)
                sortindex(i)=sortindex(j)
                i=i+1
                do while(i.lt.j.and.r(i).lt.X)
                   i=i+1
                enddo
                if(i.lt.j)then
                   r(j)=r(i)
                   sortindex(j)=sortindex(i)
                   j=j-1
                endif
             endif
          enddo
          r(i)=X
          sortindex(i)=iX
          lpoint=lpoint+1
          stack(1,lpoint)=i+1
          stack(2,lpoint)=l2
          l2=i-1
       enddo
       if(lpoint.eq.0)EXiT
       l1=stack(1,lpoint)
       l2=stack(2,lpoint)
       lpoint=lpoint-1
    enddo
    deallocate(stack)
  end subroutine Coop_quicksort_dACC

  subroutine Coop_quicksort_d_dACC(r,sortindex)
    COOP_REAL,dimension(:),intent(inOUT)::r
    COOP_REAL,dimension(:),intent(OUT)::sortindex
    COOP_INT,dimension(:,:),allocatable::stack  !!STACK
    COOP_REAL X
    COOP_INT N,l1,l2,lpoint,i,j
    COOP_REAL iX
    N=Coop_getdim("Coop_quicksortAcc",SiZE(r),Size(Sortindex))
    allocate(stack(2,N))
    l1=1
    l2=N
    lpoint=0
    do
       do while(l1.lt.l2)
          i=l1
          j=l2
          X=r(l1)
          iX=sortindex(l1)
          do while(i.lt.j)
             do while(i.lt.j.and.X.lE.r(j))
                j=j-1
             enddo
             if(i.lt.j)then
                r(i)=r(j)
                sortindex(i)=sortindex(j)
                i=i+1
                do while(i.lt.j.and.r(i).lt.X)
                   i=i+1
                enddo
                if(i.lt.j)then
                   r(j)=r(i)
                   sortindex(j)=sortindex(i)
                   j=j-1
                endif
             endif
          enddo
          r(i)=X
          sortindex(i)=iX
          lpoint=lpoint+1
          stack(1,lpoint)=i+1
          stack(2,lpoint)=l2
          l2=i-1
       enddo
       if(lpoint.eq.0)EXiT
       l1=stack(1,lpoint)
       l2=stack(2,lpoint)
       lpoint=lpoint-1
    enddo
    deallocate(stack)
  end subroutine Coop_quicksort_d_dACC


  subroutine coop_quicksort_intACC(r,sortindex,ini)
    COOP_INT,dimension(:),intent(inOUT)::r
    lOGiCAl,OPTiONAl::ini
    COOP_INT,dimension(:),intent(OUT)::sortindex
    COOP_INT,dimension(:,:),allocatable::stack  !!STACK
    COOP_INT X
    COOP_INT N,l1,l2,lpoint,i,j
    COOP_INT iX
    N=SiZE(r)
    allocate(stack(2,N))
    if(.not. PrESENT(ini))then
       do i=1,N
          sortindex(i)=i
       enddo
    ElSEif(ini)then
       do i=1,N
          sortindex(i)=i
       enddo
    endif
    l1=1
    l2=N
    lpoint=0
    do
       do while(l1.lt.l2)
          i=l1
          j=l2
          X=r(l1)
          iX=sortindex(l1)
          do while(i.lt.j)
             do while(i.lt.j.and.X.lE.r(j))
                j=j-1
             enddo
             if(i.lt.j)then
                r(i)=r(j)
                sortindex(i)=sortindex(j)
                i=i+1
                do while(i.lt.j.and.r(i).lt.X)
                   i=i+1
                enddo
                if(i.lt.j)then
                   r(j)=r(i)
                   sortindex(j)=sortindex(i)
                   j=j-1
                endif
             endif
          enddo
          r(i)=X
          sortindex(i)=iX
          lpoint=lpoint+1
          stack(1,lpoint)=i+1
          stack(2,lpoint)=l2
          l2=i-1
       enddo
       if(lpoint.eq.0)EXiT
       l1=stack(1,lpoint)
       l2=stack(2,lpoint)
       lpoint=lpoint-1
    enddo
    deallocate(stack)
  end subroutine coop_quicksort_intACC

  Subroutine Coop_quicksort_index(r, ind)
    COOP_REAL,dimension(:),intent(in)::r
    COOP_INT,dimension(:),intent(out)::ind
    COOP_INT,dimension(:,:),allocatable::stack  !!STACK
    COOP_INT ix
    COOP_INT n,l1,l2,lpoint,i,j
    n = Coop_getdim("coop_quicksort_index", size(r), size(ind))
    ind = (/ (i , i=1,n) /)
    allocate(stack(2, n))
    l1=1
    l2=N
    lpoint=0
    do
       do while(l1.lt.l2)
          i=l1
          j=l2
          ix = ind(l1)
          do while(i.lt.j)
             do while(i.lt.j.and. r(ix) .lE. r(ind(j)))
                j=j-1
             enddo
             if(i.lt.j)then
                ind(i)=ind(j)
                i=i+1
                do while(i.lt.j.and. r(ind(i)).lt. r(ix))
                   i=i+1
                enddo
                if(i.lt.j)then
                   ind(j) = ind(i)
                   j=j-1
                endif
             endif
          enddo
          ind(i) = ix
          lpoint=lpoint+1
          stack(1,lpoint)=i+1
          stack(2,lpoint)=l2
          l2=i-1
       enddo
       if(lpoint.eq.0)exit
       l1=stack(1,lpoint)
       l2=stack(2,lpoint)
       lpoint=lpoint-1
    enddo
    deallocate(stack)
  End Subroutine Coop_quicksort_index



end module coop_sort_mod
