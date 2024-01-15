module tabfunc_utils
  use basic_utils
  use general_utils
  use ode_utils
  implicit none

#ifdef HAS_FFTW
  interface tabfunc_spherical_fft
     module procedure tabfunc_spherical_fft_from_function, tabfunc_spherical_fft_from_id
  end interface tabfunc_spherical_fft
#endif

  interface tabfunc_eval
     module procedure tabfunc_eval_s, tabfunc_eval_v, tabfunc_eval_linearlog
  end interface tabfunc_eval

  interface tabfunc_derv
     module procedure tabfunc_derv_s, tabfunc_derv_v
  end interface tabfunc_derv

  interface tabfunc_int
     module procedure tabfunc_int_s, tabfunc_int_v
  end interface tabfunc_int

  
#define IND_L  tabfunc_space(id)%il
#define IND_R  tabfunc_space(id)%ir
#define IND_MAX tabfunc_space(id)%maxind
#define TAB  tabfunc_space(id)%table
#define BOUND_L tabfunc_space(id)%left
#define BOUND_R tabfunc_space(id)%right
#define DELTAX tabfunc_space(id)%dx
#define PPV tabfunc_space(id)%pp
#define CHEB_WEIGHT(i,n) ((2.*(i-1.)/(n-1.)-1.)**8+1.)

  Type tabularized_function
     real(dl), allocatable:: table(:) !!table of  f(x_i), i= 1,2,..
     real(dl), allocatable:: pp(:) !! table of -f''(x_i+dx/2)/2._dl dx^2, i= 1,2..
     real(dl) left, right, dx 
     integer maxind, il, ir
  end type tabularized_function

  Integer::TabFunc_Current_Id = 1
  Integer,Parameter::Num_TabFunc = 8192

  Integer,parameter::use_cheb_fit_limit = 30 !!when maxind is less than or equal to this number, use Chebyshev Polynomial to fit the function
  logical,dimension(Num_TabFunc):: TabFunc_used = .false.
  Type(tabularized_function),dimension(Num_TabFunc):: tabfunc_space

contains

  subroutine tabfunc_tab(id, f, xmin, xmax, tsize, sample_xmin, sample_xmax)
    external f
    real(dl) f
    real(dl) xmin, xmax, x, x1, xn
    real(dl),optional::sample_xmin, sample_xmax
    real(dl),dimension(:),allocatable:: xsample, ysample, w
    integer i, tsize, id, ns
    if(tsize .lt. 2 .or. xmin .ge. xmax) stop "invalid argument in TabFunc_tab"
    call tabfunc_activate(id, tsize, xmin, xmax)
    if(IND_MAX .le. use_cheb_fit_limit)then
       ns = (IND_MAX -1 ) * use_cheb_fit_limit + 1
       allocate(xsample(ns), ysample(ns), w(ns))
       if(present(sample_xmin))then
          x1 = sample_xmin
       else
          x1 = xmin
       endif
       if(present(sample_xmax))then
          xn = sample_xmax
       else
          xn = xmax
       endif
       call findgen(xsample, x1, xn)
       do i = 1,ns
          w(i) = CHEB_WEIGHT(i,ns)
          ysample(i) = f(xsample(i))
       enddo
       call Chebfit(xsample, ysample, w, xmin, xmax, PPV)
       TAB = ysample(1:ns:use_cheb_fit_limit) !!still tabularize, sometime useful for quick estimations
       deallocate(xsample, ysample,w)
       return
    endif
    x = xmin
    do i=1, tsize
       TAB(i) = f(x)
       x = x + DELTAX
    enddo
    x = xmin + DELTAX/2._dl
    do i=1, tsize - 1 
       PPV(i) = (f(x) - (TAB(i) + TAB(i+1))/2._dl)*4._dl 
       x = x + DELTAX
    enddo
  end subroutine tabfunc_tab

!!input id
!!output invid
  subroutine tabfunc_inv(id, invid) !!inverse of a monotonic function
    integer id, invid, ns
    real(dl) ymin, ymax, y, x
    integer i
    if(id.le.0 .or. id.gt. Num_TabFunc .or. invid .eq. id)stop "invalid argument in tabfunc_inv"
    if( .not. tabfunc_used(id)) stop "function is not initialized before calling tabfunc_inv"
    if(invid .eq. id) stop "invalid argument in tabfunc_inv"
    if(IND_MAX .gt. use_cheb_fit_limit)then
       ns = IND_MAX
    else
       ns = (IND_MAX-1)*use_cheb_fit_limit + 1
    endif
    ymin = min(TAB(1), TAB(IND_MAX))
    ymax = max(TAB(1), TAB(IND_MAX))
    if(ymax - ymin .lt. 1.e-20_dl) stop "invid failed: the function is almost a constant."
    call tabfunc_activate(invid, ns, ymin, ymax) !!initialize 
    y = ymin
    do i = 1, ns
       if(.not. tabfunc_solve_monotonic(id, x, y)) stop "error in tabfunc_inv: cannot invert a non-monotonic function"
       tabfunc_space(invid)%table(i) = x
       y = y + tabfunc_space(invid)%dx
    enddo
    y = ymin + tabfunc_space(invid)%dx/2._dl
    do i=1, ns - 1
       if(.not. tabfunc_solve_monotonic(id, x, y)) stop "error in tabfunc_inv: cannot invert a non-monotonic function"
       tabfunc_space(invid)%pp(i)= (x - (tabfunc_space(invid)%table(i) + tabfunc_space(invid)%table(i+1))/2._dl)*4._dl 
       y = y + tabfunc_space(invid)%dx
    enddo
  end subroutine tabfunc_inv

  Subroutine TabFunc_Free(id)
    integer id
    if(allocated(TAB))deallocate(TAB)
    if(allocated(PPV)) deallocate(PPV)
    TabFunc_used(id) = .false.
  end Subroutine TabFunc_Free

  subroutine tabfunc_freeall()
    integer id
    do id =  1,Num_TabFunc
       if(tabfunc_used(id)) call tabfunc_free(id)
    enddo
  end subroutine tabfunc_freeall
  
  subroutine tabfunc_activate(id, tsize, left, right)
    integer id, tsize
    real(dl),optional::left, right
    if(id.le. 0 .or. id.gt. Num_tabfunc)then
       call tabfunc_get_new_id(id)
       call TabFunc_Allocate(id, tsize)
    elseif(.not. tabfunc_used(id))then
       call TabFunc_Allocate(id, tsize)
    elseif(IND_MAX .ne. tsize)then
       call TabFunc_Allocate(id, tsize)
    else
       goto 100
    endif
    IND_MAX = tsize
    IND_L = 1
    IND_R = IND_MAX
100 continue
    if(present(left) .and. present(right))then
       BOUND_L = min(left,right)
       BOUND_R = max(left,right)
       DELTAX = (BOUND_R - BOUND_L)/(tsize - 1)
       if(DELTAX .lt. 1.e-99_dl) call ReturnError("tabfunc_active", "lower bound=upper bound?", "stop")
    endif
  end subroutine tabfunc_activate

  subroutine tabfunc_get_new_id(id)
    integer i, id
    id = -1
    do i=1, Num_TabFunc
       if(.not. tabfunc_used(i))then
          id = i 
          exit
       endif
    enddo
    if(id .eq. -1) stop "TabFunc overflow, please free some functions."
  end subroutine tabfunc_get_new_id

  Subroutine TabFunc_Allocate(id, tsize)
    integer id, tsize
    TabFunc_used(id) = .true.
    if(allocated(TAB))deallocate(TAB)
    if(allocated(PPV)) deallocate(PPV)
    allocate(TAB(tsize))
    allocate(PPV(tsize))
  end Subroutine TabFunc_Allocate

  function Tabfunc_eval_s(id, x, extrapolate)
    real(dl) tabfunc_eval_s, x, r
    integer id, pos
    logical,optional::extrapolate
    if(x .lt. BOUND_L .or. x .gt.  BOUND_R)then
       if(abs(BOUND_L - x) .lt. DELTAX * 1.e-5_dl)then
          Tabfunc_eval_s = TAB(1)
       elseif(abs(BOUND_R - x) .lt. DELTAX*1.e-5_dl)then
          Tabfunc_eval_s = TAB(IND_MAX)
       else
          if(present(extrapolate))then
             if(extrapolate)then
                if(x .lt. BOUND_L)then
                   Tabfunc_eval_s = TAB(1)
                else
                   Tabfunc_eval_s = TAB(IND_MAX)
                endif
                return
             endif
          endif
          tabfunc_eval_s = 0._dl
       endif
       return
    endif
    if(IND_MAX .le. use_cheb_fit_limit)then
       Tabfunc_eval_s = chebfit_value(BOUND_L, BOUND_R, PPV, x)
       return
    endif
    r = (x - BOUND_L)/DELTAX + 1._dl
    pos = min(max(floor(r), 1), IND_MAX-1)
    r = r-pos
    Tabfunc_eval_s = TAB(pos)*(1._dl - r) &
         + r * (TAB(pos+1) + (1._dl-r)*PPV(pos))
  end function Tabfunc_eval_s

  function tabfunc_eval_v(id, x, extrapolate)
    integer id
    real(dl),dimension(:),intent(in)::x
    real(dl) tabfunc_eval_v(size(x))
    integer i
    logical,optional::extrapolate
    if(present(extrapolate))then
       do i=1, size(x)
          tabfunc_eval_v(i) = tabfunc_eval_s(id, x(i), extrapolate)
       enddo
    else
       do i=1, size(x)
          tabfunc_eval_v(i) = tabfunc_eval_s(id, x(i))
       enddo
    endif
  end function tabfunc_eval_v
  
  function tabfunc_eval_linearlog(idlinear, idlog, x) result(y)
    integer idlinear, idlog
    real(dl) x,y
    if(x .ge. tabfunc_space(idlinear)%left .and. x .le. tabfunc_space(idlinear)%right)then
       y = tabfunc_eval_s(idlinear, x)
    else
       y = tabfunc_eval_s(idlog, log(x))
    endif
  end function tabfunc_eval_linearlog

  !!********* warning: this only works for monotonic function ******************
  !!********* root finding for a general function is more complicated, to be done ****
  Function Tabfunc_solve_monotonic(id, x, y, xmin, xmax)
    !!solve f(x) = y in xmin<x<xmax, if f is a monotonic function in xmin<x<xmax
    logical Tabfunc_solve_monotonic
    integer, intent(in)::id
    real(dl),intent(in):: y
    real(dl),intent(out)::x
    real(dl),optional::xmin, xmax
    real(dl) xl, xr, yl, yr, ym
    integer  im
    if(id.le. 0 .or. id .gt. Num_TabFunc)then
       call ReturnError("Tabfunc_solve_monotonic","function not actived","return")
       tabfunc_solve_monotonic = .false.
       return
    elseif(.not. tabfunc_used(id))then
       call ReturnError("Tabfunc_solve_monotonic","function not actived","return")
       tabfunc_solve_monotonic = .false.
       return
    endif
    if(present(xmin) .and. present(xmax))then
       xl = min(xmin, xmax)
       xr = max(xmin, xmax)
       if(xl .lt. BOUND_L-DELTAX*1.d-5 .or. xr .gt. BOUND_R+DELTAX*1.d-5) then
          write(*,*) xl, xr, BOUND_L, BOUND_R
          stop "invalid argument in tabfunc_solve_monotonic"
       endif
       IND_L = max(1, floor((xl-BOUND_L)/DELTAX + 0.999999999)) 
       IND_R = min(IND_MAX, ceiling((xr - BOUND_L)/DELTAX+1.0000001))
    endif
    if(abs(TAB(1)-y).lt. 1.e-12_dl)then
       x = BOUND_L
       goto 100
    elseif(abs(TAB(IND_MAX)-y) .lt. 1.e-12_dl)then
       x = BOUND_R
       goto 100
    endif
    do
       if(TAB(IND_L) .le. y .and. TAB(IND_R) .ge. y)then
          do while(IND_R - IND_L .gt. 1)
             im = (IND_L + IND_R)/2
             if(TAB(im) .le. y)then
                IND_L = im
             else
                IND_R = im
             endif
          enddo
          exit
       elseif(TAB(IND_L) .ge. y .and. TAB(IND_R) .le. y)then
          do while(IND_R - IND_L .gt. 1)
             im = (IND_L + IND_R)/2
             if(TAB(im) .ge. y)then
                IND_L = im
             else
                IND_R = im
             endif
          enddo
          exit
       else
          if(present(xmin) .and. present(xmax))then
             tabfunc_solve_monotonic = .false.
             return
          endif
          if((TAB(IND_L) .le. y .and.  TAB(1) .ge. y  )  &
               .or. ( TAB(IND_L) .ge. y  .and.  TAB(1) .le. y ))then
             IND_R = IND_L
             IND_L = 1
          elseif((TAB(IND_R) .le. y .and.  TAB(IND_MAX) .ge. y  )  &
               .or. ( TAB(IND_R) .ge. y  .and.  TAB(IND_MAX) .le. y ))then
             IND_L = IND_R 
             IND_R = IND_MAX
          else
             tabfunc_solve_monotonic = .false.
             return
          endif
       endif
       if(IND_R - IND_L .le. 1) exit
    enddo

    if(IND_MAX .le. use_cheb_fit_limit)then
       xl = BOUND_L + (IND_L - 1._dl) * DELTAX
       xr = BOUND_L + (IND_R - 1._dl) * DELTAX
       x = (xl+xr)/2._dl
       if(TAB(IND_L) .le. y .and. TAB(IND_R) .ge. y)then
          do while(xr - xl .ge. 1.e-12_dl)
             if(Tabfunc_eval_s(id, x) .ge. y)then
                xr = x
             else
                xl = x
             endif
             x = (xl+xr)/2._dl
          enddo
       elseif(TAB(IND_L) .ge. y .and. TAB(IND_R).le.y)then
          do while(xr - xl .ge. 1.e-12_dl)
             if(Tabfunc_eval_s(id, x) .ge. y)then
                xl = x
             else
                xr = x
             endif
             x = (xl+xr)/2._dl
          enddo
       endif
    else
       ym = TAB(IND_R) - TAB(IND_L)
       yl = y - TAB(IND_L)
       if(abs(ym) .lt. 1.e-12_dl)then
          x = 0._dl
       else
          if(PPV(IND_L) .gt. ym* 1.e-2_dl)then
             yl = 1. + yl/PPV(IND_L)
             yr = ym/PPV(IND_L)
             ym = max(0._dl, yl**2 - 4._dl*yr)
             x  = min(1._dl, (yl -  sign(sqrt(ym), yr) )/2._dl)
          else
             yr = yl/ym
             x = yr*(1._dl - (1._dl - yr)*PPV(IND_L)/ym)
             x = yr - x*(1._dl - x)*PPV(IND_L)/ym
             x = yr - x*(1._dl - x)*PPV(IND_L)/ym
          endif
       endif
       x = BOUND_L + (IND_L - 1._dl + x)*DELTAX 
    endif
100 continue
    if(present(xmin))then
       if(x.lt.xmin-1.e-12_dl)then
          tabfunc_solve_monotonic = .false.
          return
       endif
    endif
    if(present(xmax))then
       if(x.gt.xmax+1.e-12_dl)then
          tabfunc_solve_monotonic = .false.
          return
       endif
    endif
    Tabfunc_solve_monotonic = .true.
  end function Tabfunc_solve_monotonic

  !!find the first maximum
  !!this function is not very accurate
  function TabFunc_Maximum_location(id) result(x)
    real(dl) x, px, ppx, r, xi
    integer id, i, loc(1:1),j
    loc = maxloc(TAB)
    i = loc(1)
    x = BOUND_L + DELTAX * (i-1)
    xi = x
    if(IND_MAX .le. use_cheb_fit_limit)then
       do j = 1, 5
          px = tabfunc_derv_s(id, x)
          ppx = (tabfunc_derv_s(id, x+DELTAX/100._dl) - tabfunc_derv_s(id, x-DELTAX/100._dl))/DELTAX * 50._dl
          if(abs(px) .lt. abs(ppx) * DELTAX )then
             x = x - px/ppx
          endif
       enddo
       x = min(max(BOUND_L, x, xi - DELTAX),BOUND_R, xi + DELTAX) 
    else
       if(i.gt.1)then
          if( PPV(i-1) .gt. 0._dl) then
             r = (PPV(i-1) + TAB(i) - TAB(i-1))/2._dl/PPV(i-1)
             if(r.ge.0._dl .and. r .le. 1._dl)then
                x = x + (r-1) * DELTAX
                return
             endif
          endif
       endif
       if(i.lt. IND_MAX)then
          if( PPV(i) .gt. 0._dl) then
            r = (PPV(i) + TAB(i+1) - TAB(i))/2._dl/PPV(i)
             if(r.ge.0._dl .and. r .le. 1._dl)then
                x = x + r * DELTAX
                return
             endif
          endif
       end if
    endif
  End function TabFunc_Maximum_location

  !!find the first maximum
  !!this function is not very accurate
  function TabFunc_Minimum_location(id) result(x)
    real(dl) x, px, ppx, r, xi
    integer id, i, loc(1:1),j
    loc = minloc(TAB)
    i = loc(1)
    x = BOUND_L + DELTAX * (i-1)
    xi = x
    if(IND_MAX .le. use_cheb_fit_limit)then
       do j = 1, 5
          px = tabfunc_derv_s(id, x)
          ppx = (tabfunc_derv_s(id, x+DELTAX/100._dl) - tabfunc_derv_s(id, x-DELTAX/100._dl))/DELTAX * 50._dl
          if(abs(px) .lt. abs(ppx) * DELTAX )then
             x = x - px/ppx
          endif
       enddo
       x = min(max(BOUND_L, x, xi - DELTAX),BOUND_R, xi + DELTAX) 
    else
       if(i.gt.1)then
          if( PPV(i-1) .lt. 0._dl) then
             r = (PPV(i-1) + TAB(i) - TAB(i-1))/2._dl/PPV(i-1)
             if(r.ge.0._dl .and. r .le. 1._dl)then
                x = x + (r-1) * DELTAX
                return
             endif
          endif
       endif
       if(i.lt. IND_MAX)then
          if( PPV(i) .lt. 0._dl) then
            r = (PPV(i) + TAB(i+1) - TAB(i))/2._dl/PPV(i)
             if(r.ge.0._dl .and. r .le. 1._dl)then
                x = x + r * DELTAX
                return
             endif
          endif
       end if
    endif
  End function TabFunc_Minimum_location

  function Tabfunc_derv_s(id, x)
    integer id, pos
    real(dl) Tabfunc_derv_s, x, r
    if(x .lt. BOUND_L .or. x .gt.  BOUND_R)then
       tabfunc_derv_s = 0._dl
       return
    endif
    if(IND_MAX .le. use_cheb_fit_limit)then
       tabfunc_derv_s = chebfit_derv(BOUND_L, BOUND_R, PPV, x)
       return
    endif
    r = (x - BOUND_L)/DELTAX + 1._dl
    pos = min(max(floor(r), 1), IND_MAX-1)
    r = r-pos
    if(r.lt.0.5_dl)then
       if(pos .eq. 1)then
          Tabfunc_derv_s  = (TAB(pos+1)-TAB(pos) &
               + (1._dl-2._dl*r)*PPV(1)) / DELTAX
       else
          tabfunc_derv_s = ((TAB(pos+1) - TAB(pos))*(0.5_dl+r) + (TAB(pos) - TAB(pos-1))*(0.5-r))/DELTAX
       endif
    else
       if(pos .eq. IND_MAX -1 )then
          Tabfunc_derv_s  = (TAB(pos+1)-TAB(pos) &
               + (1._dl - 2._dl*r)*PPV(IND_MAX-1)) / DELTAX
       else
          tabfunc_derv_s = ((TAB(pos+1) - TAB(pos))*(1.5-r) + (TAB(pos+2) - TAB(pos+1))*(r-0.5))/DELTAX
       endif
    endif
  end function Tabfunc_derv_s

  function tabfunc_derv_v(id, x)
    integer id
    real(dl),dimension(:),intent(in)::x
    real(dl) tabfunc_derv_v(size(x))
    integer i
    do i=1, size(x)
       tabfunc_derv_v(i) =tabfunc_derv_s(id,x(i))
    end do
  end function tabfunc_derv_v

  function TabFunc_int_s(id, a, b)
    real(dl) TabFunc_int_s, xmin, xmax, rmin, rmax, a, b
    integer id, posmin, posmax
    if(a.eq.b)then
       TabFunc_int_s = 0._dl
       return
    elseif(a.lt.b)then
       xmin = a
       xmax = b
    else
       xmin = b
       xmax = a
    endif
    if(xmin .ge. BOUND_R .or. xmax .le. BOUND_L)then
       tabfunc_int_s = 0._dl
       return
    endif
    if(IND_MAX .le. use_cheb_fit_limit)then
       xmin = max(xmin, BOUND_L)
       xmax = min(xmax, BOUND_R)
       tabfunc_int_s = chebfit_int(BOUND_L, BOUND_R, PPV, xmin, xmax)
    else
       rmin = min(max(1._dl, (xmin - BOUND_L)/DELTAX + 1._dl), IND_MAX- 1.e-12_dl)
       posmin = floor(rmin)
       rmin = rmin-posmin
       rmax = min(max((xmax - BOUND_L)/DELTAX + 1._dl, 1._dl), IND_MAX  - 1.e-12_dl)
       posmax = floor(rmax)
       rmax = rmax-posmax
       if(posmin .lt. posmax)then
          TabFunc_int_s = sum(TAB(posmin:posmax)) - (TAB(posmin)+ TAB(posmax))/2._dl + sum(PPV(posmin:posmax-1))/6._dl
       else
          TabFunc_int_s = 0._dl
       endif
       TabFunc_int_s = TabFunc_int_s + rmax*(TAB(posmax) + rmax*((TAB(posmax+1)-TAB(posmax)+PPV(posmax))/2._dl  - rmax/3._dl*PPV(posmax))) 
       TabFunc_int_s = TabFunc_int_s - rmin*(TAB(posmin) + rmin*((TAB(posmin+1)-TAB(posmin)+PPV(posmin))/2._dl  - rmin/3._dl*PPV(posmin))) 
       TabFunc_int_s = TabFunc_int_s *  DELTAX
    endif
    if(a.gt. b) TabFunc_int_s = - TabFunc_int_s
  end function TabFunc_int_s

  function tabfunc_int_v(id, lbd, x)
    integer id
    real(dl) lbd
    real(dl),dimension(:),intent(in)::x
    real(dl) tabfunc_int_v(size(x))
    integer i
    do i=1,size(x)
       tabfunc_int_v(i) = tabfunc_int_s(id, lbd, x(i))
    enddo
  end function tabfunc_int_v

  function Tabfunc_eval_current(x)
    real(dl) x, Tabfunc_eval_current
    Tabfunc_eval_current = Tabfunc_eval_s(TabFunc_current_id, x)
  end function Tabfunc_eval_current

  function tabfunc_derv_current(x)
    real(dl) x, tabfunc_derv_current
    tabfunc_derv_current = tabfunc_derv_s(tabfunc_current_id, x)
  end function tabfunc_derv_current

  function tabfunc_int_current(a, b)
    real(dl) tabfunc_int_current, a, b
    tabfunc_int_current = tabfunc_int(tabfunc_current_id, a, b)
  end function tabfunc_int_current


!!******************** functionals ******************************!!
#ifdef HAS_FFTW
  !!This subroutine returns a function id for the (tabularized) Fourier transformation
  !!  \int_0^{\infty} f(x) sin(k x) / kx * 4\pi x^2 dx
  !! I am assuming f(x)=0 for x>xmax and f is a smooth (finite df/dx) function in 0<x<xmax.
  subroutine TabFunc_Spherical_FFT_from_function(id, f, xmax, tsize)
    external f
    integer tsize
    real(dl) kunit, xmax
    real(dl) f, fx(tsize), x, dx, r4int, tophat
    integer i, id
    if(xmax.le.0._dl .or. tsize .le. use_cheb_fit_limit) stop "invalid argument in TabFunc_Spherical_FFT"
    tophat = f(xmax*0.999999999999999_dl)
    kunit = const_pi / xmax /(tsize*1._dl)**(1._dl/3._dl)
    dx = const_pi/kunit/tsize
    call TabFunc_activate(id, tsize+1, 0._dl, const_pi/dx)
    TAB(1) = 0._dl 
    x = dx
    r4int = 0._dl
    do i=1, tsize
       if(x.le.xmax)then
          fx(i) = x*(f(x)-tophat) 
       else
          fx(i:tsize) = 0._dl
          exit
       endif
       TAB(1) = TAB(1) + x*fx(i)
       r4int = r4int + x**3 * fx(i)
       x = x +dx
    enddo
    r4int = r4int + tophat*xmax**5/5._dl/dx
    TAB(1) = TAB(1) * 4._dl 
    fx(2:tsize) = fx(1:tsize-1) + fx(2:tsize)
    call fft_dstii(tsize, fx, TAB(2:tsize+1))
    do i = 1, tsize
       TAB(i+1) = TAB(i+1) / (kunit * i) 
    enddo
    TAB =  TAB * (const_pi * dx)
    x = const_4pi * xmax**3 * tophat
    if(x.ne. 0._dl)then
       do i=0, tsize
          TAB(i+1) = TAB(i+1) +  x*FT_spherical_tophat(kunit * i*xmax) 
       enddo
    endif

    PPV(1) = r4int*kunit**2 * const_4pi * dx / 6._dl 

    do i = 2, tsize-1
       PPV(i) = (TAB(i) +  TAB(i+1) -TAB(i+2) -  TAB(i-1))/4._dl 
    enddo
    PPV(tsize) = (2._dl*TAB(tsize) - TAB(tsize-1) - TAB(tsize+1))/2._dl
  end subroutine TabFunc_Spherical_FFT_from_function

  subroutine tabfunc_spherical_fft_from_id( id_from, id_to, xmax, tsize)
    integer id_to, id_from, id_tmp, n
    real(dl),optional::xmax
    integer,optional:: tsize
    if(id_to .eq. id_from .or. id_from .le.0 .or. id_from .gt. Num_TabFunc) stop "invalid argument in tabfunc_spherical_fft"
    if(.not. tabfunc_used(id_from)) stop "tabfunc_spherical_fft: the input id is not activated"
    id_tmp = tabfunc_current_id
    tabfunc_current_id = id_from
    if(present(tsize))then
       n = tsize
    else
       if(  tabfunc_space(id_from)%maxind .gt. use_cheb_fit_limit)then
          n = tabfunc_space(id_from)%maxind
       else
          n = use_cheb_fit_limit * tabfunc_space(id_from)%maxind
       endif
    endif
    if(present(xmax))then
       call tabfunc_spherical_fft_from_function(id_to, tabfunc_eval_current, xmax, n)
    else
       call tabfunc_spherical_fft_from_function(id_to, tabfunc_eval_current, tabfunc_space(id_from)%right, n)
    endif
    tabfunc_current_id = id_tmp
  end subroutine tabfunc_spherical_fft_from_id

#endif

  subroutine tabfunc_rescale(id, factor)
    real(dl) factor
    integer id
    PPV = PPV * factor
    TAB = TAB * factor
  end subroutine tabfunc_rescale

  subroutine tabfunc_addconst(id, const)
    integer id
    real(dl) const
    TAB = TAB + const
    if(IND_MAX .le. use_cheb_fit_limit) PPV(1) = PPV(1) + const
  end subroutine tabfunc_addconst

!! get the tabularized function \int_ymin^ymax f(x,y) dy   
  subroutine tabfunc_tab2d(id, f, xmin, xmax, tsize, ymin, ymax)
    real(dl) f, xmin, xmax, ymin, ymax,x
    external f
    real(dl),dimension(:),allocatable:: xsample, ysample,w
    integer i, tsize, id, ns
    if(tsize .lt. 3 .or. xmin .ge. xmax) stop "invalid argument in TabFunc_tab2d"
    call tabfunc_activate(id, tsize, xmin, xmax)
    if(IND_MAX .le. use_cheb_fit_limit)then
       ns = (IND_MAX -1 ) * use_cheb_fit_limit + 1
       allocate(xsample(ns), ysample(ns),w(ns))
       call findgen(xsample, xmin, xmax)
       do i=1,ns
          w(i) = CHEB_WEIGHT(i,ns)
          ysample(i) = qromb_with_secondvar(xsample(i),f, ymin, ymax, 1.e-6_dl)
       enddo
       call Chebfit(xsample, ysample, w, xmin, xmax, PPV)
       TAB = ysample(1:ns:use_cheb_fit_limit) !!still tabularize, sometime useful for quick estimations
       deallocate(xsample, ysample,w)
       return
    endif
    x = xmin
    do i=1, tsize
       TAB(i) = qromb_with_secondvar(x, f, ymin, ymax, 1.e-6_dl)
       x = x + DELTAX
    enddo
    x = xmin + DELTAX/2._dl
    do i=1, tsize - 1 
       PPV(i) = (qromb_with_secondvar(x, f, ymin, ymax, 1.e-6_dl) - (TAB(i) + TAB(i+1))/2._dl)*4._dl 
       x = x + DELTAX
    enddo
  end subroutine tabfunc_tab2d


  !!get a tabularized function F(x) = \int_xmin^x f(t) dt 
  subroutine tabfunc_tabint(id, f, xmin, xmax, tsize)
    real(dl) f, xmin, xmax, x, y0, y1, y2, y3, y4, dxby4, dxby12, dxby6, dx, dxby2, tmp
    external f
    real(dl),dimension(:),allocatable:: xsample, ysample,w
    integer i, tsize, id, ns
    if(tsize .lt. 3 .or. xmin .ge. xmax) stop "invalid argument in TabFunc_tabint"
    call tabfunc_activate(id, tsize, xmin, xmax)
    if(IND_MAX .le. use_cheb_fit_limit)then
       ns = (IND_MAX -1 ) * use_cheb_fit_limit + 1
       allocate(xsample(ns), ysample(ns),w(ns))
       call findgen(xsample, xmin, xmax)
       do i=1,ns
          w(i) = CHEB_WEIGHT(i,ns)
       enddo
       ysample(1) = 0._dl
       x = BOUND_L
       dx = DELTAX / use_cheb_fit_limit
       dxby2 = dx/2._dl
       dxby6 = dx/6._dl
       do i=2,ns
          ysample(i) = ysample(i-1) + (f(x) + 4._dl*f(x+dxby2) + f(x+dx))*dxby6
          x = x + dx
       enddo
       call Chebfit(xsample, ysample, w, xmin, xmax, PPV)
       TAB = ysample(1:ns:use_cheb_fit_limit) !!still tabularize, sometime useful for quick estimations
       deallocate(xsample, ysample,w)
       return
    endif
    x = xmin
    dxby4 =  DELTAX / 4._dl
    dxby6 = DELTAX / 6._dl
    dxby12 = DELTAX / 12._dl
    TAB(1) = 0._dl
    y4  = f(x)
    do i=1, tsize-1
       y0 = y4
       x = x + dxby4
       y1 = f(x)
       x = x + dxby4
       y2 = f(x)
       x = x + dxby4
       y3 = f(x)
       x = x + dxby4
       y4 = f(x)
       tmp = (y0+y1*4._dl+y2)*dxby12
       TAB(i+1) = TAB(i) + (y2 + 4._dl*y3 + y4)*dxby12 + tmp
       PPV(i) = (tmp  + (TAB(i) - TAB(i+1))/2._dl)*4._dl 
    enddo
  end subroutine tabfunc_tabint

  subroutine tabfunc_tabarray(id, y, xmin, xmax, tsize, sample_xmin, sample_xmax)
    integer id
    integer,optional::tsize
    integer i, ysize
    real(dl),intent(in):: xmin, xmax, y(:)
    real(dl) x1, xn
    real(dl),optional::sample_xmin, sample_xmax
    real(dl),dimension(:),allocatable::x
    ysize = size(y)
    if(present(tsize))then
       if(tsize .ne. ysize)then
          if(tsize .gt. use_cheb_fit_limit .or. tsize .le. 0)&
               call returnError("tabfunc_tabarray", "wrong array size", "stop")
       endif
       call tabfunc_activate(id, tsize)
    else
       call tabfunc_activate(id, ysize)
    endif
    BOUND_L = min(xmin, xmax)
    BOUND_R = max(xmin, xmax)
    DELTAX = (BOUND_R - BOUND_L)/(IND_MAX - 1)
    if(IND_MAX .le. use_cheb_fit_limit)then
       allocate(x(ysize))
       if(present(sample_xmin))then
          x1 = sample_xmin
       else
          x1 = xmin
       endif
       if(present(sample_xmax))then
          xn = sample_xmax
       else
          xn = xmax
       endif
       call findgen(x, x1, xn)
       call chebfit(x, y, BOUND_L, BOUND_R, PPV)
       do i = 1, IND_MAX
          TAB(i) = chebfit_value(BOUND_L, BOUND_R, PPV, BOUND_L+(i-1)*DELTAX)
       enddo
       deallocate(x)
       return
    endif
    if(xmin .lt. xmax)then
       TAB = y(1:ysize)
    else
       TAB = y(ysize:1:-1)
    endif
    do i=2, IND_MAX-2
       PPV(i) = (TAB(i) - TAB(i+2) + TAB(i+1) - TAB(i-1))/4._dl
    enddo
    PPV(1) = TAB(2) - (TAB(1)+TAB(3))/2._dl
    PPV(IND_MAX - 1) = TAB(IND_MAX - 1) - (TAB(IND_MAX - 2)+TAB(IND_MAX))/2._dl
  end subroutine tabfunc_tabarray

  subroutine tabfunc_tabode(odert, nvar, yini, xstart, xend, tsize, id_list)
    external odert
    integer,intent(in)::nvar, tsize
    integer,intent(inout)::id_list(nvar)
    real(dl),intent(in)::xstart, xend, yini(nvar)
    real(dl),dimension(:),allocatable::xtmp, y
    real(dl),dimension(:,:),allocatable::ytmp
    integer i,ns
    if(nvar.lt.1. .or. tsize.lt.2) call returnError("tabfunc_ode","n<1?","stop")
    if(tsize .gt. use_cheb_fit_limit) then
       ns = tsize
    else
       ns = tsize*use_cheb_fit_limit
    endif
    allocate(xtmp(ns))
    allocate(y(ns))
    allocate(ytmp(nvar,ns))
    call findgen(xtmp, xstart, xend)
    ytmp(:,1) = yini
    call GenericDverk(odert, xtmp, ytmp, 1.e-6_dl)
    do i=1, nvar
       y = ytmp(i,:)  
       call tabfunc_tabarray(id_list(i), y, xstart, xend, tsize)
    enddo
    deallocate(xtmp,y,ytmp)
  end subroutine tabfunc_tabode

  subroutine tabfunc_tab_with_params(id, f, params, xmin, xmax, tsize, sample_xmin, sample_xmax)
    !!real(dl) f(type(params_pass) params, real(dl) x)
    Interface 
       function f(p, xx)
         use basic_utils
         type(params_pass)p
         real(dl) f,xx
       end function f
    End Interface
    type(params_pass) params
    real(dl) xmin, xmax, x, x1, xn
    real(dl),optional::sample_xmin, sample_xmax
    real(dl),dimension(:),allocatable:: xsample, ysample, w
    integer i, tsize, id, ns
    if(tsize .lt. 3 .or. xmin .ge. xmax)then
       print*,"tsize=",tsize
       print*,"xmin=",xmin
       print*,"xmax=",xmax
       stop "invalid argument in TabFunc_tab_params"
    endif
    call tabfunc_activate(id, tsize, xmin, xmax)
    if(IND_MAX .le. use_cheb_fit_limit)then
       ns = (IND_MAX -1 ) * use_cheb_fit_limit + 1
       allocate(xsample(ns), ysample(ns),w(ns))
       if(present(sample_xmin))then
          x1 = sample_xmin
       else
          x1 = xmin
       endif
       if(present(sample_xmax))then
          xn = sample_xmax
       else
          xn = xmax
       endif
       call findgen(xsample, x1, xn)
       do i=1,ns
          w(i) = CHEB_WEIGHT(i,ns)
          ysample(i) = f(params, xsample(i))
       enddo
       call Chebfit(xsample, ysample, w, xmin, xmax, PPV)
       TAB = ysample(1:ns:use_cheb_fit_limit) !!still tabularize, sometime useful for quick estimations
       deallocate(xsample, ysample,w)
       return
    endif
    x = xmin
    do i=1, tsize
       TAB(i) = f(params, x)
       x = x + DELTAX
    enddo
    x = xmin + DELTAX/2._dl
    do i=1, tsize - 1 
       PPV(i) = (f(params, x) - (TAB(i) + TAB(i+1))/2._dl)*4._dl 
       x = x + DELTAX
    enddo
  end subroutine tabfunc_tab_with_params

  function tabfunc_next_zero(id, xstart, zero) 
    integer id, ysign
    real(dl) xstart, x, zero, dx, y, ddx, xm
    logical tabfunc_next_zero
    if(xstart .gt. BOUND_R+1.d-12*DELTAX)then
       tabfunc_next_zero = .false.
       return
    endif
    dx = DELTAX/20.d0
    zero = max(xstart, BOUND_L) + dx/10._dl
    y = tabfunc_eval_s(id,zero)
    if(y .eq.0._dl)return
    if(y.gt.0._dl)then
       ysign = 1
    else
       ysign = -1
    endif
    do while( zero .le. BOUND_R)
       zero = zero+dx
       y = tabfunc_eval_s(id, zero)
       if( y*ysign .le. 0.d0)then
          x = zero - dx
          ddx = dx/1.e8_dl
          do while(abs(zero-x).gt. ddx)
             xm = (x+zero)/2._dl
             if(tabfunc_eval_s(id, xm)*ysign .le. 0.d0)then
                zero = xm
             else
                x=xm
             endif
          enddo
          zero = (zero+x)/2.d0
          tabfunc_next_zero = .true.
          return
       endif
    enddo
    tabfunc_next_zero = .false.
    return
  end function tabfunc_next_zero


  function tabfunc_next_solution(id, val, xstart, solution) 
    integer id, ysign
    real(dl) xstart, x, solution, dx, y, ddx, xm, val
    logical tabfunc_next_solution
    if(xstart .gt. BOUND_R+1.d-12*DELTAX)then
       tabfunc_next_solution = .false.
       return
    endif
    dx = DELTAX/20.d0
    solution = max(xstart, BOUND_L) + dx/10._dl
    y = tabfunc_eval_s(id,solution) - val
    if(y .eq.0._dl)return
    if(y.gt.0._dl)then
       ysign = 1
    else
       ysign = -1
    endif
    do while( solution .le. BOUND_R)
       solution = solution+dx
       y = tabfunc_eval_s(id, solution) - val
       if( y*ysign .le. 0.d0)then
          x = solution - dx
          ddx = dx/1.e8_dl
          do while(abs(solution-x).gt. ddx)
             xm = (x+solution)/2._dl
             if(tabfunc_eval_s(id, xm - val)*ysign.le. 0.d0)then
                solution = xm
             else
                x=xm
             endif
          enddo
          solution = (solution+x)/2.d0
          tabfunc_next_solution = .true.
          return
       endif
    enddo
    tabfunc_next_solution = .false.
    return
  end function tabfunc_next_solution

end module tabfunc_utils
