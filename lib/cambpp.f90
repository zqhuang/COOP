module camb_mypp
  implicit none
  integer,parameter::mypp_n = 1024
  real*8,parameter::mypp_EulerC = 0.5772156649d0, mypp_pi = 3.14159265359d0
  real*8::mypp_As = 2.3d-9
  real*8::mypp_ns = 0.967
  real*8::mypp_At = 0.d0
  real*8::mypp_nt = 0.d0
  integer::mypp_nknots = 0
  real*8,parameter::mypp_kpiv = 0.05
  integer::mypp_ipivot = 0  
  real*8,parameter::mypp_lnkpiv = log(mypp_kpiv)
  real*8,parameter::mypp_lnkmin = -9.22d0
  real*8::mypp_lnk_per_knot = 0.d0
  real*8::mypp_lnkmax = -0.3d0
  integer::mypp_model = 1, mypp_nleft = 0, mypp_nright = 0
  real*8,dimension(mypp_n)::mypp_lnk, mypp_lnps, mypp_lnps2, mypp_lnpt, mypp_lnpt2, mypp_lneps, mypp_lnV, mypp_phi, mypp_lnH
  real*8::mypp_dlnk

  
!!  character(LEN=1024)::cosmomc_paramnames = ''

contains

  function Wavelet_lnPk(As, ns, w, lnk)
    real*8, intent(in) :: lnk, As, ns, w(:)
    real*8  wavelet_lnpk
    real*8  lnrat
    real*8::wsum(4)
    integer i
    !lnk = log(k)
    lnrat = lnk - log(mypp_kpiv)
    wsum = 0.d0
    do i=1, 5
       wsum(1) = wsum(1) + w(i) * daubechies4_psi(lnrat + i - 3)
    enddo
    do i=6,14
       wsum(2) = wsum(2) + w(i) * daubechies4_psi(2*lnrat + i - 10)
    enddo
    do i=15, 31
       wsum(3) = wsum(3) + w(i) * daubechies4_psi(4*lnrat + i - 23)
    enddo
    !!well this part is because of mx's ''stupid setting''
    do i=32, 48
       wsum(4) = wsum(4) + w(i) * daubechies4_psi(8*lnrat + i - 48)
    enddo
    do i=49, 64
       wsum(4) = wsum(4) + w(i) * daubechies4_psi(8*lnrat + 65-i)
    enddo
    wsum(2) = wsum(2)*sqrt(2.d0)
    wsum(3) = wsum(3)*2.d0
    wsum(4) = wsum(4)*sqrt(8.d0)
    Wavelet_lnPk = log(As) + lnrat * (ns-1.0) + sum(wsum)

  contains
    
    function free_file_unit() result(lun)
      integer lun
      logical used
      do lun = 20, 99
         Inquire(Unit=lun, Opened=used)
         If (.not. used) return
      enddo
      lun = 0
    End function free_file_unit


    function daubechies2_psi(t) result(psi)
      real(kind=8)::t, psi, r
      integer,parameter::n=20000 !!totally n+1 samples
      real(kind=8),parameter::lbd = -5.d0 !!left boundary
      real(kind=8),parameter::rbd = 5.d0  !!right boundary
      real(kind=8),parameter::step = (rbd-lbd)/n  !!sample time step
      real(kind=8),save::table(0:n)  !!table to save computed psi samples
      logical,save::loaded = .false.
      integer::funit, i, ir
      if(.not. loaded)then
         funit = free_file_unit()
         open(funit, FILE = 'psi2.csv')
         do i=0, n
            read(funit, *) table(i)
         enddo
         close(funit)
         loaded = .true.
      endif
      r = (t-lbd)/step
      ir = floor(r)
      if(r .lt. 0 .or. r.ge. n)then !!out of boundary return 0
         psi = 0.d0
         return
      endif
      r = r-ir
      psi = table(ir)*(1.d0-r)+table(ir+1)*r
    end function daubechies2_psi


    function daubechies4_psi(t) result(psi)
      real(kind=8)::t, psi, r
      integer,parameter::n=20000 !!totally n+1 samples
      real(kind=8),parameter::lbd = -5.d0 !!left boundary
      real(kind=8),parameter::rbd = 5.d0  !!right boundary
      real(kind=8),parameter::step = (rbd-lbd)/n  !!sample time step
      real(kind=8),save::table(0:n)  !!table to save computed psi samples
      logical,save::loaded = .false.
      integer::funit, i, ir
      if(.not. loaded)then
         funit = free_file_unit()
         open(funit, FILE = 'psi4.csv')
         do i=0, n
            read(funit, *) table(i)
         enddo
         close(funit)
         loaded = .true.
      endif
      r = (t-lbd)/step
      ir = floor(r)
      if(r .lt. 0 .or. r.ge. n)then !!out of boundary return 0
         psi = 0.d0
         return
      endif
      r = r-ir
      psi = table(ir)*(1.d0-r)+table(ir+1)*r
    end function daubechies4_psi

    
  end function Wavelet_LnPk



  subroutine mypp_set_uniform(n, x, lower, upper, logscale)
    integer n, i
    real*8 x(n), lower, upper, rlow, dx
    logical,optional::logscale
    if(n.eq. 0)return
    x(1) = lower
    if(n.eq.1)return
    x(n) = upper
    if(present(logscale))then
       if(logscale)then
          dx = (log(upper) - log(lower))/(n-1)
          rlow = log(lower) - dx
          !$omp parallel do
          do i = 2, n-1
             x(i) = exp(rlow + dx*i)
          enddo
          !$omp end parallel do
          return
       endif
    endif
    dx = (upper-lower)/(n-1)
    rlow = lower-dx
    !$omp parallel do
    do i = 2, n-1
       x(i) =rlow + dx*i
    enddo
    !$omp end parallel do
  end subroutine mypp_set_uniform


  subroutine mypp_locate(n, x, needle, loc, res)
    integer n, loc, imin, imax
    real*8 x(n), needle, res
    if(x(1).le. x(n))then
       if(needle .lt. x(1))then
          loc = 0
          return
       endif
       if(needle .gt. x(n))then
          loc = n
          return
       endif
       imin = 1
       imax = n
       do while(imax - imin .gt. 1)
          loc = (imax + imin)/2
          if(x(loc).le. needle)then
             imin = loc
          else
             imax = loc
          endif
       enddo
       loc = imin
       if(needle .gt. x(imin))then
          res = (needle - x(imin))/(x(imax)-x(imin))
       else
          res = 0
       endif
       return
    else
       if(needle .gt. x(1))then
          loc = 0
          return
       endif
       if(needle .lt. x(n))then
          loc = n
          return
       endif
       imin = 1
       imax = n
       do while(imax - imin .gt. 1)
          loc = (imax + imin)/2
          if(x(loc).ge. needle)then
             imin = loc
          else
             imax = loc
          endif
       enddo
       loc = imin
       if(needle .lt. x(imin))then
          res = (needle - x(imin))/(x(imax)-x(imin))
       else
          res = 0
       endif
       return       
    endif
  end subroutine mypp_locate
  


  subroutine mypp_spline(n, x, y, y2, ypl, ypr)
    integer n, i
    real*8 x(n), y(n), y2(n)
    real*8, optional::ypl,ypr
    real*8 yil, yir, bet, dxr, dxl
    real*8 gam(n-1)
    if(n.le.2)then
       y2 = 0.
       return
    endif
    dxr = x(2) - x(1)
    yir=(y(2)-y(1))/dxr
    if(present(ypl))then
       y2(1)=(yir-ypl)/dxr*3.
       gam(1)= 0.5
    else
       y2(1)=0.
       gam(1)=0.
    endif
    dxr = dxr/6.
    do i=2, n-1
       dxl = dxr
       dxr=x(i+1)-x(i)
       bet=(x(i+1)-x(i-1))/3.-dxl*gam(i-1)
       if(abs(bet) .lt. 1.d-50) stop 'Error in Spline.'
       yil=yir
       yir=(y(i+1)-y(i))/dxr
       y2(i)=(yir-yil-dxl*y2(i-1))/bet
       dxr=dxr/6.
       gam(i)=dxr/bet
    enddo
    if(present(ypr))then
       bet=(x(n)-x(n-1))/3.-dxr*gam(n-1)
       if(abs(bet) .lt. 1.d-50) stop 'Error in Spline.'
       y2(n)=(ypr-yir-dxr*y2(n-1))/bet
    else
       y2(n)=0.
    endif
    do i=n-1, 1 , -1
       y2(i)=y2(i)-gam(i)*y2(i+1)
    enddo
  end subroutine mypp_spline


  subroutine mypp_splint(n, x, y, y2, xs, ys)
    integer n, l, r
    real*8 x(n), y(n), y2(n)
    real*8 xs, ys, a, b
    call mypp_locate(n, x, xs, l, b)
    if(l .lt. 1)then
       ys=y(1)
       return
    endif
    if( l .ge. n)then
       ys=y(n)
       return
    endif
    r = l + 1
    a =  1.d0 - b
    ys=y(l)*a+y(r)*b+  &
         (y2(l)*(a*a-1.)*a+y2(r)*(b*b-1.)*b)/6.*(x(r)-x(l))**2
  end subroutine mypp_splint
  

  


  subroutine mypp_setup_pp(As, ns, nknots, dlnps, r, nt)
    real*8::As, ns
    real*8,optional::r, nt
    integer::nknots
    real*8::dlnps(nknots)
    real*8,dimension(:),allocatable::lnk, lnps, lnps2
    integer  i
    real*8 dlnk
    select case(mypp_model)
    case(1)  !!spline
       if(nknots .lt. 5) stop "You need at least 5 knots for scan_spline mode"
       mypp_nknots = nknots
       mypp_nleft = nint(nknots* (mypp_lnkpiv-mypp_lnkmin) / (-mypp_lnkmin))
       mypp_nright = nknots - mypp_nleft 
       dlnk = (mypp_lnkpiv-mypp_lnkmin)/mypp_nleft
       mypp_lnk_per_knot = dlnk
       mypp_lnkmax = mypp_lnkmin + nknots * dlnk
       allocate(lnk(0:nknots), lnps(0:nknots), lnps2(0:nknots))
       call mypp_set_uniform(nknots+1, lnk, mypp_lnkmin, mypp_lnkmax)
       lnps(0:mypp_nleft-1) = dlnps(1:mypp_nleft)
       lnps(mypp_nleft) = 0.d0
       lnps(mypp_nleft+1:nknots) = dlnps(mypp_nleft+1:nknots)
       call mypp_spline(nknots+1, lnk, lnps, lnps2)

       call mypp_set_uniform(mypp_n, mypp_lnk, mypp_lnkmin, mypp_lnkmax)        
       !$omp parallel do
       do i=1, mypp_n
          call mypp_splint(nknots+1, lnk, lnps, lnps2, mypp_lnk(i), mypp_lnps(i))
       enddo
       !$omp end parallel do
       deallocate(lnk, lnps, lnps2)
       mypp_As = As
       mypp_ns = ns
       mypp_lnps = mypp_lnps + log(mypp_As) &
            + (mypp_ns - 1.d0 ) * (mypp_lnk - mypp_lnkpiv)
       call mypp_spline(mypp_n, mypp_lnk, mypp_lnps, mypp_lnps2)
       if(present(r))then
          mypp_At = r*mypp_As
          if(present(nt))then
             mypp_nt = nt
          else
             mypp_nt = -r/8.d0
          endif
          mypp_lnpt = log(mypp_At) + mypp_nt*(mypp_lnk - mypp_lnkpiv)
          mypp_lnpt2 = 0.d0
       else
          mypp_At = 0.d0
          mypp_nt = 0.d0
          mypp_lnpt = -50.d0
          mypp_lnpt2 = 0.d0
       endif
    case(2) !!wavelet
       mypp_As = As
       mypp_ns = ns       
       call mypp_set_uniform(mypp_n, mypp_lnk, mypp_lnkmin, mypp_lnkmax)
       do i=1, mypp_n
          mypp_lnps(i) = wavelet_lnPk(mypp_As, mypp_ns, dlnps, mypp_lnk(i))
       enddo
       call mypp_spline(mypp_n, mypp_lnk, mypp_lnps, mypp_lnps2)
       if(present(r))then
          mypp_At = r*mypp_As
          if(present(nt))then
             mypp_nt = nt
          else
             mypp_nt = -r/8.d0
          endif
          mypp_lnpt = log(mypp_At) + mypp_nt*(mypp_lnk - mypp_lnkpiv)
          mypp_lnpt2 = 0.d0
       else
          mypp_At = 0.d0
          mypp_nt = 0.d0
          mypp_lnpt = -50.d0
          mypp_lnpt2 = 0.d0
       endif
    end select
  end subroutine mypp_setup_pp


  function mypp_primordial_lnps(kMpc)  result(lnps)
    real*8 kMpc, lnps
    call mypp_splint(mypp_n, mypp_lnk, mypp_lnps, mypp_lnps2, log(kMpc), lnps)
  end function mypp_primordial_lnps


  function mypp_primordial_ps(kMpc)  result(ps)
    real*8 kMpc, ps
    ps = exp(mypp_primordial_lnps(kMpc))
  end function mypp_primordial_ps

  function mypp_primordial_lnpt(kMpc)  result(lnpt)
    real*8 kMpc, lnpt
    if(mypp_At .eq. 0.d0)then
       lnpt = -50.d0
    else
       call mypp_splint(mypp_n, mypp_lnk, mypp_lnpt, mypp_lnpt2, log(kMpc), lnpt)
    endif
  end function mypp_primordial_lnpt

  function mypp_primordial_pt(kMpc)  result(pt)
    real*8 kMpc, pt
    if(mypp_At .eq. 0.d0)then
       pt = 0.d0
    else
       pt = exp(mypp_primordial_lnpt(kMpc))
    endif
  end function mypp_primordial_pt

  subroutine mypp_get_eps_phi_lnV(kMpc, eps, phi, lnV)
    real*8 kMpc, eps, phi, lnV
    real*8 lnk, r
    integer il, ir
    lnk = log(kMpc)
    r= (lnk - mypp_lnkmin)/mypp_dlnk+1.d0
    il = min(max(floor(r), 1), mypp_n -1)
    ir = il + 1
    eps = exp(mypp_lneps(il)*(ir - r) + mypp_lneps(ir)*(r-il))
    phi = mypp_phi(il)*(ir-r) + mypp_phi(ir)*(r-il)
    lnV = mypp_lnV(il)*(ir-r) + mypp_lnV(ir)*(r-il)
  end subroutine mypp_get_eps_phi_lnV
  
  subroutine mypp_get_potential()
    integer::iloc(1:1)
    integer i
    real*8, parameter::max_delta = 0.1
    real*8, parameter::max_lneps = log(0.1)
    real*8 dlnk, fourdlnk, dphiby2(mypp_n), eps(mypp_n), delta(mypp_n)
    iloc = minloc(abs(mypp_lnk - mypp_lnkpiv))
    mypp_ipivot = iloc(1)
    dlnk = mypp_lnk(2)-mypp_lnk(1)
    mypp_dlnk = dlnk
    fourdlnk = dlnk*4.d0
    do i=2, mypp_n - 1
       delta(i) = (mypp_lnps(i-1)-mypp_lnps(i+1))/fourdlnk
       if(abs(delta(i)) .gt. max_delta)then
          delta(i) = sign(max_delta, delta(i))
       endif
    enddo
    delta(1) = (mypp_lnps(1)-mypp_lnps(2))/dlnk/2.d0 !delta(2)
    delta(mypp_n) = (mypp_lnps(mypp_n-1)-mypp_lnps(mypp_n))/dlnk/2.d0 !delta(mypp_n  - 1)
    i=1
    if(abs(delta(i)) .gt. max_delta)then
       delta(i) = sign(max_delta, delta(i))
    endif
    i=mypp_n
    if(abs(delta(i)) .gt. max_delta)then
       delta(i) = sign(max_delta, delta(i))
    endif    
    mypp_lneps = min(mypp_lnpt - mypp_lnps - log(16.d0), max_lneps) 
    eps = exp(mypp_lneps)
    mypp_lneps = min(max_lneps,  mypp_lneps - log( (1.d0 - (2.d0*(log(2.d0)+mypp_EulerC-1.d0))*eps ) /(1.d0-2.d0*eps+(2.d0*(2.d0-mypp_EulerC-log(2.d0)))*delta))) !!slow-roll correction
    eps = exp(mypp_lneps)
    mypp_lnH = log(mypp_pi**2/2.d0*exp(mypp_lnpt)/( 1.d0 - (2.d0*(log(2.d0)+mypp_EulerC-1.d0))*eps ))/2.d0
    mypp_lnV = 2.d0*mypp_lnH + log(3.d0*(1.d0-eps/3.d0))
    mypp_phi(1) = 0.d0
    dphiby2 = sqrt(2.d0*eps)/(1.d0-eps)*dlnk/2.d0
    do i=2, mypp_n 
       mypp_phi(i) = mypp_phi(i-1) + (dphiby2(i-1)+dphiby2(i))
    enddo
    mypp_phi = mypp_phi - mypp_phi(mypp_ipivot)
  end subroutine mypp_get_potential


end module camb_mypp
