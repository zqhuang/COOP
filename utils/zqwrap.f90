module cosmolibwrap
  use wrap_utils
  use solvep
  implicit none

#include "utils.h"

  real(dl)::pp_lnk_min = 1.e30
  real(dl)::pp_lnk_max = -1.e30
  integer,parameter::pp_n = 2048+1
  integer::pp_ipivot
  integer::pp_dip_lmin = 20
  integer::pp_dip_lmax = 30
  real(dl)::pp_dip_Amp = 0.
  real(dl)::pp_As, pp_At, pp_nsm1, pp_nrunby2, pp_nt, pp_dlnk
  real(dl),dimension(pp_n)::pp_lnk, pp_lnps, pp_lnpt, pp_lnps2, pp_lnpt2, pp_lnV, pp_phi, pp_epsilon
  SHORT_STRING :: pp_mode = "standard"
  logical::pp_scan_tensor = .false.
  logical::pp_scan_scalar = .false.
  
  interface pp_scalar_power
     module procedure pp_scalar_power_s, pp_scalar_power_v
  end interface pp_scalar_power


  interface pp_tensor_power
     module procedure pp_tensor_power_s, pp_tensor_power_v
  end interface pp_tensor_power


  interface pp_init
     module procedure pp_init_d, pp_init_s
  end interface pp_init


!  public pp_scalar_power, pp_tensor_power, pp_init, pp_mode, pp_get_potential, pp_lnk, pp_phi, pp_lnV, pp_n, pp_epsilon, pp_cosmomc_init, pp_lnk_min, pp_lnk_max, pp_higgs_phi, pp_higgs_chi, pp_higgs_epsilon, pp_higgs_eta, pp_higgs_xiv, pp_higgs_potential


contains


  subroutine pp_cosmomc_init(initpower, inflation_consistency, k_pivot)
    real,optional::k_pivot
    real(dl) kpiv, lnkmin, lnkmax
    integer num_init_power
    real(dl),dimension(:),intent(IN)::initpower
    real(dl) As, ns, nrun, r, nt, Aphiphi
    logical,optional::inflation_consistency
    real(dl),dimension(:),allocatable::lnps

    num_init_power = size(initpower)
    ns = initpower(1)
    nt = initpower(2)
    nrun = initpower(3)
    As = exp(initpower(4))*1.e-10_dl
    r = initpower(5)
    if(present(inflation_consistency))then
       if(inflation_consistency)nt = -r/8.d0
    endif
    Aphiphi = initpower(6)
    if(present(k_pivot))then
       kpiv = k_pivot
    else
       kpiv = 0.05_dl
    endif
    select case(pp_mode)
    case("standard", "elldip")
       call pp_init(As = As, ns = ns, nrun  = nrun,  r = r, nt = nt, lnkmin = log(1.d-5/kpiv), lnkmax = log(1.d0/ kpiv))
       pp_dip_Amp = initpower(7)
    case("spline2", "linear2")
       lnkmin = -6.6d0
       allocate(lnps(num_init_power - 6+3))
       select case(num_init_power)
       case(10)  !!7 knots
          lnkmax = -lnkmin*2.d0/4.d0
          lnps(1:3) = Initpower(7:9)
          lnps(4:6) = 0.
          lnps(7) =  Initpower(10)              
       case(11)  !!8 knots
          lnkmax = -lnkmin*2.d0/5.d0
          lnps(1:4) = Initpower(7:10)
          lnps(5:7) = 0.
          lnps(8) =  Initpower(11)
       case(12) !!9 knots
          lnkmax = -lnkmin*2.d0/6.d0
          lnps(1:5) = Initpower(7:11)
          lnps(6:8) = 0.
          lnps(9) =  Initpower(12)
       case(13) !!10 knots
          lnkmax = -lnkmin*3.d0/6.d0
          lnps(1:5) = Initpower(7:11)
          lnps(6:8) = 0.
          lnps(9:10) =  Initpower(12:13)
       case(14)  !! 11 knots
          lnkmax = -lnkmin*3.d0/7.d0
          lnps(1:6) = Initpower(7:12)
          lnps(7:9) = 0.
          lnps(10:11) =  Initpower(13:14)
       case default
          stop "number of knots is beyond expectation (7-11)"
       end select
       call pp_init(As = As, ns = ns, nrun  = nrun ,  r = r, nt = nt, lnkmin = lnkmin, lnkmax = lnkmax, psparams =lnps )
       deallocate(lnps)
    case("spline1","linear1")
       lnkmin = -6.6d0
       allocate(lnps(num_init_power-6+2))
       select case(num_init_power)
       case(11)  !!7 knots
          lnkmax = -lnkmin*2.d0/4.d0
          lnps(1:3) = Initpower(7:9)
          lnps(4) =  Initpower(10)
          lnps(5) = 0.
          lnps(6) = lnps(4)
          lnps(7) =  Initpower(11)              
       case(12)  !!8 knots
          lnkmax = -lnkmin*2.d0/5.d0
          lnps(1:4) = Initpower(7:10)
          lnps(5) =  Initpower(11)
          lnps(6) = 0.
          lnps(7) = lnps(5)
          lnps(8) =  Initpower(12)
       case(13) !!9 knots
          lnkmax = -lnkmin*2.d0/6.d0
          lnps(1:5) = Initpower(7:11)
          lnps(6) = Initpower(12)
          lnps(7) = 0.
          lnps(8) = Initpower(12)
          lnps(9) =  Initpower(13)
       case(14) !!10 knots
          lnkmax = -lnkmin*3.d0/6.d0
          lnps(1:5) = Initpower(7:11)
          lnps(6) = Initpower(12)
          lnps(7) = 0.
          lnps(8) = Initpower(12)
          lnps(9:10) =  Initpower(13:14)
       case(15)  !! 11 knots
          lnkmax = -lnkmin*3.d0/7.d0
          lnps(1:6) = Initpower(7:12)
          lnps(7) = Initpower(13)
          lnps(8) = 0.
          lnps(9) = Initpower(13)
          lnps(10:11) =  Initpower(14:15)
       case default
          stop "number of knots is beyond expectation (7-11)"
       end select
       call pp_init(As = As, ns =  ns, nrun  =  0.d0,  r = r, nt =  nt, lnkmin = lnkmin, lnkmax = lnkmax, psparams =lnps )
       deallocate(lnps)
    case("spline0", "linear0")
       lnkmin = -6.6d0
       allocate(lnps(num_init_power-6+1))
       select case(num_init_power)
       case(12)  !!7 knots
          lnkmax = -lnkmin*2.d0/4.d0
          lnps(1:4) = Initpower(7:10)
          lnps(5) = 0.
          lnps(6:7) =  Initpower(11:12)              
       case(13)  !!8 knots
          lnkmax = -lnkmin*2.d0/5.d0
          lnps(1:5) = Initpower(7:11)
          lnps(6) = 0.
          lnps(7:8) =  Initpower(12:13)
       case(14) !!9 knots
          lnkmax = -lnkmin*2.d0/6.d0
          lnps(1:6) = Initpower(7:12)
          lnps(7) = 0.
          lnps(8:9) =  Initpower(13:14)
       case(15) !!10 knots
          lnkmax = -lnkmin*3.d0/6.d0
          lnps(1:6) = Initpower(7:12)
          lnps(7) = 0.
          lnps(8:10) =  Initpower(13:15)
       case(16)  !! 11 knots
          lnkmax = -lnkmin*3.d0/7.d0
          lnps(1:7) = Initpower(7:13)
          lnps(8) = 0.
          lnps(9:11) =  Initpower(14:16)
       case(17)  !! 12 knots
          lnkmax = -lnkmin*3.d0/8.d0
          lnps(1:8) = Initpower(7:14)
          lnps(9) = 0.
          lnps(10:12) =  Initpower(15:17)
       case(18)  !! 13 knots
          lnkmax = -lnkmin*4.d0/8.d0
          lnps(1:8) = Initpower(7:14)
          lnps(9) = 0.
          lnps(10:13) =  Initpower(15:18)
       case(19)  !! 14 knots
          lnkmax = -lnkmin*4.d0/9.d0
          lnps(1:9) = Initpower(7:15)
          lnps(10) = 0.
          lnps(11:14) =  Initpower(16:19)
       case(20)  !! 15 knots
          lnkmax = -lnkmin*4.d0/10.d0
          lnps(1:10) = Initpower(7:16)
          lnps(11) = 0.
          lnps(12:15) =  Initpower(17:20)
       case default
          stop "number of knots is beyond expectation (7-15)"
       end select
       call pp_init(As = As, ns = 0.d0, nrun  = 0.d0,  r = r, nt =  nt, lnkmin = lnkmin, lnkmax = lnkmax, psparams =lnps )
       deallocate(lnps)
    case("m2phi2")
       call pp_init(As = As, ns = ns, nrun  =  nrun,  r = r, nt = nt, lnkmin = log(1.d-5/kpiv), lnkmax = log(1.d0/kpiv), psparams=(/ 1.d-5/sqrt(Initpower(7)), Initpower(8) /) )
    case("bump")
       call pp_init(As = As, ns=ns, nrun = nrun, r= r, nt = nt,  lnkmin = log(1.d-5/kpiv), lnkmax = log(1.d0/kpiv), psparams = (/ 1.d-6*InitPower(7), InitPower(8), InitPower(9), InitPower(10), InitPower(11) /))  !!mphi (10^-6Mp), nefolds, bump_dn, bump_amp, bump_width
    case("lambdaphi4")
       call pp_init(As = As, ns = ns, nrun  =  nrun,  r = r, nt = nt, lnkmin = log(1.d-5/kpiv), lnkmax = log(1.d0/kpiv), psparams=(/ 1.d-13*Initpower(7), Initpower(8) /) )
    case("higgs")
       call pp_init(As = As, ns = ns, nrun  =  nrun,  r = r, nt = nt, lnkmin = log(1.d-5/kpiv), lnkmax = log(1.d0/kpiv), psparams=(/ Initpower(7), Initpower(8), Initpower(9), InitPower(10) /) )
    case default
       write(*,*) trim(pp_mode)
       stop "Unknown pp_mode"       
    end select
  end subroutine pp_cosmomc_init


  subroutine pp_init_d(As, ns, nrun, r, nt, lnkmin, lnkmax, psparams, ptparams, tmp_mode)
    real(dl),parameter::q(0:5) = (/ 1.d0, 1.062970488d0, 0.2092751679d0, 0.1004277031d0, -0.02231133892d0, 0.01064541727d0 /)
    real(dl),dimension(:),intent(IN),optional::psparams, ptparams
    real(dl) As, ns, nrun, r, nt, eps, eta, xi, alpha, m, phi, V, lambda, xiphi2
    real(dl)::lnkmin, lnkmax
    real(dl),dimension(:),allocatable::lnk, lnp2
    integer np, i
    UNKNOWN_STRING,optional::tmp_mode
    SHORT_STRING::pp_mode_save
    if(present(tmp_mode))then
       pp_mode_save = pp_mode
       pp_mode = tmp_mode
    endif
    if(abs(pp_lnk_min - lnkmin).gt. 1.d-5 .or. abs(pp_lnk_max - lnkmax).gt. 1.d-5)then
       pp_lnk_min = lnkmin
       pp_lnk_max = lnkmax
       call set_uniform(pp_n, pp_lnk, pp_lnk_min, pp_lnk_max)
       pp_dlnk = pp_lnk(2) - pp_lnk(1)
    endif
    pp_ipivot = nint(-lnkmin/pp_dlnk + 1.)
    if(pp_ipivot .lt. 1 .or. pp_ipivot .gt. pp_n) stop "pp_init: bad input of lnkmin and lnkmax"
    pp_As = As
    pp_nsm1 = ns - 1.
    pp_nrunby2 = nrun/2.
    pp_At = As * r
    pp_nt = nt
    select case(trim(pp_mode))
    case("standard", "elldip")
       pp_scan_scalar = .false.
       pp_scan_tensor = .false.
    case("m2phi2")  
       pp_scan_scalar = .false.
       pp_scan_tensor = .false.
       if(present(psparams))then  !!I assume the input parameters are m and N-efolds \equiv \int_0^\phi V/V' d\phi
          m = psparams(1)
          phi = sqrt(4.d0*psparams(2))
          V = 0.5*(m*phi)**2
          eps = 2./phi**2
          eta = eps
          xi = 0.
          pp_As = V/eps/(24*const_pi2)*(1. + (6.*q(1) - 7./3.)*eps - 2.*q(1)*eta- 2.*q(2)*xi) 
          pp_nsm1 = 2.*eta - 6.*eps + 2.*q(1)*xi
          pp_nrunby2 = - xi - 12.*eps**2 + 8.*eps*eta
          pp_At = (2./3./const_pi2)*V*(1.+eps/3.)
          pp_nt = - 2.*eps
       endif
    case("lambdaphi4")  
       pp_scan_scalar = .false.
       pp_scan_tensor = .false.
       if(present(psparams))then  !!I assume the input parameters are lambda and number of e-folds
          lambda = psparams(1)
          phi = sqrt(8.d0*psparams(2))
          V = 0.25 * lambda * phi **4
          eps = 8./phi**2
          eta = 12./phi**2
          xi = 96/phi**4
          pp_As = V/eps/(24*const_pi2)*(1. + (6.*q(1) - 7./3.)*eps - 2.*q(1)*eta- 2.*q(2)*xi) 
          pp_nsm1 = 2.*eta - 6.*eps + 2.*q(1)*xi
          pp_nrunby2 =  - xi - 12.*eps**2 + 8.*eps*eta
          pp_At = (2./3./const_pi2)*V*(1.+eps/3.)
          pp_nt = - 2.*eps
       endif
    case("higgs")  !!
       pp_scan_scalar = .false.
       pp_scan_tensor = .false.
       if(present(psparams))then 
          lambda = psparams(1)
          phi = pp_higgs_phi(nefolds = psparams(2), xi = psparams(3), alpha = psparams(4)) 
          V = pp_higgs_potential(phi = phi, lambda = lambda, xi = psparams(3), alpha = psparams(4))
          eps = pp_higgs_epsilon(phi = phi,xi = psparams(3), alpha = psparams(4))
          eta = pp_higgs_eta(phi = phi,xi = psparams(3), alpha = psparams(4))
          xi = pp_higgs_xiv(phi = phi,xi = psparams(3), alpha = psparams(4))
          pp_As = V/eps/(24*const_pi2)*(1. + (6.*q(1) - 7./3.)*eps - 2.*q(1)*eta- 2.*q(2)*xi) 
          pp_nsm1 = 2.*eta - 6.*eps + 2.*q(1)*xi
          pp_nrunby2 = - xi - 12.*eps**2 + 8.*eps*eta
          pp_At = (2./3./const_pi2)*V*(1.+eps/3.)
          pp_nt = - 2.*eps
       endif
    case("bump")
       pp_scan_scalar = .true.
       pp_scan_tensor = .true.
       pp_As = 1.d0
       pp_At = 1.d0
       pp_nt = 0.d0
       pp_nsm1 = 0.d0
       pp_nrunby2 = 0.d0
       call spw_obtain_power(pp_n, pp_lnps, pp_lnpt, lnkmin, lnkmax, psparams(1), psparams(2), psparams(3), psparams(4), psparams(5))
       pp_lnps = log(pp_lnps)
       pp_lnpt = log(pp_lnpt)
    case("spline0", "spline1", "spline2", "linear0", "linear1", "linear2")  !!spline is the same as spline2
       select case(trim(pp_mode))
       case("spline0", "linear0")
          pp_nrunby2 = 0.
          pp_nsm1 = 0.
       case("spline1", "linear1")
          pp_nrunby2 = 0.
       end select
       if(present(psparams))then
          pp_scan_scalar = .true.
          np = size(psparams)
          allocate(lnk(np))
          call set_uniform(np, lnk, lnkmin, lnkmax)
          select case(trim(pp_Mode))
          case("spline0", "spline1", "spline2")
             allocate(lnp2(np))
             call splines(lnk, psparams, lnp2)
             !$omp parallel do
             do i=1, pp_n
                call splints(lnk, psparams, lnp2, pp_lnk(i), pp_lnps(i))
             enddo
             !$omp end parallel do
             deallocate(lnp2)
          case default
             !$omp parallel do
             do i=1, pp_n
                call linearinterpolate(lnk, psparams, pp_lnk(i), pp_lnps(i))
             enddo
             !$omp end parallel do
          end select
          deallocate(lnk)
       else
          pp_scan_scalar = .false.
       endif
       if(present(ptparams))then
          pp_scan_tensor = .true.
          np = size(ptparams)
          allocate(lnk(np), lnp2(np))
          call set_uniform(np, lnk, lnkmin, lnkmax)
          call splines(lnk, ptparams, lnp2)
          do i=1, pp_n
             call splints(lnk, ptparams, lnp2, pp_lnk(i), pp_lnpt(i))
          enddo
          deallocate(lnk, lnp2)
       else
          pp_scan_tensor = .false.
       endif
    case default
       write(*,*) trim(pp_mode)
       stop "UNknown parametrization of primordial power spectrum"
    end select
    if(pp_scan_scalar) &
         call splines(pp_lnk, pp_lnps, pp_lnps2)
    if(pp_scan_tensor) &
         call splines(pp_lnk, pp_lnpt, pp_lnpt2)
    if(present(tmp_mode))then
       pp_mode = pp_mode_save
    endif

  end subroutine pp_init_d


 subroutine pp_init_s(As, ns, nrun, r, nt, lnkmin, lnkmax, psparams, ptparams)
    real As, ns, nrun, r, nt
    real lnkmin, lnkmax
    real,dimension(:),intent(IN),optional::psparams, ptparams
    if(present(psparams) .and. present(ptparams))then
       call pp_init_d(dble(As),dble(ns), dble(nrun), dble(r), dble(nt), dble(lnkmin), dble(lnkmax), dble(psparams), dble(ptparams))
    elseif(present(psparams))then
       call pp_init_d(dble(As),dble(ns), dble(nrun), dble(r), dble(nt), dble(lnkmin), dble(lnkmax), dble(psparams))
    elseif(present(ptparams))then
       call pp_init_d(dble(As),dble(ns), dble(nrun), dble(r), dble(nt), dble(lnkmin), dble(lnkmax), dble(ptparams))
    else
       call pp_init_d(dble(As),dble(ns), dble(nrun), dble(r), dble(nt), dble(lnkmin), dble(lnkmax))
    endif
  end subroutine pp_init_s

  function pp_scalar_power_s(lnk) result(pw)
    real(dl) lnk, pw
    if(pp_scan_scalar)then
       call splints(pp_lnk_min, pp_dlnk, pp_n, pp_lnps, pp_lnps2, lnk, pw)
       pw = pp_As*exp(pw+(pp_nsm1 + pp_nrunby2 *lnk)*lnk)
    else
       pw = pp_As*exp((pp_nsm1 + pp_nrunby2 *lnk)*lnk)
    endif
  end function pp_scalar_power_s


  function pp_scalar_power_v(lnk) result(pw)
    real(dl),dimension(:),intent(IN):: lnk
    real(dl) pw(size(lnk))
    integer i
    !$omp parallel do
    do i=1, size(lnk)
       pw(i) = pp_scalar_power_s(lnk(i))
    enddo
    !$omp end parallel do
  end function pp_scalar_power_v



  function pp_tensor_power_s(lnk) result(pw)
    real(dl) lnk, pw
    if(pp_scan_tensor)then
       call splints(pp_lnk_min, pp_dlnk, pp_n, pp_lnpt, pp_lnpt2, lnk, pw)
       pw = pp_At * exp(pw + pp_nt*lnk)
    else
       pw = pp_At * exp(pp_nt*lnk)
    endif
  end function pp_tensor_power_s


  function pp_tensor_power_v(lnk) result(pw)
    real(dl),dimension(:),intent(IN):: lnk
    real(dl) pw(size(lnk))
    integer i
    !$omp parallel do
    do i=1, size(lnk)
       pw(i) = pp_tensor_power_s(lnk(i))
    enddo
    !$omp end parallel do
  end function pp_tensor_power_v



  subroutine pp_get_potential()
    integer i
    !$omp parallel do
    do i=1, pp_n
       pp_epsilon(i) = pp_tensor_power(pp_lnk(i))/pp_scalar_power(pp_lnk(i))/16.
    enddo
    !$omp end parallel do
    pp_phi(1) = 0
    do i=2, pp_n
       pp_phi(i) = pp_phi(i-1) + (sqrt(pp_epsilon(i))+sqrt(pp_epsilon(i-1)))
    enddo
    pp_phi = (pp_phi - pp_phi(pp_ipivot))*pp_dlnk/const_sqrt2 
    pp_lnV(1) = 0
    do i=2, pp_n
       pp_lnV(i) = pp_lnV(i-1) - (pp_epsilon(i)+pp_epsilon(i-1))
    enddo
    pp_lnV = (pp_lnV - pp_lnV(pp_ipivot))*pp_dlnk
  end subroutine pp_get_potential


!!=============================== models ==========

  function pp_higgs_phi(nefolds, xi, alpha) result(phi_n)
    real(dl),parameter::step = 2.d-3
    real(dl) nefolds, xi, alpha, phi_n
    real(dl) nn, dn1, dn2, dnm, rat, nnp, phip
    nn = 0.d0
    nnp = 0.d0
    phi_n = 0.d0
    dn1 = dndphi(phi_n)
    rat = nefolds*6.d0/step
    do while(nn.lt. rat)
       phi_n = phi_n + step
       dn2 = dndphi(phi_n)
       dnm = dndphi(phi_n - (step/2.d0))
       nnp = nn
       nn = nn + (dn1+4.d0*dnm+dn2)
       dn1 = dn2
    enddo
    phip = phi_n-step
    dn1 = dndphi(phip)
    phi_n = phi_n - step*min((nn - rat)/(nn - nnp), 1.d0)
    dn2 = dndphi(phi_n)
    dnm = dndphi((phi_n +  phip)/2.d0)
    nnp = nnp*step/6.d0
    nn = nnp + (dn1+4.d0*dnm+dn2)*(phi_n - phip)/6.d0
    phi_n = phi_n + (nefolds - nn)/dn2
  contains

    function dndphi(phi)
      real(dl) phi, dndphi
      real(dl) xiphi2, logterm
      xiphi2 = xi * phi ** 2
      logterm = alpha*log(1.d0+xiphi2)
      dndphi = phi/2.d0 * ( 1.d0 + 6.d0*xi*xiphi2/(1.d0+xiphi2))*(1.d0+logterm)/(2.d0+alpha*xiphi2 + 2.d0*logterm)
    end function dndphi

  end function pp_higgs_phi

  function pp_higgs_chi(phi, xi) result(chi)
    real(dl) phi , chi, xi, t
    t = sqrt(xi*(1.d0+6.d0*xi))
    chi = (t * asinh( t*phi) - const_sqrt6*xi*atanh(const_sqrt6*phi*xi/sqrt(1.d0+(t*phi)**2)))/xi
  end function pp_higgs_chi

  function pp_higgs_potential(phi, lambda, xi, alpha) result(V)
    real(dl) V, phi, lambda, xi, alpha
    V = (lambda/4.d0) * phi**4 /(1.d0 + xi*phi**2)**2 *(1.d0+alpha*log(1.d0+xi*phi**2))
  end function pp_higgs_potential

  function pp_higgs_epsilon(phi, xi, alpha) result(eps)
    real(dl) phi, xi, alpha, eps, xiphi2, logterm
    xiphi2  = xi*phi**2
    logterm = log(1.d0+xiphi2)*alpha
    eps = 2.d0*(2.d0+alpha*xiphi2 + 2.d0*logterm)**2/phi**2/(1.d0+xiphi2+6.d0*(xi*phi)**2)/(1.d0+logterm)**2
  end function pp_higgs_epsilon

  function pp_higgs_eta(phi, xi, alpha) result(eta)
    real(dl) phi, xi, alpha, eta, xiphi2, logterm
    xiphi2  = xi*phi**2
    logterm = log(1.d0+xiphi2)*alpha
    eta = 2 *xi * (1.d0 +6.d0*xi * xiphi2 + xi* (-6+phi**2))*(2+ alpha*xiphi2 + 2 * logterm )/((1+ xiphi2 + 6*(xi*phi)**2)**2 * (1+logterm)) - (2*(-6+3*(2 -3*alpha) *xiphi2  + alpha * xiphi2**2 + 6 * (-1 + xiphi2) *logterm ))/(phi**2*(1+xiphi2 +6 *xi *xiphi2 )*(1+ logterm))
  end function pp_higgs_eta
  

  function pp_higgs_xiv(phi, xi, alpha) result(xiv)
    real(dl) phi, xi, alpha, xiv, xiphi2, logterm
    xiphi2  = xi*phi**2
    logterm = log(1.d0+xiphi2)*alpha
    xiv = (8*(2 + alpha*xiphi2 + 2 * logterm)*(6 + 8 * (-1 + 3 * alpha)*xiphi2 -144*(-1 + 3* alpha)*xi**6* phi**8 - 6*xi**5* phi**6* (72 - 72*alpha - 8* phi**2 + 25*alpha *phi**2) + xi**2*phi**2*(12 + 5*(-6 + 7*alpha)*phi**2) + xi**4*phi**6*(alpha*(54 - 13*phi**2) + 4*(-36 + phi**2)) - 2*xi**3*phi**4*(alpha*(-102 + phi**2) + 6*(15 + phi**2)) + 2*alpha*(3 - 4*xi*phi**2 + 72*xi**6*phi**8 + 2*xi**4*phi**6*(-36 + phi**2) + 24*xi**5*phi**6*(-9 + phi**2) - 6*xi**3*phi**4*(15 + phi**2) + xi**2*(6*phi**2 - 15*phi**4))* log(1 + xiphi2)))/((phi + xi*(1 + 6*xi)*phi**3)**4*(1 + logterm)**2)
  end function pp_higgs_xiv


 
end module cosmolibwrap



