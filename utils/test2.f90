!!cosmology with k-essence dark energy and massive neutrino
!!by Zhiqi Huang, huangzhq25@mail.sysu.edu.cn 
!!use c/H_0 as the length unit, 1/H_0 the time unit, and \sqrt{3}M_p as the mass unit

#define NUM_VARS 5
#define BG_PHI y(1)
#define DOT_BG_PHI y(2)
#define HASQ y(3)
#define DELTA_PHI y(4)
#define DOT_DELTA_PHI y(5)


module kessence
  use coop_wrapper_utils
  implicit none
  !! ---------------- configuration ---------------------
  real*8, parameter:: rand_width = 3.d0  
  integer, parameter::dimX = 10 !expand to X^{dimX-1}, dimX>=2 
  integer, parameter::dimV = 10  !expand to phi^{dimV-1}, dimV>=1
  real*8, parameter:: khMpc_pivot = 0.05  !! k/h * Mpc
  real*8, parameter:: fit_amin = 0.6667d0 !!z=0.5
  !cold dark matter + baryon; this does not include massive neutrinos
  real*8, parameter::Omega_m = 0.31
  real*8, parameter:: min_cs2 = 0.d0, max_w = -1.d0/3., max_growth = 100.d0, frozen_cut  = 0.01d0
  real*8,parameter::rho0_min_ini = 0.001, rho0_max_ini = 1000.
  !! ----------------end of configuration ---------------------
  logical::check_physical = .true.
  real*8, parameter:: ksq_pivot = (khMpc_pivot * 3000.d0)**2  !! k / H_0
  real*8, parameter::zero_cut = 1.d-20
  integer::verbose = 2

  real*8:: V( 0:dimV-1, 0:dimX-1)
  real*8:: XC(0:dimX-1), dXC(0:dimX-1), d2XC(0:dimX-1)

  !spatial curvature
  real*8, parameter::Omega_k = 0.
  !effective number of species 
  real*8, parameter::N_eff = 3.046 
  !masses of species of neutrinos in eV
  real*8, dimension(3):: mass_nu =  (/ 0.001d0, 0.009d0, 0.05d0 /)
  !Hubble constant in km/s/Mpc
  real*8, parameter::H0 = 70. 
  !CMB temperature
  real*8, parameter::T_CMB = 2.726
  real*8, parameter::Riemannzeta3 = 1.2020569031595942853997, Riemannzeta5  = 1.0369277551433699263313, Riemannzeta7  = 1.0083492773819228268397, Riemannzeta9 = 1.00200839282608221441785
  real*8, parameter:: pi = 3.1415926535897932384d0

  !derived parameters
  !photon density parameter
  real*8,parameter::Omega_gamma = 2.4748d-5/(H0/100.)**2*(T_CMB/2.726)**4
  !neutrino density parameter
  real*8,parameter::eV_in_Kelvin = 11604.505
  real*8,parameter::T_CNB = T_CMB*(4./11.)**(1./3.)

  real*8::Omega_nu_reference
  real*8,dimension(3)::Omega_nu
  real*8::Omega_L
  integer::index_start, n_used
  real*8,parameter::a_ini = 1.d-4
  real*8,parameter::a_early = 1./1090.d0
  integer,parameter::nsteps = 5000
  real*8,parameter::dlna = -log(a_ini)/(nsteps-1)
  real*8::lna_list(nsteps), a_list(nsteps), ysave(5, nsteps), w(nsteps)

contains

  function poly_eval(coef, phi) result(s)
    real*8 coef(:), phi
    real*8 s, f
    integer i
    s = coef(1)
    f = 1.d0
    do i = 2, size(coef)
       f = f * phi/(i-1)
       s = s + coef(i) * f
    enddo
  end function poly_eval


  function poly_derv(coef, phi) result(s)
    real*8 coef(:), phi
    real*8 s, f
    integer i
    s = 0.d0
    f = 1.d0
    do i = 2, size(coef)
       s = s+ coef(i) * f
       f = f * phi/(i-1)        
    enddo
  end function poly_derv

  function poly_derv2(coef, phi) result(s)
    real*8 coef(:), phi
    real*8 s, f
    integer i
    s = 0.d0
    f = 1.d0
    do i = 3, size(coef)
       s = s + coef(i)*f
       f = f* (phi/(i-2))
    enddo
  end function poly_derv2

  function poly_derv3(coef, phi) result(s)
    real*8 coef(:), phi
    real*8 s, f
    integer i
    s = 0.d0
    f = 1.d0
    do i = 4, size(coef)
       s = s + coef(i)*f
       f = f* (phi/(i-3))
    enddo
  end function poly_derv3


  subroutine kessence_prepare(phi)
    real*8 phi
    integer i
    do i=0, dimX-1
       XC(i) = poly_eval(V(:, i), phi)
       dXC(i) = poly_derv(V(:, i), phi)
       d2XC(i) = poly_derv2(V(:, i), phi)
    enddo
  end subroutine kessence_prepare


  function kessence_rho(phi, X) result(rho)
    real*8::phi, X, p, rho, dpdX
    integer i
    call kessence_prepare(phi)        
    p = poly_eval(XC, X)
    dpdX = poly_derv(XC, X)
    rho = 2*X*dpdX - p
  end function kessence_rho

  function kessence_w(phi, X) result(w)
    real*8::phi, X, w, p, rho, dpdX
    call kessence_prepare(phi)            
    p = poly_eval(XC, X)
    dpdX = poly_derv(XC, X)
    rho = 2*X*dpdX - p
    w = p/rho
  end function kessence_w

  subroutine kessence_vars_all(phi, X, p, rho, dpdX, drhodX, drhodphi, S, dSdX, dSdphi)
    real*8::phi, X, p, rho, dpdX, dpdphi, drhodX, d2pdX2, d2pdXdphi, drhodphi, dSdphi, dSdX, S
    real*8::d2pdphi2, d2rhodphi2, d2rhodXdphi, d2rhodX2, d3pdX2dphi, d3pdXdphi2, d3pdX3
    call kessence_prepare(phi)        
    dpdX = poly_derv(XC, X)
    d2pdX2 = poly_derv2(XC, X)
    drhodX = dpdX + 2*X*d2pdX2
    if(abs(drhodX) < zero_cut)then
       drhodX = 0.d0
       return
    else
       p = poly_eval(XC, X)       
       dpdphi = poly_eval(dXC, X)                
       d2pdXdphi = poly_derv(dXC, X)     
       rho = 2*X*dpdX - p
       drhodphi = 2*X*d2pdXdphi - dpdphi     
       S = drhodphi / drhodX
       if(check_physical)then
          d2pdphi2 = poly_eval(d2XC, X) 
          d3pdXdphi2 = poly_derv(d2XC, X) 
          d3pdX2dphi = poly_derv2(dXC, X)
          d3pdX3 = poly_derv3(XC, X)
          d2rhodphi2 = 2*X*d3pdXdphi2 - d2pdphi2 
          d2rhodXdphi = d2pdXdphi + 2*X*d3pdX2dphi
          d2rhodX2 = 3*d2pdX2 + 2*X*d3pdX3
          dSdphi = (d2rhodphi2 - d2rhodXdphi*S)/drhodX
          dSdX = (d2rhodXdphi - d2rhodX2 * S) / drhodX
       endif
    end if
  end subroutine kessence_vars_all

  function cheb_value(x, c) result(y)
    real*8 c(:), x, y
    real*8 twox, y1, y2, y3
    integer i
    twox = x+x
    y3 = 0.d0
    y2 = 0.d0
    do i = size(c), 2, -1
       y1 = c(i) + twox * y2 - y3
       y3 = y2
       y2 = y1
    enddo
    y = c(1) + x* y2 - y3    
  end function cheb_value



  function I_rho(am)
    real*8:: am, I_rho
    real*8::lnam, invam2
    if(am < 0.101)then
       I_rho = 7./240.*pi**2 + am**2/48.-am**4*9.d-3
       return
    endif
    lnam = log(am)
    if(lnam > 3.99)then
       invam2 = 1./am**2
       I_rho = 7./240.*pi**2 * exp( lnam + log((Riemannzeta3 * 180/7/pi**4) +( (1350/7*Riemannzeta5/pi**4) + (-6075/4/pi**4 * Riemannzeta7+(172125* Riemannzeta9/8/ pi**4)*invam2)*invam2)*invam2))
    else
       I_rho = 7./240.*pi**2 * exp(cheb_value((lnam - 0.75d0)/3.25d0, &
            (/ 0.89293935180142925d0, &
            1.3851588735560008d0,              0.59201638940401391d0,              5.52167768753157942d-2, &
            -6.84994561924834600d-2, -1.78898251764955940d-2, 1.31790803821667645d-2, 5.84725640687740623d-3, &
            -2.65218714293024649d-3, -1.84616659288402659d-3, 4.82404659741137831d-4, 5.58642583176972569d-4, &
            -6.26401384210129567d-5, -1.61520019162570550d-4, -3.30934261743216566d-6, 4.47917623670222794d-5, &
            6.15434233659787645d-6, -1.15506405054980366d-5, -3.87591649596476188d-6, 3.18654185509784438d-6 /) ))
    endif
  end function I_rho


  function I_p(am)
    real*8::am, I_p
    real*8::lnam, invam2
    if( am < 0.101)then
       I_p = 7./720.*pi**2 - am**2/144. + am**4*8.d-3
       return
    endif
    lnam = log(am)
    if(lnam > 3.99)then
       invam2 = 1./am**2
       I_p = 7./240.*pi**2*((900/7/pi**4*Riemannzeta5) + ((-2025*Riemannzeta7/pi**4) +(172125*Riemannzeta9/2/pi**4)*invam2)*invam2)/am
    else
       I_p = 7./240.*pi**2 * exp(cheb_value((lnam - 0.75d0)/3.25d0, &
            (/ -1.8856808011378832d0, -1.2325642746845389d0, -0.55304532292355002d0, &
            -8.04535652389012507d-2, 4.86160741219811912d-2, 2.10996354221542545d-2, &
            -5.31380099138448365d-3, -5.02522420384828548d-3, 1.30786399004496982d-4, &
            1.00515162172776272d-3, 1.70018903096916962d-4, -1.53071349645732708d-4, &
            -6.33387811490677483d-5, 1.22421763460544607d-5, 1.36045554994012240d-5, &
            1.73976481843332439d-6, -1.66578804580377724d-6, -8.83571732180485515d-7, &
            -1.15092657356077256d-7, 1.54021206255242797d-7 /) ))
    end if
  end function I_p


  !rho_nu/rho_critical  as a function of z
  function rho_nu_scaled(a) result(s)
    real*8::a, s
    integer i
    s = 0.    
    do i = 1, size(mass_nu)
       s = s +  I_rho(mass_nu(i) * eV_in_Kelvin * a / T_CNB)
    enddo
    s = s * Omega_nu_reference
  end function rho_nu_scaled


  !rho_nu/rho_critical as a function of a
  function p_nu_scaled(a) result(s)
    real*8 a, s
    integer i
    s = 0.    
    do i=1, size(mass_nu)       
       s =  s + I_p(mass_nu(i)*eV_in_Kelvin * a/T_CNB)
    enddo
    s = s * Omega_nu_reference
  end function p_nu_scaled


  function field_eq(y, lna, yp) result(fail) !y = [\phi, \dot \phi, H]
    logical fail
    real*8:: lna, y(NUM_VARS), yp(NUM_VARS)
    real*8:: cs2, X, p_DE, rho_DE, dpdX, drhodX, drhodphi, dSdX, dSdphi, ddotphi, dHasqdlna, S, ddot_delta_phi, a, delta_rho, H, asq
    X = DOT_BG_PHI **2 / 2.
    call  kessence_vars_all(phi = BG_PHI, X = X, p = p_DE, rho = rho_DE, dpdX = dpdX, drhodX = drhodX, drhodphi = drhodphi, S = S, dSdX = dSdX, dSdphi = dSdphi)
    if(abs(drhodX) <= zero_cut)then  !non-singularity
       if(verbose>2) print*, drhodX, exp(lna),": drhodX reject"                     
       fail = .true.
       return
    endif
    cs2 = dpdX/drhodX
    if(rho_DE <= zero_cut)then
       if(verbose>2) print*, rho_DE, exp(lna),": rho_DE reject"              
       fail = .true.
       return
    endif
    if(check_physical)then
       delta_rho = drhodphi * DELTA_PHI + drhodX * DOT_BG_PHI * DOT_DELTA_PHI       
       if(abs(delta_rho) > max_growth * rho_DE)then
          if(verbose>1) print*, delta_rho/rho_DE, exp(lna), ": growth reject"       
          fail = .true.
          return
       endif
       if(p_DE >= max_w * rho_DE)then  !truncate equation of state
          if(verbose>1) print*, p_DE/rho_DE, exp(lna), ": w reject"              
          fail = .true.
          return
       endif
       if(cs2 < min_cs2)then !truncate cs^2
          if(verbose>1) print*, cs2, exp(lna),": cs^2 reject"                     
          fail = .true.
          return
       endif
    endif
    a = exp(lna)
    asq = a**2
    H = HASQ / asq
    ddotphi = -3.* H * cs2 * DOT_BG_PHI - S
    dHasqdlna = ((rho_DE - 3.*p_DE)*asq + (Omega_m + (rho_nu_scaled(a)-3.d0 * p_nu_scaled(a))/a)/a) /(2.d0 * H)
    if(check_physical) then
       ddot_delta_phi = -(3.d0*H * cs2 + dSdX * DOT_BG_PHI) * DOT_DELTA_PHI - (dSdphi + cs2/a**2*ksq_pivot)*DELTA_PHI
       yp = (/ DOT_BG_PHI/H, ddotphi/H, dHasqdlna,  DOT_DELTA_PHI / H, ddot_delta_phi / H /)
    else
       yp(1:3) = (/ DOT_BG_PHI/H, ddotphi/H, dHasqdlna /)
    endif
    fail = .false.
  end function field_eq


  function evolve_step(y, lna, dlna) result(fail) !4th order Runge Kutta Integrator
    logical fail
    real*8 lna, dlna
    real*8,dimension(NUM_VARS)::y, dy1, dy2, dy3, dy4
    if(field_eq(y, lna, dy1))then
       fail = .true.
       return
    endif
    if(field_eq(y+dy1*dlna/2., lna+dlna/2., dy2))then
       fail = .true.
       return
    endif
    if(field_eq(y+dy2*dlna/2., lna+dlna/2., dy3))then
       fail = .true.
       return
    endif
    if(field_eq(y+dy3*dlna, lna+dlna, dy4))then
       fail = .true.
       return
    endif
    if(check_physical)then
       y = y + (dy1 + 2*(dy2+dy3)+dy4)*(dlna/6.)
    else
       y(1:3) = y(1:3) + (dy1(1:3) + 2*(dy2(1:3)+dy3(1:3))+dy4(1:3))*(dlna/6.)       
    endif
    fail = .false.
  end function evolve_step


  subroutine initialize()
    integer::i
    Omega_nu_reference = 0.56202d-5/(H0/100.)**2 * N_eff/size(mass_nu)/I_rho(0.d0) !this is, Omega_nu/I_rho(0) assuming massless
    do i=1, 3
       Omega_nu(i) = I_rho(mass_nu(i)*eV_in_Kelvin/T_CNB) * Omega_nu_reference
    end do
    Omega_L = 1.0 - Omega_m - Omega_k - Omega_gamma - sum(Omega_nu)
    lna_list(1) = log(a_ini)
    lna_list(nsteps) = 0.d0
    do i=2,nsteps - 1
       lna_list(i) = lna_list(1) + dlna*(i-1)
    enddo
    a_list = exp(lna_list)
    n_used = nint(-log(fit_amin)/dlna)+1  
    index_start = nsteps - n_used + 1
  end subroutine initialize

  function set_ic(y, a) result(fail)
    logical::fail
    real*8 y(NUM_VARS), a
    real*8 p_DE, rho_DE, dpdX, drhodX, drhodphi, S, dSdX, dSdphi
    BG_PHI = 0.d0
    DOT_BG_PHI = 0.d0
    call kessence_prepare(BG_PHI)
    call kessence_vars_all(BG_PHI, DOT_BG_PHI ** 2 / 2.d0, p_DE, rho_DE, dpdX, drhodX, drhodphi, S, dSdX, dSdphi)
    if( abs(drhodX) < zero_cut .or. rho_DE < zero_cut)then
       fail = .true.
       return
    endif
    HASQ = sqrt( (Omega_k*a + Omega_m)*a + Omega_gamma + rho_nu_scaled(a) + rho_DE * a**4 )
    if(abs(drhodphi) > zero_cut )then   !!set initial delta \rho / \rho = 1
       DELTA_PHI = rho_DE / drhodphi
       DOT_DELTA_PHI  = 0.d0
    else  !!just do some random initial conditions
       DELTA_PHI = 1.d0
       DOT_DELTA_PHI = 0.d0
    endif
    fail = .false.
  end function set_ic


  function  evolve_all() result(fail)
    logical fail
    real*8 y(NUM_VARS)
    integer i
    check_physical = .true.    
    if(set_ic(y, a_ini))then
       fail = .true.
       return
    endif
    ysave(:, 1) = y
    do i = 2, nsteps
       if(evolve_step(y, lna_list(i-1), dlna))then
          if(verbose > 1)print*, "evolve all step reject"
          fail = .true.
          return
       endif
       ysave(:, i) = y
       if(ysave(3, i) < zero_cut)then
          if(verbose > 1)print*, "****************negative Hasq reject********************"
          fail = .true.
          return
       endif
    enddo
    fail = .false.
  end function evolve_all


  function  evolve_from(a, y) result(fail)
    logical fail
    real*8 y(NUM_VARS), lastH, a, step, lna
    integer i, n
    check_physical = .false.        
    if(set_ic(y, a))then
       if(verbose > 1) print*, "evolve_from set_ic reject"
       fail = .true.
       return
    endif
    lna = log(a)
    n = -lna/dlna
    step = -lna/n
    do i = 1, n
       if(evolve_step(y, lna, step))then
          if(verbose > 1) print*, "evolve step reject"
          fail = .true.
          return
       endif
       if(HASQ < zero_cut )then
          if(verbose > 1) print*, a, HASQ, "*********************negative Hasq reject***************************"
          fail = .true.
          return
       endif
       lna = lna + step        
    enddo
    fail = .false.
  end function evolve_from


  function solve_V00_scan(a, rho0_min, rho0_max, nscan) result(fail)
    logical fail
    integer nscan
    real*8 a
    real*8 y(NUM_VARS),  rho0_min, rho0_max  !!rho0 = - V00
    integer:: i, veb, best_index
    real*8::H_try(0:nscan), rho0_try(0:nscan+1)
    veb = verbose
    verbose = 0
    H_try = 0.d0
    best_index = 0
    call coop_set_uniform(nscan, rho0_try(1:nscan), rho0_min, rho0_max,  .true.)
    rho0_try(0) = rho0_min
    rho0_try(nscan+1) = rho0_max
    do i=1, nscan
       V(0, 0) =  - rho0_try(i)
       if(.not. evolve_from(a, y))then
          H_try(i) = HASQ-1.d0
          if(abs(H_try(i)) <= zero_cut)then
             rho0_min  = rho0_try(i)
             rho0_max = rho0_try(i)
             verbose = veb
             fail = .false.             
             return
          endif
          if( H_try(i-1) <= -zero_cut  .and. H_try(i) >= zero_cut)then  !!found solution
             rho0_min = rho0_try(i-1)
             rho0_max = rho0_try(i)
             V(0, 0) = -sqrt(rho0_min*rho0_max)
             verbose = veb
             fail = .false.
             return
          endif
          if(abs(H_try(i)) < H_try(best_index) ) &
               best_index = i
       endif
    enddo
    if(best_index == 0)then
       fail = .true.
       if(verbose > 1) print*, "solve_V00_scan fail"                    
    else
       V(0, 0) = -rho0_try(best_index)
       rho0_min = rho0_try(best_index-1)
       rho0_max = rho0_try(best_index+1)
    endif
    verbose = veb
  end function solve_V00_scan


  function solve_V00_binary_search(a, rho0_min, rho0_max) result(fail)
    real*8,parameter::onep = 1.d0+1.d-5
    real*8 a
    logical fail
    real*8 rho0_min, rho0_max, y(NUM_VARS), Hmax, Hmin
    integer:: i
    if(rho0_max / rho0_min < onep)then
       fail = .false.
       V(0, 0) = -sqrt(rho0_max * rho0_min)
       return
    endif
    V(0, 0) = -rho0_min
    do while(evolve_from(a, y))
       if(rho0_max / rho0_min < onep)then
          if(verbose > 1) print*, "no viable initial condition"                              
          fail = .true.
          return
       else
          rho0_min = rho0_min ** 0.98 * rho0_max ** 0.02
          V(0, 0) = -rho0_min
       endif
    enddo
    Hmin = HASQ
    if(abs(Hmin-1.d0) < zero_cut)then
       fail = .false.
       rho0_max = rho0_min
       return
    endif
    
    V(0, 0) = -rho0_max
    do while(evolve_from(a, y))
       if(rho0_max / rho0_min < onep)then
          if(verbose > 1) print*, "no viable initial condition"                              
          fail = .true.
          return
       else
          rho0_max = rho0_max ** 0.98 * rho0_min ** 0.02
          V(0, 0) = -rho0_max
       endif
    enddo
    Hmax = HASQ
    if(abs(Hmax-1.d0) < zero_cut)then
       fail = .false.
       rho0_min = rho0_max
       return
    endif
    
    if(Hmax > 1.d0 .and. Hmin < 1.d0)then
       do while(rho0_max/rho0_min > onep)
          V(0, 0) = -sqrt(rho0_max * rho0_min)
          if(evolve_from(a, y))then
             if(verbose > 1) print*, rho0_min, rho0_max, V(0, 0), "************evolve_from binary search reject********"
             fail = .true.
             return
          endif
          if(HASQ > 1.d0)then
             rho0_max = -V(0, 0)
             Hmax =  HASQ
          else
             rho0_min = -V(0, 0)
             Hmin = HASQ
          endif
       enddo
       if(abs(HASQ - 1.d0) > 2.d-3)then
          if(verbose>1) print*, HASQ, "***************binary search H reject******************"        
          fail = .true.
          return
       endif
    else
       fail = solve_V00_scan(a, rho0_min, rho0_max, 100)
       return
    endif
    fail = .false.
  end function solve_V00_binary_search


  function get_w() result(fail)
    real*8,parameter::onep = 1.d0+1.d-5    
    logical fail
    real*8::rho0, a, rho0_min, rho0_max
    integer i
    rho0_min = rho0_min_ini
    rho0_max = rho0_max_ini
    fail =  solve_V00_scan(a_ini, rho0_min, rho0_max, 100)
    if(fail) return
    do while(rho0_max / rho0_min > onep)
       fail = solve_V00_binary_search(a_ini, rho0_min, rho0_max)
       if(fail) return
    enddo
    V(0, 0) = -sqrt(rho0_min * rho0_max)
    if(evolve_all())then
       if(verbose > 1) print*, "evolve_all reject"        
       fail = .true.
       return
    endif
    rho0 = kessence_rho(ysave(1, nsteps),  ysave(2, nsteps)**2/2.d0)
    if( abs(rho0 - Omega_L) > 1.d-4)then !!
       print*, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"          
       print*, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
       print*,rho0, Omega_L, "!!!!Warning: energy conservation failed!!!!"
       print*, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
       print*, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"          
       fail = .true.
       return
    endif
    do i=1, nsteps
       w(i) = kessence_w(ysave(1, i), ysave(2, i)**2/2.d0)
       if(a_list(i) < a_early .and. abs(1.+w(i)) > frozen_cut)then
          if(verbose > 1) print*, "thawing assumption reject"        
          fail = .true.
          return          
       endif
    enddo
    fail = .false.
  end function get_w


  subroutine linearfit_w0wa(n, a, w, w0, wa)
    integer::n, i
    real*8 a(n), w(n), w0, wa, weights(n)
    real*8 mean_a, mean_w, mean_aw, mean_a2
    do i=2, n-1
       weights(i) = (a(i+1)-a(i-1))/2.d0
    enddo
    weights(1) = a(2)-a(1)
    weights(n) = a(n)-a(n-1)
    weights = weights / sum(weights)
    mean_a = sum(a*weights)
    mean_w = sum(w*weights)
    mean_aw = sum((a-mean_a)*(w-mean_w)*weights)
    mean_a2 = sum((a - mean_a)**2*weights)
    wa = -mean_aw/mean_a2
    w0 = mean_w - wa*(1.-mean_a)
  end subroutine linearfit_w0wa

  subroutine get_Gaussian(r)
    real*8 r, v(2)
    call random_number(v)
    v = 2.d0*v - 1.d0
    r=v(1)**2 + v(2)**2
    do while (R >= 1.d0)
       call random_number(v)
       v = 2.d0*v - 1.d0
       r=v(1)**2 + v(2)**2
    end do
    r = dsqrt(-2.d0*dlog(r)/r)*v(1)
  end subroutine get_Gaussian


  subroutine random_Gaussian_matrix()
    integer i, j
    do i=0, dimX-1
       do j=0, dimV-1
          call get_Gaussian(V(j, i))
       enddo
    enddo
  end subroutine random_Gaussian_matrix

  function kessence_wa(w0, Omega_m)result(wa)
    real*8::w0, Omega_m, wa
    wa = (1.d0+w0)* (-1.42*(Omega_m/0.3)**0.64) !wrat(((1.d0-Omega_m)/Omega_m)**(1.d0/3.d0))
  end function kessence_wa

  function kessence_F(x) result(F)
    real*8::x, F, y
    y = x**3
    if(x > 0.03d0)then
       F = sqrt(1.d0+1.d0/y) - log(sqrt(y) + sqrt(1.d0+y))/y
    else
       F = ((2.d0/3.d0) + y*(- 0.2d0+y*(3.d0/28.d0)))*sqrt(y)
    endif
  end function kessence_F


end module kessence


program main
  use kessence
  use coop_wrapper_utils
  implicit none
#include "constants.h"  
  integer i
  integer,parameter::want_samples = 5000
  real*8::V_save(0:dimV-1, 0:dimX-1, want_samples), w0(want_samples), wa(want_samples)
  type(coop_file)::Vdata
  type(coop_asy)::fig
  integer,parameter::plot_nx = 21, plot_ny = 21
  real*8::z(plot_nx, plot_ny)
  integer ncount
  type(coop_file)::output
  call initialize()
  if(.not. coop_file_exists("Vdata.dat"))then
     call Vdata%open("Vdata.dat", 'u')     
     ncount = 0
     do while(ncount < want_samples)
        call random_Gaussian_matrix()
        V = V*rand_width
        if(get_w())then
           if(verbose>0) print*, "...sample rejected!"
        else
           ncount = ncount + 1
           call linearfit_w0wa(n_used, a_list(index_start:nsteps), w(index_start:nsteps), w0(ncount), wa(ncount))
           V_save(:, :, ncount) = V
           if(verbose>0) print*, ncount, " samples done"
        endif
     enddo
     write(Vdata%unit) V_save
     write(Vdata%unit) w0
     write(Vdata%unit) wa      
     call Vdata%close()
  endif
  call Vdata%open("Vdata.dat", "u")
  read(Vdata%unit)V_save
  read(Vdata%unit)w0
  read(Vdata%unit)wa
  call Vdata%close()
  call fig%open("wfigs/Vstat_"//COOP_STR_OF(nint(Omega_m*100))//"_"//COOP_STR_OF(nint(rand_width))//"_"//COOP_STR_OF(dimX)//"_"//COOP_STR_OF(dimV)//".txt")  
  call fig%init(xlabel="$V_{00}$", ylabel="$V_{35}$")
  call fig%dots(x=V_save(0, 0, :), y=V_save(3, 5, :), color="skyblue")
  call fig%close()
  call output%open("stat.txt", "w")
  do i=1, want_samples
     write(output%unit, "(15G14.5)") w0(i), wa(i), V_save(0, 0, i), V_save(1, 0, i), V_save(0, 1, i), V_save(1, 1, i), V_save(3, 5, i)
  enddo
  call output%close()
end program main

