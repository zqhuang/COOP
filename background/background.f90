module coop_background_mod
  use coop_wrapper_utils
  implicit none
#include "constants.h"

  private

#if DO_EFT_DE
  !!screening: true; non-screening: false
  logical,parameter::coop_eft_de_normalize_early = .true.
#endif  
  
  public::coop_baryon, coop_cdm, coop_DE_lambda, coop_DE_w0, coop_DE_w0wa, coop_DE_quintessence, coop_radiation, coop_neutrinos_massless, coop_neutrinos_massive, coop_de_w_quintessence, coop_de_wp1_quintessence, coop_de_wp1_coupled_quintessence, coop_background_add_coupled_DE,  coop_background_add_EFT_DE,  coop_background_add_EFT_DE_with_effective_w, coop_de_aeq_fitting, coop_de_alpha_invh2, coop_de_alpha_instant, coop_de_general, coop_de_alpha_constructor

contains


  function coop_baryon(Omega_b, fcs2b) result(this)
    type(coop_species) this
    type(coop_function),optional::fcs2b
    COOP_REAL Omega_b
    if(present(fcs2b))then
       call this%init(genre = COOP_SPECIES_FLUID, name = "Baryon", id = 1, Omega=Omega_b, w = COOP_REAL_OF(0.), fcs2 = fcs2b )  !!you might want to replace the cs^2(baryon)
    else
       call this%init(genre = COOP_SPECIES_FLUID, name = "Baryon", id = 1, Omega=Omega_b, w = COOP_REAL_OF(0.), cs2 = COOP_REAL_OF(0.))  !!you might want to replace the cs^2(baryon)
    endif
  end function coop_baryon

  function coop_cdm(Omega_c) result(this)
    type(coop_species) this
    COOP_REAL Omega_c
    call this%init(genre = COOP_SPECIES_CDM, name = "CDM", id = 2, Omega=Omega_c, w = COOP_REAL_OF(0.), cs2 = COOP_REAL_OF(0.))
  end function coop_cdm

  function coop_radiation(Omega_r) result(this)
    type(coop_species) this
    COOP_REAL:: Omega_r
    call this%init(genre = COOP_SPECIES_MASSLESS, name = "Radiation", id = 3, Omega=Omega_r) 
  end function coop_radiation

  function coop_neutrinos_massless(Omega_nu) result(this)
    type(coop_species) this
    COOP_REAL:: Omega_nu
    call this%init(genre = COOP_SPECIES_MASSLESS, name = "Massless Neutrinos", id = 4, Omega=Omega_nu)
  end function coop_neutrinos_massless

  function coop_neutrinos_massive(Omega_nu, Omega_massless) result(this)
    type(coop_species) this
    COOP_REAL:: Omega_nu
    COOP_REAL:: Omega_massless
    call this%init(genre = COOP_SPECIES_MASSIVE_FERMION, name = "Massive Neutrinos", id = 4, Omega=Omega_nu, Omega_massless = Omega_massless)
  end function coop_neutrinos_massive

  
  function coop_de_lambda(Omega_Lambda) result(this)
    type(coop_species) this
    COOP_REAL Omega_Lambda
    call this%init(genre = COOP_SPECIES_LAMBDA, name = "Dark Energy",id=5, Omega = Omega_Lambda, w = COOP_REAL_OF(-1.d0), cs2 = COOP_REAL_OF(1.d0))
  end function coop_de_lambda

  function coop_de_general(Omega_Lambda, fwp1) result(this)
    type(coop_species) this
    COOP_REAL Omega_Lambda
    type(coop_function)::fwp1
    call this%init(genre = COOP_SPECIES_FLUID, name = "Dark Energy", id = 5, Omega = Omega_Lambda, cs2 = COOP_REAL_OF(1.d0), fwp1 = fwp1)
  end function coop_de_general
  
  function coop_de_w0wa(Omega_Lambda, w0, wa) result(this)
    type(coop_species) this
    COOP_REAL Omega_Lambda, w0, wa
    type(coop_function) fw0wa
    type(coop_arguments) w0wa
    if(w0 .eq. -1.d0 .and. wa.eq.0.d0)then
       this = coop_de_lambda(Omega_Lambda)
       return
    endif
    w0wa = coop_arguments_constructor(r =  (/ w0, wa /))
    fw0wa = coop_function_constructor(coop_de_wp1_w0wa, xmin = coop_min_scale_factor, xmax = COOP_REAL_OF(1.), xlog = .true., args = w0wa, name = "w0wa model 1+w")
    call this%init(genre = COOP_SPECIES_FLUID, name = "Dark Energy", id=5, Omega = Omega_Lambda, cs2 = COOP_REAL_OF(1.d0), fwp1 = fw0wa )
    call w0wa%free()
    call fw0wa%free()
  end function coop_de_w0wa

  function coop_de_w0(Omega_Lambda, w0) result(this)
    type(coop_species) this
    COOP_REAL w0, Omega_Lambda
    if(w0.eq.-1.d0)then
       this = coop_de_lambda(Omega_Lambda)
       return
    endif
    call this%init(genre = COOP_SPECIES_FLUID, name = "Dark Energy", id = 5, Omega = Omega_Lambda, w = w0, cs2 = COOP_REAL_OF(1.))
  end function coop_de_w0

  function coop_de_wp1_w0wa(a, w0wa) result(w)
    COOP_REAL a, w
    type(coop_arguments) w0wa
    w = 1.d0 + w0wa%r(1) + w0wa%r(2)*(1.d0 - a)
  end function coop_de_wp1_w0wa

  function coop_de_quintessence(Omega_Lambda, epsilon_s, epsilon_inf, zeta_s) result(this)
    type(coop_species) this
    COOP_REAL Omega_Lambda, epsilon_s, epsilon_inf, zeta_s
    type(coop_function) fq
    type(coop_arguments) arg
    if(epsilon_s .eq. 0.d0 .and. epsilon_inf .eq. 0.d0 .and. zeta_s .eq. 0.d0)then
       this = coop_de_lambda(Omega_Lambda)
       return
    endif
    arg = coop_arguments_constructor(r =  (/ Omega_Lambda, epsilon_s, epsilon_inf, zeta_s /))
    fq = coop_function_constructor(coop_de_wp1_quintessence, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args = arg, name = "quintessence 1+w")
    call this%init(genre = COOP_SPECIES_FLUID, name = "Dark Energy", id=5, Omega = Omega_Lambda, cs2 = COOP_REAL_OF(1.d0), fwp1 = fq )
    call arg%free
    call fq%free
  end function coop_de_quintessence

  function coop_de_w_quintessence(a, arg) result(w)
    COOP_REAL a, w
    type(coop_arguments) arg
    w = coop_de_wp1_quintessence(a, arg) - 1.d0
  end function coop_de_w_quintessence

  function coop_de_aeq_fitting(omega_lambda, epsilon_s, epsilon_inf, zeta_s) result(aeq)
    COOP_REAL, parameter::omega_r = 9.d-5 !!approximate value for h = 0.68    
    COOP_REAL::omega_lambda, epsilon_s, epsilon_inf, zeta_s, aeq, delta, qpsign, epss, sqrtepss, sqrtepsinf, Omega_m, diff

    if(epsilon_s .ge. 0.d0)then
       qpsign = 1.d0
       epss = epsilon_s
    else
       qpsign = -1.d0
       epss = - EPSILON_S
    endif
    Omega_m  = 1.d0 - Omega_lambda - omega_r
    sqrtepss = sqrt(epss)
    sqrtepsinf = sqrt(epsilon_inf)
    diff = sqrtepss - coop_sqrt2 * sqrtepsinf
    delta = (sqrtepsinf + (0.91-0.78*Omega_m+(0.236-0.76*Omega_m)*zeta_s)*diff)**2 &
         + (sqrtepsinf + (0.533-0.1*zeta_s)*diff)**2 
    aeq = (Omega_m/Omega_lambda)**(1.d0/(3.d0-qpsign*delta))
  end function coop_de_aeq_fitting

  function coop_de_wp1_quintessence(a, arg) result(wp1)
    COOP_REAL, parameter::omega_r = 9.d-5 !!approximate value for h = 0.68
    COOP_REAL a, wp1
    type(coop_arguments) arg
    COOP_REAL Omega_m, a_eq, epss
#define OMEGA_LAMBDA arg%r(1)
#define EPSILON_S arg%r(2)
#define EPSILON_INF arg%r(3)
#define ZETA_S arg%r(4)
    COOP_REAL::mu, mu3, sqrtepss, sqrtepsinf, diff, delta, f, f2, s0, s1, qpsign
    if(EPSILON_S .ge. 0.d0)then
       qpsign = 1.d0
       epss = EPSILON_S
    else
       qpsign = -1.d0
       epss = - EPSILON_S
    endif
    Omega_m  = 1.d0 - OMEGA_LAMBDA - omega_r
    sqrtepss = sqrt(epss)
    sqrtepsinf = sqrt(EPSILON_INF)
    diff = sqrtepss - coop_sqrt2 * sqrtepsinf
    delta = (sqrtepsinf + (0.91-0.78*Omega_m+(0.236-0.76*Omega_m)*ZETA_S)*diff)**2 &
         + (sqrtepsinf + (0.533-0.1*ZETA_S)*diff)**2 
    a_eq = (Omega_m/OMEGA_LAMBDA)**(1.d0/(3.d0-qpsign*delta))
    mu = a/a_eq
    mu3 = mu**3    
    s0 = sqrt(mu3)
    if(mu .gt. 0.05d0)then
       s1 = sqrt(1.d0+mu3)
       f = s1/s0 - log(s0+s1)/mu3
       f2 = coop_sqrt2*(1.d0-log(1.d0+mu3)/mu3) - f
    else  !!asymptotic
       f = s0*((2.d0/3.d0)-0.2d0*mu3)
       f2 = mu3*(1.d0/coop_sqrt2 - (coop_sqrt2/3.d0)*mu3)  - f
    endif
    s0 = sqrtepsinf*sqrt(((4.d0/3.d0*omega_r) + omega_m*a)/(omega_r+omega_m*a))
    wp1 = (2.d0/3.d0)*qpsign*(s0 + (sqrtepss  - coop_sqrt2*s0)*(f + ZETA_S * f2))**2
#undef EPSILON_S
#undef EPSILON_INF
#undef ZETA_S
#undef OMEGA_LAMBDA
    
  end function coop_de_wp1_quintessence


  function coop_de_wp1_coupled_quintessence(a, arg) result(wp1)
    COOP_REAL, parameter::omega_r = 9.d-5 !!approximate value for h = 0.68
    COOP_REAL a, wp1
    type(coop_arguments) arg
    COOP_REAL Omega_m, a_eq, epss
#define OMEGA_LAMBDA arg%r(1)
#define EPSILON_S arg%r(2)
#define EPSILON_INF arg%r(3)
#define ZETA_S arg%r(4)
#define BETA_S arg%r(5)    
    COOP_REAL::mu, mu3, sqrtepss, sqrtepsinf, diff, delta, f, f2, s0, s1, qpsign
    if(EPSILON_S .ge. 0.d0)then
       qpsign = 1.d0
       epss = EPSILON_S
    else
       qpsign = -1.d0
       epss = - EPSILON_S
    endif
    Omega_m  = 1.d0 - OMEGA_LAMBDA - omega_r
    sqrtepss = sqrt(epss)
    sqrtepsinf = sqrt(EPSILON_INF*a**BETA_S)
    diff = sqrtepss - coop_sqrt2 * sqrtepsinf
    delta = (sqrtepsinf + (0.91-0.78*Omega_m+(0.236-0.76*Omega_m)*ZETA_S)*diff)**2 &
         + (sqrtepsinf + (0.533-0.1*ZETA_S)*diff)**2 
    a_eq = (Omega_m/OMEGA_LAMBDA)**(1.d0/(3.d0-qpsign*delta))
    mu = a/a_eq
    mu3 = mu**3    
    s0 = sqrt(mu3)
    if(mu .gt. 0.05d0)then
       s1 = sqrt(1.d0+mu3)
       f = s1/s0 - log(s0+s1)/mu3
       f2 = coop_sqrt2*(1.d0-log(1.d0+mu3)/mu3) - f
    else  !!asymptotic
       f = s0*((2.d0/3.d0)-0.2d0*mu3)
       f2 = mu3*(1.d0/coop_sqrt2 - (coop_sqrt2/3.d0)*mu3)  - f
    endif
    s0 = sqrtepsinf*sqrt(((4.d0/3.d0*omega_r) + omega_m*a)/(omega_r+omega_m*a))
    wp1 = (2.d0/3.d0)*qpsign*(s0 + (sqrtepss  - coop_sqrt2*s0)*(f + ZETA_S * f2))**2
#undef EPSILON_S
#undef EPSILON_INF
#undef ZETA_S
#undef OMEGA_LAMBDA
    
  end function coop_de_wp1_coupled_quintessence

!!$  subroutine coop_background_add_coupled_DE_negative_powerlaw(this, Omega_c, alpha, Q, err)  !!simple model: constant Q and negative powerlaw potential V = V0 phi^{-alpha}.
!!$    class(coop_cosmology_background)::this
!!$    COOP_REAL::Omega_c
!!$    COOP_INT::err
!!$    COOP_REAL::alpha, Q
!!$#if DO_COUPLED_DE    
!!$    COOP_INT, parameter::ns = 12000
!!$    type(coop_function)::fwp1effcdm, fwp1de, fwp1effde
!!$    type(coop_ode)::ode
!!$    type(coop_species)::de, cdm, deeff
!!$    COOP_REAL::phi_ini
!!$    COOP_INT::i, i_tc_off
!!$    COOP_REAL::rho_ce, a(ns), lna(ns), y(2, ns), yp(2, ns), lnV(ns), dlnVdphi(ns), dVdphibyH2(ns), wp1de(ns), wp1effde(ns), wp1effcdm(ns), H2a4(ns), m2byH2(ns), dQdphi(ns), omde, omc, dlna, tc_w, wp1eff(ns)
!!$    err = 0
!!$    i_tc_off = 1
!!$    if(this%Omega_k() .le. Omega_c)then
!!$       err = -1
!!$       return
!!$    endif
!!$    rho_ce = omega_c * 3.d0    
!!$    call coop_set_uniform(ns, lna, log(coop_min_scale_factor), log(coop_scale_factor_today))
!!$    a = exp(lna)
!!$    dlna = lna(2) - lna(1)
!!$
!!$    
!!$    call ode%init(n=2, method = COOP_ODE_RK4)
!!$    
!!$    call ode%set_initial_conditions(xini = lna(1), yini=y(:, 1))
!!$    i = 1
!!$    call cpl_eq_get_potential(3, lna(i), y(:, i), yp(:, i), H2a4 = H2a4(i), lnV = lnV(i), wp1de = wp1de(i), wp1effde = wp1effde(i), wp1effcdm = wp1effcdm(i), omde  = omde, omc = omc)
!!$    if(err .ne. 0)goto 100
!!$    do i=2, ns
!!$       call ode%evolve(cpl_eq, lna(i))
!!$       y(:, i) = ode%y
!!$       if(err .ne. 0) goto 100
!!$       call cpl_eq_get_potential(3, lna(i), y(:, i), yp(:, i), H2a4 = H2a4(i), lnV = lnV(i), wp1de = wp1de(i), wp1effde = wp1effde(i), wp1effcdm = wp1effcdm(i), omde  = omde, omc = omc)
!!$       if(err .ne. 0)goto 100       
!!$    enddo
!!$    call fwp1effcdm%init(n = ns, xmin=coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., f = wp1effcdm, method = COOP_INTERPOLATE_LINEAR, name = "CDM 1+w", check_boundary = .false.)
!!$    call cdm%init(genre=COOP_SPECIES_FLUID, name = "CDM", id = 1, Omega = omc, w = 0.d0, fwp1eff = fwp1effcdm)
!!$
!!$    
!!$    call fwp1de%init(n = ns, xmin=coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., f = wp1de, method = COOP_INTERPOLATE_LINEAR, name = "Dark Energy 1+w", check_boundary = .false.)
!!$    call fwp1effde%init(n = ns, xmin=coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., f = wp1effde, method = COOP_INTERPOLATE_LINEAR, name= "Dark Energy effective 1+w", check_boundary = .false.)
!!$    call de%init(name = "Dark Energy", id = 5, Omega = omde, genre = COOP_SPECIES_FLUID, fwp1 = fwp1de, fwp1eff =fwp1effde)
!!$    
!!$    de%cplde_wp1 = fwp1
!!$    de%cplde_Q = fQ
!!$    !$omp parallel do
!!$    do i=1, ns
!!$       dQdphi(i)=fQ%derivative(a(i))*a(i)/max(yp(3, i), min_phi_prime)
!!$    enddo
!!$    !$omp end parallel do
!!$    call de%cplde_dQdphi_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = dQdphi, method = COOP_INTERPOLATE_LINEAR, name = "dQ / d phi")
!!$    call de%cplde_lnV_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = lnV, method = COOP_INTERPOLATE_LINEAR, name = "ln V")
!!$    call de%cplde_phi_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = y(3,:), method = COOP_INTERPOLATE_LINEAR, name = "phi")
!!$    call de%cplde_phi_prime_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = yp(3,:), method = COOP_INTERPOLATE_LINEAR, name="d phi / d lna")
!!$    
!!$    
!!$    if(i_tc_off .gt. 1)then
!!$       !$omp parallel do private(tc_w)
!!$       do i= 2, ns-1
!!$          if(fQ%eval(exp(lna(i)))*yp(3, i).gt. 0.d0)then
!!$             tc_w = (1.d0-tanh(20.d0*(lna(i)-lna(i_tc_off))))/2.d0
!!$          else
!!$             tc_w = 0.d0
!!$          endif
!!$          !!weighted sum of numeric derivative and tight-coupling approximation
!!$          dlnVdphi(i) =  (lnV(i+1)-lnV(i-1))/max(2.d0*yp(3, i)*dlna, 1.d-20)*(1.d0 - tc_w) &
!!$          - cdm%density(exp(lna(i)))*fQ%eval(exp(lna(i)))/exp(lnV(i))*tc_w
!!$       enddo
!!$       !$omp end parallel do
!!$    else
!!$       !$omp parallel do 
!!$       do i= 2, ns-1
!!$          dlnVdphi(i) =  (lnV(i+1)-lnV(i-1))/max(2.d0*yp(3, i)*dlna, 1.d-20)
!!$       enddo
!!$       !$omp end parallel do
!!$       
!!$    endif
!!$    dlnVdphi(ns) = dlnVdphi(ns-1)    
!!$    dlnVdphi(1) = dlnVdphi(2)
!!$    if(any(dlnVdphi .gt. 0.d0))then
!!$       err = 6
!!$       return
!!$    endif
!!$    dVdphibyH2 = dlnVdphi * exp(lnV + 4.d0*lna) / H2a4 
!!$
!!$    call de%cplde_dVdphibyH2_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = dVdphibyH2, method = COOP_INTERPOLATE_LINEAR, name = "d V / d phi / H^2")
!!$    dlnVdphi = min(dlnVdphi, -1.d-30)
!!$    !!compute m2a2, truncate m^2/H^2 (fast oscillations cannot be numerically resolved, but they are irrelevant for observables.)
!!$    do i = 3, ns-2
!!$       m2byH2(i) = ( -(log(-dlnVdphi(i+1)) - log(-dlnVdphi(i-1)))/max(2.d0*yp(3, i)*dlna, 1.d-20) - dlnVdphi(i) ) * (-dlnVdphi(i)) * exp(lnV(i) + 2.d0*lna(i))/(h2a4(i)*exp(-lna(i)*2.d0))
!!$    enddo
!!$    m2byH2(1:2) = m2byH2(3)
!!$    m2byH2(ns-1:ns) = m2byH2(ns)
!!$    call de%cplde_m2byH2_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = m2byH2, method = COOP_INTERPOLATE_LINEAR, name = " m^2/H^2 ")
!!$    
!!$    
!!$    call this%add_species(cdm)
!!$    call this%add_species(de)
!!$100 continue
!!$    call de%free()
!!$    call cdm%free()
!!$    call fwp1effcdm%free()
!!$    call fwp1de%free()
!!$    call fwp1effde%free()
!!$    
!!$  contains
!!$
!!$
!!$    function lnV_of_phi(phi)
!!$      COOP_REAL::phi, lnV
!!$      lnV = lnV0 - alpha*log(phi)
!!$    end function lnV_of_phi
!!$
!!$    subroutine cpl_eq(n, lna, y, yp)
!!$      COOP_INT::n
!!$      COOP_REAL::lna, y(n), yp(n), H2a4
!!$      call cpl_eq_get_potential(n, lna, y, yp, H2a4)
!!$    end subroutine cpl_eq
!!$
!!$    subroutine cpl_eq_get_potential(n, lna, y, yp, H2a4, lnV, wp1de, wp1effde, wp1effcdm, omde, omc)
!!$      COOP_INT::n
!!$      COOP_REAL, optional::lnV, wp1de, wp1effde, wp1effcdm, omde, omc
!!$      COOP_REAL::lna, y(n), yp(n)
!!$      COOP_REAL:: a, wp1, Q, phi_prime, rho_de, rhoa4_de_eff, rhoa4_de, rhoa4_c, delta_rhoa4_c, rhoa4_other, h2a4, delta, wder
!!$      a = exp(lna)
!!$      Q = fQ%eval(a)
!!$      LN_RHOC_RAT_PRIME = Q*PHID
!!$      PHI_PRIME = PHID
!!$      PHID_PRIME = 
!!$      
!!$      rhoa4_de_eff = exp(y(1)+4.d0*lna)
!!$      rhoa4_c = rho_ce*exp(y(2)+lna)
!!$      delta_rhoa4_c = rho_ce*a*(exp(y(2))-1.d0)
!!$      if(rhoa4_de_eff .lt. delta_rhoa4_c)then !!negative energy
!!$         err = 2
!!$         return
!!$      endif
!!$      rhoa4_de = rhoa4_de_eff - delta_rhoa4_c
!!$      rho_de = exp(y(1)) - delta_rhoa4_c/a**4
!!$      rhoa4_other = this%rhoa4(a)
!!$      h2a4 = (rhoa4_de + rhoa4_c + rhoa4_other)/3.d0
!!$
!!$      if(Q .gt. 0.d0 )then
!!$         wder = fwp1%derivative(a)*a -3.d0*(wp1-1.d0)*wp1
!!$         if(H2a4*rhoa4_de_eff*(wder/rhoa4_c/Q)**2  .lt. tc_tol*wp1)then
!!$            phi_prime = (wder * rhoa4_de_eff / rhoa4_c / Q)
!!$            if(phi_prime .lt. 0.d0)then
!!$               err = 5
!!$               return
!!$            endif
!!$            phi_prime = phi_prime * min(((wp1*rhoa4_de_eff)/( delta_rhoa4_c + phi_prime**2*H2a4))**20, 1.d0)
!!$            i_tc_off = i+1
!!$            goto 50
!!$         endif
!!$      endif
!!$      delta = (rhoa4_de_eff*wp1 - delta_rhoa4_c)/H2a4
!!$      if(delta.lt. 0.d0)then
!!$         err = 3
!!$         return
!!$      endif
!!$      phi_prime = sqrt(delta)
!!$50    yp(2) = Q*phi_prime
!!$      yp(3) = phi_prime
!!$      if(present(lnV))then
!!$         lnV = rho_de - phi_prime**2 * H2a4/2.d0/a**4 
!!$         if(lnV .le. 0.d0)then
!!$            err = 4
!!$            return
!!$         else
!!$            lnV = log(lnV)
!!$         endif
!!$         wp1de =  phi_prime**2*H2a4/rhoa4_de
!!$         wp1effcdm =  1.d0 - Q*phi_prime/3.d0
!!$         wp1effde = wp1de + Q*phi_prime/3.d0 * rhoa4_c/rhoa4_de
!!$         omde = rhoa4_de/H2a4/3.d0
!!$         omc = rhoa4_c/H2a4/3.d0
!!$      else
!!$         if(rhoa4_de .lt. phi_prime**2*H2a4/2.d0)then
!!$            err = 4
!!$            return
!!$         endif
!!$      endif
!!$    end subroutine cpl_eq_get_potential
!!$
!!$    
!!$#else
!!$    write(*,*) "Coupled Dark Energy model cannot be initialized"
!!$    stop "You need to set DARK_ENERGY_MODEL=COUPLED_DE in configure.in"
!!$#endif    
!!$  end subroutine coop_background_add_coupled_DE_negative_powerlaw
  
  subroutine coop_background_add_coupled_DE(this, Omega_c, fQ, fwp1, err)

    !!err = 0: success
    !!err = -1: Omega_DE < 0
    !!err = 1: w < -1
    !!err = 2: negative DE energy
    !!err = 3: negative DE kinetic energy
    !!err = 4: negaitve DE potential energy
    !!err = 5: non-monotonic potential
    COOP_REAL, parameter::tc_tol = 5.d-3
    COOP_REAL, parameter::phi_start = 0.d0
    COOP_REAL, parameter::min_phi_prime = 1.d-10
    class(coop_cosmology_background)::this
    COOP_REAL::Omega_c
    type(coop_function)::fQ, fwp1
    COOP_INT::err
#if DO_COUPLED_DE
    COOP_INT, parameter::ns = 12000
    type(coop_function)::fwp1effcdm, fwp1de, fwp1effde
    type(coop_ode)::ode
    type(coop_species)::de, cdm, deeff
    COOP_INT::i, i_tc_off
    COOP_REAL::rho_ce, a(ns), lna(ns), y(3, ns), yp(3, ns), lnV(ns), dlnVdphi(ns), dVdphibyH2(ns), wp1de(ns), wp1effde(ns), wp1effcdm(ns), H2a4(ns), m2byH2(ns), dQdphi(ns), omde, omc, dlna, tc_w
    err = 0
    i_tc_off = 1
    if(this%Omega_k() .le. Omega_c)then
       err = -1
       return
    endif
    rho_ce = omega_c * 3.d0
    call coop_set_uniform(ns, lna, log(coop_min_scale_factor), log(coop_scale_factor_today))
    a = exp(lna)
    dlna = lna(2) - lna(1)
    call ode%init(n=3, method = COOP_ODE_RK4)
    y(1, 1) = log(3.d0*(this%Omega_k() - Omega_c)) + 3.d0*coop_integrate(wp1_eval, lna(1), lna(ns))
    y(2, 1) = 0.d0
    y(3, 1) = phi_start
    call ode%set_initial_conditions(xini = lna(1), yini=y(:, 1))
    i = 1
    call cpl_eq_get_potential(3, lna(i), y(:, i), yp(:, i), H2a4 = H2a4(i), lnV = lnV(i), wp1de = wp1de(i), wp1effde = wp1effde(i), wp1effcdm = wp1effcdm(i), omde  = omde, omc = omc)
    if(err .ne. 0)goto 100
    do i=2, ns
       call ode%evolve(cpl_eq, lna(i))
       y(:, i) = ode%y
       if(err .ne. 0) goto 100
       call cpl_eq_get_potential(3, lna(i), y(:, i), yp(:, i), H2a4 = H2a4(i), lnV = lnV(i), wp1de = wp1de(i), wp1effde = wp1effde(i), wp1effcdm = wp1effcdm(i), omde  = omde, omc = omc)
       if(err .ne. 0)goto 100       
    enddo
    call fwp1effcdm%init(n = ns, xmin=coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., f = wp1effcdm, method = COOP_INTERPOLATE_LINEAR, name = "CDM 1+w", check_boundary = .false.)
    call cdm%init(genre=COOP_SPECIES_FLUID, name = "CDM", id = 1, Omega = omc, w = 0.d0, fwp1eff = fwp1effcdm)

    
    call fwp1de%init(n = ns, xmin=coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., f = wp1de, method = COOP_INTERPOLATE_LINEAR, name = "Dark Energy 1+w", check_boundary = .false.)
    call fwp1effde%init(n = ns, xmin=coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., f = wp1effde, method = COOP_INTERPOLATE_LINEAR, name= "Dark Energy effective 1+w", check_boundary = .false.)
    call de%init(name = "Dark Energy", id = 5, Omega = omde, genre = COOP_SPECIES_FLUID, fwp1 = fwp1de, fwp1eff =fwp1effde)
    
    de%cplde_wp1 = fwp1
    de%cplde_Q = fQ
    !$omp parallel do
    do i=1, ns
       dQdphi(i)=fQ%derivative(a(i))*a(i)/max(yp(3, i), min_phi_prime)
    enddo
    !$omp end parallel do
    call de%cplde_dQdphi_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = dQdphi, method = COOP_INTERPOLATE_LINEAR, name = "dQ / d phi")
    call de%cplde_lnV_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = lnV, method = COOP_INTERPOLATE_LINEAR, name = "ln V")
    call de%cplde_phi_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = y(3,:), method = COOP_INTERPOLATE_LINEAR, name = "phi")
    call de%cplde_phi_prime_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = yp(3,:), method = COOP_INTERPOLATE_LINEAR, name="d phi / d lna")
    
    
    if(i_tc_off .gt. 1)then
       !$omp parallel do private(tc_w)
       do i= 2, ns-1
          if(fQ%eval(exp(lna(i)))*yp(3, i).gt. 0.d0)then
             tc_w = (1.d0-tanh(20.d0*(lna(i)-lna(i_tc_off))))/2.d0
          else
             tc_w = 0.d0
          endif
          !!weighted sum of numeric derivative and tight-coupling approximation
          dlnVdphi(i) =  (lnV(i+1)-lnV(i-1))/max(2.d0*yp(3, i)*dlna, 1.d-20)*(1.d0 - tc_w) &
          - cdm%density(exp(lna(i)))*fQ%eval(exp(lna(i)))/exp(lnV(i))*tc_w
       enddo
       !$omp end parallel do
    else
       !$omp parallel do 
       do i= 2, ns-1
          dlnVdphi(i) =  (lnV(i+1)-lnV(i-1))/max(2.d0*yp(3, i)*dlna, 1.d-20)
       enddo
       !$omp end parallel do
       
    endif
    dlnVdphi(ns) = dlnVdphi(ns-1)    
    dlnVdphi(1) = dlnVdphi(2)
    if(any(dlnVdphi .gt. 0.d0))then
       err = 6
       return
    endif
    dVdphibyH2 = dlnVdphi * exp(lnV + 4.d0*lna) / H2a4 

    call de%cplde_dVdphibyH2_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = dVdphibyH2, method = COOP_INTERPOLATE_LINEAR, name = "d V / d phi / H^2")
    dlnVdphi = min(dlnVdphi, -1.d-30)
    !!compute m2a2, truncate m^2/H^2 (fast oscillations cannot be numerically resolved, but they are irrelevant for observables.)
    do i = 3, ns-2
       m2byH2(i) = ( -(log(-dlnVdphi(i+1)) - log(-dlnVdphi(i-1)))/max(2.d0*yp(3, i)*dlna, 1.d-20) - dlnVdphi(i) ) * (-dlnVdphi(i)) * exp(lnV(i) + 2.d0*lna(i))/(h2a4(i)*exp(-lna(i)*2.d0))
    enddo
    m2byH2(1:2) = m2byH2(3)
    m2byH2(ns-1:ns) = m2byH2(ns)
    call de%cplde_m2byH2_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = m2byH2, method = COOP_INTERPOLATE_LINEAR, name = " m^2/H^2 ")
    
    
    call this%add_species(cdm)
    call this%add_species(de)
100 continue
    call de%free()
    call cdm%free()
    call fwp1effcdm%free()
    call fwp1de%free()
    call fwp1effde%free()
    
  contains

    function wp1_eval(lna) result(wp1)
      COOP_REAL::lna, wp1
      wp1 = fwp1%eval(exp(lna))      
    end function wp1_eval
    
    !! rho_ce = \lim a-> 0 rho a^3
    !! rho_eff = rho_DE + rho_c -  rho_ce a^{-3} 
    !! y(1) = ln (rho_eff)  
    !! y(2) = ln (rho_c a^3 / rho_ce ) (>0)
    !! y(3) = phi
    !!time variable lna
    subroutine cpl_eq(n, lna, y, yp)
      COOP_INT::n
      COOP_REAL::lna, y(n), yp(n), H2a4
      call cpl_eq_get_potential(n, lna, y, yp, H2a4)
    end subroutine cpl_eq

    subroutine cpl_eq_get_potential(n, lna, y, yp, H2a4, lnV, wp1de, wp1effde, wp1effcdm, omde, omc)
      COOP_INT::n
      COOP_REAL, optional::lnV, wp1de, wp1effde, wp1effcdm, omde, omc
      COOP_REAL::lna, y(n), yp(n)
      COOP_REAL:: a, wp1, Q, phi_prime, rho_de, rhoa4_de_eff, rhoa4_de, rhoa4_c, delta_rhoa4_c, rhoa4_other, h2a4, delta, wder
      a = exp(lna)
      wp1 =fwp1%eval(a)
      !! 1 + w < 0: error
      if(wp1.lt.0.d0)then
         err = 1
         return
      endif

      Q = fQ%eval(a)
      yp(1) = -3.d0*wp1
      rhoa4_de_eff = exp(y(1)+4.d0*lna)
      rhoa4_c = rho_ce*exp(y(2)+lna)
      delta_rhoa4_c = rho_ce*a*(exp(y(2))-1.d0)
      if(rhoa4_de_eff .lt. delta_rhoa4_c)then !!negative energy
         err = 2
         return
      endif
      rhoa4_de = rhoa4_de_eff - delta_rhoa4_c
      rho_de = exp(y(1)) - delta_rhoa4_c/a**4
      rhoa4_other = this%rhoa4(a)
      h2a4 = (rhoa4_de + rhoa4_c + rhoa4_other)/3.d0

      if(Q .gt. 0.d0 )then
         wder = fwp1%derivative(a)*a -3.d0*(wp1-1.d0)*wp1
         if(H2a4*rhoa4_de_eff*(wder/rhoa4_c/Q)**2  .lt. tc_tol*wp1)then
            phi_prime = (wder * rhoa4_de_eff / rhoa4_c / Q)
            if(phi_prime .lt. 0.d0)then
               err = 5
               return
            endif
            phi_prime = phi_prime * min(((wp1*rhoa4_de_eff)/( delta_rhoa4_c + phi_prime**2*H2a4))**20, 1.d0)
            i_tc_off = i+1
            goto 50
         endif
      endif
      delta = (rhoa4_de_eff*wp1 - delta_rhoa4_c)/H2a4
      if(delta.lt. 0.d0)then
         err = 3
         return
      endif
      phi_prime = sqrt(delta)
50    yp(2) = Q*phi_prime
      yp(3) = phi_prime
      if(present(lnV))then
         lnV = rho_de - phi_prime**2 * H2a4/2.d0/a**4 
         if(lnV .le. 0.d0)then
            err = 4
            return
         else
            lnV = log(lnV)
         endif
         wp1de =  phi_prime**2*H2a4/rhoa4_de
         wp1effcdm =  1.d0 - Q*phi_prime/3.d0
         wp1effde = wp1de + Q*phi_prime/3.d0 * rhoa4_c/rhoa4_de
         omde = rhoa4_de/H2a4/3.d0
         omc = rhoa4_c/H2a4/3.d0
      else
         if(rhoa4_de .lt. phi_prime**2*H2a4/2.d0)then
            err = 4
            return
         endif
      endif
    end subroutine cpl_eq_get_potential
#else
    write(*,*) "Coupled Dark Energy model cannot be initialized"
    stop "You need to set DARK_ENERGY_MODEL=COUPLED_DE in configure.in"
#endif    
  end subroutine coop_background_add_coupled_DE


  subroutine coop_background_add_EFT_DE_with_effective_w(this, effective_wp1, err)
    COOP_INT,parameter::narr = 25000
    class(coop_cosmology_background)::this
    type(coop_function)::effective_wp1
    type(coop_species)::de, deeff
    COOP_INT::i, err, j, index_de
    COOP_REAL,dimension(narr)::lnrho,  wp1, wp1eff
    COOP_REAL::lna, lnamin, lnamax, dlna, alpha_l, alpha_r, dlnaby2, a, rhoa4de,ppra4de, ppra4tot, rhoa4tot, rhoa4de_bg, M2, rhom0, wp1_bg
#if DO_EFT_DE
    err = 0    
    rhom0=this%rhoa4(1.d0)
    if(rhom0 .gt. 3.d0)then
       err = 1
       return
    endif
    call deeff%init(name = "Dark Energy", id = 5, Omega = (3.d0-rhom0)/3.d0/coop_Mpsq0, genre = COOP_SPECIES_FLUID, fwp1 = effective_wp1)

    de%Omega = this%Omega_k()
    de%name = "Dark Energy"
    de%genre = COOP_SPECIES_EFT
    
    lnamin = log(coop_min_scale_factor)
    lnamax = log(coop_scale_factor_today)
    dlna = (lnamax-lnamin)/(narr-1.d0)
    dlnaby2 = dlna/2.d0
    lnrho(narr) = log(3.d0*de%Omega*coop_Mpsq0)
    lna = lnamin
    do i=1, narr
       a = exp(lna)
       M2 = coop_Mpsq(a)
       wp1_bg = deeff%wp1ofa(a)
       rhoa4de_bg = deeff%rhoa4(a)       
       call this%get_ppra4_rhoa4(a, ppra4tot, rhoa4tot)
       rhoa4de = rhoa4tot*(M2 - 1.d0) + rhoa4de_bg*M2
       ppra4de = ppra4tot*(M2 - 1.d0) + wp1_bg * rhoa4de_bg*M2
       if(rhoa4de .le. 0.d0)then
          err = 1  !!negative energy flag
          return
       else
          lnrho(i) = log(rhoa4de) - lna * 4.d0
          wp1(i) = ppra4de/rhoa4de
       endif
       wp1eff(i) = wp1_bg*rhoa4de_bg/rhoa4de - (coop_alphaM(a)*(rhoa4tot+rhoa4de_bg)*M2/3.d0 - (M2-1.d0)*(ppra4tot + wp1_bg * rhoa4de_bg))/rhoa4de
       lna = lna + dlna
    enddo


    lnrho = lnrho - lnrho(narr)
    call de%fwp1%init(narr, coop_min_scale_factor, coop_scale_factor_today, wp1, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = "DE 1+w(a)")    
    call de%fwp1eff%init(narr, coop_min_scale_factor, coop_scale_factor_today, wp1eff, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = "DE 1+w_eff(a)")
    call de%flnrho%init(narr,coop_min_scale_factor, coop_scale_factor_today, lnrho, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = "DE ln rho_ratio")
    call de%flnrho%set_boundary(slopeleft = -3.d0*wp1eff(1), sloperight = -3.d0*wp1eff(narr))    
    de%cs2 = 0.d0
    call this%delete_species(index_de)
    call this%add_species(de)
    call de%free()
    call deeff%free()
#else
    write(*,*) "EFT Dark Energy model cannot be initialized"
    stop "You need to set DARK_ENERGY_MODEL=EFT in configure.in"
#endif    
  end subroutine coop_background_add_EFT_DE_with_effective_w




  subroutine coop_background_add_EFT_DE(this, wp1,  err)
    COOP_INT,parameter::narr = 25000
    class(coop_cosmology_background)::this
    type(coop_function)::wp1
    type(coop_species)::de
    COOP_INT::i, err, j
    COOP_REAL,dimension(narr)::lnrho,  wp1eff
    COOP_REAL::lna, lnamin, lnamax, dlna, alpha_l, a_l, a_r, alpha_r, wp1_l, wp1_r, omega_de, om_l, om_r, rhoa4de_l, rhotot_l, rhotot_r, step, rhoa4de_r
#if DO_EFT_DE        
    err = 0
    de%name = "Dark Energy"
    de%genre = COOP_SPECIES_EFT
    de%fwp1 = wp1
    omega_de = this%Omega_k()    
    de%Omega = omega_de
    lnamin = log(coop_min_scale_factor)
    lnamax = log(coop_scale_factor_today)
    dlna = (lnamax-lnamin)/(narr-1.d0)
    lnrho(narr) = log(3.d0*omega_de*coop_Mpsq0)
    wp1_r = wp1%eval(coop_scale_factor_today)
    alpha_r = coop_alphaM(coop_scale_factor_today)
    om_r = omega_de
    wp1eff(narr) =wp1_r - alpha_r/3.d0/om_r
    a_r = coop_scale_factor_today
    lna = lnamax
    rhotot_r = 3.d0*coop_Mpsq0
    rhoa4de_r = om_r*rhotot_r
    step = (1.5d0*dlna)
    do i=narr-1, 1, -1
       lna = lna - dlna              
       a_l = exp(lna)
       alpha_l = coop_alphaM(a_l)
       wp1_l = wp1%eval(a_l)
       rhotot_l =  this%rhoa4(a_l)
       om_l = rhoa4de_r/(rhotot_l + rhoa4de_r)  !!first assuming rho_de constant
       wp1eff(i) = wp1_l - alpha_l/3.d0/om_l
       if(om_l .gt. 1.d-2)then
          do j=1, 5
             lnrho(i) = lnrho(i+1) + (wp1eff(i)+wp1eff(i+1))*step
             rhoa4de_l = exp(lnrho(i)+4.d0*lna)
             om_l = rhoa4de_l/(rhoa4de_l + rhotot_l)
             wp1eff(i) = wp1_l - alpha_l/3.d0/om_l
             if(om_l .lt. 1.d-50)then
                lnrho(1:i-1) = lnrho(i)
                wp1eff(1:i-1) = wp1eff(i)
                goto 100
             endif
          enddo          
       elseif(om_l .gt. 1.d-4)then
          do j=1, 2
             lnrho(i) = lnrho(i+1) + (wp1eff(i)+wp1eff(i+1))*step
             rhoa4de_l = exp(lnrho(i)+4.d0*lna)
             om_l = rhoa4de_l/(rhoa4de_l + rhotot_l)
             wp1eff(i) = wp1_l - alpha_l/3.d0/om_l
             if(om_l .lt. 1.d-50)then
                lnrho(1:i-1) = lnrho(i)
                wp1eff(1:i-1) = wp1eff(i)
                goto 100
             endif
          enddo
       else
          lnrho(i) = lnrho(i+1) + (wp1eff(i)+wp1eff(i+1))*step
          rhoa4de_l = exp(lnrho(i)+4.d0*lna)
          om_l = rhoa4de_l/(rhoa4de_l + rhotot_l)
          if(om_l .lt. 1.d-50)then
             lnrho(1:i-1) = lnrho(i)
             wp1eff(1:i-1) = wp1eff(i)
             goto 100
          endif
          wp1eff(i) = wp1_l - alpha_l/3.d0/om_l          
       endif
       rhotot_l = (rhotot_l + rhoa4de_l)/a_l**4
       if(rhotot_l .lt. rhotot_r)then  !!
          err = 1
          return
       endif
       rhotot_r = rhotot_l
       rhoa4de_r = rhoa4de_l
       alpha_r = alpha_l
    enddo
100 continue
    lnrho = lnrho - lnrho(narr)
    call de%fwp1eff%init(narr, coop_min_scale_factor, coop_scale_factor_today, wp1eff, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = "DE 1+w_eff(a)")
    call de%flnrho%init(narr,coop_min_scale_factor, coop_scale_factor_today, lnrho, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = "DE ln rho_ratio")
    call de%flnrho%set_boundary(slopeleft = -3.d0*wp1eff(1), sloperight = -3.d0*wp1eff(narr))    
    de%cs2 = 0.d0
    call this%add_species(de)
    call de%free()
#else
    write(*,*) "EFT Dark Energy model cannot be initialized"
    stop "You need to set DARK_ENERGY_MODEL=EFT in configure.in"
#endif    
  end subroutine coop_background_add_EFT_DE
  

  function coop_de_alpha_invh2(a, arg) result(alpha)
    type(coop_arguments)::arg  !! arg%r = (/ alpha0, Omega_m, Omega_r /) 
    COOP_REAL::alpha
    COOP_REAL::a
    alpha = arg%r(1)*a**4/(arg%r(2)*a + arg%r(3) + (1.d0 - arg%r(2) - arg%r(3))*a**4)
  end function coop_de_alpha_invh2

  function coop_de_alpha_instant(a, arg) result(alpha)
    type(coop_arguments)::arg !! arg%r = (/ alpha0, Omega_m, lambda, delta_z /)
    COOP_REAL::alpha, a, zp1, zp1cr
    zp1 = 1.d0/max(a, coop_min_scale_factor) 
    zp1cr = (((1.d0 - arg%r(2))/arg%r(2))/arg%r(3))**(1.d0/3.d0) 
    alpha = arg%r(1) * (1.d0-tanh((zp1 - zp1cr)/arg%r(4)))/2.d0
  end function coop_de_alpha_instant


  function coop_de_alpha_bump(a, arg) result(alpha)
    type(coop_arguments)::arg !! arg%r = (/ alpha0, Omega_m, lambda, delta_z /)
    COOP_REAL::alpha, a, zp1, zp1cr
    zp1 = 1.d0/max(a, coop_min_scale_factor) 
    zp1cr = (((1.d0 - arg%r(2))/arg%r(2))/arg%r(3))**(1.d0/3.d0) 
    alpha = arg%r(1) * exp(- ((zp1 - zp1cr)/arg%r(4))**2)
  end function coop_de_alpha_bump

  function coop_de_alpha_constructor(alpha0, genre) result(alpha)
    COOP_REAL, parameter::Omega_m = 0.3d0
    COOP_REAL, parameter::Omega_r = 8.d-5
    COOP_REAL, parameter::lambda = 0.001d0
    COOP_REAL, parameter::delta_z = 0.1d0
    COOP_REAL::alpha0, omegam, omegar
    COOP_UNKNOWN_STRING::genre
    type(coop_function)::alpha
    select case(COOP_LOWER_STR(genre))
    case("omega")
       alpha = coop_function_constructor( coop_de_alpha_invh2, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args= coop_arguments_constructor ( r = (/ alpha0, Omega_m, Omega_r /) ), name = "alpha")
    case("linear")
       call alpha%init_polynomial( (/ 0.d0, alpha0 /) )
    case("instant")
       alpha = coop_function_constructor( coop_de_alpha_instant, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args= coop_arguments_constructor ( r = (/ alpha0, Omega_m, lambda, delta_z /) ) , name = "alpha" )
    case("bump")
       alpha = coop_function_constructor( coop_de_alpha_bump, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args= coop_arguments_constructor ( r = (/ alpha0, Omega_m, lambda, delta_z /) ) , name = "alpha" )
    case default
       write(*,*) trim(genre)
       stop "unknown alpha function type"
    end select
  end function coop_de_alpha_constructor
  
end module coop_background_mod
