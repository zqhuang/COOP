module coop_background_mod
  use coop_wrapper_utils
  implicit none
#include "constants.h"

  private
  COOP_REAL,parameter::coop_min_de_tracking_n = 0.01d0

  public::coop_baryon, coop_cdm, coop_DE_lambda, coop_DE_w0, coop_DE_w0wa, coop_DE_quintessence, coop_radiation, coop_neutrinos_massless, coop_neutrinos_massive, coop_de_w_quintessence, coop_de_wp1_quintessence, coop_background_add_coupled_DE,  coop_background_add_coupled_DE_with_w, coop_background_add_EFT_DE, coop_de_aeq_fitting

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
  
  function coop_de_w0wa(Omega_Lambda, w0, wa) result(this)
    type(coop_species) this
    COOP_REAL Omega_Lambda, w0, wa
    type(coop_function) fw0wa
    type(coop_arguments) w0wa
    if(w0 .eq. -1.d0 .and. wa.eq.0.d0)then
       this = coop_de_lambda(Omega_Lambda)
       return
    endif
    w0wa = coop_arguments(r =  (/ w0, wa /))
    fw0wa = coop_function(coop_de_wp1_w0wa, xmin = coop_min_scale_factor, xmax = COOP_REAL_OF(1.), xlog = .true., args = w0wa, name = "w0wa model 1+w")
    call this%init(genre = COOP_SPECIES_FLUID, name = "Dark Energy", id=5, Omega = Omega_Lambda, cs2 = COOP_REAL_OF(1.d0), fwp1 = fw0wa )
    call w0wa%free
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
    arg = coop_arguments(r =  (/ Omega_Lambda, epsilon_s, epsilon_inf, zeta_s /))
    fq = coop_function(coop_de_wp1_quintessence, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args = arg, name = "quintessence 1+w")
    call this%init(genre = COOP_SPECIES_FLUID, name = "Dark Energy", id=5, Omega = Omega_Lambda, cs2 = COOP_REAL_OF(1.d0), fwp1 = fq )
    call arg%free
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


  subroutine coop_background_add_coupled_DE(this, Omega_c, Q, tracking_n, dlnQdphi, dUdphi, d2Udphi2)
    !!This Omega_c is the effective Omega_cdm (rho_m a^3 (a<<1) -> Omega_c rho_0) 
    !!Q is the coupling between DE and CDM
    !!The potential V(phi) = V0 / phi^n
    !!Tracking_n is the index n
    !!V0 will be automatically determined by the constraint Omega_k = 0
    class(coop_cosmology_background)::this
    COOP_REAL::omega_c, Q, tracking_n
    type(coop_species)::cdm, de
    type(coop_ode)::ode, ode_tc
    COOP_INT i, istart,  j
    COOP_REAL  beta, lnrhom0
    COOP_REAL::lnc_lower, lnc_upper, lnc_mid
    COOP_REAL::Ot_lower, Ot_upper, Ot_mid, lna_ini
    COOP_REAL_ARRAY::lna, lnrho_cdm, lnrho_de, wp1_de, phi_de, phidot_de, wp1eff_de, wp1eff_cdm
    COOP_REAL,parameter::a_ini = 1.d-8
    COOP_INT,parameter::nsteps = 500
    COOP_REAL::lna_coarse(nsteps)
    COOP_REAL, optional::dlnQdphi, dUdphi, d2Udphi2
    COOP_REAL::U1, U2
    logical::tight_coupling
    cdm%name = "CDM"
    de%name = "Dark Energy"
    cdm%genre = COOP_SPECIES_COUPLED
    de%genre = COOP_SPECIES_COUPLED
    lnrhom0 = log(3.d0*Omega_c)
    call coop_set_uniform(coop_default_array_size, lna, log(coop_min_scale_factor), log(coop_scale_factor_today))
    do i = 1, coop_default_array_size
       if(lna(i) .gt. log(a_ini))then
          istart = i
          lna_ini = lna(i)
          exit
       endif
    enddo
    call coop_set_uniform(nsteps, lna_coarse, lna(istart+1), log(coop_scale_factor_today))
    
    if(Q .gt. 1.d-6)then
       if(present(dlnQdphi))then
          call de%fDE_Q_of_phi%init_polynomial( (/ log(Q), dlnQdphi /), ylog = .true. )
       else
          call de%fDE_Q_of_phi%init_polynomial( (/ Q /) )
       endif
       de%DE_tracking_n = max(tracking_n, coop_min_de_tracking_n)  
    else
       de%DE_tracking_n = max(tracking_n, 0.d0)  
    endif


    if(de%DE_tracking_n .ge. 2.d0)then
       beta = 2.d0/(de%DE_tracking_n+2.d0)
    else
       if(de%fDE_Q_of_phi%initialized)then
          beta = 1.5d0/(de%DE_tracking_n+1.d0)
       else
          beta = 2.d0/(de%DE_tracking_n+2.d0)
       endif
    endif
    !!change this for dynamic U(phi) and dU/d\phi
    if(present(dUdphi))then
       U1 = dUdphi
    else
       U1 = 0.d0
    endif
    if(present(d2Udphi2))then
       U2 = d2Udphi2
    else
       U2 = 0.d0
    endif
    call de%fDE_U_of_phi%init_polynomial( (/ 0.d0, U1, U2 /), name = "U(phi)")  
    call de%fDE_dUdphi%init_polynomial( (/ U1, U2*2.d0 /), name = "dU/dphi" ) 


    call ode%init(3)
    call ode_tc%init(1)
    
#define PHI y(1)
#define PHIDOT y(2)
#define LN_RHOM y(3)
#define PHI_PRIME yp(1)
#define PHIDOT_PRIME yp(2)
#define LN_RHOM_PRIME yp(3)
    lnc_upper = 0.5d0
100 call get_ini(lna_ini, lnc_upper)
    do i = 1, nsteps
       call evolve_to(lna_coarse(i))
    enddo
    Ot_upper = (exp(ode%LN_RHOM) + de%DE_V(ode%PHI) + ode%PHIDOT**2/2.d0)/3.d0 - this%Omega_k()
    if(Ot_upper .le. 0.d0)then
       lnc_upper = lnc_upper + 0.5d0
       goto 100
    endif

    
    lnc_lower = -1.d0    
200 call get_ini(lna_ini, lnc_lower)
    do i = 1, nsteps
       call evolve_to(lna_coarse(i))
    enddo
    Ot_lower = (exp(ode%LN_RHOM) + de%DE_V(ode%PHI) + ode%PHIDOT**2/2.d0)/3.d0 - this%Omega_k()
    if(Ot_lower .ge. 0.d0)then
       lnc_lower = lnc_lower - 0.5d0
       goto 200
    endif

    do
       if(Ot_upper - Ot_lower .lt. 1.d-3)then
          lnc_mid = (lnc_upper*(-Ot_lower) + lnc_lower*Ot_upper)/(Ot_upper - Ot_lower)                 
          do i=1, istart 
             call get_ini(lna(i), lnc_mid)
             lnrho_cdm(i) = ode%LN_RHOM
             phi_de(i) = ode%PHI
             phidot_de(i) = ode%PHIDOT
             call getweff(3, lna(i), ode%y, wp1_de(i), wp1eff_cdm(i), wp1eff_de(i), lnrho_de(i))                          
          enddo
500       continue
          do i = istart + 1, coop_default_array_size
             call evolve_to(lna(i))
             lnrho_cdm(i) = ode%LN_RHOM
             phi_de(i) = ode%PHI
             phidot_de(i) = ode%PHIDOT
             call getweff(3, lna(i), ode%y, wp1_de(i), wp1eff_cdm(i), wp1eff_de(i), lnrho_de(i))                                       
          enddo
          de%Omega = (de%DE_V(ode%PHI) + ode%PHIDOT**2/2.d0)/3.d0
          cdm%Omega = this%Omega_k() - de%Omega
          lnrho_cdm = lnrho_cdm - lnrho_cdm(coop_default_array_size)
          lnrho_de = lnrho_de - lnrho_de(coop_default_array_size)
          call cdm%flnrho%init(coop_default_array_size, coop_min_scale_factor, xmax = coop_scale_factor_today, f = lnrho_cdm, fleft = lnrho_cdm(1), fright = 0.d0, slopeleft = -3.d0, sloperight = -3.d0*wp1eff_cdm(coop_default_array_size), xlog = .true.,  method = COOP_INTERPOLATE_LINEAR, check_boundary = .false., name = "CDM \ln\rho(a) ratio")
          call cdm%fwp1eff%init(coop_default_array_size, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, f = wp1eff_cdm, fleft = wp1eff_cdm(1), fright = wp1eff_cdm(coop_default_array_size), slopeleft = 0.d0, sloperight = 0.d0, xlog = .true., method = COOP_INTERPOLATE_LINEAR, check_boundary = .false., name = "CDM 1+w_eff(a)")          
          
          call de%flnrho%init(coop_default_array_size, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, f = lnrho_de, fleft = lnrho_de(1), fright = 0.d0, slopeleft = 0.d0, sloperight = -3.d0*wp1eff_de(coop_default_array_size), xlog = .true., method = COOP_INTERPOLATE_LINEAR, check_boundary = .false., name = "DE \ln\rho(a) ratio")
          call de%fwp1%init(coop_default_array_size, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, f = wp1_de, fleft = wp1_de(1), fright = wp1_de(coop_default_array_size), slopeleft = 0.d0, sloperight = 0.d0, xlog = .true., method = COOP_INTERPOLATE_LINEAR, check_boundary = .false., name = "DE 1+w(a)")
          call de%fwp1eff%init(coop_default_array_size, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, f = wp1eff_de, fleft = wp1eff_de(1), fright = wp1eff_de(coop_default_array_size), slopeleft = 0.d0, sloperight = 0.d0, xlog = .true., method = COOP_INTERPOLATE_LINEAR, check_boundary = .false., name = "DE 1+w_eff(a)")          
          call de%fDE_phi%init(coop_default_array_size, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, f = phi_de, fleft = phi_de(1), fright = phi_de(coop_default_array_size), slopeleft = 0.d0, sloperight = 0.d0, xlog = .true., method = COOP_INTERPOLATE_LINEAR, check_boundary = .false., name = "DE phi(a)")
          call de%fDE_phidot%init(coop_default_array_size, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, f = phidot_de, fleft = phidot_de(1), fright = phidot_de(coop_default_array_size), slopeleft = 0.d0, sloperight = 0.d0, xlog = .true., method = COOP_INTERPOLATE_LINEAR, check_boundary = .false., name = "DE \dot\phi(a)")          
          exit
       else
          lnc_mid = (lnc_upper + lnc_lower)/2.d0
300       call get_ini(lna_ini, lnc_mid)
          do i = 1, nsteps
             call evolve_to(lna_coarse(i))
          enddo
          Ot_mid = (exp(ode%LN_RHOM) + de%DE_V(ode%PHI) + ode%PHIDOT**2/2.d0)/3.d0 - this%Omega_k()
       endif

       if(Ot_mid .gt. 0.d0)then
          lnc_upper = lnc_mid
          Ot_upper = Ot_mid
       else
          lnc_lower = lnc_mid
          Ot_lower = Ot_mid
       endif
    enddo
    call this%add_species(cdm)
    call this%add_species(de)    
    call ode%free()
    call ode_tc%free()
  contains

    
    subroutine get_ini(lna, lnc)
      COOP_REAL lna, a, t, ddphi, hubble, norm, qnow, lnc, y(3)
      a = exp(lna)
      LN_RHOM = lnrhom0 - 3.d0*lna
      hubble = sqrt((this%rhoa4(a)/a**4+exp(LN_RHOM))/3.d0)            
      if(de%DE_tracking_n .gt. 0.d0)then
         t = 1.d0/hubble/2.d0
         PHI = exp(lnc)*t**beta
         PHIDOT = beta*PHI/t
         !!set the normalization of V
         qnow = de%DE_Q(PHI)      
         ddphi = (beta-1.d0)*PHIDOT/t
         norm = ((ddphi + 3.d0*hubble*PHIDOT) + qnow * exp(LN_RHOM)) / de%DE_tracking_n * PHI**(de%DE_tracking_n + 1.d0)
         de%DE_lnV0 = log(norm)     
         call ode%set_initial_conditions(lna, y)
         tight_coupling = ( (hubble*PHIDOT) .lt. 5.d-3*qnow*exp(LN_RHOM) )
         if(tight_coupling)then
            call ode_tc%set_initial_conditions(lna, (/ y(3) /) )
         endif
      else
         PHI = 0.d0
         PHIDOT = 0.d0
         !!set the normalization of V
         de%DE_lnV0 = lnc
         call ode%set_initial_conditions(lna, y)
         tight_coupling = .false.
      endif
    end subroutine get_ini


    subroutine evolve_to(lnaend)
      COOP_REAL::lnaend
      if(tight_coupling)then
         call ode_tc%evolve(fcn_tight_coupling, lnaend)
         ode%x = ode_tc%x
         ode%LN_RHOM = ode_tc%y(1)
         call set_tight_coupling(ode%x)
      else
         call ode%evolve(fcn, lnaend)
      endif
    end subroutine evolve_to


    subroutine fcn(n, lna, y, yp)
      COOP_INT::n
      COOP_REAL::y(n), yp(n)
      COOP_REAL::H, lna, Q, a, rhoa4, a4, qnow, rhoma4, Vp, V, KE
      a = exp(lna)
      a4 = a**4
      qnow = de%DE_Q(PHI)
      rhoma4 = exp(LN_RHOM + 4.d0*lna)
      rhoa4 =  this%rhoa4(a)
      V = de%DE_V(PHI)
      KE =  PHIDOT**2/2.d0
      rhoa4 = rhoa4 + rhoma4 + (V + KE )*a4
      H = sqrt(rhoa4/3.d0/a4)
      Vp = de%DE_dlnVdphi(PHI)*V
      PHI_PRIME = PHIDOT / H
      PHIDOT_PRIME = - 3.d0*PHIDOT - Vp/H  - (rhoma4/a4/H) * qnow
      LN_RHOM_PRIME = - 3.d0 + qnow * PHI_PRIME
    end subroutine fcn


    subroutine getweff(n, lna, y, wp1_de, wp1eff_cdm, wp1eff_de, lnrhode)
      COOP_INT::n
      COOP_REAL::y(n), yp(n)
      COOP_REAL::H, lna, Q, a, pa4, rhoa4, a4, qnow, rhoma4, rhodea4, V, KE, wp1eff_cdm, wp1eff_de, lnrhode, wp1_de
      a = exp(lna)
      a4 = a**4
      qnow = de%DE_Q(PHI)
      rhoma4 = exp(LN_RHOM + 4.d0*lna)
      call this%get_pa4_rhoa4(a, pa4, rhoa4)
      V = de%DE_V(PHI)
      KE =  PHIDOT**2/2.d0
      lnrhode = log(V+KE)
      rhodea4 = (V + KE )*a4
      rhoa4 = rhoa4 + rhoma4 + rhodea4
      H = sqrt(rhoa4/3.d0/a4)
      wp1eff_cdm = -qnow*PHIDOT/H/3.d0
      wp1_de = 2.d0*KE/(KE+V) 
      wp1eff_de = wp1_de - wp1eff_cdm*rhoma4/rhodea4
      wp1eff_cdm = wp1eff_cdm + 1.d0
    end subroutine getweff
    


    subroutine fcn_tight_coupling(n,lna, y, yp)
      COOP_INT n
      COOP_REAL::lna, a, y(n), yp(n), V, rhoa4, phi_eq_dot, phi_eq, rhom, rhotot, hubble, qnow, Qpcorr
      a = exp(lna)
      phi_eq = 1.d-5 !!initial random guess, suppose Q is not sensitive to phi
      qnow = de%DE_Q(phi_eq)      
      phi_eq  =  exp((de%DE_lnV0+log(de%DE_tracking_n/qnow) - y(1))/(de%DE_tracking_n+1.d0))
      qnow = de%DE_Q(phi_eq)      
      phi_eq  =  exp((de%DE_lnV0+log(de%DE_tracking_n/qnow) - y(1))/(de%DE_tracking_n+1.d0))
      qnow = de%DE_Q(phi_eq)            
      Qpcorr = phi_eq*de%DE_dlnQdphi(phi_eq)/(de%DE_tracking_n+1.d0)
      rhoa4 =  this%rhoa4(a)
      
      V = de%DE_V(phi_eq)      
      rhom = exp(y(1))
      rhotot= rhoa4/a**4 + rhom + V

      hubble = sqrt(rhotot/3.d0)
      phi_eq_dot = (3.d0*hubble)/(de%DE_tracking_n + 1.d0) * phi_eq/(1.d0 + Qpcorr)
      hubble = sqrt((rhotot + phi_eq_dot**2/2.d0)/3.d0)      
      phi_eq_dot = (3.d0*hubble - qnow * phi_eq_dot)/(de%DE_tracking_n + 1.d0) * phi_eq/(1.d0 + Qpcorr)
      hubble = sqrt((rhotot + phi_eq_dot**2/2.d0)/3.d0)      
      yp(1) = -3.d0 + qnow*phi_eq_dot / hubble
    end subroutine fcn_tight_coupling


    subroutine set_tight_coupling(lna)
      COOP_REAL::lna, a,  V, pa4, rhoa4, phi_eq_dot, phi_eq, rhom, rhotot, hubble, ptot, HdotbyHsq, qnow, Qpcorr
      a = exp(lna)
      phi_eq = 1.d-5 !!initial random guess, suppose Q is not sensitive to phi
      qnow = de%DE_Q(phi_eq)      
      phi_eq  =  exp((de%DE_lnV0+log(de%DE_tracking_n/qnow) - ode%LN_RHOM)/(de%DE_tracking_n+1.d0))
      qnow = de%DE_Q(phi_eq)      
      phi_eq  =  exp((de%DE_lnV0+log(de%DE_tracking_n/qnow) - ode%LN_RHOM)/(de%DE_tracking_n+1.d0))
      qnow = de%DE_Q(phi_eq)            
      Qpcorr = phi_eq*de%DE_dlnQdphi(phi_eq)/(de%DE_tracking_n+1.d0)
      rhoa4 =  this%rhoa4(a)
      
      V = de%DE_V(phi_eq)      
      rhom = exp(ode%LN_RHOM)
      rhotot= rhoa4/a**4 + rhom + V

      hubble = sqrt(rhotot/3.d0)
      phi_eq_dot = (3.d0*hubble)/(de%DE_tracking_n + 1.d0) * phi_eq/(1.d0 + Qpcorr)
      hubble = sqrt((rhotot + phi_eq_dot**2/2.d0)/3.d0)      
      phi_eq_dot = (3.d0*hubble - qnow * phi_eq_dot)/(de%DE_tracking_n + 1.d0) * phi_eq/(1.d0 + Qpcorr)
      hubble = sqrt((rhotot + phi_eq_dot**2/2.d0)/3.d0)      
      ode%PHI = phi_eq - 3.d0 * hubble**2/V/de%DE_tracking_n/(de%DE_tracking_n+1.d0)**2*(3.d0/(de%DE_tracking_n +1.d0) + 3.d0 + HdotbyHsq)*phi_eq**3
      ode%PHIDOT = phi_eq_dot
      ode%ind = 1
      tight_coupling = (hubble * ode%PHIDOT .lt. 1.d-3*qnow*rhom)
    end subroutine set_tight_coupling
    
  end subroutine coop_background_add_coupled_DE

  subroutine coop_background_add_coupled_DE_with_w(this, Omega_c, epsilon_s, epsilon_inf, zeta_s, Q, dlnQdphi, err)
    COOP_REAL::a_piv, a_zeta, aeq, a_ini
    COOP_REAL,parameter::max_dUdphi = coop_sqrt3 !!set by epsilon = epsilon_CDM
    class(coop_cosmology_background)::this
    type(coop_species)::de
    COOP_REAL::Omega_c, epsilon_s, epsilon_inf, zeta_s, Q, dlnQdphi
    COOP_REAL::tracking_n, dUdphi, d2Udphi2
    COOP_INT::i, index_de, index_cdm
    COOP_REAL::rhoc0, wp1, wp1_zeta, delta_wp1, n_low, n_high, dUdphi_low, dUdphi_high, phi_eq, d2Udphi2_trial, diffwp1, mindiff, step
    COOP_INT::err
    
    err = 0    
    index_cdm = this%num_species + 1
    index_de = this%num_species+2
    rhoc0 = 3.d0*omega_c
    aeq = coop_de_aeq_fitting(omega_lambda = this%Omega_k() - omega_c, epsilon_s = epsilon_s, epsilon_inf = epsilon_inf, zeta_s = zeta_s)
    if(aeq .gt. 1.2d0)then
       err = 1
       return
    endif
    a_ini = aeq/8.d0
    a_piv = min(aeq*1.303, coop_scale_factor_today)  !!where F_2(a/a_eq) = 0
    a_zeta = aeq*0.8d0  
    de = coop_de_quintessence(Omega_Lambda = this%Omega_k() - omega_c, epsilon_s = epsilon_s, epsilon_inf = epsilon_inf, zeta_s = zeta_s)
    
    !!search for tracking_n
    wp1 = de%wp1ofa(a_ini)    
    n_high = epsilon_inf*(4.d0/3.d0)/(1.d0-(2.d0/3.d0)*epsilon_inf)
    if(n_high .le. coop_min_de_tracking_n)then
       tracking_n = coop_min_de_tracking_n
    else
       n_low = max(n_high/3.d0, coop_min_de_tracking_n)
       do 
          tracking_n = (n_high + n_low)/2.d0
          call coop_background_add_coupled_DE(this, omega_c = omega_c, tracking_n = tracking_n, Q = Q, dlnQdphi = dlnQdphi, dUdphi = 0.d0, d2Udphi2 = 0.d0 )
          if( wp1_eff(a_ini) .gt. wp1 ) then
             n_high = tracking_n
          else
             n_low = tracking_n
          endif
          if(n_high - n_low .lt. 5.d-3)exit
          call this%delete_species(index_de)
          call this%delete_species(index_cdm)       
       enddo
    endif
    !!search for dUdphi
    wp1 = de%wp1ofa(a_piv)
    delta_wp1  = wp1_eff(a_piv) - wp1
    if(delta_wp1 .gt. 0.01d0)then
       step = 0.05
       dUdphi = step
       call this%delete_species(index_de)
       call this%delete_species(index_cdm)       
       call coop_background_add_coupled_DE(this, omega_c = omega_c, tracking_n = tracking_n, Q = Q, dlnQdphi = dlnQdphi, dUdphi =dUdphi, d2Udphi2 = 0.d0 )
       diffwp1  = wp1_eff(a_piv) - wp1
       if(abs(diffwp1) .lt. 0.01d0)goto 30
       if(diffwp1 .ge. delta_wp1 .or. diffwp1 .lt. 0.d0 )then
          err = 2
          return
       endif
       step = max(min(diffwp1/(delta_wp1-diffwp1+1.d-3)*step/2.d0, 0.5d0), 0.02d0)
       dUdphi = dUdphi + step
       delta_wp1 = diffwp1       
       do 
          call this%delete_species(index_de)
          call this%delete_species(index_cdm)       
          call coop_background_add_coupled_DE(this, omega_c = omega_c, tracking_n = tracking_n, Q = Q, dlnQdphi = dlnQdphi, dUdphi =dUdphi, d2Udphi2 = 0.d0 )
          diffwp1  = wp1_eff(a_piv) - wp1
          if(abs(diffwp1) .lt. 0.01d0)goto 30
          if(diffwp1 .ge. delta_wp1 .or. diffwp1 .lt. 0.d0 )then
             err = 3
             return
          endif
          step = max(min(diffwp1/(delta_wp1-diffwp1+1.d-3)*step/2.d0, 0.5d0), 0.02d0)
          dUdphi = dUdphi + step
          delta_wp1 = diffwp1
          if(dUdphi .gt. 30.d0)then
             err = 4
             return
          endif
       enddo          
    elseif(delta_wp1 .lt. -0.01d0)then
       dUdphi_high = 0.d0
       dUdphi_low = -max_dUdphi
       call this%delete_species(index_de)
       call this%delete_species(index_cdm)       
       call coop_background_add_coupled_DE(this, omega_c = omega_c, tracking_n = tracking_n, Q = Q, dlnQdphi = dlnQdphi, dUdphi =dUdphi_low, d2Udphi2 = 0.d0 )
       if(wp1_eff(a_piv) .lt. wp1)then
          err = 5
          return
       endif
       do  while(dUdphi_high - dUdphi_low .gt. 0.01d0)
          call this%delete_species(index_de)
          call this%delete_species(index_cdm)       
          dUdphi = (dUdphi_high + dUdphi_low)/2.d0
          call coop_background_add_coupled_DE(this, omega_c = omega_c, tracking_n = tracking_n, Q = Q, dlnQdphi = dlnQdphi, dUdphi =dUdphi, d2Udphi2 = 0.d0 )
          if(wp1_eff(a_piv) .gt. wp1)then
             dUdphi_low = dUdphi
          else
             dUdphi_high = dUdphi
          endif
       enddo
    else
       dUdphi  = 0.
    endif

30  if(abs(zeta_s) .lt. 2.d-2 .or. abs(epsilon_s - 2.d0*epsilon_inf).lt. 0.02)return
    
    !!serach for d2Udphi2
    wp1_zeta = de%wp1ofa(a_zeta)
    mindiff = (wp1_eff(a_zeta) - wp1_zeta)**2 + (wp1_eff(a_piv) - wp1)**2
    phi_eq = this%species(index_de)%de_Phi(aeq)
    d2Udphi2 = 0.d0
    do i=-5, 5
       if(i.eq.0)cycle
       call this%delete_species(index_de)
       call this%delete_species(index_cdm)       
       d2Udphi2_trial = i*0.05d0
       call coop_background_add_coupled_DE(this, omega_c = omega_c, tracking_n = tracking_n, Q = Q, dlnQdphi = dlnQdphi, dUdphi =dUdphi - d2Udphi2_trial*phi_eq, d2Udphi2 = d2Udphi2_trial )
       diffwp1 = (wp1_eff(a_zeta) - wp1_zeta)**2 + (wp1_eff(a_piv) - wp1)**2
       if(diffwp1 .lt. mindiff)then
          mindiff = diffwp1
          d2Udphi2 = d2Udphi2_trial
       endif
    enddo
    dUdphi = dUdphi - d2Udphi2*phi_eq
    call this%delete_species(index_de)
    call this%delete_species(index_cdm)
    call coop_background_add_coupled_DE(this, omega_c = omega_c, tracking_n = tracking_n, Q = Q, dlnQdphi = dlnQdphi, dUdphi =dUdphi, d2Udphi2 = d2Udphi2)    

  contains

    function wp1_eff(a) result(wp1)
      COOP_REAL::a, wp1
      wp1 = (1.d0 + this%species(index_DE)%pa2(a)/(this%species(index_DE)%rhoa2(a)+ this%species(index_CDM)%rhoa2(a) - rhoc0/a))
    end function wp1_eff       
  end subroutine coop_background_add_coupled_DE_with_w




  subroutine coop_background_add_EFT_DE(this, wp1, alpha_M, err)
    class(coop_cosmology_background)::this
    type(coop_function)::wp1, alpha_M
    type(coop_species)::de
    COOP_INT::i, err, j
    COOP_REAL_ARRAY::lnrho,  wp1eff, M2
    COOP_REAL::lna, lnamin, lnamax, dlna, alpha_l, a_l, a_r, alpha_r, wp1_l, wp1_r, omega_de, om_l, om_r, rhoa4de, rhotot_l, rhotot_r, step
    err = 0
    de%name = "Dark Energy"
    de%genre = COOP_SPECIES_EFT
    de%fwp1 = wp1
    omega_de = this%Omega_k()    
    de%Omega = omega_de
    lnamin = log(coop_min_scale_factor)
    lnamax = log(coop_scale_factor_today)
    dlna = (lnamax-lnamin)/(coop_default_array_size-1)
    lnrho(coop_default_array_size) = log(3.d0*omega_de)
    wp1_r = wp1%eval(coop_scale_factor_today)
    alpha_r = alpha_M%eval(coop_scale_factor_today)
    om_r = omega_de
    wp1eff(coop_default_array_size) =wp1_r - alpha_r/3.d0/om_r
    a_r = coop_scale_factor_today
    lna = lnamax
    rhotot_r = 3.d0
    step = (1.5d0*dlna)
    M2(coop_default_array_size)=0.d0
    do i=coop_default_array_size-1, 1, -1
       lna = lna - dlna              
       a_l = exp(lna)
       alpha_l = alpha_M%eval(a_l)
       wp1_l = wp1%eval(a_l)
       wp1eff(i) = wp1_l - alpha_l/3.d0/om_r
       rhotot_l =  this%rhoa4(a_l)              
       do j=1,3
          lnrho(i) = lnrho(i+1) + (wp1eff(i)+wp1eff(i+1))*step
          rhoa4de = exp(lnrho(i)+4.d0*lna)
          om_l = rhoa4de/(rhoa4de + rhotot_l)
          if(om_l .lt. 1.d-50)then
             lnrho(1:i-1) = lnrho(i)
             wp1eff(1:i-1) = wp1eff(i)
             M2(1:i-1) = M2(i)
             goto 100
          endif
          wp1eff(i) = wp1_l - alpha_l/3.d0/om_l
       enddo
       rhotot_l = (rhotot_l + rhoa4de)/a_l**4
       if(rhotot_l .lt. rhotot_r)then  !!
          err = 1
          return
       endif
       M2(i) = M2(i+1) - (alpha_l+alpha_r)
       rhotot_r = rhotot_l
       alpha_r = alpha_l
    enddo
100 M2 = exp(M2*(dlna/2.d0))
    lnrho = lnrho - lnrho(coop_default_array_size)
    call de%fwp1eff%init(coop_default_array_size, coop_min_scale_factor, coop_scale_factor_today, wp1eff, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = "DE 1+w_eff(a)")
    call de%flnrho%init(coop_default_array_size,coop_min_scale_factor, coop_scale_factor_today, lnrho, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = "DE ln rho_ratio")
    call de%flnrho%set_boundary(slopeleft = -3.d0*wp1eff(1), sloperight = -3.d0*wp1eff(coop_default_array_size))    
    de%cs2 = 0.d0
    call this%add_species(de)
#if DO_EFT_DE    
    this%f_alpha_M = alpha_M
    call this%f_M2%init(coop_default_array_size, coop_min_scale_factor, coop_scale_factor_today, M2, method = COOP_INTERPOLATE_LINEAR, xlog = .true., check_boundary = .false., name = "EFT M^2")
#endif    
  end subroutine coop_background_add_EFT_DE  

end module coop_background_mod
