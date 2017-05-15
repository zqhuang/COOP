module coop_background_mod
  use coop_wrapper_utils
  implicit none
#include "constants.h"

  private

#if DO_EFT_DE
  !!screening: true; non-screening: false
  logical,parameter::coop_eft_de_normalize_early = .true.
#endif  


  public::coop_baryon, coop_cdm, coop_DE_lambda, coop_DE_w0, coop_DE_w0wa, coop_DE_quintessence, coop_radiation, coop_neutrinos_massless, coop_neutrinos_massive, coop_de_w_quintessence, coop_de_wp1_quintessence, coop_de_wp1_coupled_quintessence, coop_background_add_coupled_DE,  coop_background_add_EFT_DE,  coop_background_add_EFT_DE_with_effective_w, coop_de_aeq_fitting, coop_de_alpha_invh2, coop_de_alpha_instant, coop_de_general, coop_de_alpha_constructor, coop_de_construct_alpha_from_cs2, coop_background_add_coupled_DE_with_potential, coop_convert_fofR_to_Vofphi

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
    call fw0wa%init_polynomial( (/ 1.d0+w0+wa, -wa /))
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


  subroutine coop_background_add_EFT_DE_with_effective_w(this, effective_wp1, err)
    COOP_INT,parameter::narr = 10000
    class(coop_cosmology_background)::this
    type(coop_function)::effective_wp1
    type(coop_species)::de, deeff
    COOP_INT::i, err, j
    COOP_REAL,dimension(narr)::lnrho,  wp1, wp1eff
    COOP_REAL::  a, rhoa4de,ppra4de, ppra4tot, rhoa4tot, rhoa4de_bg, M2, rhom0, wp1_bg, da
#if DO_EFT_DE
    err = 0    
    rhom0=this%rhoa4(1.d0)
    if(rhom0 .gt. 3.d0)then
       err = 1
       return
    endif
    !!define effective DE
    call deeff%init(name = "Dark Energy", id = 5, Omega = (3.d0-rhom0)/3.d0/this%Mpsq0, genre = COOP_SPECIES_FLUID, fwp1 = effective_wp1)
    deeff%Mpsq0 = this%Mpsq0

    de%name = "Dark Energy"
    de%genre = COOP_SPECIES_EFT


    da = (coop_scale_factor_today-coop_min_scale_factor)/(narr-1.d0)
    a = coop_min_scale_factor
    do i=1, narr
       M2 = this%Mpsq(a)
       wp1_bg = deeff%wp1ofa(a)
       rhoa4de_bg = deeff%rhoa4(a)       
       call this%get_ppra4_rhoa4(a, ppra4tot, rhoa4tot)
       rhoa4de = rhoa4tot*(M2 - 1.d0) + rhoa4de_bg*M2
       ppra4de = ppra4tot*(M2 - 1.d0) + wp1_bg * rhoa4de_bg*M2
       if(rhoa4de .le. 0.d0)then
          err = 1  !!negative energy flag
          return
       else
          lnrho(i) = log(max(rhoa4de/a**4, 1.d-99)) 
          wp1(i) = ppra4de/rhoa4de
       endif
       wp1eff(i) = wp1_bg*rhoa4de_bg/rhoa4de - (this%alpha_M(a)*(rhoa4tot+rhoa4de_bg)*M2/3.d0 - (M2-1.d0)*(ppra4tot + wp1_bg * rhoa4de_bg))/rhoa4de
       a = a + da
    enddo
    de%Omega = exp(lnrho(narr))/(3.d0*this%Mpsq0)
    de%Mpsq0 = this%Mpsq0

    lnrho = lnrho - lnrho(narr)
    call de%fwp1%init(narr, coop_min_scale_factor, coop_scale_factor_today, wp1, method = COOP_INTERPOLATE_LINEAR, check_boundary = .false., name = "DE 1+w(a)")    
    call de%fwp1eff%init(narr, coop_min_scale_factor, coop_scale_factor_today, wp1eff, method = COOP_INTERPOLATE_LINEAR, check_boundary = .false., name = "DE 1+w_eff(a)")
    call de%flnrho%init(narr,coop_min_scale_factor, coop_scale_factor_today, lnrho, method = COOP_INTERPOLATE_LINEAR, check_boundary = .false., name = "DE ln rho_ratio")
    call de%flnrho%set_boundary(slopeleft = 0.d0, sloperight = -3.d0*wp1eff(narr))    
    de%cs2 = 0.d0
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
    COOP_INT::i, err, j, nsteps
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
    lnrho(narr) = log(3.d0*omega_de*this%Mpsq0)
    wp1_r = wp1%eval(coop_scale_factor_today)
    alpha_r = this%alpha_M(coop_scale_factor_today)
    om_r = omega_de
    wp1eff(narr) =wp1_r - alpha_r/3.d0/om_r
    a_r = coop_scale_factor_today
    lna = lnamax
    rhotot_r = 3.d0*this%Mpsq0
    rhoa4de_r = om_r*rhotot_r
    step = (1.5d0*dlna)
    do i=narr-1, 1, -1
       lna = lna - dlna              
       a_l = exp(lna)
       alpha_l = this%alpha_M(a_l)
       wp1_l = wp1%eval(a_l)
       rhotot_l =  this%rhoa4(a_l)
       om_l = rhoa4de_r/(rhotot_l + rhoa4de_r)  !!first assuming rho_de constant
       wp1eff(i) = wp1_l - alpha_l/3.d0/om_l
       nsteps = min(nint(log10(1.d9*om_l)), 5)
       do j=1, nsteps
          lnrho(i) = lnrho(i+1) + (wp1eff(i)+wp1eff(i+1))*step
          rhoa4de_l = exp(lnrho(i)+4.d0*lna)
          om_l = rhoa4de_l/(rhoa4de_l + rhotot_l)
          wp1eff(i) = wp1_l - alpha_l/3.d0/om_l
          if(om_l .lt. 1.d-8)then
             lnrho(1:i) = lnrho(i+1)
             wp1eff(1:i) = wp1eff(i+1)
             goto 100
          endif
       enddo
       rhotot_l = (rhotot_l + rhoa4de_l)/a_l**4
!!$       if(rhotot_l .lt. rhotot_r)then  !!
!!$          err = 1
!!$          return
!!$       endif
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

  function coop_de_alpha_powerlaw(a, arg) result(alpha)
    COOP_REAL::a, alpha
    type(coop_arguments)::arg  !!arg%r = (alpha_0, index)
    alpha = arg%r(1)*a**arg%r(2)
  end function coop_de_alpha_powerlaw
  

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

  function coop_de_alpha_model_omega(a, arg) result(alpha)
    type(coop_arguments)::arg !! arg%r = (/ alpha0, omegam, w /)
    COOP_REAL::alpha, a, rhodea3
    rhodea3 =  (1.d0 - arg%r(2))*a**(-3.d0*arg%r(3))
    alpha = arg%r(1) * rhodea3/(arg%r(2)  + rhodea3)/(1.d0 - arg%r(2))
  end function coop_de_alpha_model_omega

  subroutine coop_de_construct_alpha_from_cs2(omegam, w, cs2, r_B, r_H, r_M, r_T,  alpha_B, alpha_H, alpha_K, alpha_M, alpha_T, sucess)
    COOP_REAL, parameter::alpha0_min = 1.d-4
    COOP_REAL::omegam, w, cs2, r_B, r_M, r_T, r_H, alpha_l
    type(coop_function)::alpha_B, alpha_H, alpha_K, alpha_M, alpha_T
    COOP_REAL::  hdotbyhsq, p, omegal, cs2now, cs2derv
    logical::sucess
    omegal = 1.d0 - omegam
    hdotbyhsq = -1.5d0*(omegam + omegal*(1.d0+w))
    p = -3.d0*w*omegam
    alpha_l = alpha0_min
    cs2now = cs2try(alpha_l)
    if(cs2now .gt. cs2+2.d-3)then
       do while(cs2try(alpha_l*1.001d0) .gt. cs2+2.d-2)
          alpha_l = alpha_l*1.001d0
          if(alpha_l .gt. coop_de_alpha0_max) goto 100
       enddo
       do while(cs2try(alpha_l*1.0001d0) .gt. cs2+2.d-3)
          alpha_l = alpha_l*1.0001d0
          if(alpha_l .gt. coop_de_alpha0_max) goto 100
       enddo
    elseif(cs2now .lt. cs2-2.d-3)then
       do while(cs2try(alpha_l*1.001d0) .lt. cs2-2.d-2)
          alpha_l = alpha_l*1.001d0
          if(alpha_l .gt. coop_de_alpha0_max) goto 100
       enddo
       do while(cs2try(alpha_l*1.0001d0) .lt. cs2-2.d-3)
          alpha_l = alpha_l*1.0001d0
          if(alpha_l .gt. coop_de_alpha0_max) goto 100
       enddo
    endif
    cs2derv = (cs2try(alpha_l + alpha0_min/2) - cs2try(alpha_l - alpha0_min/2))/alpha0_min
    alpha_l = alpha_l + (cs2 - cs2try(alpha_l))/sign(max(abs(cs2derv),1.d-4), cs2derv)
    if(r_B .eq. 0.d0)then
       call alpha_B%init_polynomial( (/ 0.d0 /))
    else
       alpha_B = coop_function_constructor( coop_de_alpha_model_omega, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args= coop_arguments_constructor ( r = (/ alpha_l*r_B,  omegam, w /) ) , name = "alpha" )
    endif
    if(r_H .eq. 0.d0)then
       call alpha_H%init_polynomial( (/ 0.d0 /))
    else
       alpha_H = coop_function_constructor( coop_de_alpha_model_omega, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args= coop_arguments_constructor ( r = (/ alpha_l*r_H,  omegam, w /) ) , name = "alpha" )
    endif

    if(r_M .eq. 0.d0)then
       call alpha_M%init_polynomial( (/ 0.d0 /))
    else
       alpha_M = coop_function_constructor( coop_de_alpha_model_omega, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args= coop_arguments_constructor ( r = (/ alpha_l*r_M,  omegam, w /) ) , name = "alpha" )
    endif

    if(r_T .eq. 0.d0)then
       call alpha_T%init_polynomial( (/ 0.d0 /))
    else
       alpha_T = coop_function_constructor( coop_de_alpha_model_omega, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args= coop_arguments_constructor ( r = (/ alpha_l*r_T,  omegam, w /) ) , name = "alpha" )
    endif
    alpha_K = coop_function_constructor(derived_alphaK,  xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., name = "alpha")
    sucess = .true.
    return
100 sucess = .false.
    return
    
  contains

    function derived_alphaK(a)
      COOP_REAL::a, derived_alphaK
      COOP_REAL::b, h, m, t, bp, hp, rhodea3, hdotbyhsq   
      t = alpha_T%eval(a)
      b= alpha_B%eval(a)
      h = alpha_H%eval(a)
      m = alpha_M%eval(a)
      rhodea3 = omegal*a**(-3.d0*w)
      hdotbyhsq = -1.5d0*( omegam + rhodea3*(1.d0+w))/(omegam + rhodea3)
      bp = alpha_B%derivative(a) * a
      hp = alpha_H%derivative(a) * a
      derived_alphaK =  (  -2.d0*((1.d0+b)*((1.d0+b)*(1.d0+t) - (1.d0+h)*(1.d0+m - hdotbyhsq)-hp)+(1.d0+h)*bp) - 3.d0* (1.d0+h)**2*omegam/(omegam+rhodea3) )/cs2 - 6.d0*b**2
    end function derived_alphaK

    function cs2try(alpha0)
      COOP_REAL::cs2try, alpha0
      COOP_REAL::b, h, m, t
      b = r_B * alpha0
      h = r_H*alpha0
      m = r_M * alpha0
      t = r_T*alpha0
      cs2try = (-2.d0*((1.d0+b)*((1.d0+b)*(1.d0+t) - (1.d0+h)*(1.d0+m - hdotbyhsq)- h*p)+(1.d0+ h )*b*p) + (1.d0+h)**2*(-3.d0*omegam))/alpha0
    end function cs2try
  end subroutine coop_de_construct_alpha_from_cs2

  function coop_de_alpha_constructor(alpha0, genre, pow) result(alpha)
    COOP_REAL, parameter::Omega_m = 0.314d0
    COOP_REAL, parameter::Omega_r = 8.d-5
    COOP_REAL, parameter::lambda = 0.001d0
    COOP_REAL, parameter::delta_z = 0.1d0
    COOP_REAL::alpha0, omegam, omegar
    COOP_REAL, optional::pow
    COOP_UNKNOWN_STRING::genre
    type(coop_function)::alpha
    select case(COOP_LOWER_STR(genre))
    case("omega")
       alpha = coop_function_constructor( coop_de_alpha_invh2, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args= coop_arguments_constructor ( r = (/ alpha0, Omega_m, Omega_r /) ), name = "alpha")
    case("rational")
       if(present(pow))then
          call alpha%init_rational(c_up = (/ pow*alpha0 /), alpha_up = (/ pow /), c_down = (/ 1.d0, alpha0 /), alpha_down = (/ 0.d0, pow /), name = "rational alpha")
       else
          stop "for rational function you must specify the power index"
       endif
    case("powerlaw")
       if(present(pow))then
          if(abs(nint(pow) - pow) .lt. 1.d-6)then
             select case(nint(pow))
             case(0)
                call alpha%init_polynomial( (/ alpha0 /) )
             case(1)
                call alpha%init_polynomial( (/ 0.d0, alpha0 /) )
             case(2)
                call alpha%init_polynomial( (/ 0.d0, 0.d0, alpha0 /) )
             case(3)
                call alpha%init_polynomial( (/ 0.d0, 0.d0, 0.d0, alpha0 /) )
             case(4)
                call alpha%init_polynomial( (/ 0.d0, 0.d0, 0.d0, 0.d0, alpha0 /) )

             case(5)
                call alpha%init_polynomial( (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, alpha0 /) )
             case default
                alpha = coop_function_constructor( coop_de_alpha_powerlaw, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args = coop_arguments_constructor ( r = (/ alpha0, pow /)), name = "alpha function")
             end select
          else
             alpha = coop_function_constructor( coop_de_alpha_powerlaw, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args = coop_arguments_constructor ( r = (/ alpha0, pow /)), name = "alpha function")
          endif
       else
          stop "for powerlaw alpha parametrization you need to specify the power index"
       endif
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


!!the action is 1/2 \int (R - 2 \Lambda + f(R) ) \sqrt{-g} d^4x
!! f(R) -> 0  when R-> infty
!! f'(R) < 0 
!! f''(R) > 0
!! f'(R) -> 0 when R-> infty
!!in this model the coupling Q = 1/sqrt(6) is fixed
  subroutine coop_convert_fofR_to_Vofphi(fofR, Lambda, Vofphi, force_monotonic)
    type(coop_function)::fR, Vphi
    COOP_REAL::Lambda
    type(coop_function)::fofR, Vofphi
    logical,optional::force_monotonic
    COOP_INT, parameter::n = 80000
    COOP_REAL, parameter::Rmin = 0.1d0, Rmax = 1.d99   
    COOP_REAL,dimension(n)::lnR, phi, V
    COOP_REAL::R, dfdR, f
    COOP_INT::i, imax, imin

    call coop_set_uniform(n, lnR, log(Rmin), log(Rmax))
    imin = 1
    imax = n
    do i=1, n
       R = exp(lnR(i))
       dfdR = fofR%derivative(R) 
       if(dfdR .le. -1.d0) then
          imin = i+1
          cycle
       endif
       if(abs(dfdR) .gt. 1.d-5)then
          phi(i) = -sqrt(3.d0/2.d0)*log(1.d0+dfdR)
       else
          phi(i) = -sqrt(3.d0/2.d0)*(dfdR*(1.d0 - dfdR*(1.d0/2.d0 - dfdR/3.d0)))
       endif
       f = fofR%eval(R) 
       !!this is actually (-V+Lambda)
       V(i) = (f + dfdR * (2.d0*Lambda*(dfdR+2.d0) - R) )/2.d0/(1.d0+dfdR)**2 
       if(V(i) .le. 0.d0 .or. V(i) .gt. Lambda) imin = i+1
       if(phi(i) .le. 0.d0 )then
          imax = i-1
          exit
       endif
    enddo
    if(present(force_monotonic))then
       if(force_monotonic)then
          do i= imax-1, imin, -1
             V(i) = max(V(i+1), V(i))
          enddo
       endif
    endif
    if(imax .gt. imin)then
       call Vofphi%init_nonUniform(x = phi(imin:imax), f = V(imin:imax),  xlog = .true., ylog = .true., check_boundary = .true., name = "V(phi)")
       call Vofphi%mult_const(-1.d0)
       call Vofphi%add_const(Lambda)
    else
       stop "well behaved f(R) must have f' < 0 and f''>0 for all R values"
    endif
  end subroutine coop_convert_fofR_to_Vofphi





!!fQ: Q(a)
!!fwp1: effective 1+w(a) that gives the same expansion history in wCDM model
  subroutine coop_background_add_coupled_DE(this, Omega_c, Omega_b, fQ, fwp1, err)
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
    COOP_REAL::Omega_c, omega_b, Omega_cpl
    type(coop_function)::fQ, fwp1
    COOP_INT::err
#if DO_COUPLED_DE
    COOP_INT, parameter::ns = 10000
    type(coop_function)::fwp1effcdm, fwp1de, fwp1effde
    type(coop_ode)::ode
    type(coop_species)::de, cdm, baryon
    COOP_INT::i, i_tc_off
    COOP_REAL::rho_ce, a(ns), lna(ns), y(3, ns), yp(3, ns), lnV(ns), dlnVdphi(ns), dVdphibyH2(ns), wp1de(ns), wp1effde(ns), wp1effcdm(ns), H2a4(ns), m2byH2(ns), dQdphi(ns), Q(ns), intQ(ns),  omde, omc, dlna, tc_w    
    this%Omega_c_bare = Omega_c
    this%Omega_b_bare = Omega_b
    if(fQ%is_zero .or. fwp1%is_zero)then !!not coupled
       call this%add_species(coop_baryon(Omega_b))
       call this%add_species(coop_cdm(Omega_c))
       omde = this%Omega_k()
       call de%init(name = "Dark Energy", id = 5, Omega = omde, genre = COOP_SPECIES_FLUID, fwp1 = fwp1)
       call de%cplde_Q%init_polynomial( (/ 0.d0 /))
       call de%cplde_dQdphi_lna%init_polynomial( (/ 0.d0 /) )
       call de%cplde_intQofphi%init_polynomial( (/ 0.d0 /) )
       if(fwp1%is_zero)then
          call de%cplde_lnV_lna%init_polynomial( (/ log(3.d0*omde) /) )
          call de%cplde_phi_lna%init_polynomial( (/ 0.d0 /) )
          call de%cplde_Vofphi%init_polynomial( (/ 3.d0*omde /) )
          call de%cplde_phi_prime_lna%init_polynomial( (/ 0.d0 /) )
          call de%cplde_dVdphibyH2_lna%init_polynomial( (/ 0.d0 /))
          call de%cplde_m2byH2_lna%init_polynomial( (/ 1.d30 /) )
       else
          call coop_set_uniform(ns, lna, log(coop_min_scale_factor), log(coop_scale_factor_today))
          a = exp(lna)
          dlna = lna(2) - lna(1)
          !$omp parallel do
          do i=1, ns
             lnV(i) = de%density(a(i))*(1.d0-de%wp1ofa(a(i))/2.d0)
             y(2, i) = de%density(a(i))*de%wp1ofa(a(i))
             H2a4(i) = (this%rhoa4(a(i)) + de%rhoa4(a(i)))/3.d0
          enddo
          !$omp end parallel do
          if(any(y(2, :).lt. 0.d0) .or. any(lnV .le. 0.d0))then
             err = -1  !!negative potential or kinetic energy
             return
          endif
          lnV = log(lnV)
          y(2, :) = sqrt(y(2, :))         
          yp(3, :) = y(2,:)*a**2/sqrt(H2a4)
          y(3, 1) = 0.d0
          do i=2, ns
             y(3, i) = y(3, i-1)+(yp(3, i)+yp(3, i-1))
          enddo
          y(3, :) = y(3, :)*(dlna/2.d0)
          !$omp parallel do 
          do i= 2, ns-1
             dlnVdphi(i) =  (lnV(i+1)-lnV(i-1))/max(2.d0*yp(3, i)*dlna, 1.d-20)
          enddo
          !$omp end parallel do       
          dlnVdphi(ns) = dlnVdphi(ns-1)    
          dlnVdphi(1) = dlnVdphi(2)
          dVdphibyH2 = dlnVdphi * exp(lnV + 4.d0*lna) / H2a4 
          do i = 3, ns-2
             m2byH2(i) = ((dlnVdphi(i+1) - dlnVdphi(i-1))/max(2.d0*yp(3, i)*dlna, 1.d-20) + dlnVdphi(i)**2) * exp(lnV(i) + 4.d0*lna(i))/h2a4(i)
          enddo
          m2byH2(1:2) = m2byH2(3)
          m2byH2(ns-1:ns) = m2byH2(ns-2)
          call de%cplde_lnV_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = lnV, method = COOP_INTERPOLATE_LINEAR, name = "lnV(lna)")
          call de%cplde_phi_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = y(3,:), method = COOP_INTERPOLATE_LINEAR, name = "phi")
          call de%cplde_Vofphi%init_NonUniform( x=y(3,:), f = exp(lnV), ylog = .true., name = "V(phi)")
          call de%cplde_phi_prime_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = yp(3,:), method = COOP_INTERPOLATE_LINEAR, name="d phi / d lna")          
          call de%cplde_dVdphibyH2_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = dVdphibyH2, method = COOP_INTERPOLATE_LINEAR, name = "d V / d phi / H^2")
          call de%cplde_m2byH2_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = m2byH2, method = COOP_INTERPOLATE_LINEAR, name = " m^2/H^2 ")
       endif
       call this%add_species(de)
       call de%free()
       err = 0
       return
    endif


    if(this%baryon_is_coupled)then
       Omega_cpl = Omega_c+Omega_b
    else
       Omega_cpl = Omega_c
       call this%add_species(coop_baryon(Omega_b))
    endif
    err = 0
    i_tc_off = 1

    if(this%Omega_k() .le. Omega_cpl)then
       err = -1
       return
    endif
    rho_ce = omega_cpl * 3.d0
    call coop_set_uniform(ns, lna, log(coop_min_scale_factor), log(coop_scale_factor_today))
    a = exp(lna)
    dlna = lna(2) - lna(1)
    call ode%init(n=3, method = COOP_ODE_RK4)
    y(1, 1) = log(3.d0*(this%Omega_k() - Omega_cpl)) + 3.d0*coop_integrate(wp1_eval, lna(1), lna(ns))
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
    if(all(wp1effcdm.eq.0.d0))then
       call cdm%init(genre=COOP_SPECIES_FLUID, name = "CDM", id = 1, Omega = omega_c, w = 0.d0)
    else
       call fwp1effcdm%init(n = ns, xmin=coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., f = wp1effcdm, method = COOP_INTERPOLATE_LINEAR, name = "CDM 1+w", check_boundary = .false.)
       call cdm%init(genre=COOP_SPECIES_FLUID, name = "CDM", id = 1, Omega = omc*omega_c/omega_cpl, w = 0.d0, fwp1eff = fwp1effcdm)
    endif
    call fwp1de%init(n = ns, xmin=coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., f = wp1de, method = COOP_INTERPOLATE_LINEAR, name = "Dark Energy 1+w", check_boundary = .false.)
    call fwp1effde%init(n = ns, xmin=coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., f = wp1effde, method = COOP_INTERPOLATE_LINEAR, name= "Dark Energy effective 1+w", check_boundary = .false.)
    call de%init(name = "Dark Energy", id = 5, Omega = omde, genre = COOP_SPECIES_FLUID, fwp1 = fwp1de, fwp1eff =fwp1effde)
    
    de%cplde_Q = fQ
    !$omp parallel do
    do i=1, ns
       Q(i) = fQ%eval(a(i))
       dQdphi(i)=fQ%derivative(a(i))*a(i)/max(yp(3, i), min_phi_prime)
    enddo
    !$omp end parallel do
    intQ(1) = Q(1)*y(3,1)/2.d0
    do i=2, ns
       intQ(i) = intQ(i-1) + (Q(i)+Q(i-1))*(y(3,i) - y(3,i-1))/4.d0
    enddo
    
    call de%cplde_dQdphi_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = dQdphi, method = COOP_INTERPOLATE_LINEAR, name = "dQ / d phi")
    call de%cplde_lnV_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = lnV, method = COOP_INTERPOLATE_LINEAR, name = "ln V")
    call de%cplde_phi_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = y(3,:), method = COOP_INTERPOLATE_LINEAR, name = "phi")
    call de%cplde_intQofphi%init_NonUniform(x = y(3,:), f = intQ, name = "Q(phi)")
    call de%cplde_Vofphi%init_NonUniform( x=y(3,:), f = exp(lnV), ylog = .true., name = "V(phi)")

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
    m2byH2(ns-1:ns) = m2byH2(ns-2)
    call de%cplde_m2byH2_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = m2byH2, method = COOP_INTERPOLATE_LINEAR, name = " m^2/H^2 ")

    call this%add_species(cdm)    
    if(this%baryon_is_coupled )then
       cdm%Omega = omc*omega_b/omega_cpl
       cdm%name = "BARYON"
       call this%add_species(cdm)
    endif

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



!!Vofphi is V(phi) but not normalized
!!intQofphi is \int_0^phi Q(phi) d phi
!!norm is the output that normalizes Vofphi
!!err = 0 if done without a problem
  subroutine coop_background_add_coupled_DE_with_potential(this, Omega_c, Omega_b, Vofphi, intQofphi, err, normalize_V)
    !!err = 0: success
    !!err = -1:negative Omega_L
    class(coop_cosmology_background)::this
    COOP_REAL::Omega_c, Omega_b, Omega_cpl
    type(coop_function)::Vofphi, intQofphi
    COOP_INT::err
    COOP_REAL::norm, fac
    logical, optional::normalize_V
#if DO_COUPLED_DE
    COOP_INT, parameter::ns = 10000
    type(coop_species)::de, cdm
    type(coop_function)::fwp1effcdm, fwp1effde, fwp1de
    COOP_REAL::rho_ce, lna(ns), a(ns), phi(ns), phidot(ns), wp1de(ns), wp1effde(ns), wp1effcdm(ns), V(ns), Q(ns), dQdphi(ns), phi_prime(ns), dVdphibyH2(ns), m2byH2(ns), H, omc, omde, dlna, prevH
    COOP_INT::i
    COOP_REAL, parameter::minphi = 1.d-99
    COOP_REAL::Qrhoma3_at_minphi, Vp_at_minphi, V_at_minphi, Q_at_minphi, dQdphi_at_minphi, m2_at_minphi
    COOP_REAL::norm_up, norm_down, norm_mid

    this%Omega_c_bare = Omega_c
    this%Omega_b_bare = Omega_b
    if(this%baryon_is_coupled)then
       Omega_cpl = Omega_c+Omega_b
    else
       Omega_cpl = Omega_c
       call this%add_species(coop_baryon(Omega_b))
    endif

    err = 0
    if(this%Omega_k() .le. Omega_cpl)then
       err = -1
       return
    endif
    rho_ce = omega_cpl * 3.d0
    Q_at_minphi  =  intQofphi%derivative(minphi)
    dQdphi_at_minphi  =  intQofphi%derivative2(minphi)
    Qrhoma3_at_minphi = Q_at_minphi*rho_ce*exp(intQofphi%eval(minphi))

    call coop_set_uniform(ns, lna, log(coop_min_scale_factor), log(coop_scale_factor_today))
    a = exp(lna)
    dlna = lna(2) - lna(1)
    norm_up = 2.d0
10  norm = norm_up
    call solve(H)
    if(H .lt. 1.d0)then
       norm_up = norm_up*2.d0       
       if(norm_up .gt. 1.d10)then
          err = -2
          return
       endif
       goto 10
    endif
    norm_down = 0.5d0
20  norm = norm_down
    call solve(H)
    if(H .gt. 1.d0)then
       norm_down = norm_down/2.d0
       if(norm_down .lt. 1.d-10)then
          err = -3
          return
       endif
       goto 20
    endif
    do while(norm_up - norm_down .gt. 1.d-7 .and. abs(H-1.d0).gt. 1.d-7)
       norm = sqrt(norm_up * norm_down)
       call solve(H)
!       print*, norm_down, norm_up, norm, H
       if(H .gt. 1.d0)then
          norm_up = norm
       else
          norm_down = norm
       endif
    enddo
    if(present(normalize_V))then
       if(normalize_V) &
            call Vofphi%mult_const(norm)
    endif
    omc = Omega_cpl*exp(intQofphi%eval(phi(ns)))
    if(omc .gt. this%omega_k()) then
       err = -1
       return
    endif
    omde = this%Omega_k() - omc

    if(all(wp1effcdm.eq.0.d0))then
       call cdm%init(genre=COOP_SPECIES_FLUID, name = "CDM", id = 1, Omega = omega_c, w = 0.d0)
    else
       call fwp1effcdm%init(n = ns, xmin=coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., f = wp1effcdm, method = COOP_INTERPOLATE_LINEAR, name = "CDM 1+w", check_boundary = .false.)
       call cdm%init(genre=COOP_SPECIES_FLUID, name = "CDM", id = 1, Omega = omc*omega_c/omega_cpl, w = 0.d0, fwp1eff = fwp1effcdm)
    endif

    call fwp1de%init(n = ns, xmin=coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., f = wp1de, method = COOP_INTERPOLATE_LINEAR, name = "Dark Energy 1+w", check_boundary = .false.)
    call fwp1effde%init(n = ns, xmin=coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., f = wp1effde, method = COOP_INTERPOLATE_LINEAR, name= "Dark Energy effective 1+w", check_boundary = .false.)
    call de%init(name = "Dark Energy", id = 5, Omega = omde, genre = COOP_SPECIES_FLUID, fwp1 = fwp1de, fwp1eff =fwp1effde)
    call de%cplde_Q%init(n = ns, xmin = a(1), xmax = a(ns), xlog = .true., f = Q, method = COOP_INTERPOLATE_LINEAR, name = "Q(phi)")
    call de%cplde_dQdphi_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = dQdphi, method = COOP_INTERPOLATE_LINEAR, name = "dQ / d phi")
    call de%cplde_lnV_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = log(V), method = COOP_INTERPOLATE_LINEAR, name = "ln V")
    call de%cplde_phi_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = phi, method = COOP_INTERPOLATE_LINEAR, name = "phi")
    call de%cplde_phi_prime_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = phi_prime, method = COOP_INTERPOLATE_LINEAR, name="d phi / d lna")
    call de%cplde_dVdphibyH2_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = dVdphibyH2, method = COOP_INTERPOLATE_LINEAR, name = "d V / d phi / H^2")
    call de%cplde_m2byH2_lna%init(n = ns, xmin = lna(1), xmax = lna(ns), f = m2byH2, method = COOP_INTERPOLATE_LINEAR, name = " m^2/H^2 ")
    de%cplde_Vofphi = Vofphi
    if(present(normalize_V))then
       if(.not. normalize_V)  call de%cplde_Vofphi%mult_const(norm)
    else
       call de%cplde_Vofphi%mult_const(norm)
    endif

    de%cplde_intQofphi = intQofphi

    call this%add_species(cdm)    
    if(this%baryon_is_coupled )then
       cdm%Omega = omc*omega_b/omega_cpl
       cdm%name = "BARYON"
       call this%add_species(cdm)
    endif

    call this%add_species(de)
    call de%free()
    call cdm%free()
    call fwp1effcdm%free()
    call fwp1de%free()
    call fwp1effde%free()
    
  contains      


    subroutine solve(H)

      COOP_REAL::H
      COOP_INT::i, istart
      COOP_REAL::rhode, rhom, beta, dbdlna, beta_plus, beta_minus, phinorm      
      logical::good_approx
      type(coop_ode)::ode
      call ode%init(n=2, method = COOP_ODE_RK4)

      V_at_minphi = norm*Vofphi%eval(minphi)
      Vp_at_minphi = norm*Vofphi%derivative(minphi)
      m2_at_minphi = norm*Vofphi%derivative2(minphi)

      i = 1
      !!when phi -> 0 it is numerically difficult to track what exactly phi and m^2 are (usually beyond machine precision)
      do while(Qrhoma3_at_minphi .gt.  - Vp_at_minphi*a(i)**3 )
         phi(i) = minphi
         phidot(i) = 0.d0
         phi_prime(i) = 0.d0
         V(i) = V_at_minphi
         Q(i) = Q_at_minphi
         dQdphi(i) = dQdphi_at_minphi
         rhode = V_at_minphi
         rhom = rho_ce/a(i)**3*exp(intQofphi%eval(phi(i)))
         H = sqrt((this%rhoa4(a(i))/a(i)**4 + rhode + rhom )/3.d0)
         dVdphibyH2(i) = -Qrhoma3_at_minphi/a(i)**3/H**2
         m2byH2(i) = m2_at_minphi/H **2 * (Qrhoma3_at_minphi /(-Vp_at_minphi*a(i)**3))
         wp1de(i) = 0.d0
         wp1effcdm(i) =  1.d0 
         wp1effde(i) = 0.d0
         i = i + 1
         if(a(i) .gt. 1.d-4) stop "the chameleon potential is not steep enough when phi -> 0"
      enddo

      phinorm = 1.d0
      call get_tight_coupling(a(i)*exp(dlna/2.d0), phi(i), phidot(i), good_approx, phinorm, beta_plus, 0.d0,  .true.)
      beta_minus = beta_plus
      call get_tight_coupling(a(i)*exp(-dlna/2.d0), phi(i), phidot(i), good_approx, phinorm, beta_minus, 0.d0, .false.)
      dbdlna = (beta_plus-beta_minus)/dlna
      beta = (beta_plus+beta_minus)/2.d0
      phinorm = phinorm * min(max(exp(-dbdlna*lna(i)**2), 0.1d0), 10.d0)
      do 
         call get_tight_coupling(a(i)*exp(dlna/2.d0), phi(i), phidot(i), good_approx, phinorm, beta_plus, 0.d0, .false.)
         call get_tight_coupling(a(i)*exp(-dlna/2.d0), phi(i), phidot(i), good_approx, phinorm, beta_minus, 0.d0, .false.)
         dbdlna = (beta_plus-beta_minus)/dlna
         call get_tight_coupling(a(i), phi(i), phidot(i), good_approx, phinorm, beta, dbdlna, .false.)
         phinorm = phinorm * min(max(exp(-dbdlna*lna(i)**2), 0.1d0),10.d0)
         V(i) = norm*Vofphi%eval(phi(i))
         Q(i) = intQofphi%derivative(phi(i))
         dQdphi(i) = intQofphi%derivative2(phi(i))
         rhode = phidot(i)**2/2.d0 + V(i)
         rhom = rho_ce/a(i)**3*exp(intQofphi%eval(phi(i)))
         H = sqrt((this%rhoa4(a(i))/a(i)**4 + rhode + rhom )/3.d0)
         phi_prime(i) = phidot(i)/H
         dVdphibyH2(i) = norm*Vofphi%derivative(phi(i))/H**2
         m2byH2(i) = norm*Vofphi%derivative2(phi(i))/H**2
         wp1de(i) = phidot(i)**2/rhode
         wp1effcdm(i) =  1.d0 - Q(i)*phi_prime(i)/3.d0
         wp1effde(i) = wp1de(i) + Q(i)*phi_prime(i)/3.d0 * rhom/rhode
         if(i.ge. ns)then
            istart = ns+1
            exit
         endif
         if( i.gt.1 .and. (a(i) .gt. 0.2d0 .or. .not. good_approx) )then !!start exact evolution
            call ode%set_initial_conditions(xini = lna(i-1), yini = (/ phi(i-1), phidot(i-1) /) )
            istart = i
            exit
         endif
         i = i + 1
      enddo
      do i=istart, ns         
         call ode%evolve(cpl_eq, lna(i))         
         phi(i) = ode%y(1)
         phidot(i) = ode%y(2)
         V(i) = norm * Vofphi%eval(phi(i))
         Q(i) = intQofphi%derivative(phi(i))
         dQdphi(i) = intQofphi%derivative2(phi(i))
         rhode = phidot(i)**2/2.d0 + V(i)
         rhom = rho_ce/a(i)**3*exp(intQofphi%eval(phi(i)))
         H = sqrt((this%rhoa4(a(i))/a(i)**4 + rhode + rhom )/3.d0)
         phi_prime(i) = phidot(i)/H
         dVdphibyH2(i) = norm * Vofphi%derivative(phi(i))/H**2
         m2byH2(i) = norm * Vofphi%derivative2(phi(i))/H**2
         wp1de(i) = phidot(i)**2/rhode
         wp1effcdm(i) =  1.d0 - Q(i)*phi_prime(i)/3.d0
         wp1effde(i) = wp1de(i) + Q(i)*phi_prime(i)/3.d0 * rhom/rhode
      enddo
!!$      do i=istart-10, istart+10
!!$         print*, i-istart, phi(i), phidot(i)
!!$      enddo
      call ode%free()
    end subroutine solve

    !! rho_ce = \lim a-> 0 rho a^3
    !! y(1) = \phi
    !! y(2) = \dot\phi
    !! time variable ln a
    subroutine cpl_eq(n, lna, y, yp)
      COOP_INT::n
      COOP_REAL::lna, y(n), yp(n), H2a4, a,  H
      COOP_REAL::dVdphi, q, rhoma3, qrhom
      COOP_REAL, parameter::fac = 1.d4
      a = exp(lna)
#define PHI y(1)
#define PHIDOT y(2)
      if(PHI .le. 1.d-99)then
         yp(1) = 1.d-10
         yp(2) = 1.d-10
      endif
      q = intQofphi%derivative(PHI)
      dVdphi = norm * Vofphi%derivative(PHI)
      rhoma3 = rho_ce*exp(intQofphi%eval(PHI))
      qrhom = q*rhoma3/a**3
      H2a4 = (this%rhoa4(a) + (norm * Vofphi%eval(PHI)+PHIDOT**2/2.d0)*a**4 + rhoma3*a)/3.d0
      H = sqrt(H2a4)/a**2
      yp(1) = y(2)/H 
      yp(2) = (-3.d0*H*PHIDOT - dVdphi - qrhom)/H
      if(qrhom .gt. 1.d4)then
         yp(1) = sign(min(abs(y(1)/a*fac), abs(yp(1))), yp(1))
         yp(2) = sign(min(abs(y(2)/a*fac), abs(yp(2))), yp(2))
      else
         q = (qrhom/1.d4)**4
         yp(1) = sign(min(abs(y(1)/a*fac), abs(yp(1))), yp(1))*q + yp(1)*(1.d0-q)
         yp(2) = sign(min(abs(y(2)/a*fac), abs(yp(2))), yp(2))*q + yp(2)*(1.d0-q)
      endif
#undef PHI
#undef PHIDOT
    end subroutine cpl_eq



    subroutine get_tight_coupling(a, phi, phidot, good_approx, phinorm, betabest, dbdlna, search_full_range)
      COOP_REAL,parameter::beta_max = 30.d0, beta_min = 0.1d0
      COOP_REAL,intent(INOUT)::betabest
      COOP_REAL::dbdlna
      logical::search_full_range
      COOP_REAL::a
      COOP_REAL::phi, phidot
      COOP_REAL::pa4, rhoa4, p, rho, rhom_raw, H, Hdot, phidd, phinorm
      COOP_REAL::rhom, rho_phi, p_phi, V_phi
      COOP_INT, parameter::n = 101
      COOP_REAL::beta(n), df(n), reldf(n), t1, t2, t3, beta_upper, beta_lower, diffbest, diff_upper, diff_lower
      COOP_INT::i, iloop
      logical::good_approx
      logical::search_all
      rhom_raw = rho_ce/a**3
      call this%get_pa4_rhoa4(a, pa4, rhoa4)
      p = pa4/a**4
      rho  = rhoa4/a**4
      search_all = search_full_range
10    continue
      if(search_all)then
         beta_lower = beta_min
         beta_upper = beta_max
         diff_lower = 0.d0
         diff_upper = 0.d0
         iloop = 0
         do while( (iloop .lt. 2 .or. diff_lower*diff_upper .gt. 0.d0) .and. iloop .lt. 5)
            call coop_set_uniform(n, beta,  beta_lower, beta_upper, logscale = iloop .eq. 0)
            do i = 1, n
               phi = phinorm*a**beta(i)
               rhom = rhom_raw*exp(intQofphi%eval(phi))
               V_phi =  norm*Vofphi%eval(phi)
               H = sqrt((rho + rhom + V_phi)/3.d0)
               phidot = beta(i)*H*phi
               H = sqrt((rho + rhom + V_phi+phidot**2/2.d0)/3.d0)
               phidot = (beta(i)+ dbdlna*log(a))*H*phi
               Hdot = -0.5d0*(rho + rhom + phidot**2+p)
               phidd = (beta(i) + dbdlna*log(a))*(Hdot*phi+H*phidot) &
                    + 2.d0*dbdlna*H**2*phi
               t1 = phidd + 3.d0*H*phidot
               t2 = norm*Vofphi%derivative(phi)
               t3 =  intQofphi%derivative(phi)*rhom 
               df(i) =  (t1 + t2 + t3)
               reldf(i) =abs(df(i))/max(abs(t1), abs(t2), abs(t3), 1.d0)
            enddo
            i = coop_minloc(reldf)
            betabest = beta(i)
            beta_lower = beta(max(i-1, 1))
            beta_upper = beta(min(i+1, n)) 
            diff_lower = df(max(i-1, 1))
            diff_upper = df(min(i+1, n))
            iloop = iloop+1
         enddo
      else
         beta_upper = min(betabest + max(1.d-3,abs(dbdlna) * dlna), beta_max)
         beta_lower = max(betabest - max(abs(dbdlna) * dlna, 1.d-3), beta_min)
         phi = phinorm*a**beta_upper
         rhom = rhom_raw*exp(intQofphi%eval(phi))
         V_phi =  norm*Vofphi%eval(phi)
         H = sqrt((rho + rhom + V_phi)/3.d0)
         phidot = (beta_upper+ dbdlna*log(a))*H*phi
         H = sqrt((rho + rhom + V_phi+phidot**2/2.d0)/3.d0)
         Hdot = -0.5d0*(rho + rhom + phidot**2+p)
         phidot = (beta_upper+ dbdlna*log(a))*H*phi
         phidd = (beta_upper + dbdlna*log(a))*(Hdot*phi+H*phidot) &
              + 2.d0*dbdlna*H**2*phi
         t1 = phidd + 3.d0*H*phidot
         t2 = norm*Vofphi%derivative(phi)
         t3 =  intQofphi%derivative(phi)*rhom 
         diff_upper =  t1 + t2 + t3

         phi = phinorm*a**beta_lower
         rhom = rhom_raw*exp(intQofphi%eval(phi))
         V_phi =  norm*Vofphi%eval(phi)
         H = sqrt((rho + rhom + V_phi)/3.d0)
         phidot = (beta_lower+ dbdlna*log(a))*H*phi
         H = sqrt((rho + rhom + V_phi+phidot**2/2.d0)/3.d0)
         Hdot = -0.5d0*(rho + rhom + phidot**2+p)
         phidot = (beta_lower+ dbdlna*log(a))*H*phi
         phidd = (beta_lower + dbdlna*log(a))*(Hdot*phi+H*phidot) &
              + 2.d0*dbdlna*H**2*phi
         t1 = phidd + 3.d0*H*phidot
         t2 = norm*Vofphi%derivative(phi)
         t3 =  intQofphi%derivative(phi)*rhom 
         diff_lower = t1 + t2 + t3
         if(diff_upper*diff_lower .gt. 0.d0)then
            search_all = .true.
            goto 10
         endif

      endif
      do while(beta_upper - beta_lower .gt. 1.d-10 )
         betabest = (beta_upper + beta_lower)/2.d0
         phi = phinorm*a**betabest
         rhom = rhom_raw*exp(intQofphi%eval(phi))
         V_phi =  norm*Vofphi%eval(phi)
         H = sqrt((rho + rhom + V_phi)/3.d0)
         phidot = (betabest+ dbdlna*log(a))*H*phi
         H = sqrt((rho + rhom + V_phi+phidot**2/2.d0)/3.d0)
         Hdot = -0.5d0*(rho + rhom + phidot**2+p)
         phidot = (betabest+ dbdlna*log(a))*H*phi
         phidd = (betabest + dbdlna*log(a))*(Hdot*phi+H*phidot) &
              + 2.d0*dbdlna*H**2*phi
         t1 = phidd + 3.d0*H*phidot
         t2 = norm*Vofphi%derivative(phi)
         t3 =  intQofphi%derivative(phi)*rhom 
         diffbest = ( t1 + t2 + t3)
         if(diffbest*diff_upper .le. 0.d0)then
            beta_lower = betabest
         else
            beta_upper = betabest
         endif
      enddo
      betabest = (beta_upper + beta_lower)/2.d0
      phi = phinorm*a**betabest
      rhom = rhom_raw*exp(intQofphi%eval(phi))
      V_phi =  norm*Vofphi%eval(phi)
      H = sqrt((rho + rhom + V_phi)/3.d0)
      phidot = (betabest+ dbdlna*log(a))*H*phi
      good_approx = (abs(diffbest) .lt.  max(abs(t1), abs(t2), abs(t3))* 1.d-5 .and. phidot .ge. 0.d0 .and. phi.lt.1.d-3)
      if(phi.eq. 0.d0 .or. phidot .eq. 0.d0)then
         print*, a, betabest, phinorm, phi, phidot
         stop
      endif
    end subroutine get_tight_coupling

#else
    write(*,*) "Coupled Dark Energy model cannot be initialized"
    stop "You need to set DARK_ENERGY_MODEL=COUPLED_DE in configure.in"
#endif    
  end subroutine coop_background_add_coupled_DE_with_potential


  
end module coop_background_mod
