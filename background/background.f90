module coop_background_mod
  use coop_wrapper_typedef
  implicit none
#include "constants.h"

  private
  COOP_INT:: coop_coupled_de_num_iterate = 8

  public::coop_baryon, coop_cdm, coop_DE_lambda, coop_DE_w0, coop_DE_w0wa, coop_DE_quintessence, coop_de_coupled_quintessence, coop_radiation, coop_neutrinos_massless, coop_neutrinos_massive, coop_de_w_coupled_quintessence, coop_de_w_quintessence, coop_de_iterate_coupling_equations, coop_coupled_de_num_iterate

contains


  function coop_baryon(Omega_b, fcs2b) result(this)
    type(coop_species) this
    type(coop_function),optional::fcs2b
    COOP_REAL Omega_b
    if(present(fcs2b))then
       call this%init(genre = COOP_SPECIES_FLUID, name = "Baryon", id = 1, Omega=Omega_b, w = COOP_REAL_OF(0.), fcs2 = fcs2b)  !!you might want to replace the cs^2(baryon)
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
    call this%init(genre = COOP_SPECIES_LAMBDA, name = "Cosmological Constant",id=5, Omega = Omega_Lambda, w = COOP_REAL_OF(-1.d0), cs2 = COOP_REAL_OF(1.d0))
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
    fw0wa = coop_function(coop_de_wp1_w0wa, xmin = coop_min_scale_factor, xmax = COOP_REAL_OF(1.), xlog = .true., args = w0wa)
    call this%init(genre = COOP_SPECIES_FLUID, name = "w0wa Dark Energy", id=5, Omega = Omega_Lambda, cs2 = COOP_REAL_OF(1.d0), fwp1 = fw0wa )
    call w0wa%free
  end function coop_de_w0wa

  function coop_de_w0(Omega_Lambda, w0) result(this)
    type(coop_species) this
    COOP_REAL w0, Omega_Lambda
    if(w0.eq.-1.d0)then
       this = coop_de_lambda(Omega_Lambda)
       return
    endif
    call this%init(genre = COOP_SPECIES_FLUID, name = "constant w Dark Energy", id = 5, Omega = Omega_Lambda, w = w0, cs2 = COOP_REAL_OF(1.))
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
    fq = coop_function(coop_de_wp1_quintessence, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args = arg)
    call this%init(genre = COOP_SPECIES_FLUID, name = "quintessence Dark Energy", id=5, Omega = Omega_Lambda, cs2 = COOP_REAL_OF(1.d0), fwp1 = fq )
    call arg%free
  end function coop_de_quintessence

  function coop_de_w_quintessence(a, arg) result(w)
    COOP_REAL a, w
    type(coop_arguments) arg
    w = coop_de_wp1_quintessence(a, arg) - 1.d0
  end function coop_de_w_quintessence

  function coop_de_wp1_quintessence(a, arg) result(wp1)
    COOP_REAL a, wp1
    type(coop_arguments) arg
    COOP_REAL Omega_m, a_eq
    COOP_REAL mu, mu3, s0, s1, qpsign, aux1, aux2, aux3, delta
#define OMEGA_LAMBDA arg%r(1)
#define EPSILON_S arg%r(2)
#define EPSILON_INFTY arg%r(3)
#define ZETA_S arg%r(4)
    if(EPSILON_S .ge. 0.d0)then
       qpsign = 1.d0
    else
       qpsign = -1.d0
    endif
    Omega_m  = 1.d0 - OMEGA_LAMBDA
    if(Omega_m .gt. 0.65d0 .or. Omega_m .lt. 0.05d0)then
       wp1 = qpsign * 3.d0 !!set w to be a crazy value to rule out the model
       return
    endif
    aux1 = sqrt(EPSILON_INFTY/3.d0)
    aux2 = sqrt(2.d0)*(sqrt(abs(EPSILON_S)/6.d0) - aux1)*(1.-ZETA_S)
    aux3 = 2.*ZETA_S*(sqrt(abs(EPSILON_S)/6.d0) - aux1 )
    delta=(aux1 + 0.533* aux2 + 0.307* aux3)**2 &
         + (aux1 + (0.91-0.78*Omega_m)*aux2+(0.81-1.09*Omega_m)*aux3)**2
    if(delta .lt. 0.9d0)then
       a_eq = (Omega_m/OMEGA_LAMBDA)**((1.d0/3.d0)/(1.d0 - sign(delta,EPSILON_S)))
    else
       wp1 = qpsign * 3.d0  !!set w to be a crazy value to rule out the model
       return
    endif
    mu=a/a_eq
    mu3=mu**3
    s0=sqrt(mu3)
    s1=sqrt(1.+mu3)
    if(s0.lt.2.d-4)then  
       wp1 = qpsign * EPSILON_INFTY/1.5 
    elseif(s0.lt. 1.d-2)then
       wp1 = 2.d0 * qpsign * ( aux1* sqrt(1.+ a_eq/3./(a_eq +  a)) + aux2 * (2.d0/3.d0) * s0 )**2
    else
       wp1 = 2.d0 * qpsign * (aux1*(1. + a_eq/6./(a + a_eq)) + aux2*(s1/s0-log(s0+s1)/mu3) + aux3*(1.-log(1.+mu3)/mu3))**2
    endif
#undef EPSILON_S
#undef EPSILON_INF
#undef ZETA_S
  end function coop_de_wp1_quintessence



  function coop_de_coupled_quintessence(Omega_Lambda, epsilon_s, epsilon_inf, zeta_s, at_by_aeq) result(this)
    type(coop_species) this
    COOP_REAL Omega_Lambda, epsilon_s, epsilon_inf, zeta_s, at_by_aeq
    type(coop_function) fq
    type(coop_arguments) arg
    if(epsilon_s .eq. 0.d0 .and. epsilon_inf .eq. 0.d0 .and. zeta_s .eq. 0.d0)then
       this = coop_de_lambda(Omega_Lambda)
       return
    endif
    arg = coop_arguments(r =  (/ Omega_Lambda, epsilon_s, epsilon_inf, zeta_s , at_by_aeq /))
    fq = coop_function(coop_de_wp1_coupled_quintessence, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args = arg)
    call this%init(genre = COOP_SPECIES_FLUID, name = "coupled quintessence Dark Energy", id=5, Omega = Omega_Lambda, cs2 = COOP_REAL_OF(1.d0), fwp1 = fq )
    call arg%free
  end function coop_de_coupled_quintessence

  function coop_de_w_coupled_quintessence(a, arg) result(w)
    COOP_REAL a, w
    type(coop_arguments) arg
    w = coop_de_wp1_coupled_quintessence(a, arg) - 1.d0
  end function coop_de_w_coupled_quintessence

  function coop_de_wp1_coupled_quintessence(a, arg) result(wp1)
    COOP_REAL a, wp1
    type(coop_arguments) arg
    COOP_REAL Omega_m, a_eq
    COOP_REAL mu, mu3, s0, s1, qpsign, aux1, aux2, aux3, delta
#define OMEGA_LAMBDA arg%r(1)
#define EPSILON_S arg%r(2)
#define EPSILON_INFTY arg%r(3)
#define ZETA_S arg%r(4)
#define AT_BY_AEQ arg%r(5)
    if(EPSILON_S .ge. 0.d0)then
       qpsign = 1.d0
    else
       qpsign = -1.d0
    endif
    Omega_m  = 1.d0 - OMEGA_LAMBDA
    if(Omega_m .gt. 0.65d0 .or. Omega_m .lt. 0.05d0)then
       wp1 = qpsign * 3.d0 !!set w to be a crazy value to rule out the model
       return
    endif
    aux1 = sqrt(EPSILON_INFTY/3.d0)
    aux2 = sqrt(2.d0)*(sqrt(abs(EPSILON_S)/6.d0) - aux1)*(1.-ZETA_S)
    aux3 = 2.*ZETA_S*(sqrt(abs(EPSILON_S)/6.d0) - aux1 )
    delta=(aux1 + 0.533* aux2 + 0.307* aux3)**2 &
         + (aux1 + (0.91-0.78*Omega_m)*aux2+(0.81-1.09*Omega_m)*aux3)**2
    if(delta .lt. 0.9d0)then
       a_eq = (Omega_m/OMEGA_LAMBDA)**((1.d0/3.d0)/(1.d0 - sign(delta,EPSILON_S)))
    else
       wp1 = qpsign * 3.d0 !!set w to be a crazy value to rule out the model
       return
    endif
    mu=a/a_eq
    mu3=mu**3
    s0=sqrt(mu3)
    s1=sqrt(1.+mu3)
    if(s0.lt.2.d-4)then  
       wp1 = qpsign * EPSILON_INFTY/1.5 
    elseif(s0.lt. 1.d-2)then
       wp1 = 2.d0 * qpsign * ( aux1* sqrt(1.+ a_eq/3./(a_eq +  a)) + aux2 * (2.d0/3.d0) * s0 )**2
    else
       wp1 = 2.d0 * qpsign * (aux1*(1. + a_eq/6./(a + a_eq)) + aux2*(s1/s0-log(s0+s1)/mu3) + aux3*(1.-log(1.+mu3)/mu3))**2
    endif
    if(AT_BY_AEQ .gt. 0.d0) wp1 = wp1 * tanh(a/(AT_BY_AEQ * a_eq))**3
#undef EPSILON_S
#undef EPSILON_INF
#undef ZETA_S
#undef AT_BY_AEQ
  end function coop_de_wp1_coupled_quintessence

  subroutine coop_de_iterate_coupling_equations(Q, omega_radiation, de, baryon, cdm, mnu)
    COOP_REAL,intent(IN)::Q, omega_radiation
    type(coop_species)::de, baryon, cdm
    type(coop_species),optional::mnu
    type(coop_species) mod_de, mod_baryon, mod_cdm, mod_mnu
    type(coop_function)::fw
    integer,parameter::n=4096
    COOP_REAL::a(n), wp1de(n), wp1m(n), wp1nu(n), wp1de_origin(n), wnu_origin(n), pnu, lna(n), coupl, rhode, rhob, rhoc, rhonu,  H2
    integer i, it
    !!first iteration, modify rho_de
    call coop_set_uniform(n, lna, log(coop_min_scale_factor), log(coop_scale_factor_today))
    a = exp(lna)
    !$omp parallel do 
    do i=1, n
       wp1de_origin(i) = abs(de%wp1ofa(a(i)))
       if(present(mnu))then
          wnu_origin(i) = mnu%wofa(a(i))
       endif
    enddo
    !$omp end parallel do
    do it = 1, coop_coupled_de_num_iterate
       if(present(mnu))then
          !$omp parallel do private(coupl, rhode, rhob, rhoc, rhonu, pnu, H2)
          do i=1, n
             rhode =3.d0* de%Omega*de%density_ratio(a(i))
             rhob = 3.d0*baryon%Omega*baryon%density_ratio(a(i))
             rhoc = 3.d0* cdm%Omega * cdm%density_ratio(a(i)) 
             rhonu = 3.d0*mnu%Omega * mnu%density_ratio(a(i))
             pnu = rhonu * wnu_origin(i)
             H2 = (rhode +  rhob + rhoc + rhonu + 3.d0*omega_radiation/a(i)**4)/3.d0
             coupl = min(Q*sqrt(wp1de_origin(i)*rhode/H2)/3.d0, 0.99d0) !!avoid NAN error for huge Q models
             wp1de(i) = wp1de_origin(i) + coupl*(rhob+rhoc+rhonu-3.d0*pnu)/rhode
             wp1m(i) = 1.d0 - coupl
             wp1nu(i) = 1.d0 + wnu_origin(i) - coupl*(1.d0-3.d0*pnu/rhonu)
          enddo
          !$omp end parallel do
       else
          !$omp parallel do private(coupl, rhode, rhob, rhoc,  H2)
          do i = 1, n
             rhode = 3.d0* de%Omega*de%density_ratio(a(i))
             rhob = 3.d0*baryon%Omega*baryon%density_ratio(a(i))
             rhoc = 3.d0* cdm%Omega * cdm%density_ratio(a(i)) 
             H2 = (rhode +  rhob + rhoc  + 3.d0*omega_radiation/a(i)**4)/3.d0
             coupl = min(Q*sqrt(wp1de_origin(i)*rhode/H2)/3.d0, 0.99d0) !!avoid NAN error for huge Q models
             wp1de(i) = wp1de_origin(i) + coupl*(rhob+rhoc)/rhode
             wp1m(i) = 1.d0 - coupl
          enddo
          !$omp end parallel do
       endif
       call fw%init(n, coop_min_scale_factor, coop_scale_factor_today, wp1de, xlog = .true., ylog = .false., check_boundary = .false.)
       call mod_de%init(genre = COOP_SPECIES_FLUID, name= "coupled quintessence Dark Energy", id=5, Omega = de%omega, cs2 = de%cs2, fwp1 = de%fwp1, fwp1eff = fw )
       call fw%init(n, coop_min_scale_factor, coop_scale_factor_today, wp1m, xlog = .true., ylog = .false., check_boundary = .false.)    
       call mod_baryon%init(genre = COOP_SPECIES_FLUID, name = "coupled baryon", id=5, Omega = baryon%omega, cs2 = baryon%cs2, w = baryon%wp1 - COOP_REAL_OF(1.d0), fwp1eff = fw )
       call mod_cdm%init(genre = COOP_SPECIES_FLUID, name = "coupled cdm", id=5, Omega = cdm%omega, cs2 = cdm%cs2, w = cdm%wp1 - COOP_REAL_OF(1.d0), fwp1 = fw )
       if(present(mnu))then
          call fw%init(n, coop_min_scale_factor,coop_scale_factor_today, wp1nu, xlog = .true., ylog = .false., check_boundary = .false.)       
          call mod_mnu%init(genre = COOP_SPECIES_FLUID, name = "coupled massive neutrinos", id=5, Omega = mnu%omega, fwp1 = mnu%fwp1 , fcs2 = mnu%fcs2, fwp1eff = fw )
          mnu = mod_mnu
          call mod_mnu%free()
       endif
       de = mod_de
       call mod_de%free()
       baryon = mod_baryon
       call mod_baryon%free()
       cdm = mod_cdm
       call mod_cdm%free()
    enddo
  end subroutine coop_de_iterate_coupling_equations


  function coop_zrecomb_fitting(ombh2, omch2) result(zstar)
    COOP_REAL zstar, ombh2, omch2
  !!From Hu & Sugiyama
    zstar =  1048 * (1 + 0.00124 * ombh2**(-0.738))*(1+ &
         (0.0783 * ombh2 **(-0.238)/(1+39.5* ombh2 ** 0.763)) * &
         (omch2 + ombh2)**(0.560/(1+21.1* ombh2 **1.81)))
  end function coop_zrecomb_fitting


end module coop_background_mod
