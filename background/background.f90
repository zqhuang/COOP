module coop_background_mod
  use coop_wrapper_utils
  implicit none
#include "constants.h"

  private
  COOP_INT:: coop_coupled_de_num_iterate = 8

  public::coop_baryon, coop_cdm, coop_DE_lambda, coop_DE_w0, coop_DE_w0wa, coop_DE_quintessence, coop_de_coupled_quintessence, coop_radiation, coop_neutrinos_massless, coop_neutrinos_massive, coop_de_w_coupled_quintessence, coop_de_w_quintessence, coop_de_iterate_coupling_equations, coop_coupled_de_num_iterate, coop_de_wp1_coupled_quintessence, coop_de_wp1_quintessence, coop_background_add_coupled_DE

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
    fw0wa = coop_function(coop_de_wp1_w0wa, xmin = coop_min_scale_factor, xmax = COOP_REAL_OF(1.), xlog = .true., args = w0wa)
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
    fq = coop_function(coop_de_wp1_quintessence, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args = arg)
    call this%init(genre = COOP_SPECIES_FLUID, name = "Dark Energy", id=5, Omega = Omega_Lambda, cs2 = COOP_REAL_OF(1.d0), fwp1 = fq )
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
       wp1 = 1.d0 !!set w to be a crazy value to rule out the model
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
    if(s0 .lt. 1.d-6)then  
       wp1 = qpsign * EPSILON_INFTY/1.5 
    elseif(s0.lt. 5.d-3)then
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
    call  arg%init(r =  (/ Omega_Lambda, epsilon_s, epsilon_inf, zeta_s , at_by_aeq /))
    fq = coop_function(coop_de_wp1_coupled_quintessence, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, xlog = .true., args = arg, name = "quintessence 1+w(a)")
    call this%init(genre = COOP_SPECIES_FLUID, name = "Dark Energy", id=5, Omega = Omega_Lambda, cs2 = COOP_REAL_OF(1.d0), fwp1 = fq )
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
       wp1 = 1.d0 !!set w to be a crazy value to rule out the model
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
    if(s0.lt. 1.d-6)then  
       wp1 = qpsign * EPSILON_INFTY/1.5 
    elseif(s0.lt. 5.d-3)then
       wp1 = 2.d0 * qpsign * ( aux1* sqrt(1.+ a_eq/3./(a_eq +  a)) + aux2 * (2.d0/3.d0) * s0 )**2
    else
       wp1 = 2.d0 * qpsign * (aux1*(1. + a_eq/6./(a + a_eq)) + aux2*(s1/s0-log(s0+s1)/mu3) + aux3*(1.-log(1.+mu3)/mu3))**2
    endif
    if(AT_BY_AEQ .gt. 0.d0) wp1 = wp1 * tanh(a/(AT_BY_AEQ * a_eq))**(3.d0*(1.d0-a))
#undef EPSILON_S
#undef EPSILON_INF
#undef ZETA_S
#undef AT_BY_AEQ
  end function coop_de_wp1_coupled_quintessence

  subroutine coop_de_iterate_coupling_equations(Q, dlnQdphi, omega_radiation, de, baryon, cdm, mnu)
    COOP_REAL,intent(IN)::Q, omega_radiation, dlnQdphi
    type(coop_species)::de, baryon, cdm
    type(coop_species),optional::mnu
    type(coop_species) mod_de, mod_baryon, mod_cdm, mod_mnu
    type(coop_function)::fw
    integer,parameter::n=4096
    COOP_REAL::a(n), wp1de(n), wp1m(n), wp1nu(n), wp1de_origin(n), wnu_origin(n), pnu, lna(n), coupl, Qphi(n), dphi(n), rhode, rhob, rhoc, rhonu,  H2, dlna
    integer i, it, ieq
    mod_de = de
    mod_baryon = baryon
    mod_cdm = cdm
    if(present(mnu))mod_mnu = mnu
    !!first iteration, modify rho_de
    call coop_set_uniform(n, lna, log(coop_min_scale_factor), log(coop_scale_factor_today))
    a = exp(lna)
    dlna = lna(2) - lna(1)
    ieq = n
    do i = n-1, 1, -1
       if(de%density_ratio(a(i))*de%Omega .lt.  (cdm%Omega + baryon%Omega)/a(i)**3)then
          ieq = i
          exit
       endif
    enddo

    !$omp parallel do 
    do i=1, n
       wp1de_origin(i) = abs(de%wp1ofa(a(i)))
       if(present(mnu))then
          wnu_origin(i) = mnu%wofa(a(i))
       endif
    enddo
    Qphi = Q
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
             dphi(i) = sqrt(wp1de_origin(i)*rhode/H2)
             coupl = min(Qphi(i)*dphi(i)/3.d0, 0.99d0) !!avoid NAN error for huge Q models
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
             dphi(i) = sqrt(wp1de_origin(i)*rhode/H2)
             coupl = min(Qphi(i)*dphi(i)/3.d0, 0.99d0) !!avoid NAN error for huge Q models
             wp1de(i) = wp1de_origin(i) + coupl*(rhob+rhoc)/rhode
             wp1m(i) = 1.d0 - coupl
          enddo
          !$omp end parallel do
       endif
       Qphi(ieq) = 0
       do i=ieq+1, n
          Qphi(i) = Qphi(i-1) + dlnQdphi*(dphi(i)+dphi(i-1))/2.d0
       enddo
       do i=ieq-1, 1, -1
          Qphi(i) = Qphi(i+1) - dlnQdphi*(dphi(i)+dphi(i+1))/2.d0
       enddo
       Qphi = Q * exp(Qphi*dlna)
       call fw%init(n, coop_min_scale_factor, coop_scale_factor_today, wp1de, xlog = .true., ylog = .false., check_boundary = .false., name = "quintessence 1+w_eff(a)")
       call mod_de%init(genre = COOP_SPECIES_FLUID, name= "Dark Energy", id=5, Omega = de%omega, cs2 = de%cs2, fwp1 = de%fwp1, fwp1eff = fw )
       call fw%init(n, coop_min_scale_factor, coop_scale_factor_today, wp1m, xlog = .true., ylog = .false., check_boundary = .false., name = "cdm/baryon 1+w_eff(a)")    
       call mod_baryon%init(genre = COOP_SPECIES_FLUID, name = "Baryon", id=5, Omega = baryon%omega, cs2 = baryon%cs2, w = baryon%wp1 - COOP_REAL_OF(1.d0), fwp1eff = fw )
       call mod_cdm%init(genre = COOP_SPECIES_FLUID, name = "CDM", id=5, Omega = cdm%omega, cs2 = cdm%cs2, w = cdm%wp1 - COOP_REAL_OF(1.d0), fwp1 = fw )
       if(present(mnu))then
          call fw%init(n, coop_min_scale_factor,coop_scale_factor_today, wp1nu, xlog = .true., ylog = .false., check_boundary = .false., name = "massive neutrinos 1+w(a)")       
          call mod_mnu%init(genre = COOP_SPECIES_FLUID, name = "Massive Neutrinos", id=5, Omega = mnu%omega, Omega_massless = mnu%Omega_massless, fwp1 = mnu%fwp1 , fcs2 = mnu%fcs2, fwp1eff = fw )
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


  subroutine coop_background_add_coupled_DE(this, Omega_c, Q, tracking_n)
    !!this Omega_c is the effective Omega_cdm (rho_m a^3 (a<<1) -> Omega_c rho_0) 
    !!Q is the coupling between DE and CDM
    !!the potential V(phi) = V0 / phi^n
    !!tracking_n is the index n
    !!V0 will be automatically determined by the constraint Omega_k = 0
    class(coop_cosmology_background)::this
    COOP_REAL::omega_c, Q, tracking_n
    type(coop_species)::cdm, de
    type(coop_ode)::ode, ode_tc
    COOP_INT i, istart,  j
    COOP_REAL  beta, lnrhom0
    COOP_REAL::lnc_lower, lnc_upper, lnc_mid
    COOP_REAL::Ot_lower, Ot_upper, Ot_mid, lna_ini, Vnow, KEnow
    COOP_REAL_ARRAY::lna, lnrho_cdm, lnrho_de, wp1_de, phi_de, phidot_de
    COOP_REAL,parameter::a_ini = 1.d-5
    COOP_INT,parameter::nsteps = 200
    COOP_REAL::lna_coarse(nsteps)
    logical::tight_coupling
    cdm%name = "CDM"
    de%name = "Dark Energy"
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
    
    if(Q .gt. 1.d-10)then
       call de%fDE_Q_of_phi%init_polynomial( (/ Q /) )
    endif
    de%DE_tracking_n = max(tracking_n, 0.d0)
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
    call de%fDE_U_of_phi%init_polynomial( (/ 0.d0 /), name = "U(phi)" )  !!U = 1
    call de%fDE_dUdphi%init_polynomial( (/ 0.d0 /), name = "dU/dphi" ) !! dU /d phi = 0
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
             Vnow =  de%DE_V(ode%PHI)
             KEnow  = ode%PHIDOT**2/2.d0
             lnrho_de(i) = log(Vnow + KEnow)
             wp1_de(i) = 2.d0*KEnow/(Vnow + KEnow)
             phi_de(i) = ode%PHI
             phidot_de(i) = ode%PHIDOT             
          enddo
500       continue
          do i = istart + 1, coop_default_array_size
             call evolve_to(lna(i))
             lnrho_cdm(i) = ode%LN_RHOM
             Vnow =  de%DE_V(ode%PHI)
             KEnow  = ode%PHIDOT**2/2.d0
             lnrho_de(i) = log(Vnow + KEnow)
             wp1_de(i) = 2.d0*KEnow/(Vnow + KEnow)
             phi_de(i) = ode%PHI
             phidot_de(i) = ode%PHIDOT
          enddo
          de%Omega = (de%DE_V(ode%PHI) + ode%PHIDOT**2/2.d0)/3.d0
          cdm%Omega = this%Omega_k() - de%Omega
          lnrho_cdm = lnrho_cdm - lnrho_cdm(coop_default_array_size)
          lnrho_de = lnrho_de - lnrho_de(coop_default_array_size)
          call cdm%flnrho%init(coop_default_array_size, coop_min_scale_factor, xmax = coop_scale_factor_today, f = lnrho_cdm, fleft = lnrho_cdm(1), fright = 0.d0, slopeleft = -3.d0, sloperight = -3.d0, xlog = .true.,  method = COOP_INTERPOLATE_LINEAR, check_boundary = .false., name = "CDM \ln\rho(a) ratio")
          call de%flnrho%init(coop_default_array_size, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, f = lnrho_de, fleft = lnrho_de(1), fright = 0.d0, slopeleft = 0.d0, sloperight = 0.d0, xlog = .true., method = COOP_INTERPOLATE_LINEAR, check_boundary = .false., name = "DE \ln\rho(a) ratio")
          call de%fwp1%init(coop_default_array_size, xmin = coop_min_scale_factor, xmax = coop_scale_factor_today, f = wp1_de, fleft = wp1_de(1), fright = wp1_de(coop_default_array_size), slopeleft = 0.d0, sloperight = 0.d0, xlog = .true., method = COOP_INTERPOLATE_LINEAR, check_boundary = .false., name = "DE 1+w(a)")
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
      t = 1.d0/hubble/2.d0
      PHI = exp(lnc)*t**beta
      PHIDOT = beta*PHI/t
      !!set the normalization of V
      qnow = de%DE_Q(PHI)      
      ddphi = (beta-1.d0)*PHIDOT/t
      norm = ((ddphi + 3.d0*hubble*PHIDOT) + qnow * exp(LN_RHOM)) / de%DE_tracking_n * PHI**(de%DE_tracking_n + 1.d0)
      de%DE_lnV0 = log(norm)     
      call ode%set_initial_conditions(lna, y)
      tight_coupling = ( (hubble*PHIDOT) .lt. 1.d-3*qnow*exp(LN_RHOM) )
      if(tight_coupling)then
         call ode_tc%set_initial_conditions(lna, (/ y(3) /) )
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
      COOP_REAL::H, lna, Q, a, pa4, rhoa4, a4, qnow, rhoma4, Vp, V, KE
      a = exp(lna)
      a4 = a**4
      qnow = de%DE_Q(PHI)
      rhoma4 = exp(LN_RHOM + 4.d0*lna)
      call this%get_pa4_rhoa4(a, pa4, rhoa4)
      V = de%DE_V(PHI)
      KE =  PHIDOT**2/2.d0
      rhoa4 = rhoa4 + rhoma4 + (V + KE )*a4
      pa4 = pa4 + (KE - V)*a4
      H = sqrt(rhoa4/3.d0/a4)
      Vp = de%DE_dlnVdphi(PHI)*V
      PHI_PRIME = PHIDOT / H
      PHIDOT_PRIME = - 3.d0*PHIDOT - Vp/H  - (rhoma4/a4/H) * qnow
      LN_RHOM_PRIME = - 3.d0 + qnow * PHI_PRIME
    end subroutine fcn


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
  

end module coop_background_mod
