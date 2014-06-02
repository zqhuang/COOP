module coop_background
  use coop_wrapper_typedef
  implicit none
#include "constants.h"

  private

  public::coop_baryon, coop_cdm, coop_DE_lambda, coop_DE_w0, coop_DE_w0wa, coop_DE_quintessence, coop_radiation, coop_neutrinos_massless, coop_neutrinos_massive

contains


  function coop_baryon(Omega_b, fcs2b) result(this)
    type(coop_species) this
    type(coop_function),optional::fcs2b
    COOP_REAL Omega_b
    if(present(fcs2b))then
       call this%init(gengre = COOP_SPECIES_FLUID, name = "Baryon", id = 1, Omega=Omega_b, w = COOP_REAL_OF(0.), fcs2 = fcs2b)  !!you might want to replace the cs^2(baryon)
    else
       call this%init(gengre = COOP_SPECIES_FLUID, name = "Baryon", id = 1, Omega=Omega_b, w = COOP_REAL_OF(0.), cs2 = COOP_REAL_OF(0.))  !!you might want to replace the cs^2(baryon)
    endif
  end function coop_baryon

  function coop_cdm(Omega_c) result(this)
    type(coop_species) this
    COOP_REAL Omega_c
    call this%init(gengre = COOP_SPECIES_CDM, name = "CDM", id = 2, Omega=Omega_c, w = COOP_REAL_OF(0.), cs2 = COOP_REAL_OF(0.))
  end function coop_cdm

  function coop_radiation(Omega_r) result(this)
    type(coop_species) this
    COOP_REAL:: Omega_r
    call this%init(gengre = COOP_SPECIES_MASSLESS, name = "Radiation", id = 3, Omega=Omega_r) 
  end function coop_radiation

  function coop_neutrinos_massless(Omega_nu) result(this)
    type(coop_species) this
    COOP_REAL:: Omega_nu
    call this%init(gengre = COOP_SPECIES_MASSLESS, name = "Massless Neutrinos", id = 4, Omega=Omega_nu)
  end function coop_neutrinos_massless

  function coop_neutrinos_massive(Omega_nu, Omega_massless) result(this)
    type(coop_species) this
    COOP_REAL:: Omega_nu
    COOP_REAL:: Omega_massless
    call this%init(gengre = COOP_SPECIES_MASSIVE_FERMION, name = "Massive Neutrinos", id = 4, Omega=Omega_nu, Omega_massless = Omega_massless)
  end function coop_neutrinos_massive
  
  function coop_de_lambda(Omega_Lambda) result(this)
    type(coop_species) this
    COOP_REAL Omega_Lambda
    call this%init(gengre = COOP_SPECIES_LAMBDA, name = "Cosmological Constant",id=5, Omega = Omega_Lambda, w = COOP_REAL_OF(-1.), cs2 = COOP_REAL_OF(1.))
  end function coop_de_lambda
  
  function coop_de_w0wa(Omega_Lambda, w0, wa) result(this)
    type(coop_species) this
    COOP_REAL Omega_Lambda, w0, wa
    type(coop_function) fw0wa
    type(coop_arguments) w0wa
    w0wa = coop_arguments(r =  (/ w0, wa /))
    fw0wa = coop_function(coop_de_w_w0wa, xmin = coop_min_scale_factor, xmax = COOP_REAL_OF(1.), xlog = .true., arg = w0wa)
    call this%init(gengre = COOP_SPECIES_FLUID, name = "w0wa Dark Energy", id=5, Omega = Omega_Lambda, cs2 = COOP_REAL_OF(1.d0), fw = fw0wa )
    call w0wa%free
  end function coop_de_w0wa

  function coop_de_w0(Omega_Lambda, w0) result(this)
    type(coop_species) this
    COOP_REAL w0, Omega_Lambda
    call this%init(gengre = COOP_SPECIES_FLUID, name = "constant w Dark Energy", id = 5, Omega = Omega_Lambda, w = w0, cs2 = COOP_REAL_OF(1.))
  end function coop_de_w0

  function coop_de_w_w0wa(a, w0wa) result(w)
    COOP_REAL a, w
    type(coop_arguments) w0wa
    w = w0wa%r(1) + w0wa%r(2)*(1.d0 - a)
  end function coop_de_w_w0wa

  function coop_de_quintessence(Omega_Lambda, epsilon_s, epsilon_inf, zeta_s) result(this)
    type(coop_species) this
    COOP_REAL Omega_Lambda, epsilon_s, epsilon_inf, zeta_s
    type(coop_function) fq
    type(coop_arguments) arg
    arg = coop_arguments(r =  (/ Omega_Lambda, epsilon_s, epsilon_inf, zeta_s /))
    fq = coop_function(coop_de_w_quintessence, xmin = coop_min_scale_factor, xmax = COOP_REAL_OF(1.d0), xlog = .true., arg = arg)
    call this%init(gengre = COOP_SPECIES_FLUID, name = "quintessence Dark Energy", id=5, Omega = Omega_Lambda, cs2 = COOP_REAL_OF(1.d0), fw = fq )
    call arg%free
  end function coop_de_quintessence

  function coop_de_w_quintessence(a, arg) result(w)
    COOP_REAL a, w
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
    if(Omega_m .gt. 0.6d0)then
       w = 1.d0 !!set w to be a crazy value to rule out the model
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
       w = 1.d0  !!set w to be a crazy value to rule out the model
       return
    endif
    mu=a/a_eq
    mu3=mu**3
    s0=sqrt(mu3)
    s1=sqrt(1.+mu3)
    if(s0.lt.2.d-4)then  
       w = qpsign * EPSILON_INFTY/1.5 
    elseif(s0.lt. 1.d-2)then
       w = 2.d0 * qpsign * ( aux1* sqrt(1.+ a_eq/3./(a_eq +  a)) + aux2 * (2.d0/3.d0) * s0 )**2
    else
       w = 2.d0 * qpsign * (aux1*(1. + a_eq/6./(a + a_eq)) + aux2*(s1/s0-log(s0+s1)/mu3) + aux3*(1.-log(1.+mu3)/mu3))**2
    endif
    w = w - 1.d0
#undef EPSILON_S
#undef EPSILON_INF
#undef ZETA_S
  end function coop_de_w_quintessence



end module coop_background
