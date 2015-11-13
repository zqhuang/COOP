program bgtest
  use coop_wrapper_background
  implicit none
#include "constants.h"
  COOP_REAL:: omegab = 0.05d0
  COOP_REAL:: omegac = 0.25d0
  COOP_REAL::w = -0.9d0
  COOP_REAL::cs2 = 0.9d0
  COOP_REAL::r_B, r_H, r_M, r_T, omegam
  type(coop_cosmology_background)::bg
  COOP_INT, parameter::nvars = 4
  COOP_INT, parameter::n = 512
  type(coop_function)::wp1, alpha_M, alpha_B, alpha_K, alpha_T, alpha_H
  COOP_INT::err, i, na
  type(coop_ode)::ode
  COOP_REAL::y(nvars), yp(nvars), a, da
  COOP_INT::index_de
  logical success
  omegam = omegab + omegac
  r_B = 1.d0
  r_H = 0.21d0
  r_M = 0.11d0
  r_T = 0.8d0
  call coop_de_construct_alpha_from_cs2(omegam, w, cs2, r_B, r_H, r_M, r_T, bg%f_alpha_B, bg%f_alpha_H, bg%f_alpha_K, bg%f_alpha_M, bg%f_alpha_T, success)
  if(.not. success) stop "cannot find alpha0"
  call wp1%init_polynomial( (/ 1.d0+w/) )
  call coop_EFT_DE_set_Mpsq(bg%f_alpha_M)
  call bg%init(h=0.68d0)
  call bg%add_species(coop_baryon(omegab))
  call bg%add_species(coop_cdm(omegac))  
  call bg%add_species(coop_radiation(bg%Omega_radiation()))
  call bg%add_species(coop_neutrinos_massless(bg%Omega_massless_neutrinos_per_species()*(bg%Nnu())))

  call coop_background_add_EFT_DE_with_effective_w(bg, wp1, err)
  if(err .ne. 0) stop "cannot initialize EFT DE;"
  print*, bg%omega_k()
  
  call bg%setup_background()
  call ode%init(nvars)
  call ode%set_initial_conditions(0.d0, (/ bg%HdotbyHsq(1.d0), log(bg%Hratio(1.d0)), log(bg%species(5)%density(1.d0)), log(coop_Mpsq(1.d0)) /) )
  call ode%evolve(getderv, -1.d0)

  print*, ode%y(1),  ode%y(1)/bg%HdotbyHsq(exp(ode%x))-1.d0
  print*, ode%y(2), ode%y(2) - log(bg%Hratio(exp(ode%x)))
  print*, ode%y(3), ode%y(3) - log(bg%species(5)%density(exp(ode%x)))
  print*, ode%y(4), ode%y(4) - log(coop_Mpsq(exp(ode%x)))
  call coop_de_construct_alpha_from_cs2(omegam, w, cs2, r_B, r_H, r_M, r_T, bg%f_alpha_B, bg%f_alpha_H, bg%f_alpha_K, bg%f_alpha_M, bg%f_alpha_T, success)
  na = 1000
  da = 1.d0/na
  do i = 1, na
     print*, i, i*da, bg%alphacs2(i*da)/ bg%total_alpha(i*da), bg%alpha_K(i*da)
  enddo

contains

  subroutine getderv(n, lna, y, yp)
    COOP_INT::n
    COOP_REAL::lna, y(n), yp(n), a
    a = exp(lna)
    yp(1) = bg%HddbyH3(a) - bg%HdotbyHsq(a)**2 * 2.d0
    yp(2) = bg%HdotbyHsq(a)
    yp(3) = -3.d0*bg%species(5)%wp1effofa(a)
    yp(4) = coop_alphaM(a)
  end subroutine getderv

  
end program bgtest
