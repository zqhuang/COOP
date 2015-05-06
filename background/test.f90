program bgtest
  use coop_wrapper_background
  implicit none
#include "constants.h"  
  type(coop_cosmology_background)::bg
  COOP_INT, parameter::n = 256
  COOP_REAL:: a(n), wp1(n), wp1eff(n)
  COOP_INT::i, err
  COOP_REAL::omega_c, epsilon_s, epsilon_inf, zeta_s, Q, dlnQdphi
  type(coop_asy)::fig
  type(coop_species)::de
  omega_c = 0.25d0
  epsilon_inf = 0.2d0
  epsilon_s = 0.25d0
  zeta_s = 0.5d0
  Q = 0.2d0
  dlnQdphi = 0.d0
  call bg%init(h=0.68d0)
  call bg%add_species(coop_baryon(0.049d0))
  call bg%add_species(coop_radiation(bg%Omega_radiation()))
  call bg%add_species(coop_neutrinos_massless(bg%Omega_massless_neutrinos_per_species()*(bg%Nnu())))
  de = coop_de_quintessence(Omega_Lambda = bg%Omega_k() - Omega_c, epsilon_s = epsilon_s, epsilon_inf = epsilon_inf, zeta_s = zeta_s)
  call coop_background_add_coupled_DE_with_w(bg, omega_c, epsilon_s, epsilon_inf, zeta_s, Q, dlnQdphi, err)
  if(err .ne. 0) stop "error"
  call coop_set_uniform(n, a, 0.1d0, 1.d0)
  do i=1, n
     wp1(i) = de%wp1ofa(a(i))
     wp1eff(i) = wp1_eff(a(i))
  enddo
  call fig%open("wa.txt")
  call fig%init(xlabel = "$a$", ylabel = "$1+w$")
  call fig%curve(a, wp1, color="red", linetype = "solid")
  call fig%curve(a, wp1eff, color="blue", linetype = "dotted")  
  call fig%close()
contains

    function wp1_eff(a) result(wp1)
      COOP_REAL::a, wp1
      wp1 = (1.d0 + bg%species(5)%pa2(a)/(bg%species(5)%rhoa2(a)+ bg%species(4)%rhoa2(a) - (3.d0*omega_c)/a))
    end function wp1_eff

    function quint_f(x) result(f)
      COOP_REAL x, f
      f = sqrt(1.d0+x**3)/x**1.5 - log(x**1.5+sqrt(1.d0+x**3))/x**3
    end function quint_f

    function quint_f2(x) result(f2)
      COOP_REAL x, f2
      f2 = coop_sqrt2*(1.d0 - log(1.d0+x**3)/x**3)- quint_f(x)
    end function quint_f2

end program bgtest
