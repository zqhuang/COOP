program bgtest
  use coop_wrapper_background
  implicit none
#include "constants.h"  
  type(coop_cosmology_background)::bg
  type(coop_asy)::fig
  type(coop_function)::wp1, alpha_M
  COOP_REAL_ARRAY::lna, dlnHdlna, alpha
  COOP_REAL::H, a
  COOP_INT::err, i
  call bg%init(h=0.68d0)
  call bg%add_species(coop_baryon(0.049d0))
!  call bg%add_species(coop_cdm(0.26d0))  
  call bg%add_species(coop_radiation(bg%Omega_radiation()))
  call bg%add_species(coop_neutrinos_massless(bg%Omega_massless_neutrinos_per_species()*(bg%Nnu())))
  call wp1%init_polynomial( (/ 0.d0, 0.05d0 /) )
  call alpha_M%init_polynomial( (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.05d0 /) )
  call coop_background_add_coupled_DE_with_w(bg, omega_c=0.26d0, epsilon_s = 0.02d0, epsilon_inf = 0.2d0, zeta_s = 0.d0, Q = 0.d0, dlnQdphi = 0.d0, err=err)
  print*, err
  if(err .ne. 0) stop
!  call coop_background_add_EFT_DE(bg, wp1, alpha_M, err)
 ! call bg%add_species(coop_de_w0wa(bg%Omega_k(), -0.9d0, 0.d0))
!!test energy conservation    
  call coop_set_uniform(coop_default_array_size, lna, -0.2d0, 0.d0)


  a = exp(lna(1))
  do i=1, coop_default_array_size
     dlnHdlna(i)  = bg%HdotbyHsq(exp(lna(i)))
  enddo
  dlnHdlna(1) = dlnHdlna(1)/2.d0
  dlnHdlna(coop_default_array_size) =   dlnHdlna(coop_default_array_size)/2.d0
  H = exp(- sum(dlnHdlna)*(lna(2)-lna(1)))
  print*, 3.d0*H**2/(bg%rhoa4(a)/a**4)-1.d0
end program bgtest
