program test
  use coop_wrapper_utils
  use coop_background_mod
  use coop_wrapper
  implicit none
#include "constants.h"
  integer,parameter::n = 30
  integer i, j
  COOP_REAL::w(n), omch2(n), h(n, n)
  type(coop_asy)::fig
  call coop_set_uniform(n, w, -0.8d0, 0.8d0)
  call coop_set_uniform(n, omch2, 0.1d0, 0.14d0)
  do j=1, n
     do i=1,n
        COOP_COSMO_PARAMS = coop_arguments(r= (/ 0.02214d0, omch2(i), 1.041d-2,  0.088d0, 0.d0, 0.06d0, w(j), 0.d0, 0.d0 /), i = (/ COOP_DE_QUINTESSENCE, 7, 3 /) )
        call coop_setup_global_cosmology()
        h(i,j) =  COOP_COSMO%h()
     enddo
  enddo
  h = h*100.d0
  call fig%open("hepsomch2.txt")
  call fig%init(xlabel = "$\Omega_c h^2$", ylabel = "$\epsilon_s$")
  call coop_asy_density(fig, h, xmin = omch2(1), xmax = omch2(n), ymin = w(1), ymax = w(n), label = "$H_0 ({\rm km\, s^{-1}\,Mpc^{-1}})$" )
  call coop_asy_label(fig, "$\theta_{\rm CMB} = 0.01041$ (Planck constraint)",  x = omch2(n/2), y = w(n/2), color="black") 
  call fig%close()
end program test
