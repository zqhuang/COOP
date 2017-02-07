program fR1d
  use coop_wrapper_firstorder
  use fR1d_mod
  implicit none
#include "constants.h"
  type(coop_fr1d_obj)::halo
  type(coop_asy)::fig
  call halo%init(Omega_m = 0.3d0, nr = 1024, rmax = 0.2d0, ns = 2048, a_ini = 0.03d0, delta_ini = 0.01d0, r_halo = 5.d-2, bw_halo = 2.d-3)
  call fig%open("rho.txt")
  call fig%init(xlabel = "$r$", ylabel = "$\rho$")
  call fig%plot(halo%r, halo%rho)
  call fig%close()
end program fR1d
