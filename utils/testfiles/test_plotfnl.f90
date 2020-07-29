program Test
#include "constants.h"      
  use coop_wrapper_utils
  implicit none
  type(coop_asy)::fig
  COOP_INT,parameter::n=5
  COOP_REAL::x(n), y(n), dy(n)
  call fig%open("fnl_fisher.txt")
  call fig%init(xlabel="$z$", ylabel="$F_{f_{\rm NL}^{\rm local}f_{\rm NL}^{\rm local}}$", width=4.5, height=3.6, xlog = .false., ylog = .true., xmin=0., xmax=1.5, ymin = 1.e-5, ymax = 0.1)
  x = (/ 0.15, 0.45, 0.75, 1.05, 1.35 /)
  y = (/ 9.2e-5, 7.8e-3, 1.9e-2, 4.9e-3, 5.1e-5 /)
  call fig%dots(x, y, "black", "$\Delta$")
  call fig%label("$\sigma(f_{\rm NL}^{\rm local}) = 5.6$", 0.5, 0.3)
  call fig%close()
end program Test



