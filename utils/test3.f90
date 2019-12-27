program Daubechies
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT,parameter::n=100
  COOP_REAL,dimension(n)::t, v
  COOP_INT::i
  type(coop_asy)::fig
  call coop_set_uniform(n, t, 0.0001d0, 4.d0)
  call fig%open("fig_homotopy.txt")
  call fig%init(xlabel="$t$", ylabel="$\upsilon$", width=6., height=4.5, xmin = 0., xmax = 4.2, ymin = -0.1, ymax = 1.2)
  v = tanh(t)
  call fig%plot(t, v, legend="Exact")
  v = 1.d0-exp(-t)
  call fig%plot(t, v, color="red", linetype="dotted", legend="Approx0")
  v = 0.d0
  do i = 1, 200
     v = v+exp(-i*t)/(i*i)
  enddo
  v = 1.d0 +exp(-t)*(t-coop_pi2/6.d0-1.d0+v)+(1.d0-exp(-t))*log(1.d0-exp(-t))
  call fig%plot(t, v, color="blue", linetype="dotted", legend="Approx1")
  print*, v(1)
  call fig%legend()
  call fig%close()
end program Daubechies
