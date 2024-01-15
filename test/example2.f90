program test
  use coop_wrapper
  implicit none
#include "constants.h"
  COOP_INT, parameter:: n = 128
  COOP_REAL:: x(n), y(n)
  type(coop_asy)::figure
  call coop_set_uniform(n, x, 0.d0, coop_2pi)
  y = sin(x)
  call figure%open("testfig.txt")
  call figure%init(xlabel = "$x$", ylabel = "$y$", width = 5., height = 4.)
  call figure%plot(x, y, color="red", linetype="solid", linewidth=1.5, legend = "$y = \sin x$")
  y = cos(x)
  call figure%plot(x, y, color="blue", linetype="dashed", linewidth=1., legend = "$y = \cos x$")
  call figure%legend(xratio = 0.45, yratio = 0.9)
  call figure%close()
end program test
