program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_asy)::fig
  COOP_INT, parameter::n=100
  COOP_REAL::x(n), yb(n), yt(n)
  call coop_set_uniform(n, x, 0.d0, coop_2pi)
  yb = sin(x)-0.1d0
  yt = sin(x) + cos(x)**2/2.d0
  call fig%open("testband.txt")
  call fig%init(xlabel = "x", ylabel="y")
  call coop_asy_band(fp = fig, x = x, ylower = yb, yupper = yt, colorfill = "lightgray")
  yb = sin(x)-0.05d0
  yt = sin(x) + cos(x)**2/10.d0
  call coop_asy_band(fp = fig, x = x, ylower = yb, yupper = yt, colorfill = "gray")
  call fig%close()
end program Test
