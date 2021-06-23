program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::cosmology
  type(coop_real_table)::paramtable
  COOP_INT,parameter::n = 500
  COOP_REAL,dimension(n)::a, xe, z
  type(coop_asy)::fig
  COOP_INT::i, j
  logical success
  call paramtable%insert("ombh2", 0.02225d0)
  call paramtable%insert("omch2", 0.1194d0)
  call paramtable%insert("h", 0.6748d0)
  call paramtable%insert("tau", 0.058d0)
  call paramtable%insert("As", 2.1d-9)
  call paramtable%insert("ns", 0.965d0)
  call cosmology%set_up(paramtable, success)
  if(.not. success)then
     stop "initialization failed"
  endif
  call coop_set_uniform(n, z,  0.d0, 4.d0)
  z = 10.d0**z
  do i=1, n
     xe(i) = cosmology%xeofa(1./(1.+z(i)))
  enddo
  call fig%open("xeofz.txt")
  call fig%init(width=5., height = 3.6, xlabel="$z$", ylabel = "$X_e$", xlog=.true., xmin=1., xmax=1.e4, ymin = -0.02, ymax = 1.2) !, ylog=.true.)
  call fig%plot(z, xe, linewidth=2., color="blue")
  call fig%close()

end program test
