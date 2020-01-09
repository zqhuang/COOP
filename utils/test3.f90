program Daubechies
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT,parameter::n=119
  COOP_INT::i, k
  COOP_REAL::s(n), cs(0:50), cm(0:50), x(0:50), midterm(n)
  type(coop_asy)::fig
  call coop_set_uniform(51, x, 0.d0, 5000.d0)
  cm = 0.d0
  cs = 0.d0
  open(11, file="shua.txt")
  do i=1, n
     read(11,*) s(i), midterm(i)
     k = nint(midterm(i)/100.d0)
     cm(k) = cm(k) + 1
     k = nint(s(i)/100.d0)
     cs(k) = cs(k) + 1     
  enddo
  close(11)
  call fig%open("count_shuati.txt")
  call fig%init(xlabel="number of exercises done", ylabel = "midterm score", xmin=0., xmax = 4500., ymin = 10., ymax = 105.)
  call fig%dots( s, midterm)
  call fig%close()
end program Daubechies
