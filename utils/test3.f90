module scalar
  use coop_wrapper_utils
  implicit none
#include "constants.h"
contains

  function Lagrangian(X, phi)
    
  end function Lagrangian
  

end module

program Test
  use coop_wrapper_utils
  use scalar
  implicit none
#include "constants.h"

  COOP_INT,parameter::n=1024
  COOP_REAL::k(n), P(n)
  COOP_INT::i
  type(coop_asy)::fig
  
  call fig%open("Ps.txt")
  call fig%init(xlabel = "$k$", ylabel = "$P$", xlog=.true., ylog=.true., xmin=1.e-4, xmax = 1.0, ymin = 1., ymax = 10.)  
  call coop_set_uniform(n, k, log(1.d-4), log(1.d0))
  k = exp(k)
  do i=1, n

  enddo
  call fig%plot(k, P)
  call fig%close()

contains

  
end program Test  
