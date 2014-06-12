program derivative
  use coop_wrapper
  implicit none
#include "constants.h"
  COOP_INT, parameter::n  = 2000
  type(coop_function)::func_V
  COOP_REAL phi(n), V(n)
  COOP_REAL, parameter::phi_min  = 1.d-1, phi_max = 100.d0
  COOP_INT i
  call random_number(phi) !!generate random phi between 0 and 1
  phi = exp( log(phi_min) + phi * log(phi_max/phi_min) )
  do i=1, n
     V(i) = potential(phi(i))
  enddo
  call func_V%init_NonUniform(phi, V, xlog=.true., ylog=.true.)
  print*, potential(1.5d0), func_V%eval(1.5d0)
  print*,dVdphi(99.d0), func_V%derivative(99.d0)
  print*,d2Vdphi2(0.12d0), func_V%derivative2(0.12d0)

contains

  function potential(phi)
    COOP_REAL phi, potential
    potential = phi ** 2 * exp(1./phi**2)
  end function potential


  function dVdphi(phi)
    COOP_REAL phi,dVdphi
    dVdphi = (2.d0/phi-2.d0/phi**3)*potential(phi)
  end function dVdphi

  function d2Vdphi2(phi)
    COOP_REAL phi,d2Vdphi2
    d2Vdphi2 = ((2.d0/phi-2.d0/phi**3)**2-2.d0/phi**2+6.d0/phi**4)*potential(phi)
  end function d2Vdphi2


end program derivative
