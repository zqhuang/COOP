program derivative
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_REAL, parameter::Mpl = 1.d0
  COOP_REAL, parameter::H0 = 1.d0
  type(coop_function)::V
  type(coop_arguments)::args
  COOP_REAL phi
  integer i
  complex(8)::c
#define MASS args%r(1) 
#define LAMBDA args%r(2)
  print*, kind( (1.d0, 1.d0 ))
  c = cmplx(1.d0, 1.d0)
  print*, c
  stop
  args = coop_arguments( r = (/ 1.e-2*Mpl, 3.*Mpl**2*H0**2*0.7 /) )
  V = coop_function( f = potential, xmin = 0.d0, xmax =20.d0, xlog = .false., ylog = .true., args = args, method = COOP_INTERPOLATE_CHEBYSHEV)
  do i=1, 1000
     phi = 19.9*i/1000
     print*, phi, V%derivative2(phi), d2Vdphi2(phi, args)
  enddo
contains

  function potential(phi, args) result(V)
    type(coop_arguments) args
    COOP_REAL V, phi
    V = cos(phi) + 3. 
  end function potential

  function dVdphi(phi, args) result(Vp)
    type(coop_arguments) args
    COOP_REAL Vp, phi
    Vp = -sin(phi) 
  end function dVdphi

  function d2Vdphi2(phi,args) result(Vpp)
    COOP_REAL Vpp, phi
    type(coop_arguments) args
    Vpp = -cos(phi)

  end function d2Vdphi2

end program derivative
