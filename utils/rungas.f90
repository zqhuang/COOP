program Test
  use coop_gas_mod
  implicit none
#include "constants.h"
  COOP_INT, parameter::n1 = 1, n2 = 10000
  type(coop_gas_box)::box
  COOP_INT::i, j, ir
  COOP_REAL::dt, pressure, pave, Tave, r, Veff

!  box%do_collision = .false.
  
  do ir = 0, 9
     r = 0.02d0*exp(-ir/3.d0)
     call box%init(n = 16384, L = 1.d0, r = r)
     call random_seed()
     call random_number(box%x)
     box%x = box%r + (box%L - 2.d0*box%r)*box%x
     do i=1, box%n
        box%v(1, i) =  coop_random_Gaussian()
        box%v(2, i) =  coop_random_Gaussian()
        if(coop_gas_dim .eq. 3)then
           box%v(3, i) = coop_random_Gaussian()
        endif
     enddo
     dt = 0.1d0/box%omega
     write(*,"(54A14)") "# r", " p ", " nRT ", " nRT/p   "
     do i = 1, n1
        pave = 0.d0
        Tave = 0.d0
        do j=1, n2
           call box%evolve_x(dt)
           call box%evolve_v(dt, pressure)
           pave = pressure + pave
           Tave = Tave + sum(box%v**2)
        enddo
        pave = pave/n2
        Tave = Tave/n2/coop_gas_dim
        write(*,"(5G14.5)") r, pave, Tave, Tave/pave
     enddo
  enddo
end program Test  
