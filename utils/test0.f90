program test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_COMPLEX:: z, y, Hm
  COOP_REAL,dimension(10),parameter::zeta_zeros = (/ 14.134725142d0, 21.022039639d0, 25.010857580d0, 30.424876126d0, 32.935061588d0,     37.586178159d0,     40.918719012d0,      43.327073281d0,      48.005150881d0,    49.773832478d0 /)
  COOP_COMPLEX:: eps = (0.d0, -1.d-3)
  COOP_INT::i, n
  read(*,*) n
  z =dcmplx(0.5d0,  zeta_zeros(n) )
  y = dcmplx(1.d0, 0.d0)

  print*, coop_HmZeta(y, z)
  open(11, file="trajzero_"//COOP_STR_OF(n)//".txt")
  do i=1, 50000
     z = z + eps     
     Hm = coop_Hmzeta(y, z)     
     y = y - Hm/coop_Hmzeta_dy(y, z)
     Hm = coop_Hmzeta(y, z)     
     y = y - Hm/coop_Hmzeta_dy(y, z)
     write(11, "(2E16.7)") y     
     if(mod(i, 20) .eq. 0)then
        Hm = coop_Hmzeta(y, z)
        y = y - Hm/coop_Hmzeta_dy(y, z)        
        print*, y, z, Hm
     endif
     if(real(y)<-0.99999d0 .and. abs(y) > 0.999999) exit     
  enddo
  close(11)
end program test
