program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_REAL::rec3js(10000)
  COOP_INT::twoj1, twoj2, twoj3, twom1, twom2, twom3, J, l, num_j1s
  COOP_REAL::j1, j2, j3, m1, m2, m3, w3j, j1min, j1max
  write(*,*) "j2, m2 = "
  read(*,*) j2, m2
  write(*,*) "j3, m3 = "
  read(*,*) j3, m3
  m1 = -m2 - m3
  j1min = min(abs(j2-j3), abs(m1))
  j1max = j2+j3
  twoj2 = nint(j2*2.d0)
  twoj3 = nint(j3*2.d0)
  twom1 = nint(m1*2.d0)
  twom2 = nint(m2*2.d0)
  twom3 = nint(m3*2.d0)
  call coop_ThreeJ_Array(rec3js, twoj2, twoj3, twom2, twom3, j1min, num_j1s)
  do l = 1, num_j1s
     j1 = j1min + l-1
     twoj1 = nint(j1*2.d0)
     w3j = coop_ThreeJSymbol(twoj1, twoj2, twoj3, twom1, twom2, twom3)
     write(*,*) j1, w3j, rec3js(l)
  enddo

 
end program Test
