program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT,parameter::lmax = 250
  type(coop_file)::fp
  COOP_INT l, bounds(2), lstart, lend
  COOP_REAL:: ell(2:lmax), cl(2:lmax), dcl_up(2:lmax), dcl_down(2:lmax), l_sum, cl_sum, dcl_up_sum, dcl_down_sum, w , center
  call fp%open("cls_commander_2014_rc_TT.dat","r")
  do l=2, lmax
     read(fp%unit, *) center, bounds, cl(l), dcl_up(l), dcl_down(l)
     if(nint(center).ne.l) stop "read error"
     ell(l) = l
  enddo
  call fp%close()
  call fp%open("planck14_cltt.txt")
  l = 1
  do while(l.lt.lmax)
     l_sum = 0.
     cl_sum = 0.d0
     dcl_up_sum = 0.d0
     dcl_down_sum = 0.d0
     w = 0.d0
     lstart = l+1
     do while(w.lt. 15.d0 .and. w .lt. max((l-20.d0)/10.d0, 0.99d0) .and. l.lt.lmax )
        l = l + 1
        l_sum = l_sum + l
        cl_sum = cl_sum + cl(l)
        dcl_up_sum = dcl_up_sum + 1./dcl_up(l)**2
        dcl_down_sum = dcl_down_sum + 1./dcl_down(l)**2
        w = w + 1.d0
     enddo
     lend = l
     l_sum = l_sum/w
     cl_sum = cl_sum/w
     dcl_up_sum = sqrt(1./dcl_up_sum)
     dcl_down_sum = sqrt(1./dcl_down_sum)
     write(fp%unit, "(E16.7, 2I5, 3E16.7)") l_sum, lstart, lend, cl_sum, dcl_up_sum, dcl_down_sum
  enddo
  call fp%close()
end program Test
