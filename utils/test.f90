program TestNpeak
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_nd_prob)::bao
  COOP_REAL::x(2), y
  type(coop_file)::fp
  integer i
  call fp%open("testfile.txt", "w")
  do i=1, 10000
     call random_number(x)
     x = (x*2.d0-1.d0)*3.d0
     y=exp(-sum(x**2)/2.d0 - x(1)**2*x(2)**2*0.1)
     write(fp%unit, "(10E16.7)") y, -log(y), x
  enddo
  call fp%close()
  call bao%load("testfile.txt", form = "mcmc", nvars = 2, name = "DR11 BAO Like")

  do
     write(*,*) "enter x"
     read(*,*) x
     y = bao%eval(x)
     print*, "f(x) = ", y, sum(x**2)/2.d0 + x(1)**2*x(2)**2*0.1
  enddo
  
  
end program TestNpeak
