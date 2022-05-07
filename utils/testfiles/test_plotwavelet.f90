program Daubechies
  use coop_wrapper_utils
  implicit none
#include "constants.h"
#define ORDER 4
#if ORDER > 5
  COOP_INT,parameter::n=40001
  COOP_REAL::lbd = -7.5d0
  COOP_REAL::rbd = 7.5d0  
#else  
  COOP_INT,parameter::n=20001
  COOP_REAL::lbd = -5.d0
  COOP_REAL::rbd = 5.d0  
#endif
  integer,parameter::m=512
  COOP_INT::i  
  COOP_REAL::s(n), t(n), x(m), y(m)
  type(coop_asy)::fig
  call coop_set_uniform(n, t, lbd, rbd)
  call fig%open("fig_all.txt")
  call fig%init(xlabel="$t$", ylabel="$\psi_{0,0}(t)$", width=6., height=4.5, xmin = -3.5, xmax = 3.5, ymin = -1.3, ymax = 1.9)
  call coop_set_uniform(m,x, -3.d0, 3.d0)
  do i=1, m
     if(x(i) .ge. 0.d0 .and. x(i) .lt. 0.5d0)then
        y(i) = 1.d0
     elseif(x(i).ge.0.5d0 .and. x(i).lt. 1.d0)then
        y(i) = -1.d0
     else
        y(i) = 0.d0
     endif
  enddo
  call fig%plot(x= x, y=y, linewidth=2., color="black", legend="Daubechies-1; Haar")
  
  open(10, FILE="psiwavelet/psi2.csv")
  do i=1,n
     read(10,*) s(i)
  enddo
  close(10)
  do i=1, m
     y(i) = interpsi(x(i))
  enddo
  call fig%plot(x= x, y=y, linewidth=2., color="red", legend="Daubechies-2")

  open(10, FILE="psiwavelet/psi3.csv")
  do i=1,n
     read(10,*) s(i)
  enddo
  close(10)
  do i=1, m
     y(i) = interpsi(x(i))
  enddo
  call fig%plot(x= x, y=y, linewidth=2., color="orange", legend="Daubechies-3")

  open(10, FILE="psiwavelet/psi4.csv")
  do i=1,n
     read(10,*) s(i)
  enddo
  close(10)
  do i=1, m
     y(i) = interpsi(x(i))
  enddo
  call fig%plot(x= x, y=y, linewidth=2., color="skyblue", legend="Daubechies-4")
  
  open(10, FILE="psiwavelet/psi5.csv")
  do i=1,n
     read(10,*) s(i)
  enddo
  close(10)
  do i=1, m
     y(i) = interpsi(x(i))
  enddo
  call fig%plot(x= x, y=y, linewidth=2., color="blue", legend="Daubechies-5")

  call fig%legend(xratio = 0.05, yratio=0.91)
  call fig%close()

contains

  function interpsi(tin)
    COOP_REAL::tin, interpsi
    call coop_linear_interp(n, t, s, tin, interpsi)
  end function interpsi
  
end program Daubechies
