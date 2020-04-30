program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::cosmology
  type(coop_real_table)::paramtable
  COOP_INT,parameter::n = 20000
  COOP_REAL,dimension(n)::chi, g, a, vis, ga
  COOP_REAL::chi0, maxg
  COOP_INT::i, j
  logical success
  call paramtable%insert("ombh2", 0.02225d0)
  call paramtable%insert("omch2", 0.1194d0)
  call paramtable%insert("h", 0.6748d0)
  call paramtable%insert("tau", 0.06d0)
  call paramtable%insert("As", 2.1d-9)
  call paramtable%insert("ns", 0.965d0)
  call cosmology%set_up(paramtable, success)
  print*, "distlss = ", cosmology%distlss/cosmology%H0Mpc(), lcdm_chiofa(0.313d0, 0.d0, 1.d0/1090.d0)/cosmology%H0Mpc()
  stop
  if(.not. success)then
     stop "initialization failed"
  endif
  call coop_set_uniform(n, a, -4.d0, -1.d-4)
  a = 10.d0**a
  chi0 = cosmology%distlss
  do i=1, n
     chi(i) = cosmology%tau0 - cosmology%tauofa(a(i))
     vis(i) = cosmology%visofa(a(i))
  enddo
  g = 0.d0
  do i=1, n
     do j=2, i-1
        g(i) = g(i) + (chi(j) - chi(i))/ chi(j) *vis(j)*(chi(j-1)-chi(j+1))
     enddo
     g(i) = g(i)/2.d0/chi(i)
     ga(i) = (chi0-chi(i))/chi0/chi(i)
  enddo

  maxg = maxval(g)
  do i=1, n
     if(g(i).gt. maxg*1.d-5) &
          write(*,*) a(i), g(i), g(i)/ga(i)-1.d0
  enddo
  

end program test
