program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT, parameter::n = 100
  COOP_REAL,parameter::t0 = 300.d0, p0 = 1.01325d5
  COOP_REAL,parameter::nuR = 8.31d0
  COOP_REAL,parameter::V0 = nuR*t0/p0
  COOP_REAL,parameter::nub = V0*0.01d0
  COOP_REAL,parameter::anu2 = p0*0.01d0*v0**2
  COOP_REAL:: p(n), T(n), Hp(n, n)
  COOP_INT::i, j
  type(coop_asy)::fig
  type(coop_file)::fp
  
!!$  p(1:7) = (/ 1,2,3,4,5,6,7 /)
!!$  call fp%open("numbers.dat", 'u')
!!$  write(fp%unit) p(1:7)
!!$  call fp%close()
  call coop_binfile_encrypt("numbers.dat", "ncopy.dat")
  call coop_binfile_decrypt("ncopy.dat", "ncopy2.dat")  
  stop
  call coop_set_uniform(n, p, 0.05d0*p0, 50.d0*p0)
  call coop_set_uniform(n, T, 0.02d0*t0, 2.d0*t0)
  do j=1, n
     do i=1, n
        Hp(i, j) = dHdp(p(i), T(j))*p0
     enddo
  enddo
  call fig%open("dHdp.txt")
  call fig%init(xlabel="$p/\mathrm{atm}$", ylabel = "$T/K$",width=4., height=3.2)
  call fig%density(z=Hp, xmin=p(1)/p0, xmax = p(n)/p0, ymin = T(1), ymax=T(n), label = "$\left(\frac{\partial H}{\partial p}\right)_T \left[J/\mathrm{atm}\right]$",color_table="Planck", zmin=-100.d0, zmax=100.d0)
  call fig%close()

contains

  function dHdp(p, T)

    COOP_REAL::p, T, dHdp, V
    COOP_INT::i
    V = nuR*T/max(p, anu2) + nub
    do while(abs((p+anu2/V**2)*(V-nub)/nuR/T -1.d0).gt. 1.d-6)
       V= nuR*T/(p+anu2/V**2)+nub
    enddo
    
    dHdp = V-nuR*T/(p+anu2/V**2*(1 - 2*(1-nub/V))) 
  end function dHdp

end program Test  
