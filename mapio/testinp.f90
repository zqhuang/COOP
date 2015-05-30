program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask, m2
  type(coop_healpix_inpaint)::inp
  COOP_INT,parameter ::lmax = 100
  COOP_INT,parameter ::nrun = 150
  COOP_REAL::cls(0:lmax), Cls_ave(0:lmax), ells(0:lmax)
  COOP_INT::l, ell, i, irun
  type(coop_file)::fp
  type(coop_asy)::fig
  call coop_MPI_init()
  call coop_random_init()
  call coop_set_uniform(lmax+1, ells, 0.d0, dble(lmax))
  call fp%open("planck14best_lensedCls.dat", "r")
  do l=2, lmax
     read(fp%unit, *) i, Cls(l)
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)
  enddo  
  call fp%close()
  Cls(0:1) = 0.d0  
  call mask%read("planck14/dx11_v2_commander_int_mask_040a_0256.fits")
  Cls_ave = 0.d0
  call map%init(nside = 256, nmaps = 1, genre = "TEMPERATURE", lmax = lmax)
  do irun = 1, nrun
     map%cl(0:lmax,1) = cls            
     call map%simulate()  
     write(*,*) "*************** run #", irun , " ********************"
     call inp%init(map, mask, lmax, cls)
     do i=1, 4
        call inp%upgrade()
     enddo
     inp%lMT%map = inp%lCT%map + inp%lMT%map
     call inp%lMT%map2alm(lmax = lmax)
     Cls_ave(2:lmax) = Cls_ave(2:lmax)+inp%lMT%Cl(2:lmax, 1)     
     write(*,*) "deviation: ", sum(abs((Cls_ave(2:20)/irun/Cls(2:20)-1.d0)/sqrt(1.d0/(ells(2:20)+0.5d0))))/19.d0, sum(abs((Cls_ave(2:30)/irun/Cls(2:30)-1.d0)/sqrt(1.d0/(ells(2:30)+0.5d0))))/29.d0,  sum(abs((Cls_ave(2:40)/irun/Cls(2:40)-1.d0)/sqrt(1.d0/(ells(2:40)+0.5d0))))/39.d0
     if(mod(irun,50).eq.0)then
        call fp%open("clsout/clsout"//COOP_STR_OF(irun)//".dat", "w")
        do l=2, 40
           write(fp%unit, "(I8, 2E16.7)") l, Cls(l)*l*(l+1.d0)/coop_2pi, Cls_ave(l)*l*(l+1.d0)/coop_2pi/irun
        enddo
        call fp%close()
     endif
  enddo
  call coop_MPI_finalize()  
end program test
