program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask, m2
  type(coop_healpix_inpaint)::inp
  COOP_INT,parameter ::lmax = 100
  COOP_INT,parameter ::nrun = 1000 
  COOP_REAL::cls(0:lmax), Cls_ave(0:lmax), Cls_pseudo(0:lmax)
  COOP_INT::l, ell, i, irun
  type(coop_file)::fp
  call coop_MPI_init()
  call coop_random_init()
  call fp%open("planck14best_lensedCls.dat", "r")  !!fiducial Cl's
  do l=2, lmax
     read(fp%unit, *) i, Cls(l)
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)
  enddo  
  call fp%close()
  Cls(0:1) = 0.d0
  Cls_ave = 0.d0
  call map%read("planck14/dx11_v2_smica_int_cmb_010a_1024.fits")
  call mask%read("planck14/dx11_v2_smica_int_mask_010a_1024.fits")
  map%map = map%map*mask%map
  call map%map2alm(lmax = lmax)
  Cls_pseudo(0:lmax) = map%Cl(0:lmax, 1)
  call inp%init(map, mask, lmax, cls)
  do irun = 1, nrun
     write(*,*) "***** ensemble # ", irun, " ***********"
     call inp%upgrade(reset = .true.)  !!nside -> 16
     do i=1, 2
        call inp%upgrade()    !!nside -> 32 -> 64
     enddo
     inp%lMT%map = inp%lCT%map + inp%lMT%map
!!$!======== if you want to see how the inpainted map looks like"
!!$call inp%lMT%write("inpainted_map.fits")
!!$stop     
     call inp%lMT%map2alm(lmax = lmax)
     Cls_ave(2:lmax) = Cls_ave(2:lmax)+inp%lMT%Cl(2:lmax, 1)
     if(mod(irun, 20).eq. 0)then
        call fp%open("clsout/clsout"//COOP_STR_OF(irun)//".dat", "w")
        write(fp%unit, "(A8, 3A16)") "# ell ",  "  model C_l  ", "  pseudo C_l ", " ave Cl  "
        do l=2, lmax/2
           write(fp%unit, "(I8, 3E16.7)") l, Cls(l)*l*(l+1.d0)/coop_2pi, Cls_pseudo(l)*l*(l+1.d0)/coop_2pi, Cls_ave(l)*l*(l+1.d0)/coop_2pi/irun
        enddo
        call fp%close()
     endif
  enddo
  call coop_MPI_finalize()  
end program test
