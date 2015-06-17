program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask, m2, ltmp
  type(coop_healpix_inpaint)::inp
  COOP_INT,parameter ::lmax = 128, nside = 256, nrun = 100
  COOP_REAL::cls(0:lmax), sqrtCls(0:lmax), Cls_ave(0:lmax), ells(0:lmax), Cls_sim(0:lmax)
  COOP_INT::l, ell, i, irun
  type(coop_file)::fp
  type(coop_asy)::fig
  call coop_MPI_init()
  call coop_random_init()
  !call coop_healpix_latitude_cut_mask(nside = nside, latitude_degree = 15.d0, filename = "planck14/lat15_mask_n"//COOP_STR_OF(nside)//".fits")
  call coop_set_uniform(lmax+1, ells, 0.d0, dble(lmax))
  call fp%open("planck14best_lensedCls.dat", "r")
  do l=2, lmax
     read(fp%unit, *) i, Cls(l)
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)*coop_gaussian_filter(60.d0, l)**2
  enddo  
  call fp%close()
  Cls(0:1) = 0.d0
  sqrtCls  = sqrt(Cls)
  !  call mask%read("lowl/mask_tiny.fits")
  ! call mask%read("lowl/mask_hot_bar_n0256.fits")              
  ! call mask%read("lowl/commander_dx11d2_mask_temp_n0256_likelihood_v1.fits")  
 ! call mask%read("planck14/lat15_mask_n"//COOP_STR_OF(nside)//".fits")
   call mask%read("planck14/dx11_v2_commander_int_mask_040a_0256.fits")
  Cls_ave = 0.d0
  Cls_sim = 0.d0
  call map%init(nside = nside, nmaps = 1, genre = "I", lmax = lmax)
  call m2%init(nside = nside, nmaps = 1, genre = "I", nested = .true., lmax = lmax)  
  !  call map%read("lowl/commander_dx11d2_extdata_temp_cmb_n0256_60arc_v1_cr.fits")
  do irun = 1, nrun
     print*, "sim #", irun
     call map%simulate_Tmaps(nside = nside, lmax = lmax, sqrtCls = sqrtCls)
     call inp%init(map, mask, lmax, cls)
     call inp%upgrade(nside_want = 8)
     inp%lMT%map = inp%lMT%map+inp%lCT%map
     call coop_healpix_maps_ave_udgrade(inp%lMT, m2)
     call m2%map2alm()
     cls_ave = cls_ave + m2%cl(0:lmax,1)
     call inp%upgrade(nside_want = 16)
     inp%lMT%map = inp%lMT%map+inp%lCT%map
     call coop_healpix_maps_ave_udgrade(inp%lMT, m2)
     call m2%map2alm()
     cls_sim = cls_sim + m2%cl(0:lmax,1)
  enddo
  cls_ave = cls_ave/nrun
  cls_sim=Cls_sim/nrun
  call fp%open("clr.dat","w")
  do l=2, 50
     write(fp%unit,*) l, cls_ave(l)/cls(l), Cls_sim(l)/Cls(l)
  enddo
  call fp%close()
  call coop_MPI_finalize()  
end program test
