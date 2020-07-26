program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
#ifdef HAS_HEALPIX
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  use udgrade_nr
  use coord_v_convert,only:coordsys2euler_zyz
#endif  
  implicit none
#include "constants.h"
  COOP_INT,parameter::nrun = 1000  
  type(coop_healpix_maps)::mask, map, hm
  type(coop_healpix_inpaint)::inp
  COOP_REAL,dimension(:),allocatable::cls, cls_ave, Cls_sim
  COOP_REAL::cl1, cl2, cross    
  type(coop_file)::fp, fp2
  COOP_INT::i, j, l, irun, nside, lmax
  call coop_get_command_line_argument(key = "nside", arg = nside)
  call coop_healpix_nside2lmax(nside,lmax)
  allocate(cls(0:lmax), Cls_ave(0:lmax))
  call fp%open("planck14best_lensedCls.dat", "r")
  do l=2, lmax
     read(fp%unit, *) i, Cls(l)
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)*coop_gaussian_filter(60.d0, l)**2
  enddo  
  call fp%close()
  Cls(0:1) = 0.d0
  call mask%init(nside = nside, nmaps = 1, genre = "MASK")
  mask%map = 0.
  call map%init(nside = nside, nmaps = 1, genre = "I")
  map%map  = 0.
  call hm%init(nside = nside, nmaps = 1, genre = "I")
  call inp%init(map = map, mask = mask, lmax = lmax, Cls = Cls)
  Cls_ave = 0.d0
  do irun = 1, nrun
     write(*,*) irun
     call inp%upgrade(reset = .true.)
     call coop_healpix_maps_ave_udgrade(inp%lCT, hm)
     call hm%map2alm()
     Cls_ave = Cls_ave + hm%Cl(0:lmax,1)
  enddo
  Cls_ave = Cls_ave/nrun/Cls
  call fp2%open("healpix_filters/FLH"//COOP_STR_OF(nside)//".dat", "r")  
  call fp%open("healpix_filters/InpF"//COOP_STR_OF(nside)//".dat","w")
  do l = 2, lmax
     read(fp2%unit, *) i, cl1, cl2, cross
     if(i.ne. l) stop "error"
     write(fp%unit, "(I8, 3E16.7)") l, cl1, Cls_ave(l), cross
  enddo
  call fp%close()
  call fp2%close()
end program test
