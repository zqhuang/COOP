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
  COOP_INT,parameter::nrun = 100  
  type(coop_healpix_maps)::mask, map, hm
  type(coop_healpix_inpaint)::inp
  COOP_REAL,dimension(:),allocatable::cls, cls_ave, Cls_sim, sqrtCls
  COOP_REAL::cl1, cl2, cross    
  type(coop_file)::fp, fp2
  COOP_INT::i, j, l, irun, nside, lmax
  call coop_MPI_init()
  call coop_random_init()
  call coop_get_command_line_argument(key = "nside", arg = nside)
  call coop_healpix_nside2lmax(nside,lmax)
  allocate(cls(0:lmax), Cls_ave(0:lmax), Cls_sim(0:lmax), sqrtCls(0:lmax))
  call fp%open("planck14best_lensedCls.dat", "r")
  do l=2, lmax
     read(fp%unit, *) i, Cls(l)
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)*coop_gaussian_filter(60.d0, l)**2
  enddo  
  call fp%close()
  Cls(0:1) = 0.d0
  sqrtCls = sqrt(Cls)
  !  call mask%read("lowl/mask_tiny.fits")
  ! call mask%read("lowl/mask_hot_bar_n0256.fits")              
  ! call mask%read("lowl/commander_dx11d2_mask_temp_n0256_likelihood_v1.fits")
  call mask%read("planck14/dx11_v2_commander_int_mask_040a_0256.fits")  
  !call mask%read("planck14/lat30_mask_n256.fits")

  call map%init(nside = nside, nmaps = 1, genre = "I")
  call hm%init(nside = nside, nmaps = 1, genre = "I")
  Cls_ave = 0.d0
  Cls_sim = 0.d0
  do irun = 1, nrun
     write(*,*) irun
     call map%simulate_Tmaps(nside = map%nside, lmax = lmax, sqrtCls = sqrtCls)
     call inp%init(map = map, mask = mask, lmax = lmax, Cls = Cls)
     call inp%upgrade()
    ! call inp%upgrade(reset = .true.)
     inp%lCT%map = inp%lCT%map + inp%lMT%map
     call coop_healpix_maps_ave_udgrade(inp%lCT, hm)
     call hm%map2alm()
     call inp%sim%get_cls()
     Cls_sim = Cls_sim + inp%sim%Cl(0:lmax,1)
     Cls_ave = Cls_ave + hm%Cl(0:lmax, 1)
  enddo
  Cls_ave = Cls_ave/nrun
  Cls_sim = Cls_sim/nrun
  call fp%open("check"//COOP_STR_OF(nside)//".dat","w")
  do l = 2, lmax
     write(fp%unit, "(I5, 10E14.5)") l, Cls_sim(l)/Cls(l), Cls_ave(l)/Cls(l), (inp%filter_mean(l)**2*cls_ave(l)+inp%filter_fluc(l)**2)/Cls(l)
  enddo
  call fp%close()
  call coop_MPI_Finalize()
end program test
