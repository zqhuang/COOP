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
  type(coop_healpix_maps)::hmap, lmap
  COOP_INT,parameter::lmax = 500
  COOP_REAL::Cls(0:lmax), sqrtCls(0:lmax), Cl_ll(0:lmax), Cl_lh(0:lmax), Cl_hh(0:lmax), norm
  type(coop_file)::fp
  COOP_INT::l, basenside, i, irun
  COOP_INT,parameter::nrun = 2500
  Cls(0:1) = 0.d0
  if(iargc().ge.1)then
     basenside = coop_str2int(coop_InputArgs(1))
  else
     print*, "enter nside"
     read(*,*) basenside
  endif
  call fp%open("planck14best_lensedCls.dat", "r") 
  do l=2, lmax
     read(fp%unit, *) i, Cls(l)
     if(i.ne.l) stop "error in cl file"
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)*coop_gaussian_filter(10.d0, l)**2
  enddo
  sqrtCls = sqrt(Cls)
  call fp%close()
  call lmap%init(nside = 8,  nmaps = 1, genre = "TEMPERATURE")
  call hmap%init(nside = basenside*8, nmaps = 2, genre = "TEMPERATURE", lmax = lmax)
  Cl_hh = 0
  Cl_ll = 0
  Cl_lh = 0
  coop_healpix_alm_check_done = .false.       
  do irun = 1, nrun
     write(*,*) irun
     call hmap%simulate_Tmaps(hmap%nside, lmax, sqrtCls)
     call coop_healpix_maps_ave_udgrade(from = hmap, to = lmap, imap_from = 1, imap_to = 1)
     call coop_healpix_maps_ave_udgrade(from = lmap, to = hmap, imap_from = 1, imap_to = 2)
     hmap%map(:,1) = hmap%map(:,1) - hmap%map(:,2)
     call hmap%map2alm(lmax = lmax)
     Cl_hh(2:hmap%lmax)  = Cl_hh(2:hmap%lmax)  + hmap%cl(2:hmap%lmax,coop_matsym_index(2, 1, 1))
     Cl_ll(2:hmap%lmax)  = Cl_ll(2:hmap%lmax)  + hmap%cl(2:hmap%lmax,coop_matsym_index(2, 2, 2))     
     Cl_lh(2:hmap%lmax)  = Cl_lh(2:hmap%lmax)  + hmap%cl(2:hmap%lmax,coop_matsym_index(2, 1, 2))
  enddo
  Cl_hh = Cl_hh/nrun
  Cl_ll = Cl_ll/nrun
  Cl_lh = Cl_lh/nrun
  call fp%open("filter"//COOP_STR_OF(basenside)//".dat",'w')
  do l=2, hmap%lmax
     write(fp%unit,"(I8, 4E16.7)") l, Cl_hh(l)/Cls(l), Cl_ll(l)/Cls(l), Cl_lh(l)/sqrt(Cl_hh(l)*Cl_ll(l))
  enddo
  call fp%close()
end program test
