program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools

  implicit none
#include "constants.h"
  COOP_INT,parameter::lmax = 2500
  COOP_INT,parameter::lmin = 250
  COOP_INT,parameter::nside = 2048
!  COOP_INT,parameter::nside_output = 16  
  COOP_REAL,parameter::beam_fwhm = 5.
  COOP_INT,parameter::hp_lowl = lmin-10
  COOP_INT,parameter::hp_highl = lmax+10
!!  logical,parameter::do_constrained = .false.
  type(coop_healpix_maps)::map,  lmap, mask, hmap
  integer l, m, il, i
  type(coop_file)::fp
  COOP_UNKNOWN_STRING, parameter::prefix = "act16/sim"
  COOP_REAL:: Cls(0:lmax)
!  type(coop_healpix_inpaint)::inp
  COOP_REAL::sigmasq
  call coop_MPI_Init()  
  call coop_random_init()
  !!load Cls
  Cls = 0.d0
  sigmasq = 0.
  call fp%open("planck14best_lensedCls.dat", "r")
  do l=2, lmax
     read(fp%unit, *) il, Cls(l)
     if(il.ne. l) stop "Cl file error"     
     Cls(l) = Cls(l)*(coop_2pi/l/(l+1.d0))*coop_gaussian_filter(beam_fwhm, l)**2*coop_highpass_filter(hp_lowl, hp_highl, l)**2
     sigmasq = sigmasq + (2.*l+1.)*Cls(l)
  enddo
  sigmasq = sigmasq/coop_4pi
  call fp%close()
  call map%init(nside = nside, nmaps=1, genre="TEMPERATURE", lmax = lmax)
  map%Cl = 0.
  map%Cl(2:lmax,1) = Cls(2:lmax)        
  call map%simulate()
  call map%write(trim(prefix)//"_n"//COOP_STR_OF(map%nside)//"_"//COOP_STR_OF(nint(beam_fwhm))//"a_l"//COOP_STR_OF(lmin)//"-"//COOP_STR_OF(lmax)//".fits")
  write(*,*) "map rms expected to be about: ", sqrt(sigmasq)
  write(*,*) "actual map rms: ", sqrt(sum(map%map**2)/map%npix)
!!$  if(do_constrained)then
!!$     call map%read("lowl/commander_dx11d2_extdata_temp_cmb_n0256_60arc_v1_cr.fits")
!!$     call lmap%init(nside = nside_output, nmaps = 1, genre = "TEMPERATURE", nested = .true.)
!!$     call hmap%init(nside = 64, nmaps = 1, genre = "TEMPERATURE", nested = .true.)          
!!$     call mask%read("lowl/mask_hot_bar_n0256.fits")
!!$     call map%apply_mask(mask)
!!$     call inp%init(map = map, mask = mask, lmax = lmax, cls = cls)
!!$     do i=1, 1000
!!$        write(*,*) "simulating map #"//COOP_STR_OF(i)
!!$        call inp%upgrade(reset = .true.)
!!$        call inp%upgrade()
!!$        call inp%upgrade()
!!$        call inp%upgrade()
!!$        hmap%map = inp%lMT%map + inp%lCT%map
!!$        call hmap%smooth(coop_SI_arcmin*beam_fwhm)
!!$        call coop_healpix_maps_ave_udgrade(hmap, lmap)
!!$        call lmap%write(trim(prefix)//"_inp_"//COOP_STR_OF(lmap%nside)//"_"//COOP_STR_OF(nint(beam_fwhm))//"a_"//COOP_STR_OF(i)//".fits")        
!!$     enddo

end program test
