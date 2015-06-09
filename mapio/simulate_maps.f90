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
  COOP_INT,parameter::lmax = 100
  COOP_INT,parameter::nside = 256
  COOP_INT,parameter::nside_output = 16  
  COOP_REAL,parameter::beam_fwhm = 440.
  logical,parameter::do_constrained = .true.
  type(coop_healpix_maps)::map,  lmap, mask, hmap
  integer l, m, il, i
  type(coop_file)::fp
  COOP_UNKNOWN_STRING, parameter::prefix = "mocklowl/sim"
  COOP_REAL:: sigma, Cls(0:lmax)
  type(coop_healpix_inpaint)::inp
  call coop_MPI_Init()  
  call coop_random_init()
  !!load Cls
  Cls = 0.d0
  sigma = coop_sigma_by_fwhm * beam_fwhm * coop_SI_arcmin  
  call fp%open("planck14best_lensedCls.dat", "r")
  do l=2, lmax
     read(fp%unit, *) il, Cls(l)
     if(il.ne. l) stop "Cl file error"     
     Cls(l) = Cls(l)*(coop_2pi/l/(l+1.d0))*coop_gaussian_filter(beam_fwhm, l)**2
  enddo
  call fp%close()
  if(do_constrained)then
     call map%read("lowl/commander_dx11d2_extdata_temp_cmb_n0256_60arc_v1_cr.fits")
     call lmap%init(nside = nside_output, nmaps = 1, genre = "TEMPERATURE", nested = .true.)
     call hmap%init(nside = 64, nmaps = 1, genre = "TEMPERATURE", nested = .true.)          
     call mask%read("lowl/mask_hot_bar_n0256.fits")
     call map%apply_mask(mask)
     call inp%init(map = map, mask = mask, lmax = lmax, cls = cls)
     do i=1, 1000
        write(*,*) "simulating map #"//COOP_STR_OF(i)
        call inp%upgrade(reset = .true.)
        call inp%upgrade()
        call inp%upgrade()
        call inp%upgrade()
        hmap%map = inp%lMT%map + inp%lCT%map
        call hmap%smooth(coop_SI_arcmin*beam_fwhm)
        call coop_healpix_maps_ave_udgrade(hmap, lmap)
        call lmap%write(trim(prefix)//"_inp_"//COOP_STR_OF(lmap%nside)//"_"//COOP_STR_OF(nint(beam_fwhm))//"a_"//COOP_STR_OF(i)//".fits")        
     enddo
  else
    ! Cls(2) =  254.90171756233437     
     call map%init(nside = nside, nmaps=1, genre="TEMPERATURE", lmax = lmax)
     call lmap%init(nside = nside_output, nmaps = 1, genre = "TEMPERATURE")
     do i=1, 1000
        map%Cl(:,1) = Cls        
        write(*,*) "simulating map #"//COOP_STR_OF(i)        
        call map%simulate()
        call coop_healpix_maps_ave_udgrade(map, lmap)
        call lmap%write(trim(prefix)//"_full_"//COOP_STR_OF(lmap%nside)//"_"//COOP_STR_OF(nint(beam_fwhm))//"a_"//COOP_STR_OF(i)//".fits")
     enddo
  endif
  call coop_MPI_Finalize()
end program test
