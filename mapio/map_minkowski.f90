program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask, mapcopy, smoothmask, map_count, map_p, map_a
  COOP_INT,parameter::nside = 16
  COOP_SINGLE, parameter::r = 5.d0*coop_SI_degree, Imax = 200.
  COOP_SINGLE::rmscut
  COOP_STRING::dustorcmb, prefix
  COOP_INT::resolution
  type(coop_list_realarr)::palist
  type(coop_asy)::fig
  COOP_INT::i, pix
  COOP_SINGLE::sum_pa(2)
  COOP_REAL::theta, phi
  call coop_MPI_init()
  call coop_random_init()
  call coop_get_command_line_argument(key = "r", arg = resolution)
  call coop_get_command_line_argument(key = "map", arg = dustorcmb)
  call coop_get_command_line_argument(key = "nu", arg = rmscut)
  call system("rm -f contours/*.*")  
  call map%read("dust/"//trim(dustorcmb)//"_i_n1024_"//COOP_STR_OF(resolution)//"a.fits")
  
  !!regularize point sources  
  where(map%map.gt. Imax)  
     map%map = Imax + log((map%map - Imax)*10./Imax+1.)*Imax/10.
  end where
  where(map%map .lt. -Imax)
     map%map = -Imax - log((-Imax - map%map)*10./Imax+1.)*Imax/10.
  end where
  call mask%read("dust/lat30_mask_n1024.fits")
  call smoothmask%read("dust/lat25_mask_n1024_smooth.fits")
  call map_count%init(nside = nside, nmaps = 1, genre = "T")
  map_p = map_count
  map_a = map_count
  prefix = "mm_output/"//trim(dustorcmb)//"_"//COOP_STR_OF(resolution)//"a_nu"//COOP_FILESTR_OF(dble(rmscut))
  do i= 0, map_count%npix-1
     call map_count%pix2ang(i, theta, phi)
     call map%ang2pix(theta, phi, pix)
     if(mask%map(pix, 1) .lt. 0.5)then
        map_count%map(i, 1) = map_count%bad_data
        map_a%map(i, 1) = map_a%bad_data
        map_p%map(i, 1) = map_p%bad_data
     else
        call palist%init()        
        call map%perimeter_area_list(pix, r, palist, rmscut = rmscut, sum_pa = sum_pa)
        map_count%map(i, 1) = palist%n
        map_a%map(i, 1) = sum_pa(2)/coop_SI_degree**2
        map_p%map(i, 1) = sum_pa(1)/coop_SI_degree
     endif
  enddo
  call map_count%write(trim(prefix)//"_count.fits")
  call map_p%write(trim(prefix)//"_perimeter.fits")
  call map_a%write(trim(prefix)//"_area.fits")
  where(map_p%map .ge. 0.)
     map_a%map = map_a%map/map_p%map
  end where
  call map_a%write(trim(prefix)//"_abyp.fits")       
  
  mapcopy = map  
  where(mask%map .lt. 0.5)
     mapcopy%map = mapcopy%bad_data
  end where
  call mapcopy%write(trim(prefix)//".fits")

  
  map%map = map%map*smoothmask%map

  call map%map2alm(lmax = min(ceiling(2./max(resolution*coop_SI_arcmin*coop_sigma_by_fwhm, 1.d-6)), floor(map%nside*coop_healpix_lmax_by_nside), coop_healpix_default_lmax))
  map%cl = map%cl/(sum(smoothmask%map)/smoothmask%npix)
  call map%simulate()
  
  do i= 0, map_count%npix-1
     call map_count%pix2ang(i, theta, phi)
     call map%ang2pix(theta, phi, pix)
     if(mask%map(pix, 1) .lt. 0.5)then
        map_count%map(i, 1) = map_count%bad_data
        map_a%map(i, 1) = map_a%bad_data
        map_p%map(i, 1) = map_p%bad_data
     else
        call palist%init()        
        call map%perimeter_area_list(pix, r, palist, rmscut = rmscut, sum_pa = sum_pa)
        map_count%map(i, 1) = palist%n
        map_a%map(i, 1) = sum_pa(2)/coop_SI_degree**2
        map_p%map(i, 1) = sum_pa(1)/coop_SI_degree
     endif
  enddo


  call map_count%write(trim(prefix)//"_Gaussianized_count.fits")
  call map_p%write(trim(prefix)//"_Gaussianized_perimeter.fits")
  call map_a%write(trim(prefix)//"_Gaussianized_area.fits")  
  where(map_p%map .ge. 0.)
     map_a%map = map_a%map/map_p%map
  end where
  call map_a%write(trim(prefix)//"_Gaussianized_abyp.fits")       
  where(mask%map .lt. 0.5)
     map%map = map%bad_data
  end where
  call map%write(trim(prefix)//"_Gaussianized.fits") 


  call coop_MPI_finalize()
end program test
