program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask, mapcopy, smoothmask
  COOP_INT,parameter::n = 3000, m = 128
  COOP_SINGLE, parameter::r = 5.d0*coop_SI_degree
  COOP_SINGLE::rmscut
  COOP_STRING::dustorcmb , filename
  COOP_INT::resolution
  type(coop_list_realarr)::palist
  type(coop_asy)::fig
  COOP_INT::i, pix
  COOP_SINGLE::pa(2)
  COOP_REAL::p(m), a(m)
  call coop_MPI_init()
  call coop_random_init()
  call coop_get_command_line_argument(key = "r", arg = resolution)
  call coop_get_command_line_argument(key = "map", arg = dustorcmb)
  call coop_get_command_line_argument(key = "nu", arg = rmscut)
  call system("rm -f contours/*.*")  
  call map%read("dust/"//trim(dustorcmb)//"_i_n1024_"//COOP_STR_OF(resolution)//"a.fits")
  call mask%read("dust/lat30_mask_n1024.fits")
  call smoothmask%read("dust/lat25_mask_n1024_smooth.fits")  
  call palist%init()
  do i=1, n
     pix = coop_random_index(map%npix)-1
     do while(mask%map(pix, 1) .lt. 0.5)
        pix = coop_random_index(map%npix)-1
     enddo
     call map%filament_perimeter_area_list(pix, r, palist, rmscut = rmscut, plot = "contours/contours_"//COOP_STR_OF(i)//".txt")        
  enddo
  filename = trim(dustorcmb)//"_pa_"//COOP_STR_OF(resolution)//"a_nu"//COOP_FILESTR_OF(dble(rmscut))//".txt"
  call fig%open(filename)
  call fig%init(xlabel = "perimeter (deg)", ylabel = "area (deg$^2$)", caption="black dots: original map; red dots: Gaussian simulation")
  call coop_asy_dots(fig, palist, xunit = real(coop_SI_degree), yunit = real(coop_SI_degree**2))

  map%map = map%map*smoothmask%map
  mapcopy = map
  call map%map2alm(lmax = min(ceiling(2./max(resolution*coop_SI_arcmin*coop_sigma_by_fwhm, 1.d-6)), floor(map%nside*coop_healpix_lmax_by_nside), coop_healpix_default_lmax))
  map%cl = map%cl/(sum(smoothmask%map)/smoothmask%npix)
  call map%simulate()
  call palist%init()
  do i=1, n
     pix = coop_random_index(map%npix)-1
     do while(mask%map(pix, 1) .lt. 0.5)
        pix = coop_random_index(map%npix)-1
     enddo
     call map%filament_perimeter_area_list(pix, r, palist, rmscut = rmscut, plot = "contours/contours_"//COOP_STR_OF(i)//".txt")             
  enddo
  call coop_asy_dots(fig, palist, xunit = real(coop_SI_degree), yunit = real(coop_SI_degree**2), color = "red")
  
  call coop_set_uniform(m, p, 0.d0, dble(fig%xmax))
  a = p**2/coop_4pi
 
  call coop_asy_curve(fig, p, a, color = "blue", linewidth = 1.5, legend = "$A = \frac{p^2}{4\pi}$")
  a = p*(resolution/60./2.)
  call coop_asy_curve(fig, p, a, color = "green", linewidth = 1.5, legend = "$A = {\rm FWHM} \times \frac{p}{2}$")

  a = p*(1./2.)
  call coop_asy_curve(fig, p, a, color = "violet", linewidth = 1.5, legend = "$A = 1^\circ \times \frac{p}{2}$" )  
  call fig%legend(0.1, 0.88)

  call fig%label("FWHM = $"//COOP_STR_OF(resolution)//"'$", 0.1, 0.48)
  call fig%label("$\nu = "//COOP_STR_OF(rmscut)//"$", 0.1, 0.4)  
  call fig%close()  

  call system("../utils/fasy.sh "//trim(filename))  

  call coop_MPI_finalize()
end program test
