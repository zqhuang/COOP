program test
  use coop_wrapper_utils
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
  type(coop_healpix_maps)::map, mask
  COOP_STRING::mask_file, map_file
  COOP_INT:: isim, nsim, lmax, l, i, lmin_vary
  COOP_REAL,dimension(:),allocatable::S_sim, Cls, sqrtCls
  COOP_REAL::fwhm_deg, S_data, pvalue
  type(coop_file)::fp
  type(coop_healpix_correlation_function)::c_map, c_mask, c_sim
  call coop_MPI_init()
  call coop_random_init()
  if(iargc() .lt. 2)then
     write(*,*) "./LLCORR -map MAP_FILE -mask MASK_FILE -lmax LMAX -nsim NUM_SIMULATIONS -fwhm_deg FWHM_DEG -lmin_vary LMIN_VARY"
     stop
  endif
  call coop_get_command_line_argument(key = "mask", arg = mask_file, default = "NONE")
  call coop_get_command_line_argument(key = "map", arg = map_file)
  call coop_get_command_line_argument(key = "lmax", arg = lmax, default = 60)
  call coop_get_command_line_argument(key = "lmin_vary", arg = lmin_vary, default = 2)  
  call coop_get_command_line_argument(key ="nsim", arg = nsim, default = 1000)
  call coop_get_command_line_argument(key ="fwhm_deg", arg =fwhm_deg, default = 5.d0)  
  allocate(S_sim(nsim), cls(0:lmax), sqrtCls(0:lmax))

  !!read in Cl's for fiducial LCDM model
  Cls(0:1) = 0.d0  
  call fp%open("planck14best_lensedCls.dat", "r") 
  do l=2, lmax
     read(fp%unit, *) i, Cls(l)
     if(i.ne.l) stop "error in Cl file"
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)*coop_gaussian_filter(60.d0, l)**2
  enddo  
  call fp%close()
  call c_sim%init(lmax = lmax, cls = cls)
  
  call map%read(map_file)

  if(trim(mask_file).ne. "NONE")then
     call mask%read(trim(mask_file))
     call mask%map2alm(lmax = lmax)
     call c_mask%init(lmax = lmax, cls = dble(mask%cl(0:lmax,1)))
     call c_mask%set_beam(fwhm = fwhm_deg*coop_SI_degree)    
     call map%apply_mask(mask)
  endif
  
  call map%map2alm(lmax = lmax)
  call c_map%init(lmax = lmax, cls = dble(map%cl(0:lmax,1)))
  call c_map%set_beam(fwhm = fwhm_deg*coop_SI_degree)
  if(lmin_vary .gt.2) &
       cls(2:lmin_vary-1) = map%cl(2:lmin_vary-1,1)
  if(trim(mask_file).ne. "NONE")then
     S_data = c_map%corr2int(-1.d0, 0.5d0, mask_corr = c_mask)
     sqrtCls = sqrt(cls)
     do isim = 1, nsim
        call map%simulate_Tmaps(nside = map%nside, lmax = lmax, sqrtCls = sqrtCls)
        call map%apply_mask(mask)
        call map%map2alm(lmax = lmax)
        call c_map%init(lmax = lmax, cls = dble(map%cl(0:lmax,1)))
        call c_map%set_beam(fwhm = fwhm_deg*coop_SI_degree)
        S_sim(isim) = c_map%corr2int(-1.d0, 0.5d0, mask_corr = c_mask)
        if(mod(isim, 50).eq. 0) write(*,*) "simulation #"//COOP_STR_OF(isim)
     enddo
  else
     S_data = c_map%corr2int(-1.d0, 0.5d0)
     do isim = 1, nsim
        call c_sim%simulate(fwhm = fwhm_deg*coop_SI_degree)
        S_sim(isim) = c_sim%corr2int(-1.d0, 0.5d0)
     enddo
  endif
  
  call coop_quicksort(S_sim)
  write(*,"(A, G14.5, A, G14.5)") "simulation: 95% CL:", S_sim(ceiling(nsim*0.023))," <  S_{1/2} <  ", S_sim(ceiling(nsim*0.977))
  pvalue = count(S_sim .le. S_data)/dble(nsim)
  print*, "data S = ", S_data
  print*, "p value = ", pvalue
  if(nsim .ge. 120) call coop_asy_histogram(x = log(S_sim), nbins = 15, filename = "Shalf/Shalf.txt", xlabel="$\ln S_{1/2}$", ylabel = "Probability")
  call coop_MPI_finalize()

end program test
