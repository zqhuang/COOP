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
  logical::do_constrained, fullcorr
  type(coop_healpix_maps)::map, mask, lmap, lmask
  type(coop_healpix_inpaint)::inp
  COOP_STRING::mask_file
  COOP_INT::lmax = 500
  COOP_INT::nrun, nsim
  COOP_INT,parameter::n = 30
  COOP_INT,parameter::nside_base = 16
  COOP_INT,parameter::npix_base = nside_base**2*12
  COOP_REAL,dimension(:),allocatable::cls
  COOP_REAL,dimension(:),allocatable::S_sim
  COOP_REAL,dimension(:,:),allocatable::Corr_sim
  COOP_REAL::Corr(0:n), weight(0:n), xstep, vec(3, 0:npix_base-1), pvalue, mean_Sinp,  S_inp,  xupper, xlower
  COOP_INT::l, irun, isim, i
  type(coop_file)::fp
  call coop_MPI_init()
  call coop_random_init()
  if(iargc() .lt. 2)then
     write(*,*) "./LLCORR -mask MASK_FILE -fullcorr [T/F] -lmax LMAX -constrained [T/F] -nsim NUM_SIMULATIONS"
     stop
  endif
  call coop_get_command_line_argument(key = "mask", arg = mask_file)
  call coop_get_command_line_argument(key = "fullcorr", arg = fullcorr, default = .false.)
  call coop_get_command_line_argument(key = "lmax", arg = lmax, default = 500)
  call coop_get_command_line_argument(key = "constrained", arg = do_constrained, default = .false.)    
  call coop_get_command_line_argument(key ="nsim", arg = nsim, default = 1000)
  if(do_constrained) write(*,*) "doing constrained realizations"
  allocate(Corr_sim(0:n, nsim), S_sim(nsim), cls(0:lmax))
  if(fullcorr)then
     nrun = 200
  else
     nrun = 20
  endif
  xstep = 1.5d0/(n-1.d0)
  xupper = 0.5d0+xstep/2.d0-1.d-5
  xlower = -1.d0-xstep/2.d0

  !!read in Cl's for fiducial LCDM model
  Cls(0:1) = 0.d0  
  call fp%open("planck14best_lensedCls.dat", "r") 
  do l=2, lmax
     read(fp%unit, *) i, Cls(l)
     if(i.ne.l) stop "error in Cl file"
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)*coop_gaussian_filter(60.d0, l)**2
  enddo  
  call fp%close()

  
  call map%read("lowl/commander_dx11d2_extdata_temp_cmb_n0256_60arc_v1_cr.fits")!"lowl/commander_dx11d2_extdata_temp_cmb_n0256_60arc_v1_cr.fits")
  call mask%read(trim(mask_file))
  call lmap%init(nside = nside_base, nmaps = 1, genre = "T")
  call lmask%init(nside = nside_base, nmaps = 1, genre = "MASK")    
  write(*, "(A, G12.2)") "sky coverage", sum(mask%map)/mask%npix
  call coop_healpix_maps_ave_udgrade(mask, lmask)
  do i=0, npix_base-1
     call pix2vec_nest(nside_base, i, vec(:, i))
  enddo
  call map%apply_mask(mask)
  
  do isim = 1, nsim
     weight = 0.d0
     corr = 0.d0
     if(do_constrained)then
        call lmap%read("mocklowl/sim_inp_16_440a_"//COOP_STR_OF(isim)//".fits", nested = .true.)
     else
        call lmap%read("mocklowl/sim_full_16_440a_"//COOP_STR_OF(isim)//".fits", nested = .true.)
     endif
     call get_corr()
     Corr_sim(:, isim) = Corr
     if(mod(isim, 100).eq.0) print*, "simulation #", isim
  enddo
  do isim = 1, nsim
     S_sim(isim) = sum(corr_sim(:, isim)**2)*xstep
  enddo
  call coop_quicksort(S_sim)
  write(*,"(A, G14.5, A, G14.5)") "simulation: 95% CL:", S_sim(ceiling(nsim*0.023))," <  S_{1/2} <  ", S_sim(ceiling(nsim*0.977))
  pvalue = 0.d0
  mean_sinp = 0.d0
  call inp%init(map, mask, lmax, cls)
  do irun = 1, nrun
     call inp%upgrade(reset = .true.)  !!nside -> 8
     call inp%upgrade()    !!nside ->16
     call inp%upgrade()    !!nside ->32
     call inp%upgrade()    !!nside ->64
     inp%lMT%map = inp%lCT%map + inp%lMT%map
     call inp%lMT%smooth(coop_SI_arcmin*440.d0)
     call coop_healpix_maps_ave_udgrade(inp%lMT, lmap)
     call get_corr()
     S_inp = sum(corr**2)*xstep
     Mean_sinp = Mean_sinp + S_inp
     if(mod(irun, 20).eq.0)write(*,*) "inpainting #", irun, " S = ", S_inp
     pvalue = pvalue + count(S_sim .le. S_inp)/dble(nsim)
  enddo
  pvalue = pvalue/nrun
  mean_Sinp = Mean_Sinp/nrun
  call fp%open("Shalf/Shalf.log", "a")
  write(fp%unit, "(A, 2I8, 2E16.7)") trim(mask_file), lmax, n, pvalue, mean_Sinp
  call fp%close()
  print*, "data <S> = ", mean_Sinp
  print*, "<p value> = ", pvalue
  if(nsim .ge. 120) call coop_asy_histogram(x = log(S_sim), nbins = 15, filename = "Shalf/Shalf.txt", xlabel="$\ln S_{1/2}$", ylabel = "Probability")
  call coop_MPI_finalize()

contains

  subroutine get_corr()
    COOP_REAL::x, rpos
    COOP_INT::i, j, pos
    corr = 0.d0
    weight = 0.d0
    if(fullcorr)then
       do i=0, lmap%npix-1
          do j = 0, i-1
             x = dot_product(vec(:, i), vec(:, j))
             if(x.lt.xupper)then
                rpos = (x-xlower)/xstep
                pos = floor(rpos)
                rpos = rpos - pos
                corr(pos) = corr(pos)+lmap%map(i, 1)*lmap%map(j, 1)*(1.d0-rpos)
                weight(pos) = weight(pos)+1.d0-rpos
                corr(pos+1) = corr(pos+1)+lmap%map(i, 1)*lmap%map(j, 1)*rpos
                weight(pos+1) = weight(pos+1)+rpos                 
             endif
          enddo
       enddo
    else
       do i=0, lmap%npix-1
          if(lmask%map(i, 1) .lt.0.5)cycle
          do j = 0, i-1
             if(lmask%map(j, 1) .lt.0.5)cycle              
             x = dot_product(vec(:, i), vec(:, j))
             if(x.lt.xupper)then
                rpos = (x-xlower)/xstep
                pos = floor(rpos)
                rpos = rpos - pos
                corr(pos) = corr(pos)+lmap%map(i, 1)*lmap%map(j, 1)*(1.d0-rpos)
                weight(pos) = weight(pos)+1.d0-rpos
                corr(pos+1) = corr(pos+1)+lmap%map(i, 1)*lmap%map(j, 1)*rpos
                weight(pos+1) = weight(pos+1)+rpos                 
             endif
          enddo
       enddo
    endif
    where (weight .ne. 0.d0)
       Corr = Corr/weight
    elsewhere
       corr = 0.d0
    end where
    Corr(0) = 0.d0
    Corr(n) = 0.d0  !!zero out the boundary 
  end subroutine get_corr
end program test
