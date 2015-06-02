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
  type(coop_healpix_maps)::map, mask, m2, lmap, lmask
  type(coop_healpix_inpaint)::inp
  COOP_INT,parameter ::lmax = 10
  COOP_INT::nrun
  COOP_INT,parameter ::nsim = 5000  
  COOP_REAL,parameter::radius_deg = 10.d0
  COOP_INT,parameter::n = 250
  COOP_INT,parameter::nside_base = 16
  COOP_INT,parameter::npix_base = nside_base**2*12
  logical::masked_S = .false.
  logical::remove_mean = .false.
  logical::use_commander_mask = .false.
  COOP_STRING::mask_spot = ""
  COOP_REAL,dimension(:),allocatable::cls, sqrtcls
  COOP_REAL::S_sim(nsim), S_inp
  COOP_REAL::Corr(0:n), weight(0:n), xstep, x, vec(3, 0:npix_base-1), Corr_sim(0:n, nsim), Corr_data(0:n), Corr_mean(0:n), pvalue, S_mean, rpos, mean_Sinp
  COOP_INT::l, i, irun, j, pos, iCount, isim, nside_sim
  type(coop_file)::fp
  call coop_MPI_init()
  call coop_random_init()
  if(iargc() .lt. 4)then
     write(*,*) "./LLCORR masked_sky subtract_mean use_commander_mask lmax"
     stop
  endif
  call coop_get_Input(1, masked_S)
  call coop_get_Input(2, remove_mean)
  call coop_get_Input(3, use_commander_mask)
  call coop_get_Input(4, lmax)
  allocate(cls(0:lmax), sqrtcls(0:lmax))
  if(masked_S)then
     nrun = 100
  else
     if(use_commander_mask)then
        nrun = 300        
     else
        nrun = 1000
     endif
  endif
  xstep = 2.d0/n
  call lmap%init(nside = nside_base, nmaps = 1, genre = "T")
  call lmask%init(nside = nside_base, nmaps = 1, genre = "MASK")  
  !!read in Cl's for fiducial LCDM model
  call fp%open("planck14best_lensedCls.dat", "r") 
  do l=2, lmax
     read(fp%unit, *) i, Cls(l)
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)*coop_gaussian_filter(60.d0, l)**2
  enddo  
  call fp%close()
  Cls(0:1) = 0.d0
  sqrtcls = sqrt(cls)
  call map%read("lowl/commander_dx11d2_extdata_temp_cmb_n0256_60arc_v1_cr.fits")
  if(use_commander_mask)then
     call mask%read("lowl/commander_dx11d2_mask_temp_n0256_likelihood_v1.fits")
  else
     call mask%read("planck14/dx11_v2_commander_int_mask_040a_0256.fits")     
  endif
  print*,"sky coverage", sum(mask%map)/mask%npix
  call coop_healpix_maps_ave_udgrade(mask, lmask)


!!$100 write(*,*) "enter the mask spot [NEP, SEP, NCP, SCP, NGP, SGP, COLDSPOT, NONE, NASYM, SASYM]"
  mask_spot = "NONE"
  mask_spot = adjustl(mask_spot)
  !!mask out a disk
  select case(trim(mask_spot))
  case("COLDSPOT")
     call mask%mask_disc(l_deg = 207.8d0, b_deg = -56.3d0, r_deg = radius_deg)
  case("NGP")
     call mask%mask_disc(l_deg = 0.d0, b_deg = 90.d0, r_deg = radius_deg)
  case("SGP")
     call mask%mask_disc(l_deg = 0.d0, b_deg = -90.d0, r_deg = radius_deg)
  case("NEP") 
     call mask%mask_disc(l_deg = 98.d0, b_deg = 31.d0, r_deg = radius_deg)
  case("SEP")
     call mask%mask_disc(l_deg = 278.d0, b_deg = -31.d0, r_deg = radius_deg)
  case("NCP")
     call mask%mask_disc(l_deg = 123.d0, b_deg = 28.d0, r_deg = radius_deg)     
  case("SCP")
     call mask%mask_disc(l_deg = 303.d0, b_deg = -28.d0, r_deg = radius_deg)
  case("NASYM")
     call mask%mask_disc(l_deg = 212.d0, b_deg = -13.d0, r_deg = 90.d0)
  case("SASYM")
     call mask%mask_disc(l_deg = 32.d0, b_deg = 13.d0, r_deg = 90.d0)
  case("NONE")
     !!do nothing
  case default
     write(*,*) trim(mask_spot)//": unknown mask option"
!!$     goto 100
  end select
  do i=0, npix_base-1
     call pix2vec_nest(nside_base, i, vec(:, i))
  enddo
  !!compute pseudo Cl's
  map%map = map%map*mask%map
  call inp%init(map, mask, lmax, cls)
  nside_sim = nside_base*4
  do while(nside_sim .lt. lmax .and. nside_sim .lt. map%nside)
     nside_sim  = nside_sim * 2
  enddo
  call m2%init(nside = nside_sim, nmaps = 1, genre = "T")
  do isim = 1, nsim
     weight = 0.d0
     corr = 0.d0
     call m2%simulate_Tmaps(map%nside, lmax, sqrtCls)
     call coop_healpix_maps_ave_udgrade(m2, lmap)
     if(masked_S)then
        do i=0, lmap%npix-1
           if(lmask%map(i, 1) .lt.0.5)cycle
           do j = 0, i-1
              if(lmask%map(j, 1) .lt.0.5)cycle              
              x = dot_product(vec(:, i), vec(:, j))
              if(x.lt.0.5d0)then
                 rpos = (1.d0+x)/xstep
                 pos = min(max(floor(rpos), 0), n-1)
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
           do j = 0, i-1
              x = dot_product(vec(:, i), vec(:, j))
              if(x.lt.0.5d0)then
                 rpos = (1.d0+x)/xstep
                 pos = min(max(floor(rpos), 0), n-1)
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
     Corr_sim(:, isim) = Corr
     if(mod(isim, 100).eq.0)print*, "simulation #", isim
  enddo
  if(remove_mean)then
     write(*,*) "Mean Correlation: "
     do i= 0, n
        corr_mean(i) = sum(corr_sim(i,:))/nsim  !!this gives the LCDM mean
        write(*,*) i, corr_mean(i)
     enddo
  else
     corr_mean = 0.d0
  endif
  do isim = 1, nsim
     S_sim(isim) = (sum((corr_sim(:, isim)-corr_mean)**2) - (corr_sim(0, isim)-corr_mean(0))**2/2.d0)*xstep
  enddo
  S_mean = sum(S_sim)/nsim
  print*, "simulation: S = ", S_mean, " +/- ", sqrt(sum((S_sim -S_mean)**2)/nsim)
  pvalue = 0.d0
  mean_sinp = 0.d0
  do irun = 1, nrun
     Corr = 0.d0
     weight = 0.d0

     call inp%upgrade(reset = .true.)  !!nside -> 8
     call inp%upgrade()    !!nside ->16
     inp%lMT%map = inp%lCT%map + inp%lMT%map
     if(masked_S)then
        do i=0, inp%lMT%npix-1
           if(lmask%map(i,1).lt.0.5)cycle           
           do j = 0, i-1
              if(lmask%map(j,1).lt.0.5)cycle
              x = dot_product(vec(:, i), vec(:, j))
              if(x.lt.0.5d0)then
                 rpos = (1.d0+x)/xstep
                 pos = min(max(floor(rpos), 0), n-1)
                 rpos = rpos - pos
                 corr(pos) = corr(pos)+inp%lMT%map(i, 1)*inp%lMT%map(j, 1)*(1.d0-rpos)
                 weight(pos) = weight(pos)+1.d0-rpos
                 corr(pos+1) = corr(pos+1)+inp%lMT%map(i, 1)*inp%lMT%map(j, 1)*rpos
                 weight(pos+1) = weight(pos+1)+rpos                 
              endif
           enddo
        enddo
     else
        do i=0, inp%lMT%npix-1
           do j = 0, i-1
              x = dot_product(vec(:, i), vec(:, j))
              if(x.lt.0.5d0)then
                 rpos = (1.d0+x)/xstep
                 pos = min(max(floor(rpos), 0), n-1)
                 rpos = rpos - pos
                 corr(pos) = corr(pos)+inp%lMT%map(i, 1)*inp%lMT%map(j, 1)*(1.d0-rpos)
                 weight(pos) = weight(pos)+1.d0-rpos
                 corr(pos+1) = corr(pos+1)+inp%lMT%map(i, 1)*inp%lMT%map(j, 1)*rpos
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
     S_inp = (sum((corr-corr_mean)**2) - (corr(0)-corr_mean(0))**2/2.d0)*xstep
     Mean_sinp = Mean_sinp + S_inp
     if(mod(irun, 20).eq.0)write(*,*) "inpainting #", irun, " S = ", S_inp
     pvalue = pvalue + count(S_sim .gt. S_inp)/dble(nsim)
  enddo
  pvalue = pvalue/nrun
  mean_Sinp = Mean_Sinp/nrun
  call fp%open("Shalf/S_"//COOP_STR_OF(masked_S)//COOP_STR_OF(remove_mean)//COOP_STR_OF(use_commander_mask)//".log", "a")
  write(fp%unit, "(2I8, 2E16.7)") lmax, n, pvalue, mean_Sinp
  call fp%close()
  print*, "p value = ", pvalue
  call coop_asy_histogram(x = log(S_sim), nbins = 20, filename = "Shalf/S_"//COOP_STR_OF(masked_S)//COOP_STR_OF(remove_mean)//COOP_STR_OF(use_commander_mask)//"_"//COOP_STR_OF(lmax)//"_"//COOP_STR_OF(n)//".txt", xlabel="$\ln S_{1/2}$", ylabel = "Probability")
  call coop_asy_histogram(x = log(S_sim(1:nsim/20)), nbins = 20, filename = "Shalf/S_"//COOP_STR_OF(masked_S)//COOP_STR_OF(remove_mean)//COOP_STR_OF(use_commander_mask)//"_"//COOP_STR_OF(lmax)//"_"//COOP_STR_OF(n)//".txt", xlabel="$\ln S_{1/2}$", ylabel = "Probability")  
  
  call coop_MPI_finalize()  
end program test
