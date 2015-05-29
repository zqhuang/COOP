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
  type(coop_healpix_maps)::map, mask, m2
  type(coop_healpix_inpaint)::inp
  COOP_INT,parameter ::lmax = 32
  COOP_INT,parameter ::nrun = 100
  COOP_REAL,parameter::radius_deg = 10.d0
  COOP_INT,parameter::n = 80  
  COOP_STRING::mask_spot = ""
  COOP_REAL::cls(0:lmax), Shalf, thisShalf
  COOP_REAL::Corr(0:n), weight(0:n), xstep, x, vec(3, 12*16**2-1)
  COOP_INT::l, i, irun, j, pos, iCount
  type(coop_file)::fp
  call coop_MPI_init()
  call coop_random_init()
  xstep = 2.d0/n
  
  !!read in Cl's for fiducial LCDM model
  call fp%open("planck14best_lensedCls.dat", "r") 
  do l=2, lmax
     read(fp%unit, *) i, Cls(l)
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)*coop_gaussian_filter(60.d0, l)**2
  enddo  
  call fp%close()
  Cls(0:1) = 0.d0
  call map%read("lowl/commander_dx11d2_extdata_temp_cmb_n0256_60arc_v1_cr.fits")
  !call mask%read("lowl/commander_dx11d2_mask_temp_n0256_likelihood_v1.fits")
  call mask%read("planck14/dx11_v2_commander_int_mask_040a_0256.fits")  


100 write(*,*) "enter the mask spot [NEP, SEP, NCP, SCP, NGP, SGP, COLDSPOT, NONE, NASYM, SASYM]"
  read(*,*) mask_spot
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
  case("NONE", "OUTPUT")
     !do nothing
  case default
     write(*,*) trim(mask_spot)//": unknown mask option"
     goto 100
  end select
  do i=0, 12*16**2-1
     call pix2vec_nest(16, i, vec(:, i))
  enddo
  !!compute pseudo Cl's
  map%map = map%map*mask%map  
  !!initialize
  !!Compute pixel-space covariance matrix from fiducial Cl's
  !!For more serious applications you should use a full covariance matrix from FFP9 simulations
  Corr = 0.d0
  weight = 0.d0
  call inp%init(map, mask, lmax, cls)
  do irun = 1, nrun
     write(*,*) "***** inpainting # ", irun, " ***********"
     call inp%upgrade(reset = .true.)  !!nside -> 8
     call inp%upgrade()    !!nside ->16
     inp%lMT%map = inp%lCT%map + inp%lMT%map
     do i=0, inp%lMT%npix-1
        do j = 0, i-1
           x = dot_product(vec(:, i), vec(:, j))
           if(x.lt.0.5d0)then
              pos = nint((1.d0+x)/xstep)
              corr(pos) = corr(pos)+inp%lMT%map(i, 1)*inp%lMT%map(j, 1)
              weight(pos) = weight(pos)+1.d0
           endif
        enddo
     enddo
  enddo
  where (weight .ne. 0.d0)
     Corr = Corr/weight
  elsewhere
     corr = 0.d0
  end where
  Shalf = sum(corr**2)*xstep
  print*, "S_{1/2} =", Shalf
  print*, "now estimating the p value"
  do i=0, 12*16**2-1
     call pix2vec_ring(16, i, vec(:, i))
  enddo
  call inp%lMT%convert2ring()
  do irun = 1, 1000
     weight = 0.d0
     corr = 0.d0
     call inp%lMT%simulate_Tmaps(inp%lMT%nside, inp%lmax, inp%sqrtCls)
     do i=0, inp%lMT%npix-1
        do j = 0, i-1
           x = dot_product(vec(:, i), vec(:, j))
           if(x.lt.0.5d0)then
              pos = nint((1.d0+x)/xstep)
              corr(pos) = corr(pos)+inp%lMT%map(i, 1)*inp%lMT%map(j, 1)
              weight(pos) = weight(pos)+1.d0
           endif
        enddo
     enddo
     where (weight .ne. 0.d0)
        Corr = Corr/weight
     elsewhere
        corr = 0.d0
     end where
     thisShalf =  sum(corr**2)*xstep
     if(thisShalf .le. Shalf)then
        icount = icount + 1
        print* , "***** simulation #", irun, " ***** S_{1/2} = ", thisShalf
     endif
  enddo
  print*, "p value = ", icount/1000.d0
  call coop_MPI_finalize()  
end program test
