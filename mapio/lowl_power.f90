program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask, m2
  type(coop_healpix_inpaint)::inp
  COOP_INT,parameter ::lmax = 1500
  COOP_INT,parameter ::nrun = 300
  COOP_REAL,parameter::radius_deg = 10.d0
  logical,parameter::force_output = .false.
  COOP_STRING::mask_spot = ""
  COOP_REAL::cls(0:lmax), Cls_sim(0:lmax, nrun), Cls_ave(0:lmax), delta_Cls(0:lmax), sigma
  COOP_INT::l, ell, i, irun, nside
  COOP_REAL::theta, phi
  type(coop_file)::fp
  call coop_MPI_init()
  call coop_random_init()
  !!read in Cl's for fiducial LCDM model
  call fp%open("planck14best_lensedCls.dat", "r") 
  do l=2, lmax
     read(fp%unit, *) i, Cls(l)
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)*coop_gaussian_filter(60.d0, l)**2
  enddo  
  call fp%close()
  Cls(0:1) = 0.d0

  nside = 128
  !call coop_healpix_latitude_cut_mask(nside = nside, latitude_degree = 15.d0, filename = "lowl/lat15_mask_n"//COOP_STR_OF(nside)//".fits")  
  call map%read("lowl/commander_I_n0128_60a.fits")
  call mask%read("lowl/commander_mask_n0128_60a.fits")
 ! call mask%read("planck14/dx11_v2_commander_int_mask_040a_0256.fits")  
 ! call mask%read("lowl/commander_dx11d2_mask_temp_n0256_likelihood_v1.fits")
 ! call mask%read("lowl/lat15_mask_n128.fits")
 ! call mask%read("lowl/mask_hot_bar_n0256.fits")
 ! call mask%read("lowl/mask_tiny.fits")
  call coop_get_command_line_argument(key = "mask",  arg = mask_spot, default = "OUTPUT")
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
     stop
  end select

  !!initialize
  !!Compute pixel-space covariance matrix from fiducial Cl's
  !!For more serious applications you should use a full covariance matrix from FFP9 simulations  
  call inp%init(map, mask, lmax, cls)

  !!if you want to see how the inpainted maps look like
  if(trim(mask_spot).eq."OUTPUT" .or. force_output)then
     call coop_prtsystime(.true.)
     call inp%upgrade(reset = .true., nside_want = inp%base_nside)
     m2 = inp%lMT
     m2%map = inp%lCT%map + inp%lMT%map
     call m2%write("inpainted_map_"//COOP_STR_OF(m2%nside)//".fits")            
     print*, "nside = ", inp%lMT%nside
     call coop_prtsystime()

     do while(inp%lMT%nside .lt. inp%map%nside)
        call inp%upgrade(nside_want = inp%lMT%nside*2)    !!nside -> 16-> 32 -> 64-> 128 -> 256
        m2 = inp%lMT
        m2%map = inp%lCT%map + inp%lMT%map
        call m2%write("inpainted_map_"//COOP_STR_OF(m2%nside)//".fits")            
        print*, "nside = ", inp%lMT%nside
        call coop_prtsystime()
     enddo

     
     where(inp%mask%map(:,1).eq.0.)
        inp%map%map(:,1) = map%bad_data
     end where
     call inp%map%write("original_map.fits")
     stop
  endif
  !!
  Cls_sim = 0.d0
  Cls_ave = 0.d0
  do irun = 1, nrun
     write(*,*) "***** inpainting # ", irun, " ***********"
     call inp%upgrade(reset = .true., nside_want = inp%map%nside)
     inp%lMT%map = inp%lMT%map  + inp%lCT%map !!measured map + inpainted map
     call inp%lMT%map2alm(lmax = lmax)
     Cls_sim(2:lmax, irun) = inp%lMT%Cl(2:lmax, 1)
  enddo
  do l = 2, lmax
     Cls_ave(l) = sum(Cls_sim(l, 1:nrun))/nrun
     delta_Cls(l)  = sqrt(sum((Cls_sim(l, 1:nrun)-Cls_ave(l))**2)/nrun + Cls(l)**2*2./(2.d0*l+1.d0)) !!only compute the diagonal
  enddo

  call fp%open("clsout/clsout_"//trim(mask_spot)//".dat", "w")
  write(fp%unit, "(A8, 3A16)") "# ell ",  "  model C_l  ", " ave Cl  ", "  delta C_l "
  do l=2, lmax/2
     write(fp%unit, "(I8, 3E16.7)") l, Cls(l)*l*(l+1.d0)/coop_2pi, Cls_ave(l)*l*(l+1.d0)/coop_2pi,  delta_Cls(l)*l*(l+1.d0)/coop_2pi
  enddo
  call fp%close()
  call coop_MPI_finalize()  
end program test
