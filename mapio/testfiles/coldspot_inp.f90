program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask, lowmask, lowmap, invlowmask, bandmap
  COOP_INT, parameter::lsmooth = 30
  COOP_INT, parameter::lmax = lsmooth*3, nsims = 200
  COOP_REAL, parameter::nucut = 1.d0
  COOP_REAL::q1, q2, filter, sigma, Cls(0:lmax), nu, thetaphi(2), l_deg, b_deg, SqrtCls(0:lmax), window(0:lmax), fwhm
  COOP_INT::l1, l2, l, i,ic, isim
  type(coop_file)::fp
  type(coop_stacking_options)::sto
  type(coop_healpix_inpaint)::inp
  logical::hot
  COOP_REAL::rmslow
  call coop_MPI_init()
  call coop_random_init()
  l1 = 8
  l2 = 20
  hot = .false.
  nu = 3.d0

  fwhm = 1.d0/(lsmooth+0.5d0)/coop_sigma_by_fwhm
  
  call fp%open("planck14best_lensedCls.dat", "r") 
  do l=2, lmax
     read(fp%unit, *) i, Cls(l)
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)*coop_gaussian_filter(60.d0, l)**2*exp(-l*(l+1.d0)/(lsmooth+0.5)**2)
  enddo  
  call fp%close()
  Cls(0:1) = 0.d0
  SqrtCls = sqrt(Cls)
  
  call map%read("lowl/commander_I_n0128_60a.fits") !commander_dx11d2_extdata_temp_cmb_n0256_60arc_v1_cr.fits")
  call mask%read("lowl/commander_mask_n0128_60a.fits") !commander_dx11d2_mask_temp_n0256_likelihood_v1.fits")
  call map%rotate_coor(l_deg = 208.82d0, b_deg = 32.98d0)
  call map%smooth(fwhm = fwhm)  
  call mask%rotate_coor(l_deg = 208.82d0, b_deg = 32.98d0)
  where (mask%map .gt. 0.5)
     mask%map = 1.
  elsewhere
     mask%map = 0.
  end where
  
  q1 =  (l1+0.5d0)
  q2 =  (l2+0.5d0)
  call map%map2alm(lmax = lmax)
  map%alm(0:1, :, :) = 0.  
  lowmap = map
  bandmap = map
  window(0:1) = 0.d0
  do l = 2, lmax
     window(l) = (exp(-l*(l+1.d0)/q2**2/2.d0) - exp(-l*(l+1.d0)/q1**2/2.d0))
     bandmap%alm(l, 0:l, 1) = bandmap%alm(l, 0:l, 1)* window(l)
     lowmap%alm(l, 0:l, 1) = lowmap%alm(l, 0:l, 1)*exp(-l*(l+1.d0)/q1**2/2.d0)
     sqrtCls(l) = sqrtCls(l) * window(l)
  enddo
  call bandmap%alm2map()
  call lowmap%alm2map()

  lowmask = mask
  rmslow = sqrt(sum(lowmap%map**2)/lowmap%npix)
  lowmask%map = 0.   
  where (abs(lowmap%map) .gt. rmslow*nucut .and. mask%map .gt. 0.5)
     lowmask%map = 1.
  end where
  call lowmask%write("coldspot_constrained_mask.fits")
!!$  call lowmask%savefig("coldspot_constrained_mask.gif")
!!$  stop
  invlowmask = lowmask
  invlowmask%map = (1. - lowmask%map)*mask%map

!!$  call lowmask%mask_disc(l_deg = 0.d0, b_deg = 0.d0, r_deg = 5.d0, fillvalue = 0.5d0)
!!$  call lowmask%savefig("coldspot_constrained_mask.gif", title = "green: cold spot; red: constrained region")
!!$  call map%apply_mask(mask, bad_data = .true.)
!!$  call map%savefig("coldspot_rotated_l20smooth.gif", title="cold spot rotated tothe center; l_smooth = 20")
!!$  call lowmap%apply_mask(mask, bad_data = .true.)
!!$  call lowmap%savefig("coldspot_rotated_l6smooth.gif", title="cold spot rotated to the center; l_smooth = 6")
!!$  call bandmap%apply_mask(mask, bad_data = .true.)
!!$  call bandmap%savefig("coldspot_rotated_diff_l6_l20.gif", title = "cold spot rotated to the center; diff. between l_smooth = 6 and l_smooth = 20")
!!$  stop  
!!$  call invlowmask%savefig("invlm.gif")

  call sto%init(domax = hot, peak_name = "T", Orient_name = "RANDOM", nmaps = 1)
  if(hot)then
     sto%I_lower_nu = nu
  else
     sto%I_upper_nu = -nu
  endif
  call bandmap%savefig("data.gif", title="data")
  call bandmap%get_peaks(sto, mask = invlowmask)
  do i=1, sto%peak_pix%n
     thetaphi  = sto%peak_ang%element(i)
     call coop_healpix_ang2lb(thetaphi(1), thetaphi(2), l_deg, b_deg)
     write(*,*) i, sto%peak_map%element(i)/sto%sigma_I, l_deg, b_deg
     if(maxval(abs(sto%peak_map%element(i)/sto%sigma_I)) .gt. nu)nu=maxval(abs(sto%peak_map%element(i)/sto%sigma_I))
  enddo
  write(*,*) "nu =", nu
  if(hot)then
     sto%I_lower_nu = nu
  else
     sto%I_upper_nu = -nu
  endif
  ic = 0
  call inp%init(map, lowmask, lmax, cls)
  call lowmap%convert2nested()

  lowmap%map = 0.
  do isim = 1, nsims
     call inp%upgrade(reset = .true., nside_want = mask%nside)
     bandmap%map = inp%lMT%map  + inp%lCT%map !!measured map + inpainted map
     bandmap%ordering = inp%lMT%ordering
     lowmap%map = lowmap%map+bandmap%map
     call bandmap%smooth_with_window(fwhm = 0.d0,lmax = lmax,window = window)     
     call bandmap%get_peaks(sto, mask = invlowmask)
!!$     call bandmap%savefig("sim"//COOP_STR_OF(isim)//".gif", title="simulation #"//COOP_STR_OF(isim))
     print*, isim, sto%peak_map%n             
     if(sto%peak_map%n .gt. 0) then
        ic = ic +1
        do i=1, sto%peak_pix%n
           thetaphi  = sto%peak_ang%element(i)
           call coop_healpix_ang2lb(thetaphi(1), thetaphi(2), l_deg, b_deg)
           write(*,*) i, sto%peak_map%element(i)/sto%sigma_I, l_deg, b_deg
        enddo
     endif
  enddo
  lowmap%map = lowmap%map/nsims
  call lowmap%write("coldspot_constrained_meanmap.fits")
  call lowmap%savefig("coldspot_constrained_meanmap.gif", title="constrained mean map")
  lowmask%map = 1.-lowmask%map
  call lowmap%apply_mask(lowmask)
  call lowmap%savefig("coldspot_constrained_meanmap_chopped.gif")
  write(*,"(F10.1, A)") (ic*100./nsims), "%"
  call coop_MPI_Finalize()
end program test
