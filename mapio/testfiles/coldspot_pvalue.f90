program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask, bandmap
  COOP_INT, parameter::lmax = 80, nsims = 1000
  COOP_REAL::q1, q2, filter, sigma, Cls(0:lmax), nu, thetaphi(2), l_deg, b_deg, SqrtCls(0:lmax), window(0:lmax)
  COOP_INT::l1, l2, l, i,ic, isim
  type(coop_file)::fp
  type(coop_stacking_options)::sto
  type(coop_healpix_inpaint)::inp  
  logical::hot
  COOP_REAL::rmslow
  call coop_MPI_init()
  call coop_random_init()
  l1 = 2
  l2 = 20
  hot = .false.
  nu = 3.d0
  call fp%open("planck14best_lensedCls.dat", "r") 
  do l=2, lmax
     read(fp%unit, *) i, Cls(l)
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)*coop_gaussian_filter(60.d0, l)**2
  enddo  
  call fp%close()
  Cls(0:1) = 0.d0
  SqrtCls = sqrt(Cls)
  
  call map%read("lowl/commander_I_n0128_60a.fits") !commander_dx11d2_extdata_temp_cmb_n0256_60arc_v1_cr.fits")
  call mask%read("lowl/commander_mask_n0128_60a.fits") !commander_dx11d2_mask_temp_n0256_likelihood_v1.fits")
  call map%rotate_coor(l_deg = 208.82d0, b_deg = 32.98d0)
  call mask%rotate_coor(l_deg = 208.82d0, b_deg = 32.98d0)
  where (mask%map .gt. 0.5)
     mask%map = 1.
  elsewhere
     mask%map = 0.
  end where
  q1 =  (l1+0.5d0)
  q2 =  (l2+0.5d0)
  call map%map2alm(lmax = lmax)
  call map%alm2map()
  call map%savefig("cold_spot_rotated.gif")
  map%alm(0:1, :, :) = 0.  
  bandmap = map
  window(0:1) = 0.d0
  do l = 2, lmax
     window(l) = (exp(-l*(l+1.d0)/q2**2/2.d0) - exp(-l*(l+1.d0)/q1**2/2.d0))
     bandmap%alm(l, 0:l, 1) = bandmap%alm(l, 0:l, 1)* window(l)
     sqrtCls(l) = sqrtCls(l) * window(l)
  enddo
  call bandmap%alm2map()

  call sto%init(domax = hot, peak_name = "T", Orient_name = "RANDOM", nmaps = 1)
  if(hot)then
     sto%I_lower_nu = nu
  else
     sto%I_upper_nu = -nu
  endif
  call bandmap%savefig("data.gif", title="data")
  call bandmap%get_peaks(sto, mask = mask)
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

  do isim = 1, nsims
     call bandmap%simulate_Tmaps(mask%nside, lmax, sqrtCls = sqrtCls)
     call bandmap%get_peaks(sto, mask = mask)
     if(sto%peak_map%n .gt. 0) then
        print*, "----------------------------"
        print*, isim, sto%peak_map%n        
        ic = ic +1
        do i=1, sto%peak_pix%n
           thetaphi  = sto%peak_ang%element(i)
           call coop_healpix_ang2lb(thetaphi(1), thetaphi(2), l_deg, b_deg)
           write(*,*) i, sto%peak_map%element(i)/sto%sigma_I, l_deg, b_deg
        enddo
     endif
  enddo
  write(*,"(F10.1, A)") (ic*100./nsims), "%"
  call coop_MPI_Finalize()
end program test
