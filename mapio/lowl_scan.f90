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
  logical::makefig
  COOP_INT,parameter::nside_scan = 1
  COOP_INT,parameter::npix_scan = nside_scan**2*12
  COOP_INT,parameter ::lmax = 320
  COOP_INT,parameter ::nrun = 250
  COOP_REAL,parameter::radius_deg = sqrt(4.d0/npix_scan)
  logical::loaded = .false.
  COOP_STRING::mask_spot = ""
  COOP_REAL::cls(0:lmax), Cls_ave(0:lmax, 0:npix_scan-1 ), ells(0:lmax)
  COOP_INT::l, ell, i, irun, i_scan
  COOP_REAL::theta, phi, l_deg, b_deg
  type(coop_asy)::fig
  type(coop_file)::fp
  call coop_MPI_init()
  call coop_random_init()
  !!read Cl's for fiducial LCDM model
  call fp%open("planck14best_lensedCls.dat", "r") 
  do l=2, lmax
     ells(l) = l
     read(fp%unit, *) i, Cls(l)
     Cls(l) = Cls(l)*coop_2pi/l/(l+1.d0)*coop_gaussian_filter(60.d0, l)**2
  enddo  
  call fp%close()
  Cls(0:1) = 0.d0
  ells(0) = 0.d0
  ells(1) = 1.d0
  !!read map and mask
  Cls_ave = 0.d0
  call fig%open("clsmasked_"//COOP_STR_OF(nside_scan)//".txt")
  call fig%init(xlabel="$\ell$", ylabel="$\frac{\ell(\ell+1)}{2\pi}C_l$")
  if(coop_file_exists("clmasked_nside_scan"//COOP_STR_OF(nside_scan)//".dat"))then
     write(*,*) "The output file clmasked_nside_scan"//COOP_STR_OF(nside_scan)//".dat already exist. (You can rename or delete it if you really want to redo everything.)"
     call fp%open("clmasked_nside_scan"//COOP_STR_OF(nside_scan)//".dat","ru")
     read(fp%unit) Cls_ave
     call fp%close()
     loaded = .true.
  else
     loaded = .false.
     call map%read("lowl/commander_I_n0128_60a.fits")
     call mask%read("lowl/commander_mask_n0128_60a.fits")
  endif
  call fig%curve(x = ells(2:32), y = Cls(2:32)*ells(2:32)*(ells(2:32)+1.)/coop_2pi, linetype="solid", color="red", legend = "$\Lambda$CDM", linewidth = 2.)
  do i_scan = 0, npix_scan - 1
     if(.not. loaded)then
        call pix2ang_nest(nside_scan, i_scan, theta, phi)
        l_deg = phi/coop_SI_degree
        b_deg = (coop_pio2 - theta)/coop_SI_degree
        write(*,"(A,I5,A,I5,A,I5,A)") "step ", i_scan, "; disc center (l, b) = (", nint(l_deg), ",", nint(b_deg), ")"
        call mask%mask_disc(l_deg = l_deg, b_deg = b_deg, r_deg = radius_deg)
        call inp%init(map, mask, lmax, cls)
        do irun = 1, nrun
           call inp%upgrade(reset = .true., nside_want = inp%map%nside)
           inp%lMT%map = inp%lMT%map  + inp%lCT%map !!measured map + inpainted map
           call inp%lMT%map2alm(lmax = lmax)
           Cls_ave(0:lmax, i_scan) =  Cls_ave(0:lmax, i_scan) + inp%lMT%Cl(0:lmax, 1)
        enddo
        Cls_ave(0:lmax, i_scan)  =   Cls_ave(0:lmax, i_scan) / nrun
     endif
     if(i_scan .eq. 0)then
        call fig%curve(x = ells(2:32), y = Cls_ave(2:32, i_scan)*ells(2:32)*(ells(2:32)+1.)/coop_2pi, linetype="dotted", color="gray", legend = "Masked Data", linewidth = 1.5)
     else
        call fig%curve(x = ells(2:32), y = Cls_ave(2:32, i_scan)*ells(2:32)*(ells(2:32)+1.)/coop_2pi, linetype="dotted", color="gray", linewidth = 1.5)
     endif
  enddo
  call fig%close()
  if(.not. loaded)then
     call fp%open("clmasked_nside_scan"//COOP_STR_OF(nside_scan)//".dat","u")
     write(fp%unit) cls_ave
     call fp%close()
  endif
  call coop_MPI_finalize()  
end program test
