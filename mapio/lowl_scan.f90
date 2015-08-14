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
  COOP_INT,parameter::nside_scan = 2
  COOP_INT,parameter::npix_scan = nside_scan**2*12
  COOP_INT,parameter ::lmax = 120
  COOP_INT,parameter ::nrun = 200
  logical::loaded = .false.
  COOP_REAL::cls(0:lmax), Cls_ave(0:lmax, 0:npix_scan-1 ), ells(0:lmax)
  COOP_INT::l, ell, i, irun, i_scan
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
  endif
  call fig%curve(x = ells(2:32), y = Cls(2:32)*ells(2:32)*(ells(2:32)+1.)/coop_2pi, linetype="solid", color="red", legend = "$\Lambda$CDM", linewidth = 2.)
  do i_scan = 0, npix_scan - 1
     if(.not. loaded)then
        write(*,"(A,I5,A,I5)") "step ", i_scan, " / ", npix_scan
        call map%read("lowl/commander_I_n0128_60a.fits")        
        call mask%read("lowl/commander_mask_n0128_60a.fits")
        call mask%mask_pixel(nside_scan, i_scan)
        if(nside_scan .le.2) call mask%write("scaned_mask_"//COOP_STR_OF(i_scan)//".fits")
        call inp%init(map, mask, lmax, cls)
        do irun = 1, nrun
           if(mod(irun, nrun/10).eq.0) write(*, "(A$)")"."
           call inp%upgrade(reset = .true., nside_want = inp%map%nside)
           inp%lMT%map = inp%lMT%map  + inp%lCT%map !!measured map + inpainted map
           call inp%lMT%map2alm(lmax = lmax)
           Cls_ave(0:lmax, i_scan) =  Cls_ave(0:lmax, i_scan) + inp%lMT%Cl(0:lmax, 1)
        enddo
        write(*,*) 
        Cls_ave(0:lmax, i_scan)  =   Cls_ave(0:lmax, i_scan) / nrun
     endif
     if(i_scan .eq. 0)then
        call fig%curve(x = ells(2:32), y = Cls_ave(2:32, i_scan)*ells(2:32)*(ells(2:32)+1.)/coop_2pi, linetype="dotted", color="gray", legend = "Masked Data", linewidth = 1.)
     else
        call fig%curve(x = ells(2:32), y = Cls_ave(2:32, i_scan)*ells(2:32)*(ells(2:32)+1.)/coop_2pi, linetype="dotted", color="gray", linewidth = 1.)
     endif
  enddo
  call fig%legend(0.3, 0.2)
  call fig%close()
  if(.not. loaded)then
     call fp%open("clmasked_nside_scan"//COOP_STR_OF(nside_scan)//".dat","u")
     write(fp%unit) cls_ave
     call fp%close()
  endif
  call coop_MPI_finalize()  
end program test
