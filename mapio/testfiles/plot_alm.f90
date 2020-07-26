program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask, m2
  type(coop_healpix_inpaint)::inp
  COOP_INT,parameter ::lmax = 40
  COOP_INT,parameter ::lmin = 2
  COOP_INT,parameter ::nrun = 1
  complex::alms(0:lmax, 0:lmax), alms_sim(0:lmax, 0:lmax, nrun)
  COOP_REAL::alms_plot(lmin:lmax, -lmax:lmax)
  COOP_REAL::cls(0:lmax), cls_sim(0:lmax, nrun), cls_ave(0:lmax)
  type(coop_asy)::fig
  COOP_INT::l, m, i, irun
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
  
  call map%read("lowl/commander_I_n0128_60a.fits")
  call mask%read("lowl/commander_mask_n0128_60a.fits")
  call inp%init(map, mask, lmax+50, cls)
  alms_sim = 0.d0
  cls_sim = 0.d0
  alms_plot = 0.d0
  do irun = 1, nrun
     write(*,*) "***** inpainting # ", irun, " ***********"
     call inp%upgrade(reset = .true., nside_want = 256)
     inp%lMT%map = inp%lMT%map  + inp%lCT%map !!measured map + inpainted map
     call inp%lMT%map2alm(lmax = lmax)
     alms_sim(:,:, irun) = inp%lMT%alm(0:lmax, 0:lmax, 1)
     cls_sim(:, irun) = inp%lMT%cl(0:lmax, 1)
  enddo
  do l = lmin, lmax
     do m = 0, l
        alms_plot(l,m) = sum(real(alms_sim(l, m, 1:nrun)))/nrun
     enddo
     do m = -l, -1
        alms_plot(l, m) = sum(aimag(alms_sim(l, -m, 1:nrun)))/nrun
     enddo
     cls_ave(l) = sum(cls_sim(l, 1:nrun))/nrun
  enddo
  
  call fig%open("alms.txt")
  call fig%init(xlabel = "$\ell$", ylabel = "$m$")
  call fig%density(alms_plot, dble(lmin), dble(lmax), -dble(lmax), dble(lmax), label = "$a_{lm}$", zmin = -maxval(abs(alms_plot)), zmax = maxval(abs(alms_plot)), color_table = "Planck")
  call fig%close()
  do l=2, lmax
     alms_plot(l, :) = alms_plot(l, :)*coop_sqrt2/sqrt(Cls_ave(l))
     alms_plot(l, 0) = alms_plot(l, 0)/coop_sqrt2
  enddo
  call fig%open("alms_normalized.txt")
  call fig%init(xlabel = "$\ell$", ylabel = "$m$")
  call fig%density(alms_plot, dble(lmin), dble(lmax), -dble(lmax), dble(lmax), label = "normalized $a_{lm}$", zmin = -2.5d0, zmax = 2.5d0, color_table = "Grayscale")  
  call fig%close()  
  call coop_MPI_finalize()  
end program test
