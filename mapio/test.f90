program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  implicit none
#include "constants.h"
  integer, parameter::lmax=1500, n=50
  integer, parameter::index_corr = 4 !! 1 for T, 4 for E
  COOP_REAL,parameter::fwhm = 15.d0
  COOP_REAL,parameter::nu = 1.d0
  type(coop_file)::fp, fpI(0:2), fpQU(0:2)
  integer l, il, i, m
  type(coop_arguments)::args
  
  COOP_REAL::cls(4, 0:lmax), dtheta, x, als(0:lmax), bls(0:lmax),  ell(0:lmax), sigma, l2cls(4,0:lmax), sigma2, sigma0, sigma1, cosbeta, j2, j4, j0, theta, omega, weights(4), cr(0:1, 0:2), frI(0:2), frQU(0:2)
  call coop_random_init()
  cls = 0.d0
  sigma = fwhm*coop_sigma_by_fwhm*coop_SI_arcmin
  call fp%open("planckbest_lensedtotCls.dat", "r")
  do l=2, lmax
     read(fp%unit, *) il, cls(:, l)
     ell(l)  = l
     l2cls(:,l) = cls(:, l)*coop_2pi*exp(-l*(l+1.d0)*sigma**2)     
     cls(:,l) = l2cls(:,l)/(l*(l+1.d0))
     if(il.ne.l) stop "cl file broken"
  enddo
  call fp%close()
  dtheta = 2.d0*coop_SI_degree/n
  sigma0 = sqrt(sum(Cls(1,0:lmax)*(ell+0.5d0))/coop_2pi)
  sigma1 = sqrt(sum(l2Cls(1,0:lmax)*(ell+0.5d0))/coop_2pi)
  sigma2 = sqrt(sum(l2Cls(1,0:lmax)*(ell+0.5d0)*(ell*(ell+1.d0)))/coop_2pi)
  cosbeta = sigma1**2/sigma0/sigma2
  print*, "gamma = ", cosbeta
  print*, "sigma_0 = ", sigma0
  print*, "sigma_2 = ", sigma2  
  call coop_gaussian_npeak_set_args(args, 2, sigma0, sigma1, sigma2)  
  call coop_gaussian_get_oriented_stacking_weights(nu, args, weights)
  write(*,*) "Weights = ", weights
  call sphere_correlation_init(lmax, als, bls)
  do m=0,2
     call fpI(m)%open("I_fwhm"//COOP_STR_OF(nint(fwhm))//"_m"//COOP_STR_OF(m*2)//".txt", "w")
     call fpQU(m)%open("QU_fwhm"//COOP_STR_OF(nint(fwhm))//"_m"//COOP_STR_OF(m*2)//".txt", "w")
     write(fpI(m)%unit, "(A)") "CURVE"
     write(fpQU(m)%unit, "(A)") "CURVE"
     write(fpI(m)%unit, "(A)") "51"
     write(fpQU(m)%unit, "(A)") "51"
     write(fpI(m)%unit, "(A)") "NULL"
     write(fpQU(m)%unit, "(A)") "NULL"
     write(fpI(m)%unit, "(A)") "red_dashed"
     write(fpQU(m)%unit, "(A)") "red_dashed"
     write(fpI(m)%unit, "(A)") "0"
     write(fpQU(m)%unit, "(A)") "0"
  enddo
  do i = 0, n
     theta = dtheta*i
     omega = 2.d0*sin(theta/2.d0)
     x = cos(theta)
     cr = 0.d0     

     do l = 0, lmax
        j0 = coop_bessj(0, (l+0.5d0)*omega)
        j2 = coop_bessj(2, (l+0.5d0)*omega)
        j4 = coop_bessj(4, (l+0.5d0)*omega)
        cr(0, 0) = cr(0,0) + (l+0.5d0)*j0*cls(index_corr, l)
        cr(1, 0) = cr(1,0) - (l+0.5d0)*j0*l2cls(index_corr, l)
        cr(0, 1) = cr(0,1) + (l+0.5d0)*j2*cls(index_corr, l)
        cr(1, 1) = cr(1,1) - (l+0.5d0)*j2*l2cls(index_corr, l)
        cr(0, 2) = cr(0,2) + (l+0.5d0)*j4*cls(index_corr, l)
        cr(1, 2) = cr(1,2) - (l+0.5d0)*j4*l2cls(index_corr, l)        
     enddo
     cr = cr/coop_2pi
!!$     cr(0,0) = sphere_correlation(lmax, cls(index_corr, 0:lmax), als, bls, x)
!!$     cr(1,0) = -sphere_correlation(lmax, l2cls(index_corr, 0:lmax), als, bls, x)
     call coop_gaussian_radial_modes_I(weights, cr, frI)
     call coop_gaussian_radial_modes_QU(weights, cr, frQU)
     do m=0,2
        write(fpI(m)%unit, "(2E16.7)") omega, frI(m)
        write(fpQU(m)%unit, "(2E16.7)") omega, frQU(m)
     enddo
  enddo
  
end program test
