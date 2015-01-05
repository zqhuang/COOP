program stackth
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
  integer, parameter::lmax=2000, n=50, index_temp = 1, index_pol = 4
  COOP_REAL, parameter:: r = 2.d0*coop_SI_degree
  COOP_REAL, parameter:: dr = r/n

  !!settings
  integer, parameter::index_corr = index_pol  !!index_temp
  COOP_UNKNOWN_STRING,parameter::clfile = "planckbest_lensedtotCls.dat"  
  COOP_UNKNOWN_STRING, parameter::spot_type = "PTmax"  !!"Tmax_QTUTOrient"  
  COOP_REAL, parameter::nu = 1.d0  !!threshold  
  COOP_REAL, parameter::fwhm = 15.d0
  
  
  type(coop_file)::fp, fpI(0:2), fpQU(0:2)
  type(coop_asy)::figI, figQ, figU
  integer l, il, i, j, iomega, m
  type(coop_arguments)::args
  type(coop_healpix_patch)::patchI, patchQU, patchQrUr
  
  COOP_REAL::cls(4, 0:lmax), ell(0:lmax), sigma, l2cls(4,0:lmax), sigma2, sigma0, sigma1, cosbeta, j2, j4, j0, omega, weights(4), cr(0:1, 0:2), frI(0:2, 0:n*3/2), frQU(0:2, 0:n*3/2), pomega, phi,  romega
  
  call coop_random_init()
  cls = 0.d0
  ell = 0.d0
  
  sigma = fwhm*coop_sigma_by_fwhm*coop_SI_arcmin
  call fp%open(clfile, "r")
  do l=2, lmax
     read(fp%unit, *) il, cls(:, l)
     ell(l)  = l
     l2cls(:,l) = cls(:, l)*coop_2pi*exp(-l*(l+1.d0)*sigma**2)     
     cls(:,l) = l2cls(:,l)/(l*(l+1.d0))
     if(il.ne.l) stop "cl file broken"
  enddo
  call fp%close()
  sigma0 = sqrt(sum(Cls(1,0:lmax)*(ell+0.5d0))/coop_2pi)
  sigma1 = sqrt(sum(l2Cls(1,0:lmax)*(ell+0.5d0))/coop_2pi)
  sigma2 = sqrt(sum(l2Cls(1,0:lmax)*(ell+0.5d0)*(ell*(ell+1.d0)))/coop_2pi)
  cosbeta = sigma1**2/sigma0/sigma2
  print*, "gamma = ", cosbeta
  print*, "sigma_0 = ", sigma0
  print*, "sigma_2 = ", sigma2  
  call coop_gaussian_npeak_set_args(args, 2, sigma0, sigma1, sigma2)
  select case(trim(spot_type))
  case("Tmax_QTUTOrient")
     call coop_gaussian_get_oriented_stacking_weights(nu, args, weights)
  case("Tmax")
     call coop_gaussian_get_nonoriented_stacking_weights(nu, args, weights)
  case("PTmax")
     call coop_gaussian_get_pmax_stacking_weights(nu, args, weights)
  end select
  write(*,*) "Weights = ", weights
  do m=0,2
     call fpI(m)%open("I_fwhm"//COOP_STR_OF(nint(fwhm))//"_m"//COOP_STR_OF(m*2)//".txt", "w")
     call fpQU(m)%open("QU_fwhm"//COOP_STR_OF(nint(fwhm))//"_m"//COOP_STR_OF(m*2)//".txt", "w")
     write(fpI(m)%unit, "(A)") "CURVE"
     write(fpQU(m)%unit, "(A)") "CURVE"
     write(fpI(m)%unit, "(I6)") n+1
     write(fpQU(m)%unit, "(I6)") n+1
     write(fpI(m)%unit, "(A)") "Theory"
     write(fpQU(m)%unit, "(A)") "Theory"
     write(fpI(m)%unit, "(A)") "red_dashed"
     write(fpQU(m)%unit, "(A)") "red_dashed"
     write(fpI(m)%unit, "(A)") "0"
     write(fpQU(m)%unit, "(A)") "0"
  enddo
  write(*,*) "now compute c(r)"
  do i = 0, n*3/2
     omega = dr*i
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
     call coop_gaussian_radial_modes_I(weights, cr, frI(:, i))
     call coop_gaussian_radial_modes_QU(weights, cr, frQU(:, i))
     if(i.le.n)then
        do m=0,2
           write(fpI(m)%unit, "(2E16.7)") omega, frI(m, i)
           write(fpQU(m)%unit, "(2E16.7)") omega, frQU(m, i)
        enddo
     endif
  enddo
  call patchI%init("I", n, dr, 4)
  call patchQU%init("QU", n, dr, 4)
  call patchQrUr%init("QrUr", n, dr, 4)
  patchI%nstack = 1.d0
  patchQU%nstack = 1.d0
  patchQrUr%nstack = 1.d0
  patchI%nstack_raw = 1.d0
  patchQU%nstack_raw = 1.d0
  patchQrUr%nstack_raw = 1.d0  
  do i=-n, n
     do j=-n, n
        if(i.eq.0.and.j.eq.0)then
           patchI%image(i,j,1) = frI(0,0)
           patchQU%image(i,j,1) = frQU(0, 0)
           patchQU%image(i,j,2) = 0.d0
           patchQrUr%image(i,j,1) = frQU(0, 0)
           patchQrUr%image(i,j,2) = 0.d0
        else
           phi = atan2(dble(j), dble(i))
           pomega = sqrt(dble(i)**2+dble(j)**2)
           iomega = floor(pomega)
           romega = pomega - iomega
           patchI%image(i,j,1) = (frI(0, iomega)+ frI(1, iomega)*cos(2*phi) + frI(2, iomega)*cos(4*phi) ) *(1.d0-romega)+(frI(0, iomega+1)+ frI(1, iomega+1)*cos(2*phi) + frI(2, iomega+1)*cos(4*phi) ) * romega
           patchQU%image(i, j,1) = (frI(0, iomega)+ frI(1, iomega)*cos(2*phi) + frI(2, iomega)*cos(4*phi) ) *(1.d0-romega)+(frI(0, iomega+1)+ frI(1, iomega+1)*cos(2*phi) + frI(2, iomega+1)*cos(4*phi) ) * romega
           patchQU%image(i, j, 2) = (frI(1, iomega)*sin(2*phi) + frI(2, iomega)*sin(4*phi) ) *(1.d0-romega)+(frI(0, iomega+1)+ frI(1, iomega+1)*sin(2*phi) + frI(2, iomega+1)*sin(4*phi) ) * romega
           patchQrUr%image(i, j, 1) = -patchQU%image(i, j,1)*cos(2*phi) - patchQU%image(i, j, 2)*sin(2*phi)  !!urgh, the stupid minus sign...
           patchQrUr%image(i, j, 2) = patchQU%image(i, j,1)*sin(2*phi) - patchQU%image(i, j, 2)*cos(2*phi)
        endif
     enddo
  enddo
  patchI%caption = "Best-fit $\Lambda$CDM theory"
  patchQU%caption = "Best-fit $\Lambda$CDM theory"
  patchQrUr%caption = "Best-fit $\Lambda$CDM theory"
  call patchI%plot(1, "I_stack.txt")
  call patchQU%plot(1, "Q_stack.txt")
  call patchQU%plot(2, "U_stack.txt")
  call patchQrUr%plot(1, "Qr_stack.txt")
  call patchQrUr%plot(2, "Ur_stack.txt")    
end program stackth
