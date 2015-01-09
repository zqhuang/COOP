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
  integer, parameter::lmin = 2, lmax=2000, n=50, index_temp = 1, index_pol = 4
  COOP_REAL, parameter:: r_degree = 2.d0
  COOP_REAL, parameter:: width = 2.d0*sin(r_degree*coop_SI_degree/2.d0)
  COOP_REAL, parameter:: dr = width/n

  !!settings
  integer, parameter::index_corr = index_temp  !!index_temp
  COOP_UNKNOWN_STRING,parameter::clfile = "planckbest_lensedtotCls.dat"  
  COOP_UNKNOWN_STRING, parameter::spot_type = "Tmax_QTUTOrient"  
  COOP_REAL, parameter::nu = -8.d0  !!threshold  
  COOP_REAL, parameter::fwhm = 15.d0
  COOP_INT::head_level
  COOP_STRING::prefix, line
  type(coop_file)::fp
  type(coop_asy)::figCr, fpI(0:2), fpQU(0:2)
  integer l, il, i, j, iomega, m
  type(coop_arguments)::args
  type(coop_healpix_patch)::patchI, patchQU, patchQrUr
  
  COOP_REAL::cls(4, 2:lmax), ell(2:lmax), sigma, l2cls(4,2:lmax), sigma2, sigma0, sigma1, cosbeta, j2, j4, j0, omega, weights(4), cr(0:1, 0:2, 0:n*3/2), frI(0:2, 0:n*3/2), frQU(0:2, 0:n*3/2), pomega, phi, romega, r(0:n*3/2), kr
  if(trim(coop_InputArgs(1)).ne."")then
     prefix = trim(coop_InputArgs(1))//"_"
  else
     prefix = ""
  endif
  if(trim(coop_InputArgs(2)).ne."")then
     line = coop_InputArgs(2)
     read(line,*) head_level
  else
     head_level = 0
  endif
  call coop_random_init()
  sigma = fwhm*coop_sigma_by_fwhm*coop_SI_arcmin
  call fp%open(clfile, "r")
  do l=2, lmax
     read(fp%unit, *) il, cls(:, l)
     ell(l)  = l
     l2cls(:,l) = cls(:, l)*coop_2pi*exp(-l*(l+1.d0)*sigma**2)
     l2cls(1:3, l) = l2cls(1:3, l) +  l*(l+1.)*(/ coop_Planck_TNoise(l), coop_Planck_ENoise(l), coop_Planck_BNoise(l) /)    
     cls(:,l) = l2cls(:,l)/(l*(l+1.d0))
     if(il.ne.l) stop "cl file broken"
  enddo
  call fp%close()
  sigma0 = sqrt(sum(Cls(1,:)*(ell+0.5d0))/coop_2pi)
  sigma1 = sqrt(sum(l2Cls(1,:)*(ell+0.5d0))/coop_2pi)
  sigma2 = sqrt(sum(l2Cls(1,:)*(ell+0.5d0)*(ell*(ell+1.d0)))/coop_2pi)
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
  write(*,*) "Weights = ", weights(1)*sigma0, weights(2)*sigma2, weights(3)*sigma0, weights(4)*sigma2
  do m=0,2
     call fpI(m)%open(trim(prefix)//"I_fwhm"//COOP_STR_OF(nint(fwhm))//"_m"//COOP_STR_OF(m*2)//".txt")
     call fpQU(m)%open(trim(prefix)//"QU_fwhm"//COOP_STR_OF(nint(fwhm))//"_m"//COOP_STR_OF(m*2)//".txt")
     if(head_level .ge. 1)then
        if(head_level.ge.2)then
           call fpI(m)%init(xlabel="$\varpi$", ylabel="$T_"//COOP_STR_OF(m*2)//"$")
           call fpQU(m)%init(xlabel="$\varpi$", ylabel="$P_"//COOP_STR_OF(m*2)//"$")
        endif
        write(fpI(m)%unit, "(A)") "CURVE"
        write(fpQU(m)%unit, "(A)") "CURVE"
        write(fpI(m)%unit, "(I6)") n+1
        write(fpQU(m)%unit, "(I6)") n+1
        write(fpI(m)%unit, "(A)") "Theory"
        write(fpQU(m)%unit, "(A)") "Theory"
        write(fpI(m)%unit, "(A)") "blue_dashed"
        write(fpQU(m)%unit, "(A)") "blue_dashed"
        write(fpI(m)%unit, "(A)") "0"
        write(fpQU(m)%unit, "(A)") "0"
     endif
  enddo
  cr = 0.d0       
  do i = 0, n*3/2
     omega = dr*i
     r(i) = omega
     do l = lmin, lmax
        kr = (l+0.5d0)*omega
        j0 = coop_bessj(0, kr)
        j2 = coop_bessj(2, kr)
        j4 = coop_bessj(4, kr)
        cr(0, 0, i) = cr(0, 0, i) + (l+0.5d0)*j0*cls(index_corr, l)
        cr(1, 0, i) = cr(1, 0, i) - (l+0.5d0)*j0*l2cls(index_corr, l)
        cr(0, 1, i) = cr(0, 1, i) + (l+0.5d0)*j2*cls(index_corr, l)
        cr(1, 1, i) = cr(1, 1, i) - (l+0.5d0)*j2*l2cls(index_corr, l)
        cr(0, 2, i) = cr(0, 2, i) + (l+0.5d0)*j4*cls(index_corr, l)
        cr(1, 2, i) = cr(1, 2, i) - (l+0.5d0)*j4*l2cls(index_corr, l)        
     enddo
     cr(:,:,i) = cr(:,:,i)/coop_2pi
     call coop_gaussian_radial_modes_I(weights, cr(:,:,i), frI(:, i))
     call coop_gaussian_radial_modes_QU(weights, cr(:,:,i), frQU(:, i))
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
           omega = sqrt(dble(i)**2+dble(j)**2)
           iomega = floor(omega)
           romega = omega - iomega
           patchI%image(i,j,1) = (frI(0, iomega)+ frI(1, iomega)*cos(2*phi) + frI(2, iomega)*cos(4*phi) ) *(1.d0-romega)+(frI(0, iomega+1)+ frI(1, iomega+1)*cos(2*phi) + frI(2, iomega+1)*cos(4*phi) ) * romega
           patchQU%image(i, j,1) = (frQU(0, iomega)+ frQU(1, iomega)*cos(2*phi) + frQU(2, iomega)*cos(4*phi) ) *(1.d0-romega)+(frQU(0, iomega+1)+ frQU(1, iomega+1)*cos(2*phi) + frQU(2, iomega+1)*cos(4*phi) ) * romega
           patchQU%image(i, j, 2) = (frQU(1, iomega)*sin(2*phi) + frQU(2, iomega)*sin(4*phi) ) *(1.d0-romega)+(frQU(0, iomega+1)+ frQU(1, iomega+1)*sin(2*phi) + frQU(2, iomega+1)*sin(4*phi) ) * romega
           patchQrUr%image(i, j, 1) = -patchQU%image(i, j,1)*cos(2*phi) - patchQU%image(i, j, 2)*sin(2*phi)  !!urgh, the stupid minus sign...
           patchQrUr%image(i, j, 2) = patchQU%image(i, j,1)*sin(2*phi) - patchQU%image(i, j, 2)*cos(2*phi)
        endif
     enddo
  enddo
  patchI%caption = "best-fit $\Lambda$CDM theory, $\nu="//COOP_STR_OF(nint(nu))//"$"
  patchQU%caption = "best-fit $\Lambda$CDM theory, $\nu="//COOP_STR_OF(nint(nu))//"$"
  patchQrUr%caption = "best-fit $\Lambda$CDM theory, $\nu="//COOP_STR_OF(nint(nu))//"$"
  call patchI%plot(1, trim(prefix)//"I_stack.txt")
  call patchQU%plot(1, trim(prefix)//"Q_stack.txt")
  call patchQU%plot(2, trim(prefix)//"U_stack.txt")
  call patchQrUr%plot(1, trim(prefix)//"Qr_stack.txt")
  call patchQrUr%plot(2, trim(prefix)//"Ur_stack.txt")
  call figCr%open(trim(prefix)//"cr.txt")
  call figCr%init(xlabel="$\varpi$", ylabel="$c_{m,n}(\varpi)$")
  call figCr%curve(x = r, y = cr(0,0,:)/sigma0, color="red", linewidth=1.8, linetype="solid", legend="$c_{0,0}(\varpi)/\sigma_0$")
  call figCr%curve(x = r, y = cr(0,1,:)/sigma0, color="orange", linewidth=1.5, linetype="dotted", legend="$c_{0,1}(\varpi)/\sigma_0$")
  call figCr%curve(x = r, y = cr(0,2,:)/sigma0, color="violet", linewidth=1.2, linetype="dashed", legend="$c_{0,2}(\varpi)/\sigma_0$")  
  call figCr%curve(x = r, y = -cr(1,0,:)/sigma2, color="blue", linewidth=1.8, linetype="solid", legend="$-c_{1,0}(\varpi)/\sigma_2$")  
  
  call figCr%curve(x = r, y = -cr(1,1,:)/sigma2, color="black", linewidth=1.5, linetype="dotted", legend="$-c_{1,1}(\varpi)/\sigma_2$")  

  call figCr%curve(x = r, y = -cr(1,2,:)/sigma2, color="skblue", linewidth=1.2, linetype="dashed", legend="$-c_{1,2}(\varpi)/\sigma_2$")  
  
  call figCr%legend(0.5, 0.9, 1)
  call figCr%close()
end program stackth
