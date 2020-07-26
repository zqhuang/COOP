program intermittency
  use coop_wrapper_firstorder
  use coop_zeta3d_mod
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  implicit none
#include "constants.h"
  !!parameter space
  !!----------------------------------------------------------------------------------
  COOP_UNKNOWN_STRING,parameter::dir="NG"
  COOP_REAL::NG_A_zeta=1.d-6, NG_W=1.d0, NG_alpha=1.d-2, NG_lambda = 0.d0
  COOP_REAL,parameter::NG_lnHMpc = 55.d0, G_As = 2.1e-9
  !!---------------------------------------------------------------------------------
  COOP_INT,parameter::n=2048, m = 3
  COOP_REAL,dimension(m),parameter:: alpha_list = (/ 0.1, 0.001, 0. /)
  
  COOP_REAL::lnL8, LMpc, sigma, NG_A_chi
  COOP_SINGLE::zmax, fnlmax, fnlmin, xrat, yrat
  COOP_REAL,dimension(n)::chi_L_list
  COOP_REAL,dimension(n,m)::zeta_L_list, fnl_list
  COOP_INT::i, j, delta_i, ip_min, ip_max
  type(coop_asy)::fig, figfnl
  print*, "Enter: lambda="
  read(*,*) NG_lambda
  print*, "Enter: W="
  read(*,*) NG_W
  NG_A_chi = exp(2.d0*NG_lambda*NG_W)    
  print*, "Enter: L in Mpc/h "
  read(*,*) LMpc
  lnL8 = log(LMpc/8.d0)    
  call coop_set_uniform(n, chi_L_list, -sqrt(NG_A_chi)*100.d0, sqrt(NG_A_chi)*100.d0)
  delta_i = max(nint(sqrt(NG_A_chi)/(chi_L_list(2)-chi_L_list(1))), 2)
  print*,"delta_i = ", delta_i
  ip_min = delta_i+1
  ip_max = n - delta_i
  call fig%open(dir//"/zeta_L"//COOP_STR_OF(LMpc)//"_W"//COOP_STR_OF(NG_W)//"_lam"//COOP_STR_OF(NG_lambda)//".txt")
  call figfnl%open(dir//"/fnl_L"//COOP_STR_OF(LMpc)//"_W"//COOP_STR_OF(NG_W)//"_lam"//COOP_STR_OF(NG_lambda)//".txt")
  zmax=0.d0
  fnlmax = -1.d20
  fnlmin = 1.d20
  do j=1, m
     NG_alpha = alpha_list(j)
     sigma = NG_sigma_L(lnL8)
     do i=1, n
        zeta_L_list(i, j) = NG_zeta_L(chi_L_list(i), sigma)/NG_A_zeta
        if(zeta_L_list(i, j) .gt. zmax) zmax =  zeta_L_list(i, j)
     enddo
     do i=ip_min, ip_max
        fnl_list(i, j) = effective_fnl_local(chi_L_list(i-delta_i:i+delta_i), zeta_L_list(i-delta_i:i+delta_i, j))
        if(fnl_list(i,j) .gt. fnlmax) fnlmax = fnl_list(i,j)
        if(fnl_list(i,j) .lt. fnlmin) fnlmin = fnl_list(i,j)
     enddo     
  enddo
  chi_L_list = chi_L_list/sqrt(NG_A_chi)
  call fig%init(xlabel="$e^{-\lambda W}\chi_L / \chi_{\rm arm}$", ylabel = "$\zeta_L/A_\zeta$", width=6., height=5., xmin = real(chi_L_list(ip_min)*1.02), xmax = real(chi_L_list(ip_max)*1.02), ymin=0., ymax =zmax*1.05)  
  call figfnl%init(xlabel="$e^{-\lambda W} \chi_L / \chi_{\rm arm}$", ylabel = "$f_{\rm NL}^{\rm local} A_s^2/A_\zeta^3$", width=6., height=5., xmin = real(chi_L_list(ip_min)*1.02), xmax = real(chi_L_list(ip_max)*1.02), ymin = fnlmin*1.1, ymax = fnlmax*1.1)  
  do j=1, m
     call fig%plot(x = chi_L_list(ip_min:ip_max), y = zeta_L_list(ip_min:ip_max, j), color=fig%color(j), linetype=fig%linetype(j), linewidth=fig%linewidth(j), legend="$\alpha = "//COOP_STR_OF(alpha_list(j))//"$")
     call figfnl%plot(x = chi_L_list(ip_min:ip_max), y = fnl_list(ip_min:ip_max, j), color=fig%color(j), linetype=fig%linetype(j), linewidth=fig%linewidth(j), legend="$\alpha = "//COOP_STR_OF(alpha_list(j))//"$")
  enddo
  xrat = 0.42
  if(m.gt.1)then
     call fig%legend(xratio = xrat-0.05, yratio = 0.95, cols = 1, box = .false.)
     yrat = 0.65
  else
     yrat = 0.88
     call fig%label(label="$\alpha = "//COOP_STR_OF(NG_alpha)//"$", xratio = xrat, yratio = yrat, color = "black", alignment="r")
  endif
  call fig%label(label="$L_{\rm Mpc} = "//COOP_STR_OF(LMpc)//"$", xratio = xrat, yratio = yrat-0.05, color = "black", alignment="r")
  call fig%label(label="$W = "//COOP_STR_OF(NG_W)//"$", xratio = xrat, yratio = yrat - 0.1, color = "black", alignment="r")
  call fig%label(label="$\lambda = "//COOP_STR_OF(NG_lambda)//"$", xratio = xrat, yratio = yrat - 0.15, color = "black", alignment="r")  
  call fig%close()


  if(m.gt.1)then
     call figfnl%legend(xratio = xrat-0.05, yratio = 0.95, cols = 1, box = .false.)
     yrat = 0.8
  else
     yrat = 0.9          
     call figfnl%label(label="$\alpha = "//COOP_STR_OF(NG_alpha)//"$", xratio = xrat, yratio = yrat, color = "black", alignment="r")
  endif
  call figfnl%label(label="$L_{\rm Mpc} = "//COOP_STR_OF(LMpc)//"$", xratio = xrat, yratio = yrat-0.05, color = "black", alignment="r")
  call figfnl%label(label="$W = "//COOP_STR_OF(NG_W)//"$", xratio = xrat, yratio = yrat - 0.1, color = "black", alignment="r")
  call figfnl%label(label="$\lambda = "//COOP_STR_OF(NG_lambda)//"$", xratio = xrat, yratio = yrat - 0.15, color = "black", alignment="r")  
  call figfnl%close()

contains

  function NG_sigma_L(lnL8) result(sigma_L)
    COOP_REAL::lnL8
    COOP_REAL::sigma_L
    if(lnL8 + NG_lnHMpc .le. 0.) stop "LMpc too small"
    if(NG_alpha .gt. 2.d-2)then
       sigma_L = (coop_pi/2.d0/NG_alpha)**0.25d0*sqrt(NG_A_chi*erfc(-sqrt(NG_alpha/2.d0)*lnL8))
    elseif(NG_alpha .lt. 5.d-5)then
       sigma_L = sqrt(NG_A_chi * (NG_lnHMpc + lnL8))
    else
       sigma_L = (coop_pi/2.d0/NG_alpha)**0.25d0 * sqrt(NG_A_chi * (erf(sqrt(NG_alpha/2.d0)*NG_lnHMpc)-erf(sqrt(NG_alpha/2.d0)*lnL8)))
    endif
  end function NG_sigma_L

  function NG_zeta_L(chi_L, sigma_L) result(zeta_L)
    COOP_REAL::chi_L, sigma_L
    COOP_REAL::zeta_L, t, aterm
    COOP_INT::n, i
    if(chi_L .eq. 0.d0)then
       n = 0
    else
       n = log(abs(chi_L))/NG_W
    endif
    t = 2.d0*sigma_L**2
    zeta_L = exp(n*NG_W-(exp(n*NG_W)-chi_L)**2/t) + exp(n*NG_W-(exp(n*NG_W)+chi_L)**2/t)
    i = 1
    aterm = 0.d0
    do while(aterm .gt. 1.d-12*zeta_L .or. i .lt. 3)
       aterm = exp((n-i)*NG_W-(exp((n-i)*NG_W)-chi_L)**2/t) &
            + exp((n-i)*NG_W-(exp((n-i)*NG_W)+chi_L)**2/t) &
            + exp((n+i)*NG_W-(exp((n+i)*NG_W)-chi_L)**2/t) &
            + exp((n+i)*NG_W-(exp((n+i)*NG_W)+chi_L)**2/t)
       i = i + 1
       zeta_L = zeta_L + aterm
    enddo
    zeta_L = zeta_L * NG_A_zeta/sqrt(coop_2pi)/sigma_L
  end function NG_zeta_L

  !!estimate local type fnl for a list of zeta(chi), where chi is perfectly Gaussian
  function effective_fnl_local(chi_list, zeta_list) result(fnl)
    COOP_REAL, dimension(:)::zeta_list, chi_list
    COOP_REAL::fnl
    COOP_INT::i, n
    COOP_REAL::c(3), a, b
    n = coop_getdim("effective_fnl_local", size(chi_list), size(zeta_list))
    a = minval(chi_list)
    b = maxval(chi_list)
    call coop_chebfit(n, chi_list, zeta_list, 3, a, b, c)
    fnl = 2.5d0*c(3)*c(2)**2
  end function effective_fnl_local


end program intermittency
