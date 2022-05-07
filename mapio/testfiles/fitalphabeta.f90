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
  !!fit T_sim_CMB alpha,  QU_sim_CMB beta
  COOP_INT, parameter::lmax = 2000
  COOP_INT lmin, lmax_want
  integer, parameter::nsims = 10
  type(coop_file)::fp
  COOP_INT l, i, j
  COOP_REAL alpha, beta, alpha_noise, beta_noise, chi2, r(4)
  COOP_REAL alpha_new, beta_new, alpha_noise_new, beta_noise_new, chi2_new
  COOP_REAL cl_cmb_sim(6, 0:lmax, nsims), cl_noise_sim(6, 0:lmax, nsims)
  COOP_REAL cl_data(6, 0:lmax), tenorm(0:lmax)
  lmin = coop_str2int(coop_inputArgs(1))
  lmax_want = coop_str2int(coop_inputArgs(2))
  
  call fp%open("planck14/planck14_smica_cls.txt","r")
  do l=0, lmax
     read(fp%unit, *) i, cl_data(:, l)
     if(i.ne.l) stop "error in cl file"
  enddo
  call fp%close()
  do l=0, lmax
     tenorm(l) = sqrt(cl_data(1, l)*cl_data(2, l))
  enddo
  do j = 1, nsims
     call fp%open("ffp8/ffp8_smica_cmb_"//trim(coop_Ndigits(j-1, 5))//"_cls.txt", "r")
     do l=0, lmax
        read(fp%unit, *) i, cl_cmb_sim(:, l, j)
        if(i.ne.l) stop "error in cl file"
     enddo
     call fp%close()
     call fp%open("ffp8/ffp8_smica_noise_"//trim(coop_Ndigits(j-1, 5))//"_cls.txt", "r")
     do l=0, lmax
        read(fp%unit, *) i, cl_noise_sim(:, l, j)
        if(i.ne.l) stop "error in cl file"
     enddo
     call fp%close()
  enddo

  alpha = 1.0d0
  beta = 1.0d0
  alpha_noise = 1.0d0
  beta_noise = 1.0d0
  
  alpha_new = alpha
  beta_new = beta
  alpha_noise_new = alpha_noise
  beta_noise_new = beta_noise
  call coop_random_init()
  
  chi2= chisq_new()
  do i=1, 10
     do j=1, 1000
        call random_number(r)
        r = (r-0.5d0)/100.d0/i
        alpha_new = alpha*exp(r(1))
        beta_new = beta*exp(r(2))
        alpha_noise_new = alpha_noise*exp(r(3))
        beta_noise_new = beta_noise*exp(r(3))
        chi2_new = chisq_new()
        if(chi2_new .lt. chi2)then
           alpha = alpha_new
           beta = beta_new
           alpha_noise = alpha_noise_new
           beta_noise = beta_noise_new
           chi2 = chi2_new
        endif
     enddo
  enddo
  print*, lmin, lmax_want, alpha, beta, alpha_noise, beta_noise
contains

 

    function chisq_new() result(chisq)
    COOP_REAL chisq
    COOP_INT l, imap
    COOP_REAL clrat(3)
    chisq = 0.d0
    !$omp parallel do reduction(+:chisq) private(clrat, l)
    do imap = 1, nsims
       do l=lmin, lmax_want
          clrat(1) = (alpha_new**2*cl_cmb_sim(1, l, imap) + alpha_noise_new**2*cl_noise_sim(1, l, imap))/cl_data(1, l)-1.d0
          clrat(2) = (beta_new**2*cl_cmb_sim(2, l, imap) + beta_noise_new**2*cl_noise_sim(2, l, imap))/cl_data(2, l)-1.d0
          clrat(3) = (alpha_new*beta_new*cl_cmb_sim(4, l, imap) + alpha_noise_new*beta_noise_new*cl_noise_sim(4, l,imap) - cl_data(4, l))/tenorm(l)

          chisq = chisq + sum(clrat**2)
       enddo
    enddo
    !$omp end parallel do
  end function chisq_new

  
end program test
