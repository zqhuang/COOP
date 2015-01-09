program hastack_prog
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
  COOP_INT:: n_bins, n, nmaps, nsims
  COOP_INT, parameter:: n_per_bin = 16, mmax = 4
  COOP_REAL,dimension(:,:,:),allocatable::fr, frr, fr_mean, fr_sq
  COOP_REAL,dimension(:,:), allocatable:: frr_mean
  COOP_REAL,dimension(:,:,:),allocatable::frr_cov
  type(coop_file)::fp
  COOP_STRING::prefix, output
  COOP_REAL::dr_arcmin, dr
  COOP_INT isim, i, j, m
  if(trim(coop_InputArgs(2)).ne."")then
     output = trim(coop_InputArgs(2))//"_"
  else
     output = ""
  end if
  prefix = trim(adjustl(coop_inputArgs(1)))
  call fp%open(trim(prefix)//"_info.txt", 'r')
  read(fp%unit, *) n, nmaps, dr_arcmin, nsims
  call fp%close()
  dr = dr_arcmin*coop_SI_arcmin
  call fp%open(trim(prefix)//".dat", "ru")
  n_bins = n/n_per_bin
  if(nmaps .gt. 2) stop "nmaps must be 1 or 2"
  print*, "n = ", n
  print*, "nmaps = ", nmaps
  print*, "nsim = ", nsims
  allocate(fr(0:n, 0:mmax/2, nmaps), fr_mean(0:n, 0:mmax/2, nmaps), fr_sq(0:n, 0:mmax/2, nmaps), frr(n_bins, 0:mmax/2, nsims), frr_mean(n_bins, 0:mmax/2), frr_cov(n_bins, n_bins, 0:mmax/2))
  frr_mean = 0.d0
  fr_mean = 0.d0
  frr_cov = 0.d0
  fr_sq = 0.d0
  do isim = 1, nsims
     read(fp%unit) j
     read(fp%unit) fr
     fr = fr *1.e6
     fr_mean = fr_mean + fr
     fr_sq = fr_sq + fr**2
     do m=0, mmax/2
        do j=1, n_bins
           frr(j,m,isim) = sum(fr((j-1)*n_per_bin:j*n_per_bin-1, m, 1:nmaps))/nmaps
        enddo
     enddo
     frr_mean(:, :) = frr_mean(:, :) + frr(:, :, isim)             
  enddo
  call fp%close()
  fr_mean = fr_mean/nsims
  fr_sq = fr_sq/nsims  - fr_mean**2  
  frr_mean = frr_mean / nsims
  do isim = 1, nsims
     frr(:,:, isim) = frr(:,:, isim) - frr_mean
     do i=1, n_bins
        do j = 1, n_bins
           frr_cov(i, j,:) = frr_cov(i, j, :) + frr(i, :, isim) * frr(j, :, isim)
        enddo
     enddo
  enddo
  frr_cov = frr_cov / nsims
  call fp%close()
  do m = 0, mmax/2
     call fp%open(trim(output)//"mean_profile_m"//COOP_STR_OF(m*2)//".txt", "w")
     do i = 0, n
        write(fp%unit, "(3E16.7)") dr*i,  sum(fr_mean(i, m, :))/nmaps, sqrt(sum(fr_sq(i, m, :))/nmaps)
     enddo
     call fp%close()
     call fp%open(trim(output)//"rad_cov_m"//COOP_STR_OF(m*2)//".txt")
     call coop_write_matrix(fp%unit, frr_cov(:, :, m), n_bins, n_bins)
     call fp%close()
  enddo
end program hastack_prog
