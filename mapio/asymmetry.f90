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
  COOP_UNKNOWN_STRING, parameter::color_table = "Planck"
  COOP_UNKNOWN_STRING, parameter::spots_file = "spots/simu_iqu_Tmax_threshold0_fwhm15.txt"
  COOP_REAL, parameter::zmin = 2.
  COOP_REAL, parameter::zmax = 65.
  type(coop_healpix_patch)::patch_s, patch_n, patch
  COOP_INT, parameter::mmax = 0
  COOP_INT, parameter::n = 30
  COOP_REAL, parameter::dr = coop_SI_degree * 2.d0/n
  COOP_INT, parameter:: imap = 1
  COOP_INT, parameter::n_sims = 5000
  type(coop_healpix_maps)::map, mask
  COOP_REAL frmean(0:n), cov(0:n, 0:n), diff(0:n), hfr(0:n), theta, phi, angle, chisq
  COOP_REAL,dimension(:,:),allocatable::frs
  integer::nspots
  type(coop_asy)::fig
  type(coop_file)::fp
  integer i, j, isim, ipix, npix, weight, pix
  integer,parameter::nside_scan = 16
  type(coop_healpix_disc)::disc

  !!read mask and map, initialize patch
  call map%read("simu/simu_iqu_smoothed_fwhm15arcmin.fits", nmaps_wanted = 1)
  call mask%read("simu/simu_imask.fits", nmaps_wanted = 1)
  call patch%init("T", n, dr, mmax = mmax)


  call map%stack(patch, spots_file, mask)
  call patch%get_radial_profile(imap, 0)
  patch_s = patch
  patch_n = patch
  call coop_healpix_lb2ang(l_deg = 226.d0, b_deg = -17.d0, theta = theta, phi = phi)
  call map%stack(patch_s, spots_file, mask,  hemisphere_direction=(/ theta, phi /) )
  call patch_s%get_radial_profile(imap, 0)

  call map%stack(patch_n, spots_file, mask,  hemisphere_direction=(/ coop_pi - theta, coop_pi + phi /) )
  call patch_n%get_radial_profile(imap, 0)

  print*, patch%nstack_raw, patch_s%nstack_raw, patch_n%nstack_raw, patch_s%nstack_raw + patch_n%nstack_raw
  print*, sum(patch%nstack), sum(patch_s%nstack) + sum(patch_n%nstack)

  stop

  !!read all frs
  nspots = coop_file_numlines(spots_file)
  allocate(frs(0:n, nspots))
  print*, "reading the patches"
  call fp%open(spots_file, "r")
  i = 0
  do j=1, nspots
     read(fp%unit,*) theta, phi, angle
     call ang2pix_ring(map%nside, theta, phi, pix)
     call coop_healpix_get_disc(map%nside, pix, disc)  
     call coop_healpix_fetch_patch(map, disc, angle, patch, mask)
     if(sum(patch%nstack) .ge. (n*2.d0+1.d0)**2*coop_healpix_mask_tol)then
        i = i + 1
        call patch%get_radial_profile(imap, 0)
        frs(:, i) = patch%fr(:, 0, imap)
        if(mod(i, 1000).eq.0)write(*,*) trim(coop_num2str(100.*j/nspots,"(F10.2)"))//"% done"
     endif

  enddo
  call fp%close()
  nspots = i
  write(*,*) nspots, " patches"



  !!compute the mean 
  do i = 0, n
     frmean(i) = sum(frs(i,1:nspots))/nspots
     write(*,"(I5, 4G14.5)") i, frmean(i), patch%fr(i, 0, imap), patch_s%fr(i, 0, imap), patch_n%fr(i, 0, imap)
  enddo

  !$omp parallel do
  do i=1, nspots
     frs(:, i) = frs(:, i) - frmean
  enddo
  !$omp end parallel do

  write(*,*) "Computing the covariance matrix"
  !!compute the cov
  cov = 0.d0
  do isim=1, n_sims
     hfr = 0.
     weight = 0
     do i = 1, nspots
        if(coop_random_unit() .gt. 0.5d0)then
           hfr = hfr + frs(:, i)
           weight = weight + 1
        endif
     enddo
     hfr = hfr/weight - frmean
     do j = 0, n
        do i = 0, j
           cov(i, j) = cov(i, j) +  hfr(i)*hfr(j)
        enddo
     enddo
     if(mod(isim, 500).eq.0)then
        print*, trim(coop_num2str(100.*isim/n_sims, "(F10.1)"))//"% done"
     endif
  enddo
  cov = cov/n_sims
  do j=0, n
     do i= j+1, n
        cov(i, j) = cov(j, i)
     enddo
  enddo
  call fp%open("covmat.dat", "U")
  write(fp%unit) cov
  call fp%close()
  call coop_matsym_inverse(cov)






  diff = patch_s%fr(:, 0, imap) - frmean
  do i=0, n
     print*, i, patch_s%fr(i, 0, imap)-frmean(i), patch_n%fr(i, 0, imap)- frmean(i)
  enddo
  chisq = dot_product(diff, matmul(cov, diff))
  print*, coop_IncompleteGamma((n+1.d0)/2.d0, chisq/2.d0)/Gamma((n+1.d0)/2.d0)
  

  diff = patch_n%fr(:, 0, imap) - frmean
  chisq = dot_product(diff, matmul(cov, diff))
  print*, coop_IncompleteGamma((n+1.d0)/2.d0, chisq/2.d0)/Gamma((n+1.d0)/2.d0)


end program test
