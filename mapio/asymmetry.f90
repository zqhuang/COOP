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

  COOP_UNKNOWN_STRING, parameter::prefix = "predx11"
  COOP_UNKNOWN_STRING, parameter::color_table = "Planck"
  COOP_UNKNOWN_STRING, parameter::spots_file = "spots/"//prefix//"_iqu_Tmax_threshold0_fwhm15.txt"
  COOP_REAL, parameter::zmin = 2.
  COOP_REAL, parameter::zmax = 65.
  type(coop_healpix_patch)::patch_s, patch_n, patch
  COOP_INT, parameter::mmax = 0
  COOP_INT, parameter::n = 30
  COOP_REAL, parameter::dr = 10.*coop_SI_arcmin
  COOP_INT, parameter:: imap = 1
  type(coop_healpix_maps)::map, mask, sim
  COOP_REAL frmean(0:n), cov(0:n,0: n), diff(0:n), hfr(0:n), theta, phi, angle, chisq, err(0:n)
  COOP_REAL,dimension(:,:),allocatable::frs
  COOP_REAL,dimension(:),allocatable::theta_spots, phi_spots
  integer::nspots
  type(coop_asy)::fig
  type(coop_file)::fp
  integer i, j,  npix, weight, pix, step
  integer,parameter::nside_scan = 16
  type(coop_healpix_disc)::disc

  !!read mask and map, initialize patch
  call map%read(prefix//"/"//prefix//"_iqu_smoothed_fwhm15arcmin.fits", nmaps_wanted = 1)
  call mask%read(prefix//"/"//prefix//"_imask.fits", nmaps_wanted = 1)
  call patch%init("T", n, dr, mmax = mmax)
  sim = map

  patch_s = patch
  patch_n = patch

  call coop_healpix_lb2ang(l_deg = 226.d0, b_deg = -17.d0, theta = theta, phi = phi)

  call map%stack(patch_s, spots_file, mask,  hemisphere_direction=(/ theta, phi /) )
  call patch_s%get_radial_profile(imap, 0)

  call map%stack(patch_n, spots_file, mask,  hemisphere_direction=(/ coop_pi - theta, coop_pi + phi /) )
  call patch_n%get_radial_profile(imap, 0)

  call coop_random_init()
  !!read all frs
  nspots = coop_file_numlines(spots_file)
  allocate(frs(0:n, nspots))
  allocate(theta_spots(nspots), phi_spots(nspots))
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
        theta_spots(i) = theta
        phi_spots(i) = phi
        if(mod(i, 1000).eq.0)write(*,*) trim(coop_num2str(100.*j/nspots,"(F10.2)"))//"% done"
     endif

  enddo
  call fp%close()

  call sim%map2alm(lmax=2500)
  call sim%simulate()

  nspots = i
  write(*,*) nspots, " patches"
  


  !!compute the mean 
  do i = 0, n
     frmean(i) = sum(frs(i,1:nspots))/nspots
  enddo

  !$omp parallel do
  do i=1, nspots
     frs(:, i) = frs(:, i) - frmean
  enddo
  !$omp end parallel do

  write(*,*) "Computing the covariance matrix"
  !!compute the cov
  cov = 0.d0
  npix = nside_scan**2 * 12
  do pix = 0, npix-1
     call pix2ang_ring(nside_scan, pix, theta, phi)
     hfr = 0.
     weight = 0
     do i = 1, nspots
        if(cos(theta_spots(i))*cos(theta) + sin(theta_spots(i))*sin(theta)*cos(phi_spots(i)-phi) .gt. 0.d0  )then
           hfr = hfr + frs(:, i)
           weight = weight + 1
        endif
     enddo
     hfr = hfr/weight 
     do j = 0, n
        do i = 0, j
           cov(i, j) = cov(i, j) +  hfr(i)*hfr(j)
        enddo
     enddo
     if(mod(pix, 500).eq.0)then
        print*, trim(coop_num2str(100.*pix/npix, "(F10.1)"))//"% done"
     endif
  enddo
  cov = cov/npix
  do j=0, n
     do i= j+1, n
        cov(i, j) = cov(j, i)
     enddo
     err(j) = sqrt(cov(j,j))
  enddo
  call fp%open(prefix//"covmat.dat", "U")
  write(fp%unit) cov
  call fp%close()
  
  

  call coop_matsym_inverse(cov)



  diff = patch_s%fr(:, 0, imap) - frmean(:)
  chisq = dot_product(diff, matmul(cov, diff))
  print*, chisq, coop_IncompleteGamma((n)/2.d0, chisq/2.d0)/Gamma((n)/2.d0)
  

  diff = patch_n%fr(:, 0, imap) - frmean(:)
  chisq = dot_product(diff, matmul(cov, diff))
  print*, chisq, coop_IncompleteGamma((n)/2.d0, chisq/2.d0)/Gamma((n)/2.d0)

  call fig%open(prefix//"_radial_profile.txt")
  call fig%init(xlabel = "$r$", ylabel = "$T(\mu K)$")
  call coop_asy_curve(fig, patch_n%r, patch_n%fr(:, 0, imap), legend = "North", color="red", linetype = "dashed")
  call coop_asy_curve(fig, patch_s%r, patch_s%fr(:, 0, imap), legend = "South", color="blue", linetype = "dotted")
  call coop_asy_curve(fig, patch_s%r, frmean, legend = "mean", color = "green", linetype = "solid")
  step = ceiling(15.*coop_SI_arcmin/dr)
  !err = err*sqrt(real(step))
  do i=step/2, n, step
     call coop_asy_error_bar(fig, patch%r(i), frmean(i), dy_plus = err(i), dy_minus=err(i))
  enddo
  call coop_asy_legend(fig)
  call fig%close()


end program test
