module mask_pointsource
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use coop_fitsio_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools

  implicit none
#include "constants.h"

  logical::verbose = .true.
  COOP_UNKNOWN_STRING,parameter::dir = '/home/zqhuang/work/ZAP512/'
  COOP_INT,parameter::num_freqs = 7
  COOP_INT,parameter::num_sims = 50
  COOP_SINGLE,dimension(num_freqs),parameter::freqs = (/ 23., 95., 100., 143., 150., 217., 353. /)
  COOP_INT,dimension(num_freqs),parameter::freq_names = nint(freqs)
  COOP_REAL, dimension(num_freqs),parameter::fwhms = (/ 52.8, 19., 9.651, 7.248, 4.990, 11., 4.818 /)*coop_SI_arcmin
  COOP_REAL, parameter:: fwhm_common = 52.800000001d0 * coop_SI_arcmin
  COOP_REAL, dimension(num_freqs), parameter:: smooth_fwhms = sqrt(fwhm_common**2 - fwhms**2)
  COOP_INT,parameter::nside = 512
  COOP_INT,parameter::npix_tot = nside**2*12
  COOP_UNKNOWN_STRING, parameter::mask_file = dir//'maskapo512.fits'
  COOP_UNKNOWN_STRING, parameter::binary_mask_file = dir//'mask512.fits'  
  COOP_UNKNOWN_STRING, parameter::noise_std_prefix = dir//"NoiseSTD/noise_std_"
  COOP_SINGLE:: Freq_REF_D = 353.
  COOP_SINGLE:: Freq_REF_S = 23.
  COOP_SINGLE:: ell_REF = 100.
  COOP_SINGLE:: fwhm_REF_D = 4.818 * coop_SI_arcmin
  COOP_SINGLE:: fwhm_REF_S = 52.8 * coop_SI_arcmin 


contains

  subroutine apply_QUmask(hp, mask)
    type(coop_healpix_maps)::hp, mask
    hp%map(:, 2) = hp%map(:, 2)*mask%map(:, 1)
    hp%map(:, 3) = hp%map(:, 3)*mask%map(:, 1)    
  end subroutine apply_QUmask

  subroutine smooth_QU(hp, fwhm)
    type(coop_healpix_maps)::hp
    COOP_REAL::fwhm
    if(abs(fwhm) < 1.d-5) return !do nothing
    call hp%map2alm(lmax = 300, index_list = (/ 2, 3 /) )
    call hp%filter_alm(fwhm = fwhm, index_list=(/ 2, 3/) )
    call hp%alm2map(index_list = (/ 2, 3 /) )
  end subroutine smooth_QU

  function convert_i2cmb(freq, ref) result(amp)
    COOP_SINGLE::freq, ref, hGk_t0, p, p0, amp
    hGk_t0 = 0.0479924/2.726  ! h*GHz/k_B/T0
    p = hGk_t0*freq
    p0 = hGk_t0*ref
    amp = (ref/freq)**2*exp(p0-p)*(exp(p)-1.)**2/(exp(p0)-1.)**2
  end function convert_i2cmb

  function dust_brightness_ratio(freq, ref) result(amp)
    COOP_SINGLE::freq, ref, hGk_td, p, p0, amp
    hGk_td = 0.0479924/19.6 !dust effective MBB temperature 19.6K
    amp =  (freq/ref)*(exp(hGk_td*ref)-1.)/(exp(hGk_td*freq)-1.)
  end function dust_brightness_ratio

  function dust_cmb_ratio(freq, ref) result(amp)
    COOP_SINGLE::freq, ref, amp
    amp = convert_i2cmb(freq, ref)*dust_brightness_ratio(freq, ref)
  end function dust_cmb_ratio
  
  function noise_file_name(freq, i)
    COOP_INT::freq, i
    COOP_STRING::noise_file_name
    COOP_UNKNOWN_STRING, parameter::noise_prefix = dir//'Noise_BOTH/NOISE_'    
    noise_file_name = noise_prefix//COOP_STR_OF(freq)//"GHz_BOTH_"//COOP_STR_OF(i)//".fits"
  end function noise_file_name


  subroutine get_noise_std(mask)
    type(coop_healpix_maps)::hp, mask
    COOP_REAL,dimension(:,:),allocatable::mean, var
    COOP_INT::ifreq, i
    allocate(mean(0:npix_tot-1, 3), var(0:npix_tot-1, 3))
    do ifreq = 1, num_freqs       
       if(.not. coop_file_exists(noise_std_prefix//COOP_STR_OF(freq_names(ifreq))//"GHz.fits"))then
          mean = 0.d0
          var = 0.d0
          do i=0, num_sims-1
             if(verbose)then
                print*, "loading "//trim(noise_file_name(freq_names(ifreq), i))
             endif             
             call hp%read(noise_file_name(freq_names(ifreq), i), nmaps_wanted = 3)
             call apply_QUmask(hp, mask)
             call smooth_QU(hp, smooth_fwhms(ifreq))
             mean = mean + hp%map
             var = var + hp%map**2
          enddo
          mean = mean/num_sims
          hp%map =  sqrt(var/num_sims - mean**2)
          call hp%write(noise_std_prefix//COOP_STR_OF(freq_names(ifreq))//"GHz.fits")
          if(verbose)then
             print*, noise_std_prefix//COOP_STR_OF(freq_names(ifreq))//"GHz.fits is done"
          endif
       else
          if(verbose)then
             print*, "loading "//noise_std_prefix//COOP_STR_OF(freq_names(ifreq))//"GHz.fits"
          endif          
       endif
    enddo
    deallocate(mean, var)
  end subroutine get_noise_std

  function fg_file_name(freq)
    COOP_INT::freq
    COOP_STRING::fg_file_name
    COOP_UNKNOWN_STRING,parameter::fg_prefix = dir//"Foreground/FG_"
    fg_file_name = fg_prefix//COOP_STR_OF(freq)//"GHz_BOTH.fits"
  end function fg_file_name


  function fginp_file_name(freq)
    COOP_INT::freq
    COOP_STRING::fginp_file_name
    COOP_UNKNOWN_STRING,parameter::fginp_prefix = dir//"FGINPmaps/FG_"
    fginp_file_name = fginp_prefix//COOP_STR_OF(freq)//"GHz_BOTH.fits"
  end function fginp_file_name
  
  

  function sim_file_name(freq)
    COOP_INT::freq
    COOP_STRING::sim_file_name
    COOP_UNKNOWN_STRING,parameter::sim_prefix = dir//"simplesims/sim_d1s1c1_"
    sim_file_name = sim_prefix//COOP_STR_OF(freq)//"GHz.fits"
  end function sim_file_name

  function data_file_name(freq)
    COOP_INT::freq
    COOP_STRING::data_file_name
    COOP_UNKNOWN_STRING,parameter::data_prefix = dir//"DCmaps/DC1_"
    data_file_name = data_prefix//COOP_STR_OF(freq)//"GHz_BOTH.fits"
  end function data_file_name

  function inp_file_name(freq)
    COOP_INT::freq
    COOP_STRING::inp_file_name
    COOP_UNKNOWN_STRING,parameter::inp_prefix = dir//"INPmaps/DC1_"
    inp_file_name = inp_prefix//COOP_STR_OF(freq)//"GHz_BOTH.fits"
  end function inp_file_name
  
  subroutine symmat3_inv(a)
    COOP_SINGLE::a(3,3), det
    det =  a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(3,2)) + a(2,1)*(a(3,2)*a(3,1)-a(2,1)*a(3,3)) + a(3,1)*(a(2,1)*a(3,2)- a(2,2)*a(3,1)) 
    a = Reshape( (/ a(2,2)*a(3,3) - a(3,2)*a(3,2), a(3,2)*a(3,1)-a(2,1)*a(3,3), a(2,1)*a(3,2)- a(2,2)*a(3,1),  &
         a(3,2)*a(3,1)-a(2,1)*a(3,3), a(3,3)*a(1,1)-a(3,1)*a(3,1), a(2,1)*a(3,1)-a(1,1)*a(3,2), &
         a(2,1)*a(3,2)- a(2,2)*a(3,1), a(2,1)*a(3,1)-a(1,1)*a(3,2),  a(1,1)*a(2,2)-a(2,1)*a(2,1) /), (/ 3, 3 /) ) &
         / det
  end subroutine symmat3_inv


  subroutine three_vec_proj(vtot, vecs, coefs, res)
    COOP_SINGLE::vtot(num_freqs), vecs(num_freqs, 3), coefs(3), res
    COOP_SINGLE::b(3), mat(3, 3)
    COOP_INT::i, j
    do i=1, 3
       b(i) = dot_product(vtot, vecs(:, i))
       do j=1, i
          mat(i, j) = dot_product(vecs(:, i), vecs(:, j))
       enddo
    enddo
    call symmat3_inv(mat)
    coefs = matmul(mat, b)
    res = sum((vtot - matmul(vecs, coefs))**2)/num_freqs
  end subroutine three_vec_proj


  subroutine six_vec_proj(qtot, utot, qvecs, uvecs, coefs, res)
    COOP_SINGLE::qtot(num_freqs), utot(num_freqs), qvecs(num_freqs, 3), uvecs(num_freqs, 3), coefs(3), res
    COOP_SINGLE::b(3), mat(3, 3)
    COOP_INT::i, j
    do i=1, 3
       b(i) = dot_product(qtot, qvecs(:, i)) + dot_product(utot, uvecs(:, i))
       do j=1, i
          mat(i, j) = dot_product(qvecs(:, i), qvecs(:, j)) + dot_product(uvecs(:, i), uvecs(:, j))
       enddo
    enddo
    call symmat3_inv(mat)
    coefs = matmul(mat, b)
    res = ( sum((qtot - matmul(qvecs, coefs))**2) +  sum((utot - matmul(uvecs, coefs))**2) )/(num_freqs*2)
  end subroutine six_vec_proj

  
end module mask_pointsource

program test
  use mask_pointsource
  implicit none
  COOP_INT::npix_mask, i, j, ifreq, istart, nneigh, list(8), npix_inp
  COOP_INT,dimension(:),allocatable::mask_pix, inp_pix
  COOP_SINGLE,dimension(:,:),allocatable::Q_std, U_std, Q_data, U_data, coefs
  COOP_SINGLE,dimension(:),allocatable::res, tmp
  COOP_INT,dimension(:), allocatable::res_ind
  COOP_INT,dimension(:,:), allocatable::coefs_ind
  COOP_SINGLE::res_cut, coefs_cut(3)
  type(coop_healpix_maps)::hp, mask, hpd, binary_mask
  COOP_SINGLE,dimension(:),allocatable::inpmask
  COOP_SINGLE,parameter::mask_tiny = 1.e-3, threshold = 0.998
  call mask%read(mask_file, nmaps_wanted=1)
  call binary_mask%read(binary_mask_file, nmaps_wanted=1)    
  call get_noise_std(mask)
  npix_mask = count(mask%map > mask_tiny)
  print*, "before masking point sources, mask size = ", npix_mask  
  allocate(mask_pix(npix_mask))
  j = 0
  do i=0, npix_tot - 1
     if(mask%map(i, 1) > mask_tiny)then
        j=j+1
        mask_pix(j) = i
     endif
  enddo
  allocate( Q_std(num_freqs, npix_mask), U_std(num_freqs, npix_mask) )
  allocate( Q_data(num_freqs, npix_mask), U_data(num_freqs, npix_mask) )
  allocate( coefs(3, npix_mask),  coefs_ind(npix_mask, 3))
  allocate( res(npix_mask), res_ind(npix_mask) )
  do ifreq = 1, num_freqs
     call hp%read(noise_std_prefix//COOP_STR_OF(freq_names(ifreq))//"GHz.fits", nmaps_wanted=3)
     call hpd%read(data_file_name(freq_names(ifreq)), nmaps_wanted = 3)
     call apply_QUmask(hpd, mask)     
     call smooth_QU(hpd, smooth_fwhms(ifreq))
     do i=1, npix_mask
        Q_std(ifreq, i) = hp%map(mask_pix(i), 2)
        U_std(ifreq, i) = hp%map(mask_pix(i), 3)
        Q_data(ifreq, i) = hpd%map(mask_pix(i), 2)/Q_std(ifreq, i)
        U_data(ifreq, i) = hpd%map(mask_pix(i), 3)/U_std(ifreq, i)        
     enddo
  enddo
  call do_projection(1.6, -3.)
  istart = nint(npix_mask*threshold)
  tmp = res
  call coop_quicksort_index(tmp, res_ind)
  call coop_asy_histogram(dble(res), nbins=31, filename = 'histo_res.txt')
  mask%map(mask_pix(res_ind(istart:npix_mask)), 1) = -1.
  res_cut = max(tmp(istart), 2.)
  print*, "noise cut = ", res_cut
  do i=1, 3
     tmp = abs(coefs(i, :))
     call coop_quicksort_index(tmp, coefs_ind(:, i))
     call coop_asy_histogram(dble(coefs(i, :)), nbins=31, filename = 'histo_coefs'//COOP_STR_OF(i)//'.txt')
     if(i<3) mask%map(mask_pix(coefs_ind(istart:npix_mask, i)), 1) = -1.     
     coefs_cut(i) = tmp(istart)
     print*, "component #", i, " cut = ", coefs_cut(i)
  enddo
  npix_inp = count(mask%map < -mask_tiny)
  allocate(inp_pix(npix_inp))
  call mask%convert2nested()
  j = 0
  do i=0, npix_tot - 1
     if(mask%map(i, 1) < -mask_tiny)then
        j = j + 1
        inp_pix(j) = i
     endif
  enddo
  allocate(inpmask(0:mask%npix-1))
  if(.true.)then !!inpaint foreground
     do ifreq = 1, num_freqs
        inpmask = mask%map(:, 1)
        call hpd%read(fg_file_name(freq_names(ifreq)), nmaps_wanted = 3)
        call apply_QUmask(hpd, binary_mask)
        call hpd%convert2nested()
        j = 0
        do while(any(inpmask < -mask_tiny))
           j = j + 1
           print*, "doing inpaint step #", j, ", npix TBD:", count(inpmask<-mask_tiny)
           call inp_step()
        enddo
        call hpd%convert2ring()
        call hpd%write(fginp_file_name(freq_names(ifreq)))
     enddo
  else  !!inpaint maps
     do ifreq = 1, num_freqs
        inpmask = mask%map(:, 1)
        call hpd%read(data_file_name(freq_names(ifreq)), nmaps_wanted = 3)
        call apply_QUmask(hpd, binary_mask)
        call hpd%convert2nested()
        j = 0
        do while(any(inpmask < -mask_tiny))
           j = j + 1
           print*, "doing inpaint step #", j, ", npix TBD:", count(inpmask<-mask_tiny)
           call inp_step()
        enddo
        call hpd%convert2ring()
        call hpd%write(inp_file_name(freq_names(ifreq)))
     enddo
  endif
  
contains


  subroutine inp_step()
    COOP_INT::i, j
    COOP_SINGLE:: weight, sum_qu(2)
    do j=1, npix_inp
       if(inpmask(inp_pix(j)) < -mask_tiny) then
          call neighbours_nest(hpd%nside, inp_pix(j), list, nneigh)
          weight = 0.
          sum_qu = 0.
          do i=1, nneigh
             if(inpmask(list(i)) > mask_tiny )then
                weight = weight + 1.
                sum_qu = sum_qu + hpd%map(list(i), 2:3)
             endif
          enddo
          if(weight > 0.)then
             hpd%map(inp_pix(j), 2:3) = sum_qu/weight
             inpmask(inp_pix(j)) = 1.
          endif
       endif
    enddo
  end subroutine inp_step

  subroutine do_projection(beta_dust, beta_sync)
    COOP_SINGLE::beta_dust, beta_sync
    COOP_SINGLE, dimension(num_freqs) :: vec_dust, vec_sync
    COOP_SINGLE, dimension(num_freqs, 3) :: qvecs, uvecs
    COOP_INT::ifreq, i
    do ifreq = 1, num_freqs
       vec_dust(ifreq) =  (Freqs(ifreq)/Freq_REF_D)**beta_dust * dust_cmb_ratio(Freqs(ifreq), Freq_REF_D)
       vec_sync(ifreq) = (Freqs(ifreq)/Freq_REF_S)**beta_sync * convert_i2cmb(Freqs(ifreq), Freq_REF_S)
    enddo
    !$omp parallel do private(i, qvecs, uvecs)
    do i = 1,npix_mask
       qvecs(:, 1) = vec_dust / Q_std(:, i)
       qvecs(:, 2) = vec_sync / Q_std(:, i)
       qvecs(:, 3) = 1. / Q_std(:, i)
       uvecs(:, 1) = vec_dust / U_std(:, i)
       uvecs(:, 2) = vec_sync / U_std(:, i)
       uvecs(:, 3) = 1. / U_std(:, i)
       call six_vec_proj(Q_data(:, i), U_data(:, i), qvecs, uvecs, coefs(:, i), res(i))
    enddo
    !$omp end parallel do
  end subroutine do_projection

  subroutine inpaint_scale(beta_dust, beta_sync)
    COOP_SINGLE::beta_dust, beta_sync, c(3)
    COOP_SINGLE, dimension(num_freqs) :: vec_dust, vec_sync
    COOP_INT::ifreq, i, ipix, irand
    do ifreq = 1, num_freqs
       vec_dust(ifreq) =  (Freqs(ifreq)/Freq_REF_D)**beta_dust * dust_cmb_ratio(Freqs(ifreq), Freq_REF_D)
       vec_sync(ifreq) = (Freqs(ifreq)/Freq_REF_S)**beta_sync * convert_i2cmb(Freqs(ifreq), Freq_REF_S)
    enddo
    !$omp parallel do private(i, c)
    do i = 1,npix_mask
       c = coefs(:, i)
       do while(abs(c(1)) > coefs_cut(1) )
          c(1) = coefs(1, coop_random_index(npix_mask))
       enddo
       do while(abs(c(2)) > coefs_cut(2) )
          c(2) = coefs(2, coop_random_index(npix_mask))
       enddo
       if(res(i) > res_cut)then
          Q_data(:, i) = (c(1)*vec_dust + c(2)*vec_sync + c(3)) + Q_std(:, i)*coop_random_Gaussian_vector(num_freqs)
          U_data(:, i) = (c(1)*vec_dust + c(2)*vec_sync + c(3)) + U_std(:, i)*coop_random_Gaussian_vector(num_freqs)
       else
          Q_data(:, i) = Q_data(:, i)*Q_std(:, i)
          U_data(:, i) = U_data(:, i)*U_std(:, i)          
          if(abs(c(1)-coefs(1, i))>1.e-10 )then
             Q_data(:, i) = Q_data(:, i) + (c(1)-coefs(1, i))*vec_dust
             U_data(:, i) = U_data(:, i) + (c(1)-coefs(1, i))*vec_dust             
          endif
          if(abs(c(2)-coefs(2, i))>1.e-10 )then
             Q_data(:, i) = Q_data(:, i) + (c(2)-coefs(2, i))*vec_sync
             U_data(:, i) = U_data(:, i) + (c(2)-coefs(2, i))*vec_sync             
          endif
       endif
    enddo
    !$omp end parallel do
  end subroutine inpaint_scale
  
end program test
