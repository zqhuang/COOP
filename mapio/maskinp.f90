module disk_inpaint
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
  COOP_INT,parameter::nside = 512
  COOP_INT,parameter::npix_tot = nside**2*12
  COOP_REAL,parameter::di_ell = 90.d0
  COOP_INT,parameter::di_mmax = 2
  COOP_INT,parameter::di_nt = 50
  COOP_INT,parameter::di_nphi = 12
  COOP_REAL::di_phil(di_nphi), di_sin2phil(di_nphi), di_cos2phil(di_nphi)
  COOP_REAL::di_smax, di_rhomax, di_rmax 
  COOP_REAL::di_jpz(0:di_mmax*2, di_nt)
  COOP_REAL::di_mu(0:2*di_mmax)
  COOP_INT::nlist_max

contains

  !! rho flat coordinate, r curved coordinate
  function di_rho_of_r(r) result(rho)
    COOP_REAL::r, rho
    rho = 2.*sin(r/2.)
  end function di_rho_of_r

  function di_r_of_rho(rho) result(r)
    COOP_REAL::rho, r
    r = 2.*asin(rho/2.)
  end function di_r_of_rho
  

  subroutine di_initialize(smax)
    COOP_REAL,optional::smax
    logical::success
    COOP_INT::n, i
    call coop_import_matrix("jnp_zeros.txt", di_jpz, 2*di_mmax+1, di_nt, success)
    if(.not. success) stop 'Unknown error when loading jnp_zeros.txt'
    if(present(smax))then
       di_smax = max(di_jpz(2*di_mmax, 1), smax)
    else
       di_smax = di_jpz(2*di_mmax, 1)       
    endif
    di_rhomax = di_smax/di_ell
    di_rmax = di_r_of_rho(di_rhomax)
    do n = 0, 2*di_mmax
       i = 0
       do while(di_jpz(n, i)<= di_smax)
          i = i+1
          if(i >= di_nt)stop 'Error in di_initialize: nt too small'
       enddo
       di_mu(n) = di_jpz(n, i-1)
    enddo
    nlist_max = min(ceiling(di_rhomax**2*npix_tot/3.9)+1024, npix_tot)
    print*, "disk size: ", di_rmax/coop_SI_degree, " deg, "
    print*, "flat correction: ", di_rhomax/di_rmax
    print*, "disk pixels <= ", nlist_max
    do i=1, di_nphi
       di_phil(i) = coop_2pi * (i-1.d0)/di_nphi
       di_sin2phil(i) = sin(2*di_phil(i))
       di_cos2phil(i) = cos(2*di_phil(i))       
    enddo
  end subroutine di_initialize

  function di_cos(s, theta) result(c)
    COOP_REAL::s, theta, c, sig
    COOP_INT::m
    c = bessel_Jn(0, s)
    sig = 2.d0
    do m=1, di_mmax
       sig = -sig
       if(s .le. di_mu(2*m)) c = c + sig*bessel_Jn(2*m, s)*cos(2*m*theta)
    enddo
  end function di_cos

  function di_sin(s, theta) result(c)
    COOP_REAL::s, theta, c, sig
    COOP_INT::m
    c = 0.
    sig = 2.d0
    do m=1, di_mmax
       sig = -sig
       if(s .le. di_mu(2*m-1)) c = c + sig*bessel_Jn(2*m-1, s)*cos((2*m-1)*theta)
    enddo
  end function di_sin
  
end module disk_inpaint

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
  use disk_inpaint
  
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

  subroutine symmat2_inv(a)
    COOP_REAL::a(2,2), det
    det =  a(1,1)*a(2,2) - a(2,1)**2
    a = Reshape( (/ a(2,2), -a(2, 1), -a(2, 1), a(1, 1) /), (/ 2, 2 /) ) &
         / det
  end subroutine symmat2_inv

  
  
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



program MaskAndInpaint
  use disk_inpaint
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
  COOP_SINGLE,parameter::mask_tiny = 1.e-3, threshold = 0.999
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
  if(npix_inp > 0)then
     write(*, *) "Inpainting ", npix_inp, " pixels."
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
     call  di_initialize()             
     do ifreq = 1, num_freqs
        inpmask = mask%map(:, 1)
        call hpd%read(data_file_name(freq_names(ifreq)), nmaps_wanted = 3)
        write(*,*) "inpainting file: "//trim(data_file_name(freq_names(ifreq)))        
        call hpd%convert2nested()
        call do_inpaint()
        call hpd%convert2ring()
        call hpd%write(inp_file_name(freq_names(ifreq)))
        write(*,*) "file: "//trim(data_file_name(freq_names(ifreq)))//" has been inpainted and saved as "//trim(inp_file_name(freq_names(ifreq)))
     enddo
  endif
  
contains


  subroutine do_inpaint()
    COOP_INT::i, j, iphi, listpix(nlist_max), nlist
    COOP_REAL::rho, phi, theta, s, sin2theta, cos2theta, QUrot
    type(coop_healpix_disc)::disc
    COOP_REAL::B_ell_Sum(di_nphi), B_ell_clean(di_nphi), B_ell_Re(di_nphi), mat(2, 2), b(2)
    do j=1, npix_inp
       B_ell_Re = 0.d0
       call hpd%get_disc(inp_pix(j), disc)
       call hpd%query_disc(inp_pix(j), di_rmax, listpix, nlist)
       if(any(mask%map(listpix(1:nlist), 1) .eq. 0.))then  !!do not inpaint points around the boundaries
!          print*, "skipping point near the boundary"
          cycle
       endif
       do i=1, nlist
          call disc%pix2ang(listpix(i), rho, phi)
          s = di_ell*rho
          do iphi = 1, di_nphi                 
             theta = phi - di_phil(iphi)
             QUrot = hpd%map(listpix(i), 2) * di_sin2phil(iphi) - hpd%map(listpix(i), 3)*di_cos2phil(iphi)
             B_ell_Re(iphi) = B_ell_Re(iphi) + di_cos(s, theta)  * QUrot  !!real part
          enddo
       enddo
!       B_ell_Sum =  hpd%map(inp_pix(j), 2)*di_sin2phil - hpd%map(inp_pix(j), 3) * di_cos2phil + B_ell_Re
       mat(1, 1) = sum(di_sin2phil**2)
       mat(2, 2) = sum(di_cos2phil**2)
       mat(2, 1) = -sum(di_sin2phil*di_cos2phil)
       call symmat2_inv(mat)
       b(1) = -sum(di_sin2phil*B_ell_Re)
       b(2) = sum(di_cos2phil*B_ell_Re)     
       hpd%map(inp_pix(j), 2:3) = matmul(mat, b)
 !      B_ell_clean = hpd%map(inp_pix(j), 2)*di_sin2phil - hpd%map(inp_pix(j), 3) * di_cos2phil + B_ell_Re
 !      print*, sum(B_ell_clean**2)/sum(B_ell_Sum**2)
    enddo
  end subroutine do_inpaint

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

end program MaskAndInpaint
