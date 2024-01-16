program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  use coop_fitsio_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::hm, mask, hmsaved
  COOP_STRING::prefix 
  COOP_INT,parameter::lmax_mask = 200
  COOP_INT::lmin, lmax, l
  COOP_REAL,dimension(:,:),allocatable::kernel
  COOP_REAL,dimension(:),allocatable::Cl_pseudo, Cl, SqrtCls
  COOP_INT,parameter::nnu = 25, nfwhm = 25
  COOP_REAL,parameter::minfwhm = 10.d0
  COOP_INT::l_lower
  COOP_REAL::dfwhm, dnu, maxIntensity
  COOP_STRING::mapfile, maskfile
  COOP_REAL::nu, na_arr(0:nfwhm, 0:nnu), nuarr(0:nnu)
  COOP_REAL:: fwhm_arcmin
  COOP_INT::inu, ifwhm, na
  COOP_SINGLE::rms, mean
  type(coop_asy)::figure
  type(coop_function)::mapping
  logical::cold, gaussianize, lognorm, sim
  call coop_get_command_line_argument(key = 'prefix', arg = prefix, default="Commander_Dust")
  call coop_get_command_line_argument(key = 'map', arg = mapfile)
  call hm%open(mapfile, nmaps_wanted = 1)
  call coop_get_command_line_argument(key = 'mask', arg = maskfile)
  call mask%open(maskfile, nmaps_wanted = 1)
  call coop_get_command_line_argument(key = 'lognorm', arg = lognorm, default = .false.)
  call coop_get_command_line_argument(key = 'gaussianize', arg = gaussianize, default = .false.)
  call coop_get_command_line_argument(key = 'sim', arg = sim, default = .false.)
  call coop_get_command_line_argument(key = 'fwhm', arg = fwhm_arcmin, default=0.d0)
  call coop_get_command_line_argument(key = 'nu', arg = nu, default = -1.1d30)
  call coop_get_command_line_argument(key = 'cold', arg = cold, default = .false.)

  call coop_get_command_line_argument(key = 'l_lower', arg = l_lower, default = 30)

  maxIntensity = maxval(abs(hm%map(:,1)*mask%map(:,1)))
  hm%map(:,1) = min(hm%map(:,1), maxIntensity)
  if(lognorm)then
     hm%map(:,1) = log(max(hm%map(:,1), maxIntensity*1.d-8))
     hm%map(:,1) = hm%map(:,1) - sum(hm%map(:,1))/hm%npix
     call hm%write("lognorm.fits")
  elseif(gaussianize)then
     call hm%gaussianize(mask, 1, mapping)
     call hm%write("gaussianize.fits")
  endif
  if(sim)then
     write(*,*) "Doing Gaussian sim with reconstructed power spectrum."
     call coop_random_init()
     coop_healpix_want_cls = .true.
     lmin = l_lower - 10
     lmax = min(nint(4./(max(fwhm_arcmin, minfwhm)*coop_SI_arcmin)), 1500)
     allocate(kernel(lmin:lmax, lmin:lmax), Cl_pseudo(lmin:lmax), Cl(lmin:lmax), SqrtCls(0:lmax))
     call mask%convert2ring()
     call hm%convert2ring()
     hm%map(:,1) = hm%map(:,1)*mask%map(:,1)
     call mask%get_mask_kernel(lmax_mask = lmax_mask, lmin = lmin, lmax = lmax, kernel = kernel)
     call hm%map2alm(lmax = lmax)
     call hm%get_cls()

     Cl_pseudo(lmin:lmax) = hm%cl(lmin:lmax,1)
     call coop_pseudoCl2Cl(kernel = kernel, lmin = lmin, lmax = lmax, cl_pseudo = cl_pseudo, cl = cl)
     sqrtCls = 0.d0
     sqrtCls(lmin:lmax) = sqrt(max(cl(lmin:lmax), 0.d0))
!!$     hm%cl = 0.
!!$     hm%cl(lmin:lmax,1) = cl(lmin:lmax)
     call hm%simulate_Tmaps(nside = hm%nside, lmax = lmax, lmin = lmin, sqrtCls = sqrtCls)
  endif


  if(fwhm_arcmin .gt. 0.d0)then
     call hm%smooth(fwhm = fwhm_arcmin*coop_SI_arcmin, l_lower = l_lower)
     call hm%write("tmp.fits")
     call system("rm -f tmp.gif")
     call system("map2gif -inp tmp.fits -out tmp.gif -bar T")
     call system("rm -f tmp.fits")
     call hm%convert2nested()
     call mask%convert2nested()
     mean  = sum(hm%map(:,1), mask%map(:,1).gt. 0.5)/count(mask%map(:,1).gt. 0.5)
     rms = sqrt(sum(hm%map(:,1)**2, mask%map(:,1).gt. 0.5)/count(mask%map(:,1).gt. 0.5))       
     call get_num_areas()
     print*, "num areas = ", na
  else
     dnu = 3.d0/nnu
     dfwhm = nint((60.d0-minfwhm)/nfwhm)
     hmsaved = hm
     if(sim)then
        if(lognorm)then
           call figure%open(trim(prefix)//"_num_hotspots_lognorm__GaussianSim.txt")
        elseif(gaussianize)then
           call figure%open(trim(prefix)//"_num_hotspots_equalize_GaussianSim.txt")
        else
           call figure%open(trim(prefix)//"_num_hotspots_GaussianSim.txt")
        endif
     else
        if(lognorm)then
           call figure%open(trim(prefix)//"_num_hotspots_lognorm.txt")
        elseif(gaussianize)then
           call figure%open(trim(prefix)//"_num_hotspots_equalize.txt")
        else
           call figure%open(trim(prefix)//"_num_hotspots.txt")
        endif
     endif
     call figure%init(xlabel = "FWHM (arcmin)", ylabel = "$\nu$")
     do ifwhm = 0, nfwhm
        fwhm_arcmin = ifwhm * dfwhm + minfwhm
        hm = hmsaved
        call hm%smooth(fwhm = fwhm_arcmin*coop_SI_arcmin, l_lower = l_lower)
        call hm%convert2nested()
        call mask%convert2nested()
        mean  = sum(hm%map(:,1), mask%map(:,1).gt. 0.5)/count(mask%map(:,1).gt. 0.5)
        rms = sqrt(sum(hm%map(:,1)**2, mask%map(:,1).gt. 0.5)/count(mask%map(:,1).gt. 0.5))       
        do inu = 0, nnu
           nu = inu*dnu
           nuarr(i) = nu
           call get_num_areas()
           na_arr(ifwhm, inu) = log10(dble(na))
           print*, fwhm_arcmin, nu, na
        enddo
     enddo
     call figure%density(na_arr, xmin = minfwhm, xmax = minfwhm+dfwhm*nfwhm, ymin = 0.d0, ymax = dnu*nnu, label = "$\log_{10}N_{I>\nu\sigma_I}$", zmin = 1.2d0, zmax = 3.8d0)
     call figure%close()
  endif
  call hm%free()
  call mask%free()
  
contains

  subroutine get_num_areas()
    COOP_SINGLE::threshold
    threshold = rms*nu + mean
    na = hm%num_areas(imap = 1, mask = mask, threshold = threshold, cold = cold)
  end subroutine get_num_areas
end program test
