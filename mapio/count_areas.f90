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
  COOP_INT,parameter::lmax_mask = 200, highpass_l = 30
  COOP_INT,parameter::nnu = 12, nfwhm = 2, lmin = highpass_l - 10, lmax = 1500
  COOP_REAL,parameter::minfwhm = 10.d0, maxfwhm = 30.d0, minnu = 0.d0, maxnu = 3.d0
 
  COOP_INT::l
  COOP_REAL,dimension(:,:),allocatable::kernel
  COOP_REAL,dimension(:),allocatable::Cl_pseudo, Cl, SqrtCls, Cl_mask
  COOP_INT::l_lower
  COOP_REAL::dfwhm, dnu, maxIntensity, minIntensity
  COOP_STRING::mapfile, maskfile
  COOP_REAL::nu, Narr(0:nnu), nuarr(0:nnu), NGaussian(0:nnu)
  COOP_REAL:: fwhm_arcmin
  COOP_INT::inu, ifwhm, na
  COOP_SINGLE::rms, mean
  type(coop_asy)::figure
  type(coop_function)::mapping
  logical::cold, gaussianize, lognorm, sim
  if(iargc().lt.2)then
     write(*,*) "Syntax:"
     write(*,*) "./NAREAS -prefix PREFIX -map MAP -mask MASK -log T/F -eq T/F -cold T/F "
     stop
  endif
  coop_healpix_want_cls = .true.
  allocate(cl_mask(0:lmax_mask), cl(lmin:lmax), cl_pseudo(lmin:lmax), kernel(lmin:lmax, lmin:lmax), sqrtcls(0:lmax))

  call coop_get_command_line_argument(key = 'prefix', arg = prefix)
  call coop_get_command_line_argument(key = 'map', arg = mapfile)
  call hm%open(mapfile, nmaps_wanted = 1)
  call coop_get_command_line_argument(key = 'mask', arg = maskfile)
  call mask%open(maskfile, nmaps_wanted = 1)
  call coop_get_command_line_argument(key = 'log', arg = lognorm, default = .false.)
  call coop_get_command_line_argument(key = 'eq', arg = gaussianize, default = .false.)
  call coop_get_command_line_argument(key = 'cold', arg = cold, default = .false.)

  call coop_get_command_line_argument(key = 'l_lower', arg = l_lower, default = highpass_l)
  call hm%convert2ring()
  call mask%convert2ring()
  maxIntensity = maxval(abs(hm%map(:,1)*mask%map(:,1)))
  minIntensity = minval(abs(hm%map(:,1)*mask%map(:,1)))
  hm%map(:,1) = max(min(hm%map(:,1), maxIntensity), minIntensity)
  if(lognorm)then
     hm%map(:,1) = log(max(hm%map(:,1), maxIntensity*1.d-8))
  elseif(gaussianize)then
     call hm%gaussianize(mask, 1, mapping)
  endif
  hmsaved = hm


  hm%map(:,1) = hm%map(:,1)*mask%map(:,1)
  call mask%map2alm(lmax = lmax_mask)
  cl_mask = mask%cl(0:lmax_mask, 1)
  call coop_pseudoCl_get_kernel(lmax_mask = lmax_mask, Cl_mask = Cl_mask, lmin = lmin, lmax = lmax, kernel = kernel)
  call hm%map2alm(lmax = lmax)
  cl_pseudo = hm%cl(lmin:lmax,1)
  call coop_pseudoCl2Cl(kernel = kernel, lmin = lmin, lmax = lmax, cl_pseudo = cl_pseudo, cl = cl)

  if(lognorm)then
     call figure%open(trim(prefix)//"_num_hotspots_lognorm.txt")
  elseif(gaussianize)then
     call figure%open(trim(prefix)//"_num_hotspots_equalize.txt")
  else
     call figure%open(trim(prefix)//"_num_hotspots.txt")
  endif
  call figure%init(xlabel = "$\nu$", ylabel = "$N/N_{\rm Gaussian}$")

  sqrtCls = 0.d0
  dnu = (maxnu-minnu)/nnu
  dfwhm = nint((maxfwhm-minfwhm)/nfwhm)

  do ifwhm = 0, nfwhm
     fwhm_arcmin = ifwhm * dfwhm + minfwhm
     hm = hmsaved
     call hm%smooth(fwhm = fwhm_arcmin*coop_SI_arcmin, l_lower = l_lower)
     call hm%convert2nested()
     call mask%convert2nested()
     mean  = sum(hm%map(:,1), mask%map(:,1).gt. 0.5)/count(mask%map(:,1).gt. 0.5)
     rms = sqrt(sum(hm%map(:,1)**2, mask%map(:,1).gt. 0.5)/count(mask%map(:,1).gt. 0.5))       
     do inu = 0, nnu
        nu = minnu + inu*dnu
        call get_num_areas()
        Narr(inu) = na
        print*, fwhm_arcmin, nu, na
        nuarr(inu) = nu
     enddo
     do l = lmin, lmax
        sqrtCls(l) = sqrt(max(cl(l), 0.d0))*coop_gaussian_filter(fwhm_arcmin, l)
     enddo
     call hm%simulate_Tmaps(nside = hm%nside, lmax = lmax, lmin = lmin, sqrtCls = sqrtCls)

     call hm%convert2nested()
     call mask%convert2nested()
     mean  = sum(hm%map(:,1), mask%map(:,1).gt. 0.5)/count(mask%map(:,1).gt. 0.5)
     rms = sqrt(sum(hm%map(:,1)**2, mask%map(:,1).gt. 0.5)/count(mask%map(:,1).gt. 0.5))       
     do inu = 0, nnu
        nu = minnu + inu*dnu
        call get_num_areas()
        NGaussian(inu) = na
        print*, fwhm_arcmin, nu, na
     enddo
     call figure%plot(nuarr, Narr/NGaussian, color = figure%color(ifwhm+1), linetype = figure%linetype(ifwhm+1), linewidth = 2., legend = "FWHM "//COOP_STR_OF(nint(fwhm_arcmin))//"'")
  enddo
  call figure%legend(0.05, 0.95)
  call figure%close()

  call hm%free()
  call mask%free()
  
contains

  subroutine get_num_areas()
    COOP_SINGLE::threshold
    threshold = rms*nu + mean
    na = hm%num_areas(imap = 1, mask = mask, threshold = threshold, cold = cold)
  end subroutine get_num_areas
end program test
