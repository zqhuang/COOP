program shells
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
  type(coop_cosmology_firstorder)::cosmology
  COOP_STRING::fnl_option = "g"
  type(coop_function)::deltaN
  COOP_REAL,parameter::sigma_zeta = 4.6d-5
  COOP_INT,parameter::lmax = 50
  COOP_INT,parameter::nside = 256
  COOP_REAL,parameter::fwhm_arcmin = 40.d0
  COOP_UNKNOWN_STRING,parameter::dir = "zetaproj/"
  COOP_UNKNOWN_STRING,parameter::prefixtrans = dir//"testmap"
  COOP_UNKNOWN_STRING,parameter::prefix3d = dir//"sim3d"
  COOP_REAL::sigma_chi = 1.d-7
  COOP_REAL::mean_chi = 4.3d-6
  COOP_REAL::gp_width = 3.
  COOP_REAL::gp_A = 2.e-4
  COOP_REAL::gp_mean = 12.
  COOP_SINGLE,parameter::resolution_hMpc = 10.!!In unit of h^{-1}Mpc
  COOP_SINGLE,parameter::resolution = resolution_hMpc/3000. 
  type(coop_healpix_maps)::hm, zeta, hmng
  COOP_INT:: i, j, ix, ipix, itrial
  COOP_INT,parameter::n_x = 256, n_rms = 128
  COOP_SINGLE,dimension(:),allocatable::x_table, rms_table
  COOP_SINGLE,dimension(:,:),allocatable::fnl_table
  COOP_SINGLE::dx, xmin, rmsmin, drms, diff

  if(iargc().gt.0)then
     fnl_option = trim(coop_InputArgs(1))
  else
     write(*,*) "./Test i/b/g"
     stop
  endif
  call coop_file_load_function("deltaN-LUT-1.875", 1, 2, deltaN, .false.)
  call cosmology%Set_Planck_bestfit()
  call cosmology%compute_source(0)
  call coop_feedback("source done")
  if(trim(fnl_option).ne."i")then
     allocate(x_table(n_x), rms_table(n_rms), fnl_table(n_x, n_rms))
     call coop_set_uniform(n_x, x_table, -20., 20.)
     call coop_set_uniform(n_rms, rms_table, 0.,  sqrt(3.*log(0.05/resolution)))
     xmin = x_table(1)
     dx = x_table(2) - x_table(1)
     rmsmin = rms_table(1)
     drms = rms_table(2) - rms_table(1)
     do ix = 1, n_x
        do ipix = 1, n_rms
           fnl_table(ix, ipix) = fnl_smooth(x_table(ix), rms_table(ipix))
        enddo
     enddo
  endif
  call coop_random_init()
  select case(trim(fnl_option))
  case("i")
     call coop_zeta3d_generate_cmb(hm, zeta, cosmology, identity, lmax, nside, prefix3d = prefix3d, prefixtrans = prefixtrans, prefixmap=trim(dir)//trim(prefixmap), fwhm_arcmin= fwhm_arcmin, writefile = .true.)
  case("b")
     call coop_zeta3d_generate_cmb(hm, zeta, cosmology, fnl, lmax, nside, prefix3d=prefix3d, prefixtrans=prefixtrans, prefixmap=trim(dir)//"gp_meanchi"//COOP_STR_OF(nint(mean_chi*1.e7))//"_sigmachi"//COOP_STR_OF(nint(sigma_chi*1.e7)) , fwhm_arcmin = fwhm_arcmin, writefile=.true.)
  case("g")
     call coop_zeta3d_generate_cmb(hm, zeta, cosmology, fnl, lmax, nside, prefix3d=prefix3d, prefixtrans=prefixtrans, prefixmap = trim(dir)//"gs_amp"//COOP_STR_OF(nint(1.e6*gp_A))//"_mean"//COOP_STR_OF(nint(gp_mean))//"_width"//COOP_STR_OF(nint(gp_width))  , fwhm_arcmin = fwhm_arcmin, writefile=.true.)
  end select

contains

!!define the fnl function of zeta; my convention of \zeta is:  \zeta = - \delta N (this determines the relative sign between \zeta and \Phi when I set up initial conditions in the Boltzmann code).

  subroutine identity(x, pixsize)
    COOP_SINGLE x, pixsize
    return
  end subroutine identity


  function fnl_raw(x) result(y)
    COOP_SINGLE x, y
    select case(trim(fnl_option))
    case("i")
       y = x
    case("b")
       y =  deltaN%eval(abs(mean_chi+x*sigma_chi))  
    case("g")
       y = gp_A* exp(-(abs(x) - gp_mean)**2/(2.d0*gp_width**2)) 
    end select
    return
  end function fnl_raw

  function fnl_smooth(x, rms) result(y)
    COOP_INT::n
    COOP_SINGLE x, y
    COOP_SINGLE rms, dx, dw, w
    COOP_INT::i
    n = max(1, ceiling(rms*2./0.2))
    y = fnl_raw(x)
    dx = rms*2./n
    w = 1.d0
    do i = 1, n
       dw = exp(-(2.*i/n)**2/2.)
       y = y + (fnl_raw(x+i*dx)+fnl_raw(x-i*dx))*dw
       w = w + 2.d0*dw
    enddo
    y= y/w
  end function fnl_smooth

  subroutine fnl(x, pixsize) 
    COOP_SINGLE x
    COOP_SINGLE pixsize, rms
    COOP_INT::ix, irms
    COOP_REAL::rx, rrms
    x = x/sigma_zeta
    if(pixsize .lt. resolution)then
       x = fnl_raw(x)
       return
    else
       rms = sqrt(3.*log(pixsize/resolution))
    endif
    rx = (x-xmin)/dx+1.d0
    rrms = (rms - rmsmin)/drms+1.d0+1.d-10
    ix = floor(rx)
    irms =floor(rrms)
    if(ix .lt. 1 .or. ix .ge. n_x)then
       stop "x overflow"
    endif
    if(irms .lt. 1 .or. irms .gt. n_rms)then
       stop "rms overflow"
    endif
    rrms = rrms - irms
    rx = rx - ix
    x = (fnl_table(ix, irms)*(1.d0-rrms)+fnl_table(ix, irms+1)*rrms)*(1.d0-rx) +  (fnl_table(ix+1, irms)*(1.d0-rrms)+fnl_table(ix+1, irms+1)*rrms)*rx
  end subroutine fnl

end program shells
