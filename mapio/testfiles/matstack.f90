program stackth
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
  integer, parameter::lmin = 2, lmax=2000, lmax_vary = 30
  logical:: do_orient = .false.
  COOP_REAL,parameter::fsky = 0.6d0
  COOP_REAL, parameter:: r_degree = 2.d0
  COOP_REAL, parameter:: width = 2.d0*sin(r_degree*coop_SI_degree/2.d0)
  COOP_INT,parameter:: n=36
  COOP_REAL, parameter:: dr = width/n
  logical,parameter::flat = .false. !!use nonflat is actually faster
  COOP_STRING::Im0_file, Im1_file
  COOP_STRING::clfile
  COOP_REAL::nu !! threshold
  COOP_REAL::fwhm !!fwhm in arcmin
  COOP_STRING::line
  integer l, il, i, j, iomega, m, dsize
  type(coop_arguments)::args
  type(coop_healpix_patch)::patchI
  COOP_REAL::cls(lmin:lmax), ell(lmin:lmax), lsq(lmin:lmax), sigma, sigma2, sigma0, sigma1, omega, weights(4), frI(0:2, 0:n*3/2), phi, romega, r(0:n*3/2), kr, cr(0:1, 0:2), cls_begin(lmin:lmax), sql(lmin:lmax), cls_prev(lmin:lmax), chisq, chisq_prev, stepsize, res_prev, res_now, rel_weight
  type(coop_file)::fp
  COOP_REAL,dimension(:,:),allocatable::Pl0, Pl2
  COOP_REAL,dimension(:),allocatable::fr_data, fr_theory, vol
  COOP_INT::loop, big_loop
  
  Im0_file = coop_InputArgs(4)  
  if(trim(Im0_file).eq."")then
     write(*,*) "Syntax:"
     write(*,*) "./MatStack clfile nu fwhm_arcmin Im0_file [Im1file]"
     stop
  endif
  clfile = trim(coop_InputArgs(1))
  line = coop_InputArgs(2)
  read(line, *) nu
  line = coop_InputArgs(3)
  read(line, *) fwhm
  Im1_file = trim(coop_InputArgs(5))
  do_orient = (trim(Im1_file) .ne. "")
  
  allocate(Pl0(0:lmax, 0:n*3/2), Pl2(0:lmax, 0:n*3/2))
  if(flat)then  !!get Plms
     do i= 0, n*3/2
        !$omp parallel do private(l, kr)
        do l = 0, lmax
           kr = sqrt(l*(l+1.d0))*(i*dr)
           Pl0(l, i) = coop_bessj(0, kr)
           Pl2(l, i) = coop_bessj(2, kr)
        enddo
        !$omp end parallel do
     enddo
  else
     do i=0, n*3/2
        call coop_get_normalized_Plm_array(m=0, lmax=lmax, x = 1.d0-(dr*i)**2/2.d0, Plms = Pl0(0:lmax, i))
        call coop_get_normalized_Plm_array(m=2, lmax=lmax, x = 1.d0-(dr*i)**2/2.d0, Plms = Pl2(0:lmax, i))
     enddo
  endif

  if(do_orient)then
     dsize = (n+1)*2
  else
     dsize = n+1
  endif
  allocate(fr_data(dsize), fr_theory(dsize), vol(dsize))
  
  call coop_random_init()

  !!load cls
  sigma = fwhm*coop_sigma_by_fwhm*coop_SI_arcmin

  call fp%open(clfile, "r")
  do l=lmin, lmax
     read(fp%unit, *) il, cls_begin(l)
     ell(l)  = l
     lsq(l) = l*(l+1.d0)
     sql(l) = sqrt(2.d0*l+1.d0)
     cls_begin(l) = cls_begin(l)*(coop_2pi*exp(-lsq(l)*sigma**2)/lsq(l))
     if(il.ne.l) stop "cl file broken"
  enddo
  call fp%close()
  cls = cls_begin
  
  call fp%open(Im0_file, "r")
  do i=0, n
     read(fp%unit, *) r(i), fr_data(i+1)
     vol(i+1) = i * coop_2pi
  enddo
  call fp%close()
  if(do_orient)then
     call fp%open(Im1_file, "r")
     do i=0, n
        read(fp%unit, *) r(i), fr_data(i+2+n)
        vol(i+2+n) = i * coop_2pi
     enddo
     call fp%close()
  end if

  call patchI%init("I", n, dr)

  patchI%nstack = 1.d0
  patchI%nstack_raw = 1.d0
  stepsize = 0.15d0
  rel_weight = 1.d0
  do big_loop = 1, 1
     call get_chisq(chisq, res_now)
     print*, 0,  chisq
     do loop = 1, 10000
        res_prev = res_now
        chisq_prev= chisq
        il = coop_random_index(lmax_vary-lmin+1)+lmin-1
        cls_prev(il) = cls(il)
        cls(il) = cls(il) * exp(coop_random_Gaussian()/sql(il)*stepsize)
        call get_chisq(chisq, res_now)
        if(chisq .gt. chisq_prev)then
           cls(il) = cls_prev(il)
           chisq = chisq_prev
           res_now = res_prev
        else
           print*, big_loop, loop, res_now,  chisq
        endif
     enddo
     rel_weight = rel_weight*0.95d0
  enddo
  
  call fp%open("final_cls.dat", "w")
  do l = lmin, lmax
     write(fp%unit, "(I5, 4E16.7)") l, cls(l)/(coop_2pi*exp(-lsq(l)*sigma**2)/lsq(l)), cls_begin(l)/(coop_2pi*exp(-lsq(l)*sigma**2)/lsq(l)), 0. , 0.
  enddo
  call fp%close()
  
contains

  subroutine get_chisq(chi2, res)
    COOP_REAL chi2, res
    sigma0 = sqrt(sum(Cls*(ell+0.5d0))/coop_2pi)
    sigma1 = sqrt(sum(Cls*(ell+0.5d0)*(ell*(ell+1.d0)))/coop_2pi)
    sigma2 = sqrt(sum(Cls*(ell+0.5d0)*(ell*(ell+1.d0))**2)/coop_2pi)
    call coop_gaussian_npeak_set_args(args, 2, sigma0, sigma1, sigma2)

    
    if(do_orient)then
       call coop_gaussian_get_oriented_stacking_weights(nu, args, weights)
    else
       call coop_gaussian_get_nonoriented_stacking_weights(nu, args, weights)
    endif
    do i=0, n*3/2
       cr = 0.d0       
       do l= lmin, lmax
          cr(0, 0) = cr(0, 0) + (l+0.5d0)*Pl0(l,i)*cls(l)
          cr(1, 0) = cr(1, 0) - (l+0.5d0)*Pl0(l,i)*cls(l)*lsq(l)
          if(do_orient)then
             cr(0, 1) = cr(0, 1) + (l+0.5d0)*Pl2(l,i)*cls(l)
             cr(1, 1) = cr(1, 1) - (l+0.5d0)*Pl2(l,i)*cls(l)*lsq(l)
          endif
       enddo
       cr  = cr/coop_2pi
       call coop_gaussian_radial_modes_I(weights, cr, frI(:, i))
    enddo
    do i=-n, n
       do j=-n, n
          if(i.eq.0.and.j.eq.0)then
             patchI%image(i,j,1) = frI(0,0)
          else
             phi = atan2(dble(j), dble(i))
             omega = sqrt(dble(i)**2+dble(j)**2)
             iomega = floor(omega)
             romega = omega - iomega
             patchI%image(i,j,1) = (frI(0, iomega)+ frI(1, iomega)*cos(2*phi) + frI(2, iomega)*cos(4*phi) ) *(1.d0-romega)+(frI(0, iomega+1)+ frI(1, iomega+1)*cos(2*phi) + frI(2, iomega+1)*cos(4*phi) ) * romega
          endif
       enddo
    enddo  
  !!test the integrator
    call patchI%get_all_radial_profiles()
    fr_theory(1:n+1) = patchI%fr(:, 0, 1)
    if(do_orient) fr_theory(n+2:2*n+2) = patchI%fr(:, 1, 1)
    res = sum((fr_data - fr_theory)**2*vol) 
    chi2  = res + rel_weight*fsky*sum((log(cls_begin(lmin:lmax_vary)/cls(lmin:lmax_vary)) - 1.d0 + cls(lmin:lmax_vary)/cls_begin(lmin:lmax_vary))*(2.d0*ell(lmin:lmax_vary)+1.d0))
    
  end subroutine get_chisq
  
end program stackth
