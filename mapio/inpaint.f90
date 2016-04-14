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
  COOP_INT::lcut, npix, nlm, nmaps, i, j, lwork, liwork, nside, l, m, imap
  type(coop_healpix_maps)::map, mask, sim, fullmap, fullmask
  logical::want_fluc
  COOP_STRING::genre, fmap, fmask
  COOP_INT,dimension(:),allocatable::listpix
  COOP_REAL,dimension(:,:),allocatable::mat, yvec
  COOP_REAL,dimension(:), allocatable::s, work
  COOP_INT,dimension(:), allocatable::iwork
  COOP_REAL::rcond = 1.d-8
  COOP_INT::rank, info
  COOP_INT,parameter::lmax = 2500
  COOP_REAL::ClTT(0:lmax), ClEE(0:lmax), ClBB(0:lmax), ClTE(0:lmax)
  type(coop_file)::fp
  
#if HAS_LAPACK
  call coop_get_command_line_argument(key = 'lcut', arg = lcut)
  call coop_get_command_line_argument(key = 'nside', arg = nside)
  call coop_get_command_line_argument(key = 'map', arg = fmap)
  call coop_get_command_line_argument(key = 'mask', arg = fmask)
  call coop_get_command_line_argument(key = 'genre', arg = genre)
  call coop_get_command_line_argument(key = 'fluc', arg = want_fluc, default  =.false.)
  call map%read(fmap)
  call mask%read(fmask)
  fullmap = map
  fullmask = mask
  call coop_healpix_maps_diffuse(from = map, to = fullmap, mask = mask, fwhm = 0.5d0/lcut)
  call coop_healpix_mask_diffuse(from = mask, to = fullmask, fwhm = 0.5d0/lcut)
  if(fullmap%nside .ne. fullmask%nside) call fullmask%udgrade_safe(fullmap%nside)
  call fullmap%convert2ring()
  call fullmask%convert2ring()
  call map%convert2ring()
  call mask%convert2ring()
  call map%map2alm(lmax=lcut*2)  !!to avoid noise induced svd instability
  call map%alm2map()
  sim = fullmap
  call map%udgrade_safe(nside)
  call mask%udgrade_safe(nside)
  call map%convert2ring()
  call mask%convert2ring()
  ClTT = 0.
  CLEE = 0.
  ClBB = 0.
  if(want_fluc)then
     call sim%allocate_alms( max(lcut, min(lmax, floor(sim%nside*coop_healpix_lmax_by_nside))) )
  else
     call sim%allocate_alms(lcut)
  endif
  sim%alm = 0.
  if(want_fluc)then
     call fp%open("planck14best_lensedCls.dat", "r")
     do l = 2, sim%lmax
        read(fp%unit, *) i, ClTT(l), ClEE(l), ClBB(l), ClTE(l)
        ClTT(l) = ClTT(l)*coop_2pi/l/(l+1.)
        ClEE(l) = ClEE(l)*coop_2pi/l/(l+1.)
        ClBB(l) = ClBB(l)*coop_2pi/l/(l+1.)
        ClTE(l) = ClTE(l)*coop_2pi/l/(l+1.)
     enddo
     call fp%close()
  endif
  select case(trim(genre))
  case("I")
     nlm = (lcut+1)**2 
     nmaps = 1
     if(want_fluc)then
        write(*,*) "adding random fluctuations up to l = ", sim%lmax
        call coop_random_init()
        ClTT= sqrt(ClTT)
        do l = lcut+1, sim%lmax
           sim%alm(l, 0, 1) = ClTT(l)*coop_random_complex_Gaussian(.true.)
           do m = 1, l
              sim%alm(l, m, 1) = ClTT(l)*coop_random_complex_Gaussian()
           enddo
        enddo
     endif
  case("QU")
     nlm = (lcut-1)*(lcut+3)  !l = 0, 1 excluded
     nmaps = 2
     if(want_fluc)then
        write(*,*) "adding random fluctuations up to l = ", sim%lmax
        call coop_random_init()
        ClEE= sqrt(ClEE)
        ClBB= sqrt(ClBB)
        do l = lcut+1, sim%lmax
           sim%alm(l, 0, 1) = ClEE(l)*coop_random_complex_Gaussian(.true.)
           sim%alm(l, 0, 2) = ClBB(l)*coop_random_complex_Gaussian(.true.)
           do m = 1, l
              sim%alm(l, m, 1) = ClEE(l)*coop_random_complex_Gaussian()
              sim%alm(l, m, 2) = ClBB(l)*coop_random_complex_Gaussian()
           enddo
        enddo
     endif

  case default
     stop "only support genre  = I or QU"
  end select

  if(map%nmaps .ne. nmaps) stop "input map nmaps ! =  genre nmaps"
  npix = count(mask%map(:,1) .gt. 0.5)
  lwork = 200*max(npix, nlm)*nmaps + 100000
  liwork = lwork
  print*, npix*nmaps, nlm*nmaps
  allocate(listpix(npix), mat(npix*nmaps, nlm*nmaps), yvec(max(npix, nlm)*nmaps,1), s(min(npix, nlm)*nmaps), work(lwork), iwork(liwork))
  j = 0
  do i=0, mask%npix-1
     if(mask%map(i,1) .gt. 0.5)then
        j = j + 1
        listpix(j) = i
        if(j.ge.npix)exit
     endif
  enddo
  call coop_healpix_get_alm2pix(lcut, map%nside, listpix, genre, mat)
  do i = 1, nmaps
     yvec((i-1)*npix+1:i*npix,1) = map%map(listpix,i)
  enddo
  call coop_prtsystime(.true.)
  call dgelsd(npix*nmaps, nlm*nmaps, 1, mat, npix*nmaps, yvec, max(npix, nlm)*nmaps, s, rcond, rank, work, lwork, iwork, info)
  call coop_prtsystime()
  i = 1
  do l=2*(nmaps-1),lcut
     do imap = 1, nmaps
        sim%alm(l, 0, imap) = yvec(i, 1)
        i = i + 1
     enddo
     do m = 1, l
        do imap = 1, nmaps
           sim%alm(l, m, imap) = cmplx(yvec(i, 1), yvec(i+1,1))
           i = i + 2
        enddo
     enddo
  enddo
  call sim%alm2map()
  do imap = 1, nmaps
     fullmap%map(:, imap)= fullmask%map(:,1)* fullmap%map(:, imap)+(1.- fullmask%map(:,1))*sim%map(:,imap)
  enddo
  call fullmap%write(coop_file_add_postfix(trim(fmap), "_INP"))
#else
  stop "you need install COOP with lapack to do inpainting"
#endif
end program shells

