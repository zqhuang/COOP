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
  COOP_INT::lcut, npix, nlm,  i, j, lwork, liwork, nside, l, m, imap
  type(coop_healpix_maps)::map, mask, sim, fullmap, fullmask
  logical::want_fluc
  COOP_STRING::genre, fmap, fmask, fcl
  COOP_INT,dimension(:),allocatable::listpix
  COOP_REAL,dimension(:,:),allocatable::mat, yvec
  COOP_REAL,dimension(:), allocatable::s, work
  COOP_INT,dimension(:), allocatable::iwork
  COOP_REAL::rcond = 1.d-8
  COOP_INT::rank, info, ncols
  COOP_INT,parameter::lmax = 2500
  COOP_REAL::ClTT(0:lmax), ClEE(0:lmax), ClBB(0:lmax), ClTE(0:lmax), readcls(6)
  COOP_REAL::sqrtClTT(0:lmax), sqrtClEE(0:lmax), sqrtClBB(0:lmax), sqrtCluEuE(0:lmax), EbyT(0:lmax), covMat
  type(coop_file)::fp
  
#if HAS_LAPACK
  if(iargc()<1)then
     write(*,*) "Syntax:"
     write(*,*) "./INP -lcut ... -nside ... -map ... -mask ... -genre ... -clfile ..."
     write(*,*) "genre can be I, QU, IQU, TEB"
     stop
  endif
  call coop_get_command_line_argument(key = 'lcut', arg = lcut)
  call coop_get_command_line_argument(key = 'nside', arg = nside)
  call coop_get_command_line_argument(key = 'map', arg = fmap)
  call coop_get_command_line_argument(key = 'mask', arg = fmask)
  call coop_get_command_line_argument(key = 'genre', arg = genre)
  call coop_get_command_line_argument(key = 'clfile', arg = fcl, default  ="") !!initial C_l input
  want_fluc = trim(fcl).ne.""
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
  ClTE = 0.
  if(want_fluc)then
     call sim%allocate_alms( max(lcut, min(lmax, floor(sim%nside*coop_healpix_lmax_by_nside))) )
  else
     call sim%allocate_alms( lcut )
  endif
  sim%alm = 0.
  if(want_fluc)then
     ncols = coop_file_NumColumns(fcl)
     if(ncols .lt. 2) stop "cl file is broken"
     call fp%open_skip_comments(fcl)
     do l = 2, sim%lmax
        read(fp%unit, *) i, readcls(1:ncols-1)
        ClTT(l) = readcls(1) * coop_2pi/l/(l+1.)
        if(ncols.gt.2)ClEE(l) = readcls(2) * coop_2pi/l/(l+1.)
        if(ncols.gt.3)ClBB(l) = readcls(3) * coop_2pi/l/(l+1.)
        if(ncols.gt.4)ClTE(l) = readcls(4) * coop_2pi/l/(l+1.)
     enddo
     call fp%close()
     sqrtClTT = sqrt(max(ClTT, 0.d0))
     sqrtClEE= sqrt(max(ClEE, 0.d0))
     sqrtClBB= sqrt(max(ClBB, 0.d0))
     sqrtCluEuE = sqrt(max(ClEE - ClTE**2/ClTT, 0.d0)) 
     EbyT = ClTE/max(ClTT, 1.d-30)
  endif
  select case(trim(genre))
  case("I")
     if(map%nmaps .ne. 1 .or. map%spin(1) .ne. 0) stop "The input map is not a simple I map"
     nlm = (lcut+1)**2 
     if(want_fluc)then
        write(*,*) "adding random fluctuations up to l = ", sim%lmax
        call coop_random_init()
        do l = lcut+1, sim%lmax
           sim%alm(l, 0, 1) = sqrtClTT(l)*coop_random_complex_Gaussian(.true.)
           do m = 1, l
              sim%alm(l, m, 1) = sqrtClTT(l)*coop_random_complex_Gaussian()
           enddo
        enddo
     endif
  case("QU")
     if(map%nmaps .ne. 2 .or. map%spin(1) .ne. 2) stop "The input map is not a simple QU map"
     if(map%spin(2) .ne. 2) stop "The input map is not a simple QU map"
     nlm = (lcut-1)*(lcut+3)  !l = 0, 1 excluded
     if(want_fluc)then
        write(*,*) "adding random fluctuations up to l = ", sim%lmax
        call coop_random_init()
        do l = lcut+1, sim%lmax
           sim%alm(l, 0, 1) = sqrtClEE(l)*coop_random_complex_Gaussian(.true.)
           sim%alm(l, 0, 2) = sqrtClBB(l)*coop_random_complex_Gaussian(.true.)
           do m = 1, l
              sim%alm(l, m, 1) = sqrtClEE(l)*coop_random_complex_Gaussian()
              sim%alm(l, m, 2) = sqrtClBB(l)*coop_random_complex_Gaussian()
           enddo
        enddo
     endif
  case("IQU")
     if(map%nmaps .ne. 3 .or. map%spin(1) .ne. 0) stop "The input map is not a simple IQU map"
     if(map%spin(2) .ne. 2 .or. map%spin(3).ne.2) stop "The input map is not a simple IQU map"
     nlm = (lcut-1)*(lcut+3)  !l = 0, 1 excluded
     if(want_fluc)then
        write(*,*) "adding random fluctuations up to l = ", sim%lmax
        call coop_random_init()
        do l = lcut+1, sim%lmax
           sim%alm(l, 0, 1) = sqrtClTT(l)*coop_random_complex_Gaussian(.true.)
           sim%alm(l, 0, 2) = sqrtCluEuE(l)*coop_random_complex_Gaussian(.true.) + EbyT(l)*sim%alm(l, 0, 1)
           sim%alm(l, 0, 3) = sqrtClBB(l)*coop_random_complex_Gaussian(.true.) 
           do m = 1, l
              sim%alm(l, m, 1) = sqrtClTT(l)*coop_random_complex_Gaussian()
              sim%alm(l, m, 2) = sqrtCluEuE(l)*coop_random_complex_Gaussian() + EbyT(l)*sim%alm(l, m, 1)
              sim%alm(l, m, 3) = sqrtClBB(l)*coop_random_complex_Gaussian()
           enddo
        enddo
     endif
  case("TEB")
     if(map%nmaps .ne. 3 .or. map%spin(1) .ne. 0) stop "The input map is not a simple TEB map"
     if(map%spin(2) .ne. 0 .or. map%spin(3).ne.0) stop "The input map is not a simple TEB map"
     nlm = (lcut-1)*(lcut+3)  !l = 0, 1 excluded
     if(want_fluc)then
        write(*,*) "adding random fluctuations up to l = ", sim%lmax
        call coop_random_init()
        do l = lcut+1, sim%lmax
           sim%alm(l, 0, 1) = sqrtClTT(l)*coop_random_complex_Gaussian(.true.)
           sim%alm(l, 0, 2) = sqrtCluEuE(l)*coop_random_complex_Gaussian(.true.) + EbyT(l)*sim%alm(l, 0, 1)
           sim%alm(l, 0, 3) = sqrtClBB(l)*coop_random_complex_Gaussian(.true.) 
           do m = 1, l
              sim%alm(l, m, 1) = sqrtClTT(l)*coop_random_complex_Gaussian()
              sim%alm(l, m, 2) = sqrtCluEuE(l)*coop_random_complex_Gaussian() + EbyT(l)*sim%alm(l, m, 1)
              sim%alm(l, m, 3) = sqrtClBB(l)*coop_random_complex_Gaussian()
           enddo
        enddo
     endif
  case default
     stop "only support genre  = I or QU"
  end select

  npix = count(mask%map(:,1) .gt. 0.5)
  if(npix .eq. 0)goto 200
  lwork = 200*max(npix, nlm)*map%nmaps + 100000
  liwork = lwork
  print*, npix*map%nmaps, nlm*map%nmaps
  allocate(listpix(npix), mat(npix*map%nmaps, nlm*map%nmaps), yvec(max(npix, nlm)*map%nmaps,1), s(min(npix, nlm)*map%nmaps), work(lwork), iwork(liwork))
  j = 0
  do i=0, mask%npix-1
     if(mask%map(i,1) .gt. 0.5)then
        j = j + 1
        listpix(j) = i
        if(j.ge.npix)exit
     endif
  enddo
  call coop_healpix_get_alm2pix(nmaps = map%nmaps, spin = map%spin, lmax = lcut, nside = map%nside, listpix = listpix, mat = mat)
  do i = 1, map%nmaps
     yvec((i-1)*npix+1:i*npix,1) = map%map(listpix,i)
  enddo
  call coop_prtsystime(.true.)
  call dgelsd(npix*map%nmaps, nlm*map%nmaps, 1, mat, npix*map%nmaps, yvec, max(npix, nlm)*map%nmaps, s, rcond, rank, work, lwork, iwork, info)
  if(info .ne. 0) stop "SVD solver failed"
  call coop_prtsystime()
  i = 1
  do l=sim%spin(1),lcut
     do imap = 1, map%nmaps
        if(map%spin(imap) .gt. l) cycle
        sim%alm(l, 0, imap) = yvec(i, 1)
        i = i + 1
     enddo
     do m = 1, l
        do imap = 1, map%nmaps
           if(map%spin(imap) .gt. l) cycle
           sim%alm(l, m, imap) = cmplx(yvec(i, 1), yvec(i+1,1))
           i = i + 2
        enddo
     enddo
  enddo
200 call sim%alm2map()
  do imap = 1, map%nmaps
     fullmap%map(:, imap) = fullmask%map(:,1) * fullmap%map(:, imap) + (1.- fullmask%map(:,1)) * sim%map(:,imap)
  enddo
  call fullmap%write(coop_file_add_postfix(trim(fmap), "_INP"))
#else
  stop "you need install COOP with lapack to do inpainting"
#endif
end program shells

