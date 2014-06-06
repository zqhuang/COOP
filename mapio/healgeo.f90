module healpix_geometry
!!I always assume ring order
  use wrap_utils
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  implicit none
#include "firstorder.h"
#include "zeta.h"

  integer, parameter::hg_default_lmax=3500
  integer, parameter::hg_mask_tol = 0.75
  integer::hg_inpainting_lowl=5

  integer,parameter::hg_index_TT = 1
  integer,parameter::hg_index_EE = 2
  integer,parameter::hg_index_BB = 3
  integer,parameter::hg_index_TE = 4
  integer,parameter::hg_index_EB = 5
  integer,parameter::hg_index_TB = 6

  type hg_disc
     integer nside
     integer center
     real(dl) theta, phi
     real(dl),dimension(3)::nx, ny, nz
  end type hg_disc

  type hg_maps
     integer npix, nside, nmaps, ordering, lmax, iq, iu, mask_npix, maskpol_npix
     character(LEN=80),dimension(64)::header
     integer,dimension(:),allocatable::spin
     real, dimension(:,:),allocatable::map
     complex, dimension(:,:,:),allocatable::alm
     real, dimension(:,:),allocatable::Cl
     integer,dimension(:),allocatable::mask_listpix, maskpol_listpix
     real chisq, mcmc_temperature
  end type hg_maps



#define COS2RADIUS(cosx) (sqrt(2.d0*(1.d0 - (cosx))))
#define RADIUS2COS(r)  (1.d0-(r)**2/2.d0)

contains

  subroutine hg_simulate(hgm)
    type(hg_maps) hgm
    real,dimension(:),allocatable::sqrtCls
    real,dimension(:, :),allocatable::Cls_sqrteig
    real,dimension(:,:,:),allocatable::Cls_rot
    integer l
    if(hgm%nmaps.eq.1 .and. hgm%spin(1).eq.0)then
       allocate(sqrtCls(0:hgm%lmax))
       !$omp parallel do
       do l = 0, hgm%lmax
          sqrtCls(l) = sqrt(hgm%Cl(l,1))
       enddo
       !$omp end parallel do
       call hg_simulate_Tmaps(hgm, hgm%nside, hgm%lmax, sqrtCls)
       deallocate(sqrtCls)
    elseif(hgm%nmaps.eq.3 .and. hgm%iq .eq.2)then
       allocate(Cls_sqrteig(3, 0:hgm%lmax), Cls_rot(3,3,0:hgm%lmax))
       call hg_Cls2Rot(hgm%lmax, hgm%Cl, Cls_sqrteig, Cls_rot)
       call hg_simulate_TQUmaps(hgm, hgm%nside, hgm%lmax, Cls_sqrteig, Cls_rot)
       deallocate(Cls_sqrteig, Cls_rot)
    else
       stop "unknown hg_simulate mode"
    endif
  end subroutine hg_simulate


  subroutine hg_simulate_Tmaps(hgm, nside, lmax, sqrtCls)
    type(hg_maps) hgm
    integer nside
    integer lmax
    real sqrtCls(0:lmax)
    integer l,m
    call hg_ini_map(hgm, nside = nside, nmaps = 1, spin = (/ 0 /), lmax = lmax)
    !$omp parallel do private(l, m)
    do l=0, lmax
       hgm%alm(l, 0, 1) = random_complex_Gaussian(.true.)*SqrtCls(l)     
       do m = 1, l
          hgm%alm(l, m, 1) = random_complex_Gaussian()*SqrtCls(l)
       enddo
    enddo
    !$omp end parallel do
    call hg_alm2map(hgm)
  end subroutine hg_simulate_Tmaps


  subroutine hg_get_Cls(hgm) !!I assume you have already called    hgm_map2alm(hgm)
    type(hg_maps)hgm
    integer l, m, i, j, k
    if(.not.allocated(hgm%alm)) stop "hg_get_Cls: you have to call hg_map2alm before calling this subroutine"
    !$omp parallel do private(i,j,k,l)
    do i=1, hgm%nmaps
       do j=1, i
          k = matsym_index(hgm%nmaps, i, j)
          do l = 0, hgm%lmax
             hgm%Cl(l, k) = (sum(MULT_REAL(hgm%alm(l, 1:l, i), hgm%alm(l, 1:l, j))) + 0.5d0 * MULT_REAL(hgm%alm(l,0,i), hgm%alm(l,0,j)) )/(l+0.5d0)
          enddo
       enddo
    enddo
    !$omp end parallel do
  end subroutine hg_get_Cls

  subroutine hg_Cls2Rot(lmax, Cls, Cls_sqrteig, Cls_rot)
    integer lmax
    real,dimension(0:lmax, 6),intent(IN)::Cls !!ordering is TT, EE, BB, TE, EB, TB
    real, dimension(3, 0:lmax),intent(OUT)::Cls_sqrteig
    real, dimension(3, 3, 0:lmax),intent(OUT)::Cls_rot
    integer l
    real(dl) a2(2,2), a3(3,3)
    real(dl) psi2(2,2), psi3(3,3)
    if(all(Cls(:,hg_index_EB).eq.0.d0) .and. all(Cls(:, hg_index_TB).eq.0.d0))then
       Cls_sqrteig(hg_index_BB,:) = sqrt(Cls(:, hg_index_BB))
       Cls_rot(hg_index_BB, hg_index_TT, :) = 0
       Cls_rot(hg_index_BB, hg_index_EE, :) = 0
       Cls_rot(hg_index_TT, hg_index_BB, :) = 0
       Cls_rot(hg_index_EE, hg_index_BB, :) = 0
       Cls_rot(hg_index_BB, hg_index_BB, :) = 1.
       do l=0, lmax
          a2(1,1) = Cls(l, hg_index_TT)
          a2(2,2) = Cls(l, hg_index_EE)
          a2(1,2) = Cls(l, hg_index_TE)
          a2(2,1) = a2(1,2)
          call matsymdiag_small(2, a2, psi2)
          Cls_sqrteig(hg_index_TT, l) = sqrt(a2(1,1))
          Cls_sqrteig(hg_index_EE, l) = sqrt(a2(2,2))
          Cls_rot( hg_index_TT:hg_index_EE, hg_index_TT:hg_index_EE, l) = psi2
       end do
    else
       do l=0, lmax
          a3(1,1) = Cls(l,hg_index_TT)
          a3(2,2) = Cls(l,hg_index_EE)
          a3(3,3) = Cls(l,hg_index_BB)
          a3(1,2) = Cls(l,hg_index_TE)
          a3(2,3) = Cls(l,hg_index_EB)
          a3(3,1) = Cls(l,hg_index_TB)
          a3(2,1) = a3(1,2)
          a3(3,2) = a3(2,3)
          a3(1,3) = a3(3,1)
          call matsymdiag_small(3, a3, psi3)
          Cls_sqrteig(hg_index_TT, l) = sqrt(a3(1,1))
          Cls_sqrteig(hg_index_EE, l) = sqrt(a3(2,2))
          Cls_sqrteig(hg_index_BB, l) = sqrt(a3(3,3))
          Cls_rot(:, :, l) = psi3
       end do
    endif
  end subroutine hg_Cls2Rot

  subroutine hg_iqu2TQTUT(map_file)
    UNKNOWN_STRING map_file
    type(hg_maps) hgm
    call hg_read_map(map_file, hgm, nmaps_wanted = 3, spin = (/ 0, 2 , 2 /) )
    call hg_map2alm(hgm)
    hgm%alm(:,:,2) = hgm%alm(:,:,1)
    hgm%alm(:,:,3) = 0
    call hg_alm2map(hgm)
    call hg_write_map(file_add_postfix(trim(str_replace(trim(map_file), "iqu", "TQTUT")), "_transTQTUT"), hgm)
    call hg_free_map(hgm)
  end subroutine hg_iqu2TQTUT

  subroutine hg_iqu2TEB(map_file)
    UNKNOWN_STRING map_file
    type(hg_maps) hgm
    call hg_read_map(map_file, hgm, nmaps_wanted = 3, spin = (/ 0, 2, 2 /) )
    call hg_map2alm(hgm)
    hgm%spin = 0
    hgm%iq = 0
    hgm%iu = 0
    call hg_alm2map(hgm)
    call hg_write_map(file_add_postfix(trim(str_replace(trim(map_file), "iqu", "TEB")), "_transTEB"), hgm)
    call hg_free_map(hgm)
  end subroutine hg_iqu2TEB


  subroutine hg_simulate_TQUmaps(hgm, nside, lmax, Cls_sqrteig, Cls_rot)
    type(hg_maps) hgm
    integer lmax, nside
    real,dimension(3, 0:lmax)::Cls_sqrteig
    real,dimension(3, 3, 0:lmax)::Cls_rot
    integer l, m
    call hg_ini_map(hgm, nside = nside, nmaps = 3, spin = (/ 0, 2, 2 /), lmax = lmax)    
    !$omp parallel do private(l, m)
    do l=0, lmax
       hgm%alm(l, 0, :) = matmul(Cls_rot(:, :, l), Cls_sqrteig(:,l) * (/ random_complex_Gaussian(.true.), random_complex_Gaussian(.true.), random_complex_Gaussian(.true.) /) )
       do m = 1, l
          hgm%alm(l, m, :) = matmul(Cls_rot(:, :, l), Cls_sqrteig(:,l) * (/ random_complex_Gaussian(), random_complex_Gaussian(), random_complex_Gaussian() /) )
       enddo
    enddo
    !$omp end parallel do
    call hg_alm2map(hgm)
  end subroutine hg_simulate_TQUmaps


  subroutine hg_ini_map(hgm, nside, nmaps, spin, lmax)
    type(hg_maps) hgm
    integer:: nside, nmaps
    integer:: spin(nmaps)
    integer, optional::lmax
    if(allocated(hgm%map))then
       if(hgm%nside .eq. nside .and. hgm%nmaps.eq.nmaps)then
          goto 100
       endif
       deallocate(hgm%map)
    endif
    if(allocated(hgm%spin))deallocate(hgm%spin)
    hgm%nside = nside
    hgm%nmaps = nmaps
    hgm%npix = nside2npix(nside)
    allocate(hgm%map(0:hgm%npix - 1, nmaps))
    allocate(hgm%spin(hgm%nmaps))
100 hgm%spin = spin
    if(all(hgm%spin(1:hgm%nmaps-1) .ne. 2))then
       hgm%iq = 0
       hgm%iu = 0
    else
       hgm%iq = 1
       do while(hgm%spin(hgm%iq) .ne. 2)
          hgm%iq = hgm%iq + 1
       enddo
       hgm%iu = hgm%iq + 1
       if(hgm%spin(hgm%iu).ne.2)then
          hgm%iq = 0
          hgm%iu = 0
       endif
    endif
    if(present(lmax))then
       if(allocated(hgm%alm))then
          if(hgm%lmax  .eq. lmax)then
             goto 200
          endif
          deallocate(hgm%alm)
       endif
       if(allocated(hgm%cl))deallocate(hgm%cl)
       hgm%lmax = lmax
       allocate(hgm%alm(0:hgm%lmax, 0:hgm%lmax, hgm%nmaps))
       allocate(hgm%cl(0:hgm%lmax, hgm%nmaps*(hgm%nmaps+1)/2))
    endif
200 hgm%ordering = RING !!default ordering
    call write_minimal_header(hgm%header,dtype = 'MAP', nside=hgm%nside, order = hgm%ordering, creator='Zhiqi Huang', version = 'CosmoLib', units='muK', polar=any(hgm%spin.eq.2) )
    hgm%maskpol_npix = 0
  end subroutine hg_ini_map

  subroutine hg_free_map(hgm)
    type(hg_maps) hgm
    if(allocated(hgm%map))deallocate(hgm%map)
    if(allocated(hgm%alm))deallocate(hgm%alm)
    if(allocated(hgm%cl))deallocate(hgm%cl)
    if(allocated(hgm%spin))deallocate(hgm%spin)
    if(allocated(hgm%mask_listpix))deallocate(hgm%mask_listpix)
  end subroutine hg_free_map

  subroutine hg_read_map(filename, hgm, nmaps_wanted, spin)
    type(hg_maps) hgm
    UNKNOWN_STRING filename
    integer,optional::nmaps_wanted
    integer,dimension(:),optional::spin
    integer(8) npixtot
    integer nmaps_actual
    if(.not. file_exists(filename))then
       write(*,*) trim(filename)
       stop "cannot find the file"
    endif
    npixtot = getsize_fits(trim(filename), nmaps = nmaps_actual, nside = hgm%nside, ordering = hgm%ordering)
    hgm%npix =nside2npix(hgm%nside)
    if(present(nmaps_wanted))then       
       hgm%nmaps = nmaps_wanted
       if(nmaps_wanted .lt. nmaps_actual)then
          nmaps_actual = nmaps_wanted
       endif
    else
       hgm%nmaps = nmaps_actual
    endif
    if(allocated(hgm%spin))then
       if(size(hgm%spin) .ne. hgm%nmaps)then
          deallocate(hgm%spin)
          allocate(hgm%spin(hgm%nmaps))
       endif
    else
       allocate(hgm%spin(hgm%nmaps))
    endif
    if(present(spin))then
       if(size(spin).ne. hgm%nmaps)then
          stop "hg_read_map: the list of spins should have the same size as nmaps"
       else
          hgm%spin = spin
          hgm%iq = 1
          do while(hgm%spin(hgm%iq) .ne.2 .and. hgm%iq .lt. hgm%nmaps)
             hgm%iq = hgm%iq + 1
          enddo
          if(hgm%iq .ge. hgm%nmaps)then
             hgm%iq = 0
             hgm%iu = 0
          else
             hgm%iu = hgm%iq + 1
          endif
       endif
    else
       if( (index(filename, "iqu") .ne. 0 .or. index(filename, "IQU") .ne. 0 .or.  index(filename, "tqu") .ne. 0  .or.  index(filename, "TQU") .ne. 0  .or.  index(filename, "TQTUT") .ne. 0   .or.  index(filename, "IQTUT") .ne. 0) .and. hgm%nmaps .ge. 3 )then
          hgm%spin = 0
          hgm%spin(2:3) = 2
          hgm%iq = 2
          hgm%iu = 3
          write(*,*) hgm%nmaps, " for nmaps >= 3 I assume it is an IQU map"
       elseif(((index(filename, "qu") .ne. 0 .or. index(filename, "QU") .ne. 0) .or.  index(filename, "QTUT") .ne. 0) .and. hgm%nmaps .ge. 2 )then
          hgm%spin = 0
          hgm%spin(1:2) = 2
          hgm%iq = 1
          hgm%iu = 2
          write(*,*) hgm%nmaps, " for nmaps >= 2 I assume it is an QU map"
       else
          hgm%iq = 0
          hgm%iu = 0
          hgm%spin = 0
          if(hgm%nmaps .gt. 1) write(*,*) "hg_read_map: I am assuming the maps are all scalar"
       endif
    endif
    if(allocated(hgm%map))then
       if(size(hgm%map, 1).ne. hgm%npix .or. size(hgm%map, 2).ne.hgm%nmaps)then
          deallocate(hgm%map)
          allocate(hgm%map(0:hgm%npix-1, hgm%nmaps))
       endif
    else
       allocate(hgm%map(0:hgm%npix-1, hgm%nmaps))
    endif

    call input_map(trim(filename), hgm%map, hgm%npix, nmaps_actual, fmissval = 0.)
    call hg_convert_to_ring(hgm)
    call write_minimal_header(hgm%header,dtype = 'MAP', nside=hgm%nside, order = hgm%ordering, creator='Zhiqi Huang', version = 'CosmoLib', units='muK', polar=any(hgm%spin.eq.2) )

  end subroutine hg_read_map

  subroutine hg_copy_map(hgm1, hgm2)
    type(hg_maps)::hgm1, hgm2
    hgm2%nside = hgm1%nside
    hgm2%npix = hgm1%npix
    hgm2%ordering = hgm1%ordering
    hgm2%nmaps = hgm1%nmaps
    hgm2%lmax = hgm1%lmax
    hgm2%chisq = hgm1%chisq
    hgm2%mcmc_temperature = hgm1%mcmc_temperature
    hgm2%iq = hgm1%iq
    hgm2%iu = hgm1%iu
    if(allocated(hgm2%map))then
       if(size(hgm2%map, 1) .ne. hgm2%npix .or. size(hgm2%map, 2).ne. hgm2%nmaps)then
          deallocate(hgm2%map)
          allocate(hgm2%map(0:hgm2%npix-1, hgm2%nmaps))
       endif
    else
       allocate(hgm2%map(0:hgm2%npix-1, hgm2%nmaps))
    endif
    hgm2%map = hgm1%map
    if(allocated(hgm2%spin))then
       if(size(hgm2%spin).ne. hgm2%nmaps)then
          deallocate(hgm2%spin)
          allocate(hgm2%spin(hgm2%nmaps))
       endif
    else
       allocate(hgm2%spin(hgm2%nmaps))
    endif
    hgm2%spin = hgm1%spin
    if(allocated(hgm1%alm) .and. allocated(hgm1%cl))then
       if(allocated(hgm2%alm))then
          if(size(hgm2%alm, 1) .ne. hgm2%lmax+1 .or. size(hgm2%alm, 2) .ne. hgm2%lmax+1 .or. size(hgm2%alm, 3) .ne. hgm2%nmaps)then
             deallocate(hgm2%alm)
             allocate(hgm2%alm(0:hgm2%lmax,0:hgm2%lmax, hgm2%nmaps))
          endif
       else
          allocate(hgm2%alm(0:hgm2%lmax,0:hgm2%lmax, hgm2%nmaps))
       endif
       if(allocated(hgm2%cl))then
          if(size(hgm2%cl, 1) .ne. hgm2%lmax+1 .or. size(hgm2%cl, 2) .ne. hgm2%nmaps*(hgm2%nmaps+1)/2)then
             deallocate(hgm2%cl)
             allocate(hgm2%cl(0:hgm2%lmax, hgm2%nmaps*(hgm2%nmaps+1)/2))
          endif
       else
          allocate(hgm2%cl(0:hgm2%lmax, hgm2%nmaps*(hgm2%nmaps+1)/2))
       endif
       hgm2%alm = hgm1%alm
       hgm2%cl = hgm1%cl
    endif
    if(allocated(hgm1%mask_listpix))then
       if(allocated(hgm2%mask_listpix))then
          if(hgm2%mask_npix .ne. hgm1%mask_npix)then
             hgm2%mask_npix = hgm1%mask_npix
             deallocate(hgm2%mask_listpix)
             allocate(hgm2%mask_listpix(hgm2%mask_npix))
          endif
       else
          hgm2%mask_npix = hgm1%mask_npix
          allocate(hgm2%mask_listpix(hgm2%mask_npix))
       endif
       hgm2%mask_listpix = hgm1%mask_listpix
    endif
    if(allocated(hgm1%maskpol_listpix))then
       if(allocated(hgm2%maskpol_listpix))then
          if(hgm2%maskpol_npix .ne. hgm1%maskpol_npix)then
             hgm2%maskpol_npix = hgm1%maskpol_npix
             deallocate(hgm2%maskpol_listpix)
             allocate(hgm2%maskpol_listpix(hgm2%maskpol_npix))
          endif
       else
          hgm2%maskpol_npix = hgm1%maskpol_npix
          allocate(hgm2%maskpol_listpix(hgm2%maskpol_npix))
       endif
       hgm2%maskpol_listpix = hgm1%maskpol_listpix
    endif

  end subroutine hg_copy_map

  subroutine hg_write_map(filename, hgm, index_list)
    type(hg_maps)hgm
    UNKNOWN_STRING filename
    integer,dimension(:),optional::index_list
    logical pol
    if(present(index_list))then
       if(any(index_list .lt. 1 .or. index_list .gt. hgm%nmaps)) stop "hg_write_map: index out of range"
       pol = any(hgm%spin(index_list).eq.2)
    else
       pol =any(hgm%spin.eq.2)
    endif
    call delete_file(trim(filename))
    if(allocated(hgm%alm))then
       call write_minimal_header(hgm%header,dtype = 'MAP', nside=hgm%nside, order = hgm%ordering, creator='Zhiqi Huang', version = 'CosmoLib', units='muK', nlmax = hgm%lmax, nmmax = hgm%lmax, polar= pol)
    else
       call write_minimal_header(hgm%header,dtype = 'MAP', nside=hgm%nside, order = hgm%ordering, creator='Zhiqi Huang', version = 'CosmoLib', units='muK', polar= pol )
    endif
    if(present(index_list))then
       call output_map(hgm%map(:, index_list), hgm%header, trim(filename))
    else
       call output_map(hgm%map, hgm%header, trim(filename))
    endif
  end subroutine hg_write_map

  subroutine hg_convert_to_nested(hgm)
    type(hg_maps) hgm
    if(.not. allocated(hgm%map)) stop "hg_convert_to_nested: map is not allocated yet"
    if(hgm%ordering .eq. RING) then
       call convert_ring2nest(hgm%nside, hgm%map)
       hgm%ordering = NESTED
    elseif(hgm%ordering .ne. NESTED)then
       write(*,*) "ordering = ", hgm%ordering
       stop "hg_convert_to_nested: unknown ordering"
    endif
  end subroutine hg_convert_to_nested

  subroutine hg_convert_to_ring(hgm)
    type(hg_maps) hgm
    if(.not. allocated(hgm%map)) stop "hg_convert_to_ring: map is not allocated yet"
    if(hgm%ordering .eq. NESTED)then
       call convert_nest2ring(hgm%nside, hgm%map)
       hgm%ordering = RING
    elseif(hgm%ordering .ne. RING)then
       write(*,*) "ordering = ", hgm%ordering
       stop "hg_convert_to_ring: UNKNOWN ordering"
    endif
  end subroutine hg_convert_to_ring

  subroutine hg_map2alm(hgm, lmax)
    type(hg_maps) hgm
    integer,optional::lmax
    integer i, l
    complex, dimension(:,:,:),allocatable::alm
    call hg_convert_to_ring(hgm)
    if(present(lmax))then
       if(lmax .gt. hgm%nside*3)then
          write(*,*) "lmax > nside x 3 is not recommended"
          stop
       else
          hgm%lmax = lmax          
       endif
    else
       hgm%lmax =  min(hg_default_lmax, hgm%nside*2)
    endif
    if(allocated(hgm%alm))then
       if(size(hgm%alm, 1) .ne. hgm%lmax+1 .or. size(hgm%alm, 2) .ne. hgm%lmax+1 .or. size(hgm%alm, 3) .ne. hgm%nmaps)then
          deallocate(hgm%alm)
          allocate(hgm%alm(0:hgm%lmax, 0:hgm%lmax, hgm%nmaps))
       endif
    else
       allocate(hgm%alm(0:hgm%lmax, 0:hgm%lmax, hgm%nmaps))
    endif
    hgm%alm = 0.
    if(allocated(hgm%cl))then
       if(size(hgm%cl, 1) .ne. hgm%lmax+1 .or. size(hgm%cl, 2) .ne. hgm%nmaps*(hgm%nmaps+1)/2)then
          deallocate(hgm%cl)
          allocate(hgm%cl(0:hgm%lmax, hgm%nmaps*(hgm%nmaps+1)/2))
       endif
    else
       allocate(hgm%cl(0:hgm%lmax, hgm%nmaps*(hgm%nmaps+1)/2))
    endif

    i = 1
    do while(i.le. hgm%nmaps)
       if(hgm%spin(i).eq.0)then
          call map2alm(hgm%nside, hgm%lmax, hgm%lmax, hgm%map(:,i), hgm%alm(:,:,i:i))
          i = i + 1
       else
          if(.not. allocated(alm))allocate(alm(2, 0:hgm%lmax, 0:hgm%lmax))
          if(i.lt. hgm%nmaps)then
             if(hgm%spin(i+1) .eq. hgm%spin(i))then
                call map2alm_spin(hgm%nside, hgm%lmax, hgm%lmax, hgm%spin(i), hgm%map(:,i:i+1), alm)
                hgm%alm(:,:,i) = alm(1, :, :)
                hgm%alm(:,:,i+1) = alm(2, :, :)
                i = i + 2
                cycle
             endif
          endif
          write(*,*) hgm%spin
          stop "hg_map2alm: nonzero spin maps must appear in pairs"
       endif
    enddo
    if(allocated(alm))deallocate(alm)
    call hg_get_Cls(hgm)
  end subroutine hg_map2alm


  subroutine hg_alm2map(hgm)
    type(hg_maps) hgm
    integer i
    complex,dimension(:,:,:),allocatable::alm
    i = 1
    do while(i.le. hgm%nmaps)
       if(hgm%spin(i).eq.0)then
          call alm2map(hgm%nside, hgm%lmax, hgm%lmax, hgm%alm(:,:,i:i), hgm%map(:,i))
          i = i + 1
       else
          if(.not.allocated(alm))allocate(alm(2,0:hgm%lmax, 0:hgm%lmax))
          if(i.lt. hgm%nmaps)then
             alm(1,:,:) = hgm%alm(:,:,i)
             alm(2,:,:) = hgm%alm(:,:,i+1)
             if(hgm%spin(i+1) .eq. hgm%spin(i))then
                call alm2map_spin(hgm%nside, hgm%lmax, hgm%lmax, hgm%spin(i), alm, hgm%map(:,i:i+1))
                i = i + 2
                cycle
             endif
          endif
          stop "hg_alm2map: nonzero spin maps must appear in pairs"
       endif
    enddo
    if(allocated(alm))deallocate(alm)
  end subroutine hg_alm2map

  subroutine hg_filter_alm(hgm, fwhm, lpower, window)
    type(hg_maps) hgm
    real,optional::window(0:hgm%lmax)
    real,optional::fwhm
    real,optional::lpower
    integer l
    real c, w(0:hgm%lmax)
    w = 1.
    if(present(fwhm))then
       c = sign((const_sigmabyfwhm * fwhm)**2/2., dble(fwhm))
       !$omp parallel do
       do l = 0,  hgm%lmax
          w(l) = w(l)*exp(-l*(l+1.)*c)
       enddo
       !$omp end parallel do
    endif
    if(present(lpower))then
       !$omp parallel do
       do l = 0,  hgm%lmax
          w(l) = w(l)*(l*(l+1.))**(lpower/2.)
       enddo
       !$omp end parallel do
    endif
    if(present(window))then
       !$omp parallel do
       do l = 0,  hgm%lmax
          w(l) = w(l)*window(l)
       enddo
       !$omp end parallel do       
    endif
    !$omp parallel do
    do l = 0, hgm%lmax
       hgm%alm(l,:,:) = hgm%alm(l,:,:)*w(l)
       hgm%Cl(l,:) = hgm%Cl(l,:)*w(l)**2
    enddo
    !$omp end parallel do
  end subroutine hg_filter_alm


  subroutine split_angular_mode(n, qmap, umap, m, nr, fr)
    integer n, m, nr
    real qmap(-n:n,-n:n), umap(-n:n,-n:n)
    real fr(0:nr), w(0:nr), q, u, r, phi, fpoint
    integer i, j, ir
    fr = 0
    w = 0
    do i = -n, n
       do j = -n, n
          r = sqrt(real(i)**2 + real(j)**2)
          ir = floor(r)          
          if(ir .le. nr)then
             phi = m*POLAR_ANGLE(real(i), real(j))
             r = r - ir
             fpoint = qmap(i,j)*cos(phi) + umap(i,j)*sin(phi)
             fr(ir) = fr(ir) + fpoint*(1.d0-r)
             w(ir) = w(ir) + 1.d0-r
             if(ir.ne.nr)then
                fr(ir+1) = fr(ir+1)+fpoint*r
                w(ir+1) = w(ir+1)+r
             endif
          endif
       enddo
    enddo
    do ir = 0, nr
       if(w(ir).gt.0.)then
          fr(ir) = fr(ir)/w(ir)
       endif
    enddo
    do ir = n, nr !!damp the amplitude in the corners
       fr(ir) = fr(ir)*exp(-((ir-n+1.)/n*(1.414/0.414))**2)
    enddo
    if(m.ne.0) fr(0) = 0
  end subroutine split_angular_mode


  subroutine map_filter_modes(n, qmap, umap, ms)
    integer n, nm
    integer ms(:)
    real qmap(-n:n, -n:n), umap(-n:n, -n:n)
    real,dimension(:,:),allocatable::fr
    real r, phi, s1, s2
    integer i, j, ir, nr, im
    nm = size(ms)
    nr = ceiling(const_sqrt2*n)+1
    allocate(fr(0:nr, nm))
    do im = 1, nm
       call split_angular_mode(n, qmap, umap, ms(im), nr, fr(0:nr, im))
    enddo
    qmap = 0
    umap = 0
    do i=-n, n
       do j=-n,n
          r = sqrt(real(i)**2 + real(j)**2)
          ir = floor(r)
          r = r - ir
          phi = POLAR_ANGLE(real(i), real(j))
          do  im =1, nm
             qmap(i,j) = qmap(i,j) + (fr(ir,im)*(1.-r)+fr(ir+1,im)*r)*cos(ms(im)*phi)
             umap(i,j) = umap(i,j) + (fr(ir,im)*(1.-r)+fr(ir+1,im)*r)*sin(ms(im)*phi)
             
          enddo
       enddo
    enddo
    deallocate(fr)
  end subroutine map_filter_modes


  subroutine hg_get_disc(nside, pix, disc)
    integer pix, nside
    type(hg_disc) disc
    real(dl) r
    disc%nside  = nside
    disc%center = pix
    call pix2ang_ring(nside, pix, disc%theta, disc%phi)
    call ang2vec(disc%theta, disc%phi, disc%nz)
    disc%nx = (/  sin(disc%phi) , - cos(disc%phi) , 0.d0 /)
    call vector_cross_product(disc%nz, disc%nx, disc%ny)
  end subroutine hg_get_disc

  subroutine hg_disc_pix2ang(disc, pix, r, phi)
    type(hg_disc) disc
    integer pix
    real(dl) r, phi, vec(3), x, y
    if(pix .eq. disc%center)then
       r = 0
       phi = 0
       return
    endif
    call pix2vec_ring(disc%nside, pix, vec)
    r = COS2RADIUS( dot_product(vec, disc%nz) )
    x = dot_product(vec, disc%nx)
    y = dot_product(vec, disc%ny)
    phi = POLAR_ANGLE(x, y)
  end subroutine hg_disc_pix2ang

  subroutine hg_disc_ang2pix(disc, r, phi, pix)
    type(hg_disc) disc
    real(dl) r !!in unit of radian
    real(dl) phi, vec(3), cost, sint
    integer pix
    cost = RADIUS2COS(r)
    sint = sqrt(1.d0 - cost**2)
    vec = sint*cos(phi)* disc%nx + sint*sin(phi)*disc%ny + cost*disc%nz
    call vec2pix_ring(disc%nside, vec, pix)
  end subroutine hg_disc_ang2pix


  subroutine hg_disc_pix2xy(disc, pix, x, y)
    type(hg_disc) disc
    integer pix
    real(dl) r, phi, vec(3), x, y
    if(pix .eq. disc%center)then
       x = 0
       y = 0
       return
    endif
    call pix2vec_ring(disc%nside, pix, vec)
    r = COS2RADIUS( dot_product(vec, disc%nz) )
    x = dot_product(vec, disc%nx)
    y = dot_product(vec, disc%ny)
    r = r/sqrt(x**2+y**2)
    x = r * x
    y = r * y
  end subroutine hg_disc_pix2xy


  subroutine hg_disc_xy2pix(disc, x, y, pix)
    type(hg_disc) disc
    real(dl) x, y !!in unit of radian
    real(dl) vec(3), cost, sint, r
    integer pix
    r = sqrt(x**2+y**2)
    if(r.lt.1.d-8)then
       pix = disc%center
       return
    endif
    cost = RADIUS2COS(r)
    sint = sqrt(1.d0 - cost**2)
    vec = sint*(x/r)* disc%nx + sint*(y/r)*disc%ny + cost*disc%nz
    call vec2pix_ring(disc%nside, vec, pix)
  end subroutine hg_disc_xy2pix

  subroutine hg_rotate_qu(qu, phi)
    real qu(2)
    real(dl) phi, cosp, sinp
    cosp = cos(2.d0*phi)
    sinp = sin(2.d0*phi)
    qu = (/ qu(1)*cosp + qu(2)*sinp,  -qu(1)*sinp + qu(2)*cosp /)
  end subroutine hg_rotate_qu

  subroutine hg_stack(map, disc, n, rpix, angle, stack_option, image, uimage, mask, counter)
    character(LEN=*) stack_option
    type(hg_disc) disc
    real qu(2)
    integer n
    real image(-n:n, -n:n), tmpq(-n:n, -n:n), tmpu(-n:n, -n:n)
    real,optional::uimage(-n:n, -n:n)
    real(dl) rpix, angle
    integer i, j, pix
    real(dl) r, phi,  x, y
    type(hg_maps) map
    type(hg_maps),optional::mask
    real mask_count
    real counter
    tmpq = 0
    tmpu = 0
    mask_count = 0
    select case(trim(stack_option))
    case("T", "E", "B")
       do i = -n, n
          do j = -n, n
             x = rpix*i
             y = rpix*j
             r = sqrt(x**2+y**2)
             phi = POLAR_ANGLE(x, y) + angle
             call hg_disc_ang2pix(disc, r, phi, pix)
             if(present(mask))then
                if(mask%map(pix,1) .gt. 0.5)then
                   mask_count = mask_count + mask%map(pix, 1)
                   tmpq(i, j) = tmpq(i, j) + map%map(pix,1)*mask%map(pix,1)
                endif
             else 
                tmpq(i, j) = tmpq(i, j) + map%map(pix,1)
             endif
          enddo
       enddo
    case("Q", "U")
       do i = -n, n
          do j = -n, n
             x = rpix*i
             y = rpix*j
             r = sqrt(x**2+y**2)
             phi = POLAR_ANGLE(x, y) + angle
             call hg_disc_ang2pix(disc, r, phi, pix)
             if(present(mask)) then
                if(mask%map(pix, 1) .gt. 0.5)then
                   qu = map%map(pix, map%iq:map%iu)
                   call hg_rotate_qu(qu, angle)
                   mask_count = mask_count + mask%map(pix, 1) 
                   tmpq(i, j) = tmpq(i, j) + qu(1)*mask%map(pix, 1)
                   tmpu(i, j) = tmpu(i, j) + qu(2)*mask%map(pix, 1)
                endif
             else
                qu = map%map(pix, map%iq:map%iu)
                call hg_rotate_qu(qu, angle)
                tmpq(i, j) = tmpq(i, j) + qu(1)
                tmpu(i, j) = tmpu(i, j) + qu(2)
             endif
          enddo
       enddo
    case("Qr", "QR", "Ur", "UR")
       do i = -n, n
          do j = -n, n
             x = rpix*i
             y = rpix*j
             r = sqrt(x**2+y**2)
             phi = POLAR_ANGLE(x, y) + angle
             call hg_disc_ang2pix(disc, r, phi, pix)
             if(present(mask))then
                if(mask%map(pix, 1) .gt. 0.5)then
                   qu = map%map(pix, map%iq:map%iu)
                   call hg_rotate_qu(qu, phi)
                   tmpq(i, j) = tmpq(i, j) + qu(1)*mask%map(pix, 1)
                   tmpu(i, j) = tmpu(i, j) + qu(2)*mask%map(pix, 1)
                   mask_count = mask_count + mask%map(pix, 1) 
                endif
             else
                qu = map%map(pix, map%iq:map%iu)
                call hg_rotate_qu(qu, phi)
                tmpq(i, j) = tmpq(i, j) + qu(1)
                tmpu(i, j) = tmpu(i, j) + qu(2)
                
             endif
          enddo
       enddo
    case default
       write(*,*) "hg_stack: Unknown stack_option"//trim(stack_option)
       stop
    end select
    if(present(mask) .and. mask_count .lt. (2*n+1.)**2*hg_mask_tol)then
       return
    endif
    if(present(mask))then
       counter = counter + mask_count/(2*n+1.)**2
    else
       counter = counter + 1
    endif
    select case(trim(stack_option))
    case("T", "E", "B")
       image = image + tmpq
    case("Q", "Qr", "QR")
       image = image +  tmpq
       if(present(uimage)) uimage = uimage + tmpu
    case("U", "Ur", "UR")
       image = image + tmpu
       if(present(uimage)) stop "hg_stack: for U stacking this should not happen"
    end select
  end subroutine hg_stack


  subroutine hg_stack_io(map_file, mean_image_file, spots_file, rmax, r_resolution, stack_option, title, headless_vector, m_filter, caption, mask_file, smooth_fwhm)
    UNKNOWN_STRING, optional::mask_file
    UNKNOWN_STRING::stack_option, title, map_file, mean_image_file, spots_file
    logical,optional::headless_vector
    UNKNOWN_STRING, optional::caption
    STRING the_caption, the_title
    integer,optional::m_filter(:)
    integer,parameter::n_threads = 8  !!must be >=2
    integer::mhs 
    real,dimension(:,:,:),allocatable:: image, uimage
    real(dl),dimension(:),allocatable::theta, phi, angle_rotate
    real(dl) norm, thislen, x, y, xshift, yshift
    real  rot
    integer n, nblocks, pix, i,   j, k, space, ithread,  nspots, imf
    real nstack(n_threads)
    real(dl) rmax, r_resolution
    type(hg_disc) disc(n_threads)
    type(file_pointer) fp
    real,dimension(:),allocatable::xstart, ystart, xend, yend
    type(hg_maps) mask, map
    logical do_mask
    real,optional::smooth_fwhm
    real sigma
    integer nw
    real,dimension(:,:),allocatable::wrap_image, window
    STRING mpost

    call hg_read_map(map_file, map)
    call hg_convert_to_ring(map)
    do_mask = .false.
    if(present(mask_file))then
       if(trim(mask_file).ne."")then
          call hg_read_map(mask_file, mask, nmaps_wanted = 1)
          call hg_convert_to_ring(mask)
          if(map%nside .ne. mask%nside) stop "hg_stack_io: mask must have the same resolution"
          do_mask = .true.
       endif
    endif
    disc%nside = map%nside
    

    if(trim(stack_option).ne. "T" .and. trim(stack_option).ne. "E" .and. trim(stack_option).ne. "B" .and. map%iq.eq. 0)then
       if(map%nmaps .lt. 2)then
          write(*,*) "hg_stack_io: the input map for stacking does not contain Q, U component"
       else
          write(*,*) "hg_stack_io: you have to specify the Q, U components in the input map by naming the file as *IQU*.fits (the 2nd, 3rd maps are Q, U maps) or *QU*.fits (the first two maps are Q, U maps)"
       endif
       stop
    endif

    if((trim(stack_option).eq."Qr" .or. trim(stack_option).eq."QR" .or. trim(stack_option).eq."Q").and. present(headless_vector))then
       if(headless_vector)then
          nblocks = 2  !!want headless vectors, too
       else
          nblocks = 1
       endif
    else
       nblocks = 1
    endif
    n = max(nint(rmax/r_resolution), 5)
    allocate(image(-n:n, -n:n, n_threads))
    nstack = 0
    image = 0.
    if(nblocks.ne.1)then
       allocate(uimage(-n:n,-n:n,n_threads))
       uimage = 0.
    endif
    if(.not. file_exists(spots_file))then
       write(*,*) trim(spots_file)
       write(*,*) "File does not exist"
       stop
    endif
    nspots = file_numlines(spots_file)
    allocate(theta(nspots), phi(nspots), angle_rotate(nspots))
    fp = open_file(trim(spots_file), "r")
    do i=1, nspots
       read(fp%unit, *) theta(i), phi(i), angle_rotate(i)       
    enddo
    call close_file(fp)
    if(nblocks.eq.1)then
       !$omp parallel do private(i, ithread, pix) 
       do ithread = 1, n_threads
          do i = ithread, nspots, n_threads
             call ang2pix_ring(map%nside, theta(i), phi(i), pix)
             call hg_get_disc(map%nside, pix, disc(ithread))
             if(do_mask)then
                call hg_stack(map, disc(ithread), n, r_resolution, angle_rotate(i), stack_option, image(:,:,ithread), mask = mask, counter =  nstack(ithread))
             else
                call hg_stack(map, disc(ithread), n, r_resolution, angle_rotate(i), stack_option, image(:,:,ithread), counter = nstack(ithread))
             endif
          enddo
       enddo
       !$omp end parallel do
    else
       !$omp parallel do private(i, ithread, pix)
       do ithread = 1, n_threads
          do i = ithread, nspots, n_threads
             call ang2pix_ring(map%nside, theta(i), phi(i), pix)
             call hg_get_disc(map%nside, pix, disc(ithread))
             if(do_mask)then
                call hg_stack( map, disc(ithread), n, r_resolution, angle_rotate(i), stack_option, image(:,:,ithread), uimage(:,:, ithread), mask = mask, counter =  nstack(ithread))
             else
                call hg_stack( map, disc(ithread), n, r_resolution, angle_rotate(i), stack_option, image(:,:,ithread), uimage(:,:, ithread), counter = nstack(ithread))
             endif
          enddo
       enddo
       !$omp end parallel do
    endif
    deallocate(theta, phi, angle_rotate)    

    if(sum(nstack) .gt.0)then
       do i=-n, n
          do j=-n, n
             image(i, j, 1) = sum(image(i,j,:))/real(sum(nstack))
          enddo
       enddo
       if(nblocks.ne.1)then
          do i=-n, n
             do j=-n, n
                uimage(i, j, 1) = sum(uimage(i,j,:))/real(sum(nstack))
             enddo
          enddo
       endif
    else
       write(*,*) "no spots are listed in file "//trim(spots_file)
       stop
    endif

    if(present(smooth_fwhm))then
       sigma = smooth_fwhm * const_sigmabyfwhm
       nw = ceiling((sigma*3.)/r_resolution)
       if(nw.gt.1)then
          allocate(window(-nw:nw, -nw:nw), wrap_image(-n-nw:n+nw, -n-nw:n+nw))
          wrap_image = 0
          sigma = (r_resolution/sigma)**2/2.
          do i=-nw, nw
             do j=-nw, nw
                window(i,j) = exp(-(i**2+j**2)*sigma)
             enddo
          enddo
          window = window/sum(window)
          wrap_image(-n:n, -n:n) = image(-n:n, -n:n, 1)
          !$omp parallel do private(i,j)
          do i=-n, n
             do j=-n, n
                image(i, j, 1) = sum(window*wrap_image(i-nw:i+nw, j-nw:j+nw))
             enddo
          enddo
          !$omp end parallel do
          if(nblocks.ne.1)then
             wrap_image(-n:n, -n:n) = uimage(-n:n, -n:n, 1)
             !$omp parallel do private(i,j)
             do i=-n, n
                do j=-n, n
                   uimage(i, j, 1) = sum(window*wrap_image(i-nw:i+nw, j-nw:j+nw))
                enddo
             enddo
             !$omp end parallel do
          endif
          deallocate(window, wrap_image)
       endif
    endif
    
    write(*,*) "Stacked on "//trim(num2str(nint(sum(nstack))))//" spots"

    if(present(caption))then
       if(caption(1:1).eq."#")then
          the_caption = trim(num2str(nint(sum(nstack))))//" patches on "//trim(caption(2:))
       else
          the_caption = trim(caption)
       endif
    else
       the_caption =  trim(num2str(nint(sum(nstack))))//" patches stacked"
    endif

    if(nblocks .ne. 1)then
       space = nint((n - 1)/7.)
       mhs = floor(real(n)/space - 0.5)
       allocate(xstart((2*mhs+1)**2), ystart((2*mhs+1)**2), xend((2*mhs+1)**2), yend((2*mhs+1)**2))
    endif

    call do_write_file(trim(mean_image_file))
    if(stack_option.eq."Q" .and. present(m_filter) .and. nblocks.ne.1)then
       image(:,:,2) = image(:,:,1)
       uimage(:,:,2) = uimage(:,:,1)
       do imf = 1, size(m_filter)
          call map_filter_modes(n, image(:,:,1), uimage(:,:,1), m_filter(imf:imf))
          call do_write_file(trim(file_add_postfix(mean_image_file, "_m="//trim(num2str(m_filter(imf))))))
          image(:,:,1) = image(:,:,2)
          uimage(:,:,1) = uimage(:,:,2)

          if(imf .eq. 1)then
             mpost = "_m="//trim(num2str(m_filter(imf)))
          else
             mpost=trim(mpost)//"-"//trim(num2str(m_filter(imf)))
          endif
       enddo
       if(size(m_filter).gt.1)then
          call map_filter_modes(n, image(:,:,1), uimage(:,:,1), m_filter)
          call do_write_file(trim(file_add_postfix(mean_image_file, trim(mpost))))
       endif
    endif
    if(nblocks.ne.1) deallocate(uimage, xstart, ystart, xend, yend)    
    deallocate(image) 
    call hg_free_map(map)
    call hg_free_map(mask)

  contains 

    subroutine do_write_file(fname)
      UNKNOWN_STRING fname
      fp = open_file(trim(fname),"w")
      call asy_init(fp, caption = trim(the_caption),  xlabel = "$2\sin(\theta/2)\cos\phi$", ylabel = "$2\sin(\theta/2)\sin\phi$", nblocks = nblocks) 
      call asy_density(fp, image(-n:n,-n:n,1), real(-r_resolution*n), real(r_resolution*n), real(-r_resolution*n), real(r_resolution*n), trim(title))
      if(nblocks .ne. 1)then !!headless vectors
         norm = r_resolution*space/sqrt(maxval(image(:,:,1)**2+uimage(:,:,1)**2))/2. * 0.975  !!*0.95 to avoid overlap
         k = 1
         do i = -mhs*space, mhs*space, space
            do j = -mhs*space, mhs*space, space
               x = i*r_resolution
               y = j*r_resolution
               rot = 0.5*POLAR_ANGLE(image(i, j,1), uimage(i, j,1)) 
               if(stack_option .eq. "Qr" .or. stack_option .eq. "QR") rot = rot + POLAR_ANGLE(x, y)
               thislen = norm*sqrt(uimage(i,j,1)**2 + image(i,j,1)**2)
               xshift = thislen*cos(rot)
               yshift = thislen*sin(rot)
               xstart(k) =  x - xshift
               ystart(k) = y - yshift
               xend(k) = x + xshift
               yend(K) = y + yshift
               k = k + 1
            enddo
         enddo
         call asy_lines(fp, xstart, ystart, xend, yend, "black", "solid", 2.)
      endif
      call close_file(fp)
    end subroutine do_write_file

  end subroutine hg_stack_io

  subroutine hg_smooth_mapfile(mapfile, filter_fwhm)
    UNKNOWN_STRING mapfile
    type(hg_maps) map
    real(dl) filter_fwhm
    call hg_read_map(mapfile, map)
    call hg_smooth_map(map, filter_fwhm)
    call hg_write_map(trim(file_add_postfix(trim(mapfile),"_smoothed_fwhm"//trim(num2str(nint(filter_fwhm/const_arcmin)))//"arcmin")), map)
    call hg_free_map(map)
  end subroutine hg_smooth_mapfile

  subroutine hg_smooth_map(map, filter_fwhm)
    type(hg_maps) map
    real(dl) filter_fwhm
    integer lmax
    lmax = min(ceiling(3./max(abs(filter_fwhm)*const_sigmabyfwhm, 1.d-6)), map%nside*3)
    if((lmax*filter_fwhm*const_sigmabyfwhm).lt. 0.02) return
    write(*,*) "Smoothing with lmax = ", lmax
    call hg_map2alm(map, lmax)
    call hg_filter_alm(map, fwhm = real(filter_fwhm))
    call hg_alm2map(map)
  end subroutine hg_smooth_map


  


  subroutine hg_getQU(Emap_file, QUmap_file)
    UNKNOWN_STRING Emap_file, QUmap_file
    type(hg_maps) hge, hgqu
    call hg_read_map(Emap_file, hge, nmaps_wanted = 1)
    call hg_map2alm(hge)
    call hg_ini_map(hgqu, nside = hge%nside, nmaps = 2, spin =  (/ 2, 2 /), lmax=hge%lmax)
    hgqu%alm(:, :, 1) = hge%alm(:, :, 1)
    hgqu%alm(:, :, 2) = 0
    call hg_alm2map(hgqu)
    call hg_write_map(QUmap_file, hgqu)
    call hg_free_map(hgqu)
    call hg_free_map(hge)
  end subroutine hg_getQU

  subroutine hg_export_spots(map_file, spots_file, spot_type, threshold, mask_file, filter_fwhm)
    UNKNOWN_STRING map_file, spots_file, spot_type
    UNKNOWN_STRING, optional::mask_file
    real(dl),optional::threshold
    type(file_pointer) fp
    real(dl) theta, phi, rotate_angle, fcut
    type(hg_maps) mask, map
    real(dl),optional::filter_fwhm
    real(dl) total_weight
    integer i, iq, iu, j
    integer nneigh, list(8)
    logical do_mask
    call hg_read_map(trim(map_file), map)
    if(present(filter_fwhm)) call hg_smooth_map(map, filter_fwhm)
    do_mask = .false.
    if(present(mask_file))then
       if(trim(mask_file).ne."")then
          call hg_read_map(mask_file, mask, nmaps_wanted = 1)
          do i=1, map%nmaps
             map%map(:,i) = map%map(:,i)*mask%map(:,1)
          enddo
          if(mask%nside .ne. map%nside) stop "hg_export_spots: mask and map must have the same nside"
          total_weight = sum(mask%map(:,1))
          do_mask = .true.
       else
          total_weight = map%npix
       endif
    else
       total_weight = map%npix
    endif
    call hg_convert_to_nested(map)
    if(do_mask)call hg_convert_to_nested(mask)

    fp = open_file(trim(spots_file),"w")
    select case(trim(spot_type))
    case("Tmax", "Emax", "Bmax")  !!random orientation
       if(present(threshold))then
          fcut = threshold*sqrt(sum(map%map(:,1)**2)/total_weight)
       else
          fcut = -1.e30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if( map%map(i,1) .lt. fcut .or. mask%map(i, 1) .le. 0.5 ) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if ( all(map%map(list(1:nneigh),1).lt.map%map(i,1)) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                call random_number(rotate_angle)
                rotate_angle = rotate_angle*const_2pi
                call pix2ang_nest(map%nside, i, theta, phi)
                write(fp%unit, "(3E16.7)") theta, phi, rotate_angle
             endif
          else
             if(map%map(i,1).lt.fcut )cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh),1).lt. map%map(i,1)))then
                call random_number(rotate_angle)
                rotate_angle = rotate_angle*const_2pi
                call pix2ang_nest(map%nside, i, theta, phi)
                write(fp%unit, "(3E16.7)") theta, phi, rotate_angle
             endif
          endif
       enddo
    case("Tmin", "Emin", "Bmin")  !!random orientation
       if(present(threshold))then
          fcut = - threshold*sqrt(sum(map%map(:,1)**2)/total_weight)
       else
          fcut = 1.e30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if(map%map(i,1).gt.fcut .or. mask%map(i, 1) .le. 0.5)cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if( all(map%map(list(1:nneigh),1) .gt. map%map(i,1)) .and. all(mask%map(list(1:nneigh), 1) .gt. 0.5 ) ) then
                call random_number(rotate_angle)
                rotate_angle = rotate_angle*const_2pi
                call pix2ang_nest(map%nside, i, theta, phi)
                write(fp%unit, "(3E16.7)") theta, phi, rotate_angle
             endif
          else
             if(map%map(i,1).gt.fcut )cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh),1).gt.map%map(i,1)))then
                call random_number(rotate_angle)
                rotate_angle = rotate_angle*const_2pi
                call pix2ang_nest(map%nside, i, theta, phi)
                write(fp%unit, "(3E16.7)") theta, phi, rotate_angle
             endif
          endif
       enddo
    case("TQUmax") !!oriented with QU, minima of T
       if(map%nmaps .lt. 3) stop "For TQUmax mode you need 3 maps"
       if(present(threshold))then
          fcut = threshold*sqrt(sum(map%map(:,1)**2)/total_weight)
       else
          fcut = -1.e30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if( map%map(i,1) .lt. fcut .or. mask%map(i,1) .le. 0.5)cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh),1) .lt. map%map(i,1)) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5)) then
                call pix2ang_nest(map%nside, i, theta, phi)
                rotate_angle = POLAR_ANGLE(map%map(i, 2), map%map(i, 3))/2.
                write(fp%unit, "(3E16.7)") theta, phi, rotate_angle
             endif
          else
             if(map%map(i,1) .lt. fcut)cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh),1) .lt. map%map(i,1)))then
                call pix2ang_nest(map%nside, i, theta, phi)
                rotate_angle = POLAR_ANGLE(map%map(i, 2), map%map(i, 3))/2.
                write(fp%unit, "(3E16.7)") theta, phi, rotate_angle
             endif
          endif
       enddo
    case("TQUmin") !!oriented with QU, minima of T
       if(map%nmaps .lt. 3) stop "For TQUmin mode you need 3 maps"
       if(present(threshold))then
          fcut = -threshold*sqrt(sum(map%map(:,1)**2)/total_weight)
       else
          fcut = 1.e30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if( map%map(i,1) .gt. fcut .or. mask%map(i,1) .le. 0.5) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if( all( map%map(list(1:nneigh), 1) .gt. map%map(i,1) ) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5) ) then
                call pix2ang_nest( map%nside, i, theta, phi )
                rotate_angle = POLAR_ANGLE(map%map(i, 2), map%map(i, 3))/2.
                write(fp%unit, "(3E16.7)") theta, phi, rotate_angle
             endif
          else
             if(map%map(i,1) .gt. fcut)cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh),1).gt. map%map(i,1)))then
                call pix2ang_nest(map%nside, i, theta, phi)
                rotate_angle = POLAR_ANGLE(map%map(i, 2), map%map(i, 3))/2.
                write(fp%unit, "(3E16.7)") theta, phi, rotate_angle
             endif
          endif
       enddo
    case("Pmax")
       select case(map%nmaps)
       case(2:3)
          iq = map%nmaps - 1
          iu = iq + 1
       case default
          stop "For polarization you need to specify Q, U maps in the input file."
       end select
       if(present(threshold))then
          fcut = threshold**2*(sum(map%map(:,iq)**2 + map%map(:,iu)**2)/total_weight)
       else
          fcut = -1.e30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if( map%map(i,iq)**2 + map%map(i,iu)**2 .lt. fcut  .or.  mask%map(i, 1) .le. 0.5 ) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if( all( map%map(list(1:nneigh),iq)**2 + map%map(list(1:nneigh),iu)**2 .lt. map%map(i,iq)**2 + map%map(i, iu)**2 ) .and. all(mask%map(list(1:nneigh), 1) .gt. 0.5) ) then
                rotate_angle = POLAR_ANGLE(map%map(i, iq), map%map(i, iu))/2.
                call pix2ang_nest(map%nside, i, theta, phi)
                write(fp%unit, "(3E16.7)") theta, phi, rotate_angle
             endif
          else
             if( map%map(i, iq)**2 + map%map(i, iu)**2 .lt. fcut) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if( all( map%map(list(1:nneigh), iq)**2 + map%map(list(1:nneigh), iu)**2 .lt. map%map(i, iq)**2 + map%map(i, iu)**2 ) )then
                rotate_angle = POLAR_ANGLE(map%map(i, iq), map%map(i, iu))/2.
                call pix2ang_nest(map%nside, i, theta, phi)
                write(fp%unit, "(3E16.7)") theta, phi, rotate_angle
             endif
          endif
       enddo
    case("Pmin")
       select case(map%nmaps)
       case(2:3)
          iq = map%nmaps - 1
          iu = iq + 1
       case default
          stop "For polarization you need to specify Q, U maps in the input file."
       end select
       if(present(threshold))then
          fcut = threshold**2*(sum(map%map(:, iq)**2+map%map(:, iu)**2)/total_weight)
       else
          fcut = 1.d30
       endif
       do i=0, map%npix-1
          if(do_mask)then
             if(map%map(i,iq)**2 + map%map(i,iu)**2 .gt. fcut .or. mask%map(i, 1) .le. 0.5 ) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh), iq)**2 + map%map(list(1:nneigh), iu)**2 .gt. map%map(i, iq)**2 + map%map(i, iu)**2) .and. all(mask%map(list(1:nneigh), 1) .gt. 0.5) ) then
                rotate_angle = POLAR_ANGLE(map%map(i,iq), map%map(i,iu))/2.
                call pix2ang_nest(map%nside, i, theta, phi)
                write(fp%unit, "(3E16.7)") theta, phi, rotate_angle
             endif
          else
             if(map%map(i, iq)**2 + map%map(i, iu)**2 .gt. fcut) cycle
             call neighbours_nest(map%nside, i, list, nneigh)  
             if(all(map%map(list(1:nneigh), iq)**2 + map%map(list(1:nneigh), iu)**2 .gt. map%map(i, iq)**2 + map%map(i, iu)**2))then
                rotate_angle = POLAR_ANGLE(map%map(i, iq), map%map(i, iu))/2.
                call pix2ang_nest(map%nside, i, theta, phi)
                write(fp%unit, "(3E16.7)") theta, phi, rotate_angle
             endif
          endif
       enddo
    case default
       write(*,*) trim(spot_type)
       stop "unknown spot type"
    end select
    call hg_free_map(map)
    call hg_free_map(mask)
    call close_file(fp)
  end subroutine hg_export_spots


  subroutine hg_output_map(map, header, fname)
    real, dimension(:,:):: map
    character(LEN=80),dimension(:):: header
    UNKNOWN_STRING fname
    call delete_file(trim(fname))
    call output_map(map, header, fname)
  end subroutine hg_output_map
  
  subroutine hg_mask_map(mapfile, maskfile, output, index_list)
    type(hg_maps) map, mask
    integer i
    UNKNOWN_STRING mapfile, maskfile, output
    integer,dimension(:),optional::index_list
    call hg_read_map(mapfile, map)
    call hg_read_map(maskfile, mask, nmaps_wanted = 1)
    if(mask%ordering .eq. RING)then
       call hg_convert_to_ring(map)
    elseif(mask%ordering .eq. NESTED)then
       call hg_convert_to_nested(map)
    else
       write(*,*) "hg_mask_map: mask file has unknown ordering."
       stop
    endif
    call delete_file(output)
    if(present(index_list))then
       if(any(index_list .lt. 1 .or. index_list .gt. map%nmaps)) stop "hg_write_map: index out of range"
       do i=1, size(index_list)
          map%map(:,index_list(i)) = map%map(:, index_list(i))*mask%map(:,1)
       enddo
       call hg_write_map(trim(file_add_postfix(output,"_masked")), map, index_list)
       call hg_write_map(output, map)
    else
       do i=1, map%nmaps
          map%map(:,i) = map%map(:, i)*mask%map(:,1)
       enddo
       call hg_write_map(output, map)
    endif
    call hg_free_map(map)
    call hg_free_map(mask)
  end subroutine hg_mask_map




  subroutine hg_smooth_maskfile(mask_file, smoothscale, output)
    UNKNOWN_STRING mask_file
    real smoothscale
    UNKNOWN_STRING, optional::output
    type(hg_maps) hgm
    call hg_read_map(mask_file, hgm, nmaps_wanted = 1)
    call hg_smooth_mask(hgm, smoothscale)
    if(present(output))then
       call hg_write_map(trim(output), hgm)
    else
       call hg_write_map(trim(file_add_postfix(trim(mask_file), "_smoothed")), hgm)
    endif
    call hg_free_map(hgm)
  end subroutine hg_smooth_maskfile

  subroutine hg_smooth_mask(hgm, smoothscale)
    real,parameter::nefolds = 4
    type(hg_maps) hgm, hgs
    integer,dimension(:),allocatable::listpix
    integer i, j, nsteps
    real smoothscale, decay
    nsteps = ceiling(smoothscale*hgm%nside*nefolds/2.)
    if(nsteps .le. 0 .or. nsteps .gt. 200)stop "hg_smooth_mask: invalid input of smoothscale"
    decay = exp(-nefolds/nsteps/2.)
    call hg_convert_to_nested(hgm)
    hgm%mask_npix = count(hgm%map(:,1) .lt. 1.)
    allocate(hgm%mask_listpix(hgm%mask_npix))
    j = 0
    do i = 0, hgm%npix - 1
       if(hgm%map(i,1).lt. 1.)then
          j = j + 1
          hgm%mask_listpix(j) = i
       endif
    enddo
    call hg_copy_map(hgm, hgs)
    do i=1, nsteps
       call hg_iterate_mask(hgm, hgs, decay)
       call hg_iterate_mask(hgs, hgm, decay)
    enddo
    call hg_free_map(hgs)

  contains 

    subroutine hg_iterate_mask(hgm_from, hgm_to, decay)  
      type(hg_maps) hgm_from, hgm_to
      integer list(8), nneigh
      integer i
      real decay
      !$omp parallel do private(list, nneigh, i)
      do i = 1, hgm%mask_npix
         call neighbours_nest(hgm_from%nside, hgm%mask_listpix(i), list, nneigh)
         hgm_to%map(hgm_from%mask_listpix(i), 1) = max(maxval(hgm_from%map(list(1:nneigh), 1)) * decay , hgm_from%map(hgm_from%mask_listpix(i), 1))
      enddo
      !$omp end parallel do
    end subroutine hg_iterate_mask
  end subroutine hg_smooth_mask


  subroutine hg_inpainting(mode, map_file, mask_file, maskpol_file)
    UNKNOWN_STRING map_file, mask_file, mode
    UNKNOWN_STRING,optional:: maskpol_file
    integer,parameter:: output_steps = 50, total_steps = 1000, burnin = 20
    integer step, naccept, weight
    logical accept
    type(hg_maps) map, simumap, mask, maskpol, mapmean
    if(trim(mode) .eq. "I")then
       call hg_read_map(map_file, map, nmaps_wanted = 1)
    else
       call hg_read_map(map_file, map, nmaps_wanted = 3)
    endif
    call hg_copy_map(map, mapmean)
    mapmean%map = 0
    call hg_read_map(mask_file, mask, nmaps_wanted = 1)
    if(all( mask%map .eq. 1. .or. mask%map .eq. 0.))then
          write(*,*) "smoothing the temperature mask "//trim(mask_file)
       call hg_smooth_mask(mask, 2.*real(const_degree))
    endif
    map%mask_npix = floor(map%npix - sum(mask%map(:,1)))
    if(present(maskpol_file)) then
       call hg_read_map(maskpol_file, maskpol, nmaps_wanted = 1)
       if(all(maskpol%map .eq. 1. .or. maskpol%map .eq. 0.))then
          write(*,*) "smoothing the polarization mask "//trim(maskpol_file)
          call hg_smooth_mask(maskpol, 2.*real(const_degree))
       endif
       map%maskpol_npix = floor(map%npix - sum(maskpol%map(:,1)))
    else
       map%maskpol_npix = 0
    endif
    call hg_inpainting_init(map, simumap)
    naccept = 0
    weight = 0
    step = 1
    do while(step .lt. total_steps)
       if(present(maskpol_file))then
          call hg_inpainting_step(accept, map, simumap, mask, maskpol)
       else
          call hg_inpainting_step(accept, map, simumap, mask)
       endif
       if(accept) naccept = naccept + 1
       if(naccept .eq. burnin/2 .and. accept)then
          map%Cl = simumap%Cl
       endif
       if(naccept .ge. burnin)then
          mapmean%map = mapmean%map + map%map
          if(weight.eq.0)then
             write(*,*) "initial trials done: resetting step # = 1"
             step = 1
          endif
          weight = weight + 1
       endif
       write(*,*) step, accept, "temperature = ", map%mcmc_temperature, "chisq = ", map%chisq, simumap%chisq
       if(mod(step, output_steps).eq.0 )then
          call hg_write_map(trim(file_add_postfix(trim(map_file), "_inp"//trim(str_ndigits(step, 4)))), map)
          if(weight .gt. 0)then
             simumap%map =  mapmean%map/weight
             call hg_write_map(trim(file_add_postfix(trim(map_file), "_mean"//trim(str_ndigits(step, 4)))), simumap)
          endif
       endif
    enddo
    call hg_free_map(map)
    call hg_free_map(mapmean)
    call hg_free_map(simumap)
    call hg_free_map(mask)
    call hg_free_map(maskpol)
  end subroutine hg_inpainting

  subroutine hg_inpainting_init(map, simumap)
    type(hg_maps) map, simumap
    call hg_map2alm(map)
    map%cl(0:1, :) = 0.
    map%cl(:, hg_index_TT) = max(map%cl(:, hg_index_TT), 1.e-8)*(map%npix/real(map%npix - map%mask_npix))
    map%cl(:, hg_index_EE) = max(map%cl(:, hg_index_EE), 1.e-10)*(map%npix/real(map%npix-map%maskpol_npix))
    map%cl(:, hg_index_TE) = map%cl(:, hg_index_TE)*sqrt((map%npix/real(map%npix-map%maskpol_npix))*(map%npix/real(map%npix-map%mask_npix))*0.9999) !!0.9999 is to avoid 0 determinant
    map%cl(:, hg_index_BB) = max(map%cl(:, hg_index_BB), 1.e-12)*(map%npix/real(map%npix-map%maskpol_npix))
    
    map%cl(:, hg_index_EB) = 0. 
    map%Cl(:, hg_index_TB) = 0.
    map%Cl(0:1,:) = 0.
    map%mcmc_temperature = 20.     !!start with a high temperature
    hg_inpainting_lowl = 5
    call hg_copy_map(map, simumap)
    call hg_inpainting_get_chisq(map, simumap)
    map%chisq = simumap%chisq
  end subroutine hg_inpainting_init
  
  subroutine hg_inpainting_step(accept, map, simumap, mask, maskpol)
    type(hg_maps) map, simumap, mask
    type(hg_maps),optional::maskpol
    logical accept
    integer i
    simumap%Cl = map%Cl
    call hg_simulate(simumap)
    simumap%map(:,1) =  map%map(:,1) * mask%map(:,1) + simumap%map(:,1) * sqrt(1.-mask%map(:,1)**2)
    if(map%nmaps.eq.3)then
       if(present(maskpol))then
          simumap%map(:,2) =  map%map(:,2) * maskpol%map(:,1) + simumap%map(:,2) * sqrt(1.-maskpol%map(:,1)**2)
          simumap%map(:,3) =  map%map(:,3) * maskpol%map(:,1) + simumap%map(:,3) * sqrt(1.-maskpol%map(:,1)**2)
       else
          stop "hg_inpainting_step: need polarization mask"
       endif
    endif
    call hg_inpainting_accept_reject(accept, map, simumap)
  end subroutine hg_inpainting_step


  subroutine hg_inpainting_get_chisq(map, simumap)
    type(hg_maps) map, simumap
    integer l
    real chisq
    call hg_map2alm(simumap)
    if(map%nmaps .eq. 1)then
       chisq = 0.
       !$omp parallel do reduction(+:chisq)
       do l = 2, hg_inpainting_lowl 
          chisq = chisq + (2*l+1)*(simumap%Cl(l,1)/map%Cl(l,1) - log(simumap%Cl(l,1)/map%Cl(l,1))-1.d0)
       enddo
       !$omp end parallel do
       simumap%chisq = chisq 

    else
       chisq = 0.
       !$omp parallel do reduction(+:chisq)
       do l = 2, hg_inpainting_lowl 
          chisq = chisq + (2*l+1)*(simumap%Cl(l,hg_index_BB)/map%Cl(l,hg_index_BB) - log(simumap%Cl(l,hg_index_BB)/map%Cl(l,hg_index_BB))-3.d0 - log((simumap%Cl(l, hg_index_TT)*simumap%Cl(l, hg_index_EE) - simumap%Cl(l, hg_index_TE)**2)/(map%Cl(l, hg_index_TT)*map%Cl(l, hg_index_EE) - map%Cl(l, hg_index_TE)**2)) + (map%Cl(l, hg_index_TT)*simumap%Cl(l, hg_index_EE) + simumap%Cl(l, hg_index_TT)*map%Cl(l, hg_index_EE) - 2.*map%Cl(l, hg_index_TE)*simumap%Cl(l, hg_index_TE) )/(map%Cl(l, hg_index_TT)*map%Cl(l, hg_index_EE) - map%Cl(l, hg_index_TE)**2))
       enddo
       !$omp end parallel do
       simumap%chisq = chisq 
    endif
  end subroutine hg_inpainting_get_chisq

  subroutine hg_inpainting_accept_reject(accept, map, simumap)
    type(hg_maps) map, simumap
    logical accept
    call hg_inpainting_get_chisq(map, simumap)
    if(simumap%chisq - map%chisq .lt. random_exp()*2*map%mcmc_temperature)then
       accept = .true.
       map%chisq = simumap%chisq
       map%map = simumap%map
       map%alm = simumap%alm
       map%mcmc_temperature = max(1., map%mcmc_temperature*0.95)
       hg_inpainting_lowl = min(20, hg_inpainting_lowl + 1)
    else
       accept = .false.
       map%mcmc_temperature = map%mcmc_temperature*1.01
       hg_inpainting_lowl = max(5, hg_inpainting_lowl - 1)
    endif
  end subroutine hg_inpainting_accept_reject


  subroutine hg_split_map(mapfile)
    UNKNOWN_STRING mapfile
    type(hg_maps) map
    integer i
    call hg_read_map(mapfile, map)
    do i=1, map%nmaps
       call hg_write_map(trim(file_add_postfix(mapfile, "_submap"//trim(str_ndigits(i, 3)))), map, (/ i /) )
    enddo
    call hg_free_map(map)
  end subroutine hg_split_map

  subroutine hg_plot_spots(spotsfile, mapfile)
    UNKNOWN_STRING spotsfile, mapfile
    type(file_pointer) fp
    real(dl) theta, phi, angle_rotate
    integer pix
    type(hg_maps) hgm
    call hg_ini_map(hgm, 64, 1, (/ 0 /) )
    hgm%map = 0
    fp = open_file(trim(spotsfile), "r")
    do
       read(fp%unit, *, END=100, ERR=100) theta, phi, angle_rotate
       call ang2pix_ring(hgm%nside, theta, phi, pix)
       hgm%map(pix, 1) = 1.d0
    enddo
100 call close_file(fp)
    call hg_write_map(trim(mapfile), hgm)
    call hg_free_map(hgm)
  end subroutine hg_plot_spots

end module healpix_geometry
