module coop_fitsio_mod
  use coop_wrapper_typedef
  use coop_MPI_mod
  use coop_list_mod
  use coop_random_mod
  use coop_special_function_mod
  use coop_cholesky_mod  
  use coop_matrix_mod
  use coop_interpolation_mod
  use coop_integrate_mod
  use coop_ode_mod
  use coop_file_mod
  use coop_asy_mod
  use coop_fft_mod
  use coop_jl_mod
  use coop_gaussian_peak_stat_mod
  use coop_nd_prob_mod
  use coop_evalstr_mod
  use coop_sphere_mod
  implicit none

#include "constants.h"
#define COOP_FITSIO_CARD character(LEN=80)
  COOP_INT,parameter::coop_fitsio_image_hdu = 0
  COOP_INT,parameter::coop_fitsio_ascII_table_hdu = 1
  COOP_INT,parameter::coop_fitsio_binary_table_hdu = 2
  COOP_INT,parameter::coop_fitsio_datatype_short_int = 21
  COOP_INT,parameter::coop_fitsio_datatype_int = 41
  COOP_INT,parameter::coop_fitsio_datatype_single = 42
  COOP_INT,parameter::coop_fitsio_datatype_double = 82
!!not used
  COOP_INT,parameter::coop_fitsio_datatype_byte = 11
  COOP_INT,parameter::coop_fitsio_datatype_logical = 14
  COOP_INT,parameter::coop_fitsio_datatype_character = 16
  COOP_INT,parameter::coop_fitsio_datatype_complex = 83
  COOP_INT,parameter::coop_fitsio_datatype_double_complex = 163

  type coop_fits_file
     COOP_INT::unit = 0
     COOP_STRING::filename = ""
     COOP_INT::nhdus = 0
     COOP_INT::chdu = 0
     COOP_INT::status = 0
     COOP_INT::rwmode = 0 !!default readonly, set to 1 for readwrite
     COOP_INT::blocksize= 0
     COOP_INT::hdutype = 0
     type(coop_dictionary)::header
   contains
     procedure::init => coop_fits_file_init
     procedure::open => coop_fits_file_open
     procedure::open_image => coop_fits_file_open_image
     procedure::open_table => coop_fits_file_open_table
     procedure::close => coop_fits_file_close
     procedure::create_hdu => coop_fits_file_create_hdu
     procedure::create_binary_table => coop_fits_file_create_binary_table
     procedure::create_ascii_table => coop_fits_file_create_ascii_table
     procedure::create_image => coop_fits_file_create_image
     procedure::update_hdu_info => coop_fits_file_update_hdu_info
     procedure::move_to_hdu => coop_fits_file_move_to_hdu
     procedure::check_error => coop_fits_file_check_error
     procedure::report_error => coop_fits_file_report_error
     procedure::load_header => coop_fits_file_load_header
     procedure::add_key => coop_fits_file_add_key
     procedure::load_image_1d => coop_fits_file_load_image_1d
     procedure::load_image => coop_fits_file_load_image_1d
     procedure::load_image_2d => coop_fits_file_load_image_2d
     procedure::load_image_3d => coop_fits_file_load_image_3d
     procedure::load_double_column => coop_fits_file_load_double_column
     procedure::load_single_column => coop_fits_file_load_single_column
     procedure::load_int_column => coop_fits_file_load_int_column
     procedure::get_nrows_ncols => coop_fits_file_get_nrows_ncols
     procedure::get_naxes => coop_fits_file_get_naxes
     procedure::get_bitpix => coop_fits_file_get_bitpix
  end type coop_fits_file

  type coop_cls
     COOP_STRING::genre = ''
     COOP_STRING::unit='Unkown'
     COOP_INT::lmin = 0
     COOP_INT::lmax = -1
     COOP_INT::num_fields = 0
     COOP_INT::num_cls = 0
     COOP_REAL,dimension(:,:),allocatable::Cls
     COOP_INT,dimension(:),allocatable::spin
   contains
     procedure::free => coop_cls_free
     procedure::smooth => coop_Cls_smooth
     procedure::init => coop_cls_init
     procedure::load => coop_cls_load
     procedure::read => coop_cls_load
     procedure::dump => coop_cls_dump
     procedure::write => coop_cls_dump
     procedure::load_fits => coop_cls_load_fits
     procedure::dump_fits => coop_cls_dump_fits
     procedure::read_fits => coop_cls_load_fits
     procedure::write_fits => coop_cls_dump_fits
     procedure::load_txt => coop_cls_load_txt
     procedure::dump_txt => coop_cls_dump_txt
     procedure::read_txt => coop_cls_load_txt
     procedure::write_txt => coop_cls_dump_txt
     procedure::filter => coop_cls_filter
     procedure::select_maps => coop_cls_select_maps
  end type coop_cls

  type, extends(coop_cls)::coop_binned_cls
     COOP_INT::nb = 0
     COOP_REAL,dimension(:),allocatable::lb
     COOP_REAL,dimension(:,:),allocatable::Cbs, wb, wl
     COOP_INT,dimension(:),allocatable::lmin_b, lmax_b
     COOP_INT,dimension(:),allocatable::ibmin_l, ibmax_l
   contains
     procedure::alloc => coop_binned_Cls_alloc
     procedure::bin => coop_binned_Cls_bin
     procedure::unbin => coop_binned_Cls_unbin
  end type coop_binned_cls

  interface coop_fits_file_write_binary_table
     module procedure coop_fits_file_write_binary_table_d, coop_fits_file_write_binary_table_s
  end interface coop_fits_file_write_binary_table

  interface coop_fits_file_write_binary_table_with_indices
     module procedure coop_fits_file_write_binary_table_with_indices_s, coop_fits_file_write_binary_table_with_indices_d
  end interface coop_fits_file_write_binary_table_with_indices

contains

  subroutine coop_cls_smooth(this, delta_l)
    class(coop_cls)::this
    COOP_INT::l, delta_l, lp
    COOP_REAL,dimension(:,:),allocatable::cls_tmp
    COOP_REAL::w(-delta_l:delta_l)
    allocate(cls_tmp(this%lmin-delta_l:this%lmax+delta_l, this%num_cls))
    do l  = this%lmin, this%lmax
       cls_tmp(l, :) = this%cls(l, :)*(l+0.5)**2
    enddo
    do l = this%lmin-delta_l, this%lmin - 1
       cls_tmp(l, :) = cls_tmp(this%lmin, :)
    enddo
    do l = this%lmax+1, this%lmax+delta_l
       cls_tmp(l, :) = cls_tmp(this%lmax, :)
    enddo

    do l = -delta_l, delta_l
       w(l) = exp(- (2.d0*l/delta_l)**2)
    enddo
    w = w/sum(w)
    do l = this%lmin, this%lmax
       this%cls(l, :) = 0.d0
       do lp = - delta_l, delta_l
          this%cls(l,:) = this%cls(l,:) + w(lp)*cls_tmp(l+lp,:)
       enddo
       this%cls(l, :) = this%cls(l, :)/(l+0.5)**2
    enddo
    deallocate(cls_tmp)
  end subroutine coop_cls_smooth

  subroutine coop_cls_init(this, lmin, lmax, spin, cls, unit, genre)
    class(coop_cls)::this
    COOP_INT::lmin, lmax
    COOP_UNKNOWN_STRING::genre
    COOP_INT,optional::spin(:)
    COOP_REAL,optional::cls(:,:)
    COOP_UNKNOWN_STRING,optional::unit
    call this%free()
    this%lmin = lmin
    this%lmax = lmax
    this%genre = trim(adjustl(genre))
    this%num_fields = len_trim(this%genre)
    this%num_cls = this%num_fields*(this%num_fields+1)/2
    allocate( this%spin(this%num_fields), this%Cls(this%lmin:this%lmax, this%num_cls))
    if(present(cls))then
       if(size(cls, 1).ne. this%lmax-this%lmin+1 .or. size(cls,2).ne. this%num_cls) stop "cls_init: wrong size of cls array"
       this%cls = cls
    endif
    if(present(spin))then
       if(size(spin).ne. this%num_fields) stop "cls_init: wrong size of spin array"
       this%spin = spin
    endif
    if(present(unit))this%unit = trim(adjustl(unit))
  end subroutine coop_cls_init



  subroutine coop_binned_Cls_alloc(this, nb, ells)
    class(coop_binned_Cls)::this
    COOP_INT::nb
    COOP_REAL,dimension(:),optional::ells
    COOP_REAL,dimension(:,:), allocatable::mat
    COOP_REAL::delta_l,  det
    COOP_INT::i, ib, l, lp, np, nbsize
    if(this%lmax .lt. this%lmin .or. this%num_fields .eq. 0) stop "binned_Cls_alloc: you need to load cls first"
    if(nb .gt. this%lmax - this%lmin+1) stop "binned_Cls_alloc: too many bins (>lmax-lmin+1)"
    this%nb = nb
    COOP_DEALLOC(this%lb)
    COOP_DEALLOC(this%lmin_b)
    COOP_DEALLOC(this%lmax_b)
    COOP_DEALLOC(this%ibmin_l)
    COOP_DEALLOC(this%ibmax_l)
    COOP_DEALLOC(this%Cbs)
    COOP_DEALLOC(this%wb)
    COOP_DEALLOC(this%wl)
    allocate(this%lb(nb), this%Cbs(nb, this%num_cls), this%wb(this%lmin:this%lmax, this%nb), this%lmin_b(nb), this%lmax_b(nb), this%ibmin_l(this%lmin:this%lmax), this%ibmax_l(this%lmin:this%lmax), this%wl(this%nb, this%lmin:this%lmax))
    if(present(ells))then
       this%lb = ells(1:nb)
       call coop_quicksort(this%lb)
       if(this%lb(1) .lt. this%lmin .or. this%lb(this%nb) .gt. this%lmax)then
          stop "binned_Cls_alloc: input ell array exceeds the range [lmin, lmax]"
       endif
    else
       delta_l = (this%lmax - this%lmin)/this%nb
       if(delta_l .le. 3.d0 .or. nb .lt. 50)then
          call coop_set_uniform(this%nb, this%lb, this%lmin+delta_l/2.d0, this%lmax-delta_l/2.d0)
       else
          i = 1
          this%lb(i) = min(this%lmin * 1.025d0 + 0.5d0, this%lmin+delta_l/2.d0)
          do 
             delta_l = (this%lmax - this%lb(i))/(this%nb - i + 0.5)
             if(this%lb(i) .gt. delta_l*10.d0 .or. i .ge. nb*3/4)then
                call coop_set_uniform(this%nb - i, this%lb(i+1:this%nb), this%lb(i)+delta_l, this%lmax - delta_l/2.d0)
                exit
             endif
             this%lb(i+1) = max(this%lb(i)*1.1d0, this%lb(i) + 1.d0)
             i = i + 1
          enddo
       endif
    endif
    this%wb = 0.d0
    this%ibmin_l = this%nb
    this%ibmax_l = 1
    do ib = 1, this%nb 
       if(ib .eq. 1)then
          delta_l = this%lb(ib) - this%lmin+0.5d0
       elseif(ib.eq.this%nb)then
          delta_l = this%lmax - this%lb(ib)+0.5d0
       else
          delta_l = (this%lb(ib+1) - this%lb(ib-1))/2.d0
       end if
       this%lmin_b(ib) = max(this%lmin, floor(this%lb(ib) - 2.5d0*delta_l))
       this%lmax_b(ib) = min(this%lmax, floor(this%lb(ib) + 2.5d0*delta_l))
       do l = this%lmin_b(ib), this%lmax_b(ib)
          this%ibmin_l(l) = min(ib, this%ibmin_l(l))
          this%ibmax_l(l) = max(ib, this%ibmax_l(l))
          this%wb(l, ib) = exp(-((l-this%lb(ib))/delta_l)**2)
       enddo
       this%wb(this%lmin_b(ib):this%lmax_b(ib), ib) = this%wb(this%lmin_b(ib):this%lmax_b(ib), ib)/ sum(this%wb(this%lmin_b(ib):this%lmax_b(ib), ib))
    enddo
    this%wl = transpose(this%wb)
    do l = this%lmin, this%lmax
       this%wl(this%ibmin_l(l):this%ibmax_l(l), l) =   this%wl(this%ibmin_l(l):this%ibmax_l(l), l)/sum( this%wl(this%ibmin_l(l):this%ibmax_l(l), l))
    enddo
!!$    this%wl = 0.d0
!!$    do l = this%lmin, this%lmax
!!$       if(l .lt. this%lb(1))then
!!$          this%ibmin_l(l) = 1
!!$          this%ibmax_l(l) = 1
!!$          this%wl(1, l) = 1.d0
!!$          cycle
!!$       endif
!!$       if(this%ibmin_l(l) .lt. this%nb)then
!!$          do while(this%lb(this%ibmin_l(l)+1) .lt. l .and. this%ibmin_l(l) .lt. this%nb)
!!$             this%ibmin_l(l) = this%ibmin_l(l) + 1
!!$          enddo
!!$       endif
!!$       if(this%ibmin_l(l) .lt. this%nb)then
!!$          this%ibmax_l(l) = this%ibmin_l(l)+1
!!$          this%wl(this%ibmin_l(l), l) = (this%lb(this%ibmax_l(l)) - l)/(this%lb( this%ibmax_l(l)) - this%lb( this%ibmin_l(l) ))
!!$          this%wl(this%ibmax_l(l), l) = (l-this%lb(this%ibmin_l(l)))/(this%lb( this%ibmax_l(l) ) - this%lb( this%ibmin_l(l) ))
!!$       else
!!$          this%ibmax_l(l) = this%ibmin_l(l)
!!$          this%wl(this%ibmin_l(l), l) = 1.d0
!!$       endif
!!$    enddo
!!$    do l = this%lmin, this%lmax
!!$       allocate(mat(this%ibmin_l(l):this%ibmax_l(l), this%ibmin_l(l):this%ibmax_l(l)))
!!$       mat = 0.d0
!!$       delta_l = this%lb(this%ibmax_l(l)) - this%lb(this%ibmin_l(l))
!!$       do np = this%ibmin_l(l), this%ibmax_l(l)
!!$          do ib = this%ibmin_l(l), this%ibmax_l(l)
!!$             do lp = this%lmin_b(ib), this%lmax_b(ib)
!!$                mat(ib, np) = mat(ib,np) + ((lp-l)/delta_l)**(np-this%ibmin_l(l))*this%wb(lp, ib)
!!$             enddo
!!$          enddo
!!$       enddo
!!$       call coop_matrix_inverse(mat)
!!$       do ib = this%ibmin_l(l), this%ibmax_l(l)
!!$          this%wl(ib, l) = mat(this%ibmin_l(l), ib)
!!$       enddo
!!$       deallocate(mat)
!!$    enddo
!!$    do while(l.le. this%lmax)
!!$       this%ibmin_l(l) = this%nb-1
!!$       this%ibmax_l(l) = this%nb
!!$       this%wl(this%nb-1, l) = (this%lb(this%nb)-l)/(this%lb(this%nb) - this%lb(this%nb-1))
!!$       this%wl(this%nb, l) =  (l - this%lb(this%nb-1))/(this%lb(this%nb) - this%lb(this%nb-1))
!!$       l = l + 1
!!$    enddo
    
    !!put in the l(l+1) factors
    !$omp parallel do private(ib, l)
    do ib = 1, this%nb
       do l = this%lmin_b(ib), this%lmax_b(ib)
          this%wb(l, ib) = this%wb(l, ib)*((l+0.5d0)/(this%lb(ib)+0.5))**2
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(ib, l)
    do l = this%lmin, this%lmax
       do ib = this%ibmin_l(l), this%ibmax_l(l)
          this%wl(ib, l) = this%wl(ib, l)*( (this%lb(ib) + 0.5d0)/(l + 0.5d0) )**2
       enddo
    enddo
    !$omp end parallel do
  end subroutine coop_binned_Cls_alloc


  subroutine coop_binned_Cls_bin(this)
    class(coop_binned_Cls)::this
    COOP_INT::l, ib
    this%Cbs = 0.d0
    !$omp parallel do private(ib, l)
    do ib = 1, this%nb
       do l = this%lmin_b(ib), this%lmax_b(ib)
          this%Cbs(ib, :) = this%Cbs(ib,:) + this%Cls(l, :)*this%wb(l, ib)
       enddo
    enddo
    !$omp end parallel do
  end subroutine coop_binned_Cls_bin


  subroutine coop_binned_Cls_unbin(this)
    class(coop_binned_Cls)::this
    COOP_INT::l, ib
    this%Cls = 0.d0
    !$omp parallel do private(ib, l)
    do l = this%lmin, this%lmax
       do ib = this%ibmin_l(l), this%ibmax_l(l)
          this%Cls(l, :) = this%Cls(l,:) + this%Cbs(ib, :)*this%wl(ib, l)
       enddo
    enddo
    !$omp end parallel do
  end subroutine coop_binned_Cls_unbin


  function Coop_highpass_filter(l1, l2, l) result(w)
    COOP_INT l1, l2, l
    COOP_REAL w
    if(l.ge.l2)then
       w = 1.d0
       return
    endif
    if(l.le. l1)then
       w = 0.d0
       return
    endif
    w = sin(dble(l-l1)/dble(l2-l1)*coop_pio2)
  end function Coop_highpass_filter

  function coop_gaussian_filter(fwhm_arcmin, l) result(w)
    COOP_INT l
    COOP_REAL fwhm_arcmin, w
    w = exp(-((coop_sigma_by_fwhm*coop_SI_arcmin/coop_sqrt2)*fwhm_arcmin)**2*l*(l+1.d0))
  end function coop_gaussian_filter

  function coop_lowpass_filter(l1, l2, l) result(w)
    COOP_INT l1, l2, l
    COOP_REAL w
    if(l.le.l1)then
       w = 1.d0
       return
    endif
    if(l.ge. l2)then
       w = 0.d0
       return
    endif
    w = sin(dble(l2-l)/dble(l2-l1)*coop_pio2)
  end function coop_lowpass_filter


  subroutine coop_cls_select_maps(this, index_list)
    class(coop_cls)::this
    COOP_INT,dimension(:):: index_list
    COOP_REAL,dimension(:,:),allocatable::Cls_tmp, Cbs_tmp
    COOP_INT,dimension(:),allocatable::spin_tmp
    COOP_STRING::genre_tmp
    COOP_INT::n, i, j
    if(any(index_list .gt. this%num_fields)) stop "cls_select_maps: index overflow"
    n = size(index_list)
    allocate(Cls_tmp(this%lmin:this%lmax, n*(n+1)/2), spin_tmp(n))
    select type(this)
    type is(coop_binned_Cls)
       if(this%nb .gt. 0)then
          allocate(Cbs_tmp(this%nb,  n*(n+1)/2))
       endif
    end select

    genre_tmp = ''
    do i=1, n
       do j = 1, i
          Cls_tmp(this%lmin:this%lmax, COOP_MATSYM_INDEX(n, i, j)) = this%Cls(this%lmin:this%lmax, COOP_MATSYM_INDEX(this%num_fields, index_list(i), index_list(j)))
          select type(this)
          type is(coop_binned_Cls)
             if(this%nb .gt. 0)then
                Cbs_tmp(1:this%nb, COOP_MATSYM_INDEX(n, i, j)) = this%Cbs(1:this%nb, COOP_MATSYM_INDEX(this%num_fields, index_list(i), index_list(j)))
             endif
          end select
       enddo
       spin_tmp(i)  = this%spin(index_list(i))
       genre_tmp(i:i) = this%genre(index_list(i):index_list(i))
    enddo
    deallocate(this%Cls)
    deallocate(this%spin)
    this%num_fields = n
    this%num_cls = n*(n+1)/2
    allocate(this%spin(n), this%cls(this%lmin:this%lmax, this%num_cls))
    this%spin  = spin_tmp
    this%cls = cls_tmp
    this%genre = genre_tmp
    deallocate(Cls_tmp,spin_tmp)
    select type(this)
    type is(coop_binned_Cls)
       if(this%nb .gt. 0)then
          this%Cbs = Cbs_tmp
          deallocate(Cbs_tmp)
       end if
    end select
  end subroutine coop_cls_select_maps

  subroutine coop_cls_filter(this, lmin, lmax, fwhm_arcmin, highpass_l1, highpass_l2, lowpass_l1, lowpass_l2)
    class(coop_cls)::this
    COOP_REAL, dimension(:,:),allocatable::Cls_tmp
    COOP_REAL, optional::fwhm_arcmin
    COOP_INT, optional::highpass_l1, highpass_l2, lowpass_l1, lowpass_l2, lmin, lmax
    COOP_INT::l, lmax_expect, lmin_expect, lmin_final, lmax_final
    COOP_REAL::beam
    if(.not. allocated(this%Cls)) stop "cls_filter: cls not allocated"
    if(present(highpass_l1) .and. present(highpass_l2))then
       lmin_expect = max(this%lmin, highpass_l1)
    else
       lmin_expect = this%lmin
    endif
    if(present(lowpass_l1) .and. present(lowpass_l2))then
       lmax_expect = min(this%lmax, lowpass_l2)
    else
       lmax_expect = this%lmax
    endif
    if(present(fwhm_arcmin))then
       lmax_expect = min(lmax_expect, ceiling(3.d0/(coop_sigma_by_fwhm*fwhm_arcmin*coop_SI_arcmin)))
    endif
    if(present(fwhm_arcmin) .or. present(highpass_l1).and.present(highpass_l2) .or. present(lowpass_l1) .and. present(lowpass_l2))then
       !$omp parallel do private(beam)
       do l=lmin_expect, lmax_expect
          beam = 1.d0
          if(present(fwhm_arcmin)) &
               beam = beam*coop_Gaussian_filter(fwhm_arcmin=fwhm_arcmin, l = l)**2 
          if(present(highpass_l1) .and. present(highpass_l2)) &
               beam = beam*coop_highpass_filter(l1 = highpass_l1, l2 = highpass_l2, l = l)**2 
          if(present(lowpass_l1) .and. present(lowpass_l2)) &
               beam = beam*coop_lowpass_filter(l1 = lowpass_l1, l2 = lowpass_l2, l = l)**2 
          if(abs(beam-1.d0) .gt. 1.d-10) &
               this%cls(l, :) = this%cls(l, :)*beam
       enddo
       !$omp end parallel do
    endif
    if(present(lmin))then
       lmin_final = lmin
    else
       lmin_final = lmin_expect
    endif
    if(present(lmax))then
       lmax_final = lmax
    else
       lmax_final =  lmax_expect
    endif
    if(lmin_final .ne. this%lmin .or. lmax_final.ne.this%lmax)then
    !!reallocate cls if needed
       allocate(Cls_tmp(lmin_expect:lmax_expect, this%num_cls))
       Cls_tmp = this%Cls(lmin_expect:lmax_expect,:)
       deallocate(this%Cls)
       this%lmin = lmin_final
       this%lmax = lmax_final
       allocate(this%cls(this%lmin:this%lmax, this%num_cls))
       this%cls(this%lmin:lmin_expect-1,:) = 0.d0
       this%cls(lmax_expect+1:this%lmax,:) = 0.d0
       this%cls(max(lmin_expect, this%lmin):min(lmax_expect,this%lmax), :) = cls_tmp(max(lmin_expect, this%lmin):min(lmax_expect,this%lmax), :)
       deallocate(Cls_tmp)
       select type(this)
       type is(coop_binned_Cls)
          if(this%nb .gt.0)&
               call this%bin()
       end select
    elseif(present(fwhm_arcmin) .or. present(highpass_l1).and.present(highpass_l2) .or. present(lowpass_l1) .and. present(lowpass_l2))then
       select type(this)
       type is(coop_binned_Cls)
          if(this%nb .gt.0)&
               call this%bin()
       end select
    endif
  end subroutine coop_cls_filter

  subroutine coop_cls_free(this)
    class(coop_cls)::this
    COOP_DEALLOC(this%cls)
    COOP_DEALLOC(this%spin)
    this%lmin = 0
    this%lmax = 01
    this%num_fields = 0
    this%num_cls  = 0
    this%genre = ''
    this%unit = 'Unknown'
    select type(this)
    type is(coop_binned_Cls)
       COOP_DEALLOC(this%lb)
       COOP_DEALLOC(this%Cbs)
       COOP_DEALLOC(this%wb)
       COOP_DEALLOC(this%wl)
       COOP_DEALLOC(this%lmin_b)
       COOP_DEALLOC(this%lmax_b)
       COOP_DEALLOC(this%ibmin_l)
       COOP_DEALLOC(this%ibmax_l)
       this%nb = 0
    end select
  end subroutine coop_cls_free

  subroutine coop_fits_file_get_naxes(this, naxis, naxes)
    class(coop_fits_file)::this
    COOP_INT::naxis, naxes(:) 
#if HAS_CFITSIO
    select case(this%hdutype)
    case(coop_fitsio_image_hdu)
       call ftgidm(this%unit, naxis, this%status)
       call ftgisz(this%unit, size(naxes), naxes, this%status)
    case default
       call this%report_error("selected HDU is not an image or primary array")
    end select
    call this%check_error()
#else
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_get_naxes

  subroutine coop_fits_file_get_nrows_ncols(this, nrows, ncols)
    class(coop_fits_file)::this
    COOP_INT::nrows, ncols, naxes(2), naxis
#if HAS_CFITSIO
    select case(this%hdutype)
    case(coop_fitsio_image_hdu)
       call ftgidm(this%unit, naxis, this%status)
       if(naxis .ne. 2)then
          write(*,*) trim(this%filename)
          write(*,*) "CHDU = ", this%chdu
          call this%header%print()
          write(*,*) "NAXIS = ", naxis
          write(*,*) "ERROR: coop_fitsio only support nrows_ncols for 2d iamge"
          stop
       endif
       call ftgisz(this%unit, 2, naxes, this%status)
       nrows = naxes(1)
       ncols = naxes(2)
    case(coop_fitsio_binary_table_hdu, coop_fitsio_ascII_table_hdu)
       call ftgnrw(this%unit, nrows, this%status)
       call ftgncl(this%unit, ncols, this%status)
    end select
    call this%check_error()
#else
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_get_nrows_ncols

  subroutine coop_fits_file_read_binary_table(filename, table, first_row, first_col, cols)
    COOP_UNKNOWN_STRING::filename
    COOP_REAL::table(:,:)
    COOP_INT,optional::first_row, first_col
    COOP_INT,dimension(:),optional::cols
#if HAS_CFITSIO    
    type(coop_fits_file)::fp
    COOP_INT::i, nrows, ncols, frow, maxcol, maxrow
    COOP_INT::cols_want(size(table,2))
    call fp%open_table(filename)
    if(present(first_row))then
       frow = first_row
    else
       frow = 1
    endif
    if(present(cols))then
       cols_want = cols(1:size(table,2))
    else
       if(present(first_col))then
          cols_want = (/ (i, i= first_col, first_col + size(table, 2)-1) /)
       else
          cols_want = (/ (i, i=1, size(table, 2) ) /)
       endif
    endif
    maxcol = maxval(cols_want)
    maxrow = frow + size(table, 1) - 1
    call coop_dictionary_lookup(fp%header, "NAXIS2", nrows)
    call coop_dictionary_lookup(fp%header, "NAXIS1", ncols)
    if(nrows .lt. maxrow .or. ncols .lt. maxcol)then
       write(*,*) "file "//trim(filename)
       write(*,*) "number of columns = ", ncols
       write(*,*) "number of rowss = ", nrows
       write(*,*) "rows wanted: ", maxrow
       write(*,*) "columns wanted: ", maxcol
       stop "table size overflow"
    endif
    do i=1, size(table, 2)
       call fp%load_double_column(cols_want(i), table(:,i), first_row = frow)
    enddo
    call fp%close()
#else
    stop "you have to install COOP with CFITSIO"
#endif
  end subroutine coop_fits_file_read_binary_table



  subroutine coop_fits_file_load_double_column(this, col, data, bad_value, first_row)
    class(coop_fits_file)::this
    COOP_INT::col
    COOP_REAL,dimension(:)::data
    COOP_SINGLE,dimension(:),allocatable::single_data
    COOP_INT,dimension(:),allocatable::int_data
    COOP_SHORT_INT,dimension(:),allocatable::short_int_data
    COOP_REAL, optional::bad_value
    COOP_REAL::nulval
    COOP_INT::ncols, width, repeat, datacode, nrows, frow
    COOP_INT,optional::first_row
    logical anyf
#if HAS_CFITSIO
    if(present(first_row))then
       frow = first_row
    else
       frow = 1
    endif
    if(present(bad_value))then
       nulval = bad_value
    else
       nulval = 0.d0
    endif
    call this%get_nrows_ncols(nrows, ncols)
    if(col .le. 0 .or. col.gt. ncols)then
       write(*,*) trim(this%filename)
       write(*,*) "number of table columns:", ncols
       write(*,*) "Error: want column: ", col
       stop
    endif
    if(nrows .lt. frow+size(data)-1)then
       write(*,*) trim(this%filename)
       write(*,*) "number of table rows:", nrows
       write(*,*) "Error: want rows ", frow + size(data) - 1
       stop
    endif
    call ftgtcl(this%unit, col, datacode, repeat, width, this%status)
    if(repeat.gt.1)then
       call this%report_error("repeat >1 is not supported in coop_fitsio")
    endif
    select case(datacode)
    case(coop_fitsio_datatype_double)
       call ftgcvd(this%unit, col, frow, 1, size(data), nulval, data, anyf, this%status)
    case(coop_fitsio_datatype_single)
       allocate(single_data(size(data)))
       call ftgcve(this%unit, col, frow, 1, size(data), real(nulval), single_data, anyf, this%status)
       data = single_data
       deallocate(single_data)
    case(coop_fitsio_datatype_int)
       write(*,*) "Warning: integer column in "//trim(this%filename)//" is loaded as double"
       allocate(int_data(size(data)))
       call ftgcvj(this%unit, col, frow, 1, size(data), nint(nulval), int_data, anyf, this%status)
       data = int_data
       deallocate(int_data)
    case(coop_fitsio_datatype_short_int)
       write(*,*) "Warning: short integer column in "//trim(this%filename)//" is loaded as double"
       allocate(short_int_data(size(data)))
       call ftgcvi(this%unit, col, frow, 1, size(data), nint(nulval, coop_short_int_length), short_int_data, anyf, this%status)
       data = short_int_data
       deallocate(short_int_data)
    case default
       write(*,*) trim(this%filename)
       write(*,*) "data code = ", datacode
       write(*,*) "not supported in coop_fits_file_load_double_column."
       stop
    end select
    call this%check_error()
#else
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_load_double_column


  subroutine coop_fits_file_load_single_column(this, col, data, bad_value, first_row)
    class(coop_fits_file)::this
    COOP_INT::col
    COOP_SINGLE,dimension(:)::data
    COOP_REAL,dimension(:),allocatable::double_data
    COOP_INT,dimension(:),allocatable::int_data
    COOP_SHORT_INT,dimension(:),allocatable::short_int_data
    COOP_SINGLE, optional::bad_value
    COOP_SINGLE::nulval
    COOP_INT::ncols, width, repeat, datacode, nrows, frow
    COOP_INT,optional::first_row
    logical anyf
    character, dimension(:),allocatable::byte_data
#if HAS_CFITSIO
    if(present(first_row))then
       frow = first_row
    else
       frow = 1
    endif
    if(present(bad_value))then
       nulval = bad_value
    else
       nulval = 0.
    endif
    call this%get_nrows_ncols(nrows, ncols)
    if(col .le. 0 .or. col.gt. ncols)then
       write(*,*) trim(this%filename)
       write(*,*) "number of table columns:", ncols
       write(*,*) "Error: want column: ", col
       stop
    endif
    if(nrows .lt. frow+size(data)-1)then
       write(*,*) trim(this%filename)
       write(*,*) "number of table rows:", nrows
       write(*,*) "Error: want rows ", frow + size(data) - 1
       stop
    endif
    call ftgtcl(this%unit, col, datacode, repeat, width, this%status)
    if(repeat.gt.1)then
       call this%report_error("repeat >1 is not supported in coop_fitsio")
    endif
    select case(datacode)
    case(coop_fitsio_datatype_double)
       allocate(double_data(size(data)))
       call ftgcvd(this%unit, col, frow, 1, size(data), dble(nulval), double_data, anyf, this%status)
       data = double_data
       deallocate(double_data)
    case(coop_fitsio_datatype_single)
       call ftgcve(this%unit, col, frow, 1, size(data), nulval, data, anyf, this%status)
    case(coop_fitsio_datatype_int)
       write(*,*) "Warning: integer column in "//trim(this%filename)//" is loaded as single"
       allocate(int_data(size(data)))
       call ftgcvj(this%unit, col, frow, 1, size(data), nint(nulval), int_data, anyf, this%status)
       data = int_data
       deallocate(int_data)
    case(coop_fitsio_datatype_short_int)
       write(*,*) "Warning: short integer column in "//trim(this%filename)//" is loaded as single"
       allocate(short_int_data(size(data)))
       call ftgcvi(this%unit, col, frow, 1, size(data), nint(nulval, coop_short_int_length), short_int_data, anyf, this%status)
       data = short_int_data
       deallocate(short_int_data)

    case(coop_fitsio_datatype_byte)
       write(*,*) "Warning: byte column in "//trim(this%filename)//" is loaded as single"
       allocate(byte_data(size(data)))
       if(nulval .eq. 0.)then
          call ftgcvb(this%unit, col, frow, 1, size(data), '0', byte_data, anyf, this%status)
       else
          call ftgcvb(this%unit, col, frow, 1, size(data), '1', byte_data, anyf, this%status)
       endif
       data = ichar(byte_data)
       deallocate(byte_data)
    case default
       write(*,*) trim(this%filename)
       write(*,*) "data code = ", datacode
       write(*,*) "not supported in coop_fits_file_load_single_column."
       stop
    end select
    call this%check_error()
#else
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_load_single_column


  subroutine coop_fits_file_load_int_column(this, col, data, bad_value, first_row)
    class(coop_fits_file)::this
    COOP_INT::col
    COOP_INT,dimension(:)::data
    COOP_SINGLE,dimension(:),allocatable::single_data
    COOP_REAL,dimension(:),allocatable::double_data
    COOP_SHORT_INT,dimension(:),allocatable::short_int_data
    COOP_INT, optional::bad_value
    COOP_REAL::nulval
    COOP_SHORT_INT::short_int_nulval
    COOP_REAL::double_nulval
    COOP_INT::ncols, width, repeat, datacode, nrows
    COOP_INT,optional::first_row
    COOP_INT::frow
    logical anyf
    character, dimension(:),allocatable::byte_data
#if HAS_CFITSIO
    if(present(first_row))then
       frow = first_row
    else
       frow = 1
    endif
    if(present(bad_value))then
       nulval = bad_value
    else
       nulval = 0.d0
    endif
    call this%get_nrows_ncols(nrows, ncols)
    if(col .le. 0 .or. col.gt. ncols)then
       write(*,*) trim(this%filename)
       write(*,*) "number of table columns:", ncols
       write(*,*) "Error: want column: ", col
       stop
    endif
    if(nrows .lt. size(data)+frow-1)then
       write(*,*) trim(this%filename)
       write(*,*) "number of table rows:", nrows
       write(*,*) "Error: want rows ", size(data)+frow-1
       stop
    endif
    call ftgtcl(this%unit, col, datacode, repeat, width, this%status)
    if(repeat.gt.1)then
       call this%report_error("repeat >1 is not supported in coop_fitsio")
    endif
    select case(datacode)
    case(coop_fitsio_datatype_int)
       call ftgcvj(this%unit, col, frow, 1, size(data), nulval, data, anyf, this%status)
    case(coop_fitsio_datatype_single)
       write(*,*) "Warning: float column in "//trim(this%filename)//" is loaded as int"       
       allocate(single_data(size(data)))
       call ftgcve(this%unit, col, frow, 1, size(data), real(nulval), single_data, anyf, this%status)
       data = nint(single_data)
       deallocate(single_data)
    case(coop_fitsio_datatype_double)
       write(*,*) "Warning: double column in "//trim(this%filename)//" is loaded as int"
       allocate(double_data(size(data)))
       call ftgcvd(this%unit, col, frow, 1, size(data), nint(nulval), double_data, anyf, this%status)
       data = nint(double_data)
       deallocate(double_data)
    case(coop_fitsio_datatype_short_int)
       allocate(short_int_data(size(data)))
       short_int_nulval = nint(nulval)
       call ftgcvi(this%unit, col, frow, 1, size(data), short_int_nulval, short_int_data, anyf, this%status)
       data = short_int_data
       deallocate(short_int_data)
    case(coop_fitsio_datatype_byte)
       write(*,*) "Warning: byte column in "//trim(this%filename)//" is loaded as double"
       allocate(byte_data(size(data)))
       if(nulval .eq. 0.)then
          call ftgcvb(this%unit, col, frow, 1, size(data), '0', byte_data, anyf, this%status)
       else
          call ftgcvb(this%unit, col, frow, 1, size(data), '1', byte_data, anyf, this%status)
       endif
       data = ichar(byte_data)
       deallocate(byte_data)
    case default
       write(*,*) trim(this%filename)
       write(*,*) "data code = ", datacode
       write(*,*) "not supported in coop_fits_file_load_double_column."
       stop
    end select
    call this%check_error()
#else
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_load_int_column

  subroutine coop_fits_file_get_bitpix(this, bitpix)
    class(coop_fits_file)::this
    COOP_INT::bitpix
#if HAS_CFITSIO
    call ftgidt(this%unit, bitpix, this%status)
    call this%check_error()
#else
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_get_bitpix

  subroutine coop_fits_file_load_image_1d(this, image, bad_value)
    class(coop_fits_file)::this
    COOP_REAL,dimension(:)::image
    COOP_SINGLE,dimension(:),allocatable::fimage
    COOP_INT::bitpix
    logical anynul
    COOP_REAL,optional::bad_value
    COOP_INT::naxis, naxes(1024)
#if HAS_CFITSIO
    call this%get_naxes(naxis, naxes)
    if(naxis .gt. 1024)then
       call this%report_error("naxis exceeding 1024 is not supported by coop_fitsio.")
    endif
    if( product(naxes(1:naxis)) .ne. size(image))then
       write(*,*) trim(this%filename)
       write(*,*) "naxis = ", naxis
       write(*,*) "naxes = ", naxes(1:naxis)
       write(*,*) "wanted size = ", size(image)
       write(*,*) "size of wanted image does not agree with fits image size"
       stop
    endif
    call this%get_bitpix(bitpix)
    select case(bitpix)
    case(-64)
       if(present(bad_value))then
          call ftgpvd(this%unit, 1, 1, size(image), bad_value, image, anynul, this%status)
       else
          call ftgpvd(this%unit, 1, 1, size(image), 0.d0, image, anynul, this%status)
       endif
    case(-32)
       allocate(fimage(size(image)))
       if(present(bad_value))then
          call ftgpve(this%unit, 1, 1, size(image), real(bad_value), fimage, anynul, this%status)
       else
          call ftgpve(this%unit, 1, 1, size(image), 0., fimage, anynul, this%status)
       endif
       image =fimage
       deallocate(fimage)
    case default
       write(*,*) trim(this%filename)
       write(*,*) "bitpix = ", bitpix
       write(*,*) "image data are supposed to be double or float numbers"
       stop
    end select
    call this%check_error()
#else
    stop "CFITSIO is not installed"
#endif    
  end subroutine coop_fits_file_load_image_1d


  subroutine coop_fits_file_load_image_2d(this, image, bad_value)
    class(coop_fits_file)::this
    COOP_REAL,dimension(:,:)::image
    COOP_SINGLE,dimension(:,:),allocatable::fimage
    COOP_INT::bitpix
    logical anynul
    COOP_REAL,optional::bad_value
    COOP_INT::naxis, naxes(2)
#if HAS_CFITSIO
    call this%get_naxes(naxis, naxes)
    if(naxis .ne. 2)then
       write(*,*) trim(this%filename)
       write(*,*) "naxis = ", naxis
       write(*,*) "wanted image is 2d"
       stop
    endif
    if(product(naxes(1:naxis)) .ne. size(image))then
       write(*,*) trim(this%filename)
       write(*,*) "naxis = ", naxis
       write(*,*) "naxes = ", naxes(1:naxis)
       write(*,*) "wanted size = ", size(image)
       write(*,*) "size of wanted image does not agree with fits image size"
       stop
    endif
    call this%get_bitpix(bitpix)
    select case(bitpix)
    case(-64)
       if(present(bad_value))then
          call ftgpvd(this%unit, 1, 1, size(image), bad_value, image, anynul, this%status)
       else
          call ftgpvd(this%unit, 1, 1, size(image), 0.d0, image, anynul, this%status)
       endif
    case(-32)
       allocate(fimage(size(image,1), size(image,2)))
       if(present(bad_value))then
          call ftgpve(this%unit, 1, 1, size(image), real(bad_value), fimage, anynul, this%status)
       else
          call ftgpve(this%unit, 1, 1, size(image), 0., fimage, anynul, this%status)
       endif
       image =fimage
       deallocate(fimage)
    case default
       write(*,*) trim(this%filename)
       write(*,*) "bitpix = ", bitpix
       write(*,*) "image data are supposed to be double or float numbers"
       stop
    end select
    call this%check_error()
#else
    stop "CFITSIO is not installed"
#endif    
  end subroutine coop_fits_file_load_image_2d


  subroutine coop_fits_file_load_image_3d(this, image, bad_value)
    class(coop_fits_file)::this
    COOP_REAL,dimension(:,:,:)::image
    COOP_SINGLE,dimension(:,:,:),allocatable::fimage
    COOP_INT::bitpix
    logical anynul
    COOP_INT::naxis, naxes(3)
    COOP_REAL,optional::bad_value
#if HAS_CFITSIO
    call this%get_naxes(naxis, naxes)
    if(naxis .ne. 3)then
       write(*,*) trim(this%filename)
       write(*,*) "naxis = ", naxis
       write(*,*) "wanted image is 2d"
       stop
    endif
    if(product(naxes(1:naxis)) .ne. size(image))then
       write(*,*) trim(this%filename)
       write(*,*) "naxis = ", naxis
       write(*,*) "naxes = ", naxes(1:naxis)
       write(*,*) "wanted size = ", size(image)
       write(*,*) "size of wanted image does not agree with fits image size"
       stop
    endif
    call this%get_bitpix(bitpix)
    select case(bitpix)
    case(-64)
       if(present(bad_value))then
          call ftgpvd(this%unit, 1, 1, size(image), bad_value, image, anynul, this%status)
       else
          call ftgpvd(this%unit, 1, 1, size(image), 0.d0, image, anynul, this%status)
       endif
    case(-32)
       allocate(fimage(size(image,1), size(image,2), size(image,3)))
       if(present(bad_value))then
          call ftgpve(this%unit, 1, 1, size(image), real(bad_value), fimage, anynul, this%status)
       else
          call ftgpve(this%unit, 1, 1, size(image), 0., fimage, anynul, this%status)
       endif
       image =fimage
       deallocate(fimage)
    case default
       write(*,*) trim(this%filename)
       write(*,*) "bitpix = ", bitpix
       write(*,*) "image data are supposed to be double or float numbers"
       stop
    end select
    call this%check_error()
#else
    stop "CFITSIO is not installed"
#endif    
  end subroutine coop_fits_file_load_image_3d



  subroutine coop_fits_file_add_key(this, key, val)
    class(coop_fits_file)::this    
    COOP_UNKNOWN_STRING::key, val
#if HAS_CFITSIO
    if(this%header%index(key).ne.0)then
       call this%header%insert(key, val)
       call ftprec(this%unit, trim(adjustl(key))//"="//trim(adjustl(val)), this%status) 
       call this%check_error()
    else
       call this%header%insert(key, val, overwrite = .true.)
       call ftucrd(this%unit, trim(adjustl(key)), trim(adjustl(key))//"="//trim(adjustl(val)), this%status) 
       call this%check_error()
    endif
#else
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_add_key

  subroutine coop_fits_file_load_header(this)
    class(coop_fits_file)::this    
    COOP_INT::n, i, nadd, pos
    COOP_FITSIO_CARD::card
    COOP_SHORT_STRING::key
    COOP_STRING::val
    call this%header%free()
#if HAS_CFITSIO
    call ftghsp(this%unit, n,  nadd, this%status)
    do i=1, n
       call ftgrec(this%unit, i, card, this%status)
       pos = scan(card, "=")
       if(pos .ne. 0)then
          key = trim(adjustl(card(1:pos-1)))
          val = trim(adjustl(coop_string_strip_quotes(card(pos+1:))))
          if(trim(key).ne.'' .and. trim(val).ne.'') &
               call this%header%insert(key, val )
       endif
    enddo
    call this%check_error()
#else
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_load_header

  subroutine coop_fits_file_report_error(this, errtext)
    class(coop_fits_file)::this    
    COOP_UNKNOWN_STRING::errtext
    write(*,*) "Error in "//trim(this%filename)
    write(*,*) "HDU: "//COOP_STR_OF(this%chdu)
    write(*,*) trim(errtext)
    stop
  end subroutine coop_fits_file_report_error

  subroutine coop_fits_file_check_error(this)
    class(coop_fits_file)::this    
    character(LEN=30)::errtext
#if HAS_CFITSIO
    if(this%status .ne. 0)then
       write(*,*) "Error when reading/writing: "//trim(this%filename)
       call ftgerr(this%status, errtext)
       write(*,*) trim(errtext)
       stop
    endif
#else
    stop "You need to specify cfitio library in configure.in."
#endif
  end subroutine coop_fits_file_check_error

  subroutine coop_fits_file_write_image(image, filename, header)
    type(coop_dictionary),optional::header
    COOP_REAL,dimension(:)::image
    COOP_UNKNOWN_STRING::filename
    type(coop_fits_file)::fp
    COOP_INT::i
    COOP_INT::bitpix, naxis, group, fpixel, nelements, naxes(1024)
    COOP_FITSIO_CARD::card
#if HAS_CFITSIO
    call system('rm -f '//trim(filename))
    nelements = size(image)
    if(present(header))then
       call coop_dictionary_lookup(header, "BITPIX", bitpix, -64)
       call coop_dictionary_lookup(header, "NAXIS", naxis, 0)
       if(naxis .gt. 1024) stop "fits_file_write_image: huge naxis?"
       if(naxis .ge.1)then
          do i = 1, naxis
             call coop_dictionary_lookup(header, "NAXIS"//COOP_STR_OF(i), naxes(i))
          enddo
          if(product(naxes(1:naxis)).ne.nelements)then
             write(*,*) "size of image", nelements
             write(*,*) "naxes = ", naxes(1:naxis)
             write(*,*) "fits_file_write_image: size does not match"
             stop
          endif
       else
          naxis = 1
          naxes(1) = nelements
       endif
    else
       bitpix = -64
       naxis = 1
       naxes(1) = nelements
    endif
    call fp%init(filename)
    call fp%create_image(bitpix, naxis, naxes)
    select case(bitpix)
    case(-64)
       call ftpprd(fp%unit, group, fpixel, nelements, image, fp%status)
    case(-32)
       call ftppre(fp%unit, group, fpixel, nelements, real(image), fp%status)
    case default
       write(*,*) "BITPIX = ", bitpix
       write(*,*) "fits_file_write_image: BITPIX must be -32 or -64 for float image"
       stop
    end select
    if(present(header))then
       do i = 1, header%n
          if(header%key(i)(1:5).eq. "NAXIS")cycle
          select case(trim(header%key(i)))
          case("SIMPLE", "NAXIS", "EXTEND", "NAXIS1", "NAXIS2", "NAXIS3","BITPIX" , "PCOUNT", "GCOUNT")
             !!do nothing
          case default
             card(1:8) = trim(header%key(i))
             card(9:)= "="//trim(header%val(i))
             call ftprec(fp%unit, card, fp%status)
          end select
       enddo
    endif
    call fp%check_error()
    call fp%close()
#else
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_write_image

  subroutine coop_fits_file_get_cls_dimension(filename, lmin, lmax, num_fields, genre)
   type(coop_fits_file)::fp
   COOP_INT::lmin, lmax, num_fields
   COOP_UNKNOWN_STRING::filename, genre
   call fp%open_table(filename)
   call coop_dictionary_lookup(fp%header, "LMIN", lmin)
   call coop_dictionary_lookup(fp%header, "LMAX", lmax)
   call coop_dictionary_lookup(fp%header, "GENRE", genre)
   num_fields = len_trim(genre)
   call fp%close()
 end subroutine coop_fits_file_get_cls_dimension

 subroutine coop_cls_load_txt(this, filename)
   class(coop_cls)::this
   COOP_UNKNOWN_STRING::filename
   type(coop_file)::fp
   COOP_INT::l, il, numlines, numcols
   COOP_STRING::junk, firstline
   call fp%open(filename, 'r')
   read(fp%unit, "(A)", ERR=100, END=100) firstline
   if(firstline(1:2) .eq. "# ")then    !!COOP text Cl format
      this%genre = trim(adjustl(firstline(3:)))
      read(fp%unit, "(A2, I8)", ERR=100, END=100) junk, this%lmin
      read(fp%unit, "(A2, I8)", ERR=100, END=100) junk, this%lmax
      read(fp%unit, "(A2, I8)", ERR=100, END=100) junk, this%num_fields
      this%num_cls = this%num_fields*(this%num_fields+1)/2
      allocate(this%spin(this%num_fields), this%Cls(this%lmin:this%lmax, this%num_cls))
      read(fp%unit, "(A2, "//COOP_STR_OF(this%num_fields)//"I3)") junk, this%spin
      read(fp%unit, "(A2, A18)") junk, this%unit
      this%unit = trim(adjustl(this%unit))
      do l = this%lmin, this%lmax
         read(fp%unit, *, ERR=100, END=100) il, this%cls(l, :)
         if(il.ne.l) goto 100
         if(l.gt.0)then
            this%cls(l, :) = this%cls(l,:)/(l*(l+1.d0)/coop_2pi)
         else
            this%cls(l, :) = this%cls(l,:)*coop_2pi
         endif
      enddo
      call fp%close()
      return
   else
      read(firstline, *, ERR=100) this%lmin
      call fp%close()
      numlines = coop_file_NumLines(filename)
      numcols = coop_file_NumColumns(filename)
      this%lmax = this%lmin + numlines - 1
      this%unit = "muK"
      call fp%open(filename)
      select case(numcols)
      case(2, 4)  !! camb TT file
         this%genre = "I"
         this%num_fields = 1
         this%num_cls = 1
         allocate(this%spin(this%num_fields), this%Cls(this%lmin:this%lmax, this%num_cls))
         this%spin = 0
         do l = this%lmin, this%lmax
            read(fp%unit, *, ERR=100, END=100) il, this%cls(l, 1)
            if(il.ne.l) goto 100
            if(l.gt.0)then
               this%cls(l, :) = this%cls(l,:)/(l*(l+1.d0)/coop_2pi)
            else
               this%cls(l, :) = this%cls(l,:)*coop_2pi
            endif
         enddo
      case(5, 7)  !!camb TT, EE, BB, TE , ...
         this%genre = "IQU"
         this%num_fields = 3
         this%num_cls = 6
         allocate(this%spin(this%num_fields), this%Cls(this%lmin:this%lmax, this%num_cls))         
         this%spin = (/ 0, 2, 2 /)
         this%cls(:, coop_TEB_index_TB) = 0.d0
         this%cls(:, coop_TEB_index_EB) = 0.d0
         do l = this%lmin, this%lmax
            read(fp%unit, *, ERR=100) il, this%cls(l, coop_TEB_index_TT),  this%cls(l, coop_TEB_index_EE), this%cls(l, coop_TEB_index_BB), this%cls(l, coop_TEB_index_TE)
            if(il.ne.l) goto 100
            if(l.gt.0)then
               this%cls(l, :) = this%cls(l,:)/(l*(l+1.d0)/coop_2pi)
            else
               this%cls(l, :) = this%cls(l,:)*coop_2pi
            endif
         enddo
      case default
         goto 100
      end select
      call fp%close()
      return
   endif
100  write(*,*) "Error in the cl input file "//trim(filename)
     stop
 end subroutine coop_cls_load_txt

 subroutine coop_cls_dump_txt(this, filename)
   COOP_UNKNOWN_STRING::filename
   class(coop_cls)::this
   type(coop_file)::fp
   COOP_INT::l
   call fp%open(filename)
   write(fp%unit, "(A2, A18)") "# ", trim(this%genre)
   write(fp%unit, "(A2, I8)") "# ", this%lmin
   write(fp%unit, "(A2, I8)") "# ", this%lmax
   write(fp%unit, "(A2, I8)") "# ", this%num_fields
   write(fp%unit, "(A2, "//COOP_STR_OF(this%num_fields)//"I3)") "# ", this%spin
   if(trim(this%unit).ne.'')then
      write(fp%unit, "(A2, A18)") "# ", trim(this%unit)
   else
      write(fp%unit, "(A2, A18)") "# ", "Unkown"
   endif
   do l = this%lmin, this%lmax
      if(l.gt.0)then
         write(fp%unit, "(I8, "//COOP_STR_OF(this%num_cls)//"E16.7)") l, this%Cls(l, :)*l*(l+1.d0)/coop_2pi
      else
         write(fp%unit, "(I8, "//COOP_STR_OF(this%num_cls)//"E16.7)") l, this%Cls(l, :)/coop_2pi
      endif
   enddo
   call fp%close()
 end subroutine coop_cls_dump_txt

 subroutine coop_cls_load_fits(this, filename)
   class(coop_cls)::this
   COOP_UNKNOWN_STRING::filename
   COOP_STRING::genre
   COOP_INT::i
   call this%free()
   call coop_fits_file_get_cls_dimension(filename, this%lmin, this%lmax, this%num_fields, this%genre)   
   this%num_cls = this%num_fields*(this%num_fields+1)/2
   allocate( this%spin(this%num_fields), this%Cls(this%lmin:this%lmax, this%num_cls))
   call coop_fits_file_load_cls(this%lmin, this%lmax, this%cls, filename, this%genre, this%spin, this%unit)
 end subroutine coop_cls_load_fits

 subroutine  coop_cls_dump_fits(this, filename)
   class(coop_cls)::this
   COOP_UNKNOWN_STRING::filename
   call coop_fits_file_write_cls(this%lmin, this%lmax, this%cls, filename, this%genre, this%spin, this%unit)
 end subroutine coop_cls_dump_fits

 subroutine coop_cls_load(this, filename)
   class(coop_cls)::this
   COOP_UNKNOWN_STRING::filename
   select case(trim(coop_str_numLowerAlpha(coop_file_postfix_of(filename))))
   case("fits")
      call this%load_fits(filename)
   case default
      call this%load_txt(filename)
   end select
 end subroutine coop_cls_load

 subroutine coop_cls_dump(this, filename)
   class(coop_cls)::this
   COOP_UNKNOWN_STRING::filename
   select case(trim(coop_str_numLowerAlpha(coop_file_postfix_of(filename))))
   case("fits")
      call this%dump_fits(filename)
   case default
      call this%dump_txt(filename)
   end select
 end subroutine coop_cls_dump


 subroutine coop_fits_file_load_cls(lmin, lmax, cls, filename, genre, spin, unit)
   type(coop_fits_file)::fp
   COOP_INT::lmin, lmax, llmin, llmax, numcls, i, j, numfields, col, ind
   COOP_REAL::Cls(lmin:lmax, *)
   COOP_INT,dimension(:),optional::spin
   COOP_UNKNOWN_STRING::filename, genre
   COOP_UNKNOWN_STRING,optional::unit
   COOP_STRING::genre_saved
#if HAS_CFITSIO
   if(size(cls,1).ne.lmax-lmin+1) stop "fits_file_load_cls: Size of Cls does not agree with lmin, lmax inputs"
   call fp%open_table(filename)
   call coop_dictionary_lookup(fp%header, "LMIN", llmin)
   call coop_dictionary_lookup(fp%header, "LMAX", llmax)
   call coop_dictionary_lookup(fp%header, "GENRE", genre_saved)
   if(present(unit))then
      call coop_dictionary_lookup(fp%header, "UNIT", unit)
   end if
   if(lmin .lt. llmin .or. lmax .gt. llmax)then
      write(*,*) "L range in file "//trim(filename)
      write(*,*) llmin, llmax
      write(*,*) "Cannot load Cls for lmin = ", lmin, "; lmax = ", lmax
      stop
   endif
   numfields = len_trim(genre)
   numcls = (numfields+1)*numfields/2
   if(present(spin))then
      do i = 1, numfields
         ind = scan(trim(genre_saved), genre(i:i))
         if(ind .ne. 0)then
            call coop_dictionary_lookup(fp%header, "SPIN"//COOP_STR_OF(ind), spin(i))
         else
            spin(i) = 0
         endif
      enddo
   endif
   i = verify(trim(genre), trim(genre_saved))
   if(i.ne. 0)then
      write(*,*) "saved cls in "//trim(filename)//" are for fields: "//trim(genre_saved)
      write(*,*) "cannot load cls for fields: "//trim(genre)
      stop
   endif
   do i=1, numfields
      do j=i, numfields
         call coop_dictionary_lookup(fp%header, "COL"//genre(i:i)//genre(j:j), col, 0)
         if(col .eq. 0)then
            write(*,*) "CL_"//genre(i:i)//genre(j:j)//" cannot be found in "//trim(filename)
            stop
         endif
         call fp%load_double_column(col = col, data = Cls(lmin:lmax, COOP_MATSYM_INDEX(numfields, i, j)), first_row = lmin - llmin + 1) 
      enddo
   enddo
   call fp%close()
#endif
  end subroutine coop_fits_file_load_cls

  subroutine coop_fits_file_write_cls(lmin, lmax, cls, filename, genre, spin, unit)
    COOP_INT::lmin, lmax, numcls, numfields, i, j
    COOP_REAL::Cls(:,:)
    COOP_INT,dimension(:),optional::spin
    type(coop_dictionary)::header
    COOP_UNKNOWN_STRING::filename
    COOP_UNKNOWN_STRING, optional::unit
    COOP_UNKNOWN_STRING::genre
#if HAS_CFITSIO
    if(size(cls, 1) .ne. lmax-lmin+1) stop "fits_file_write_cls: Size of Cls does not agree with lmin, lmax inputs"
    call header%init()
    call header%insert("LMIN", COOP_STR_OF(lmin))
    call header%insert("LMAX", COOP_STR_OF(lmax))
    call header%insert("GENRE", trim(adjustl(genre)))
    if(present(unit))then
       call header%insert("UNIT", trim(adjustl(unit)))
    endif
    numcls = size(Cls,2)
    numfields = nint(sqrt(numcls*2+0.25d0)-0.5d0)
    if(present(spin))then
       do i=1, numfields
          call header%insert("SPIN"//COOP_STR_OF(i), COOP_STR_OF(spin(i)))
       enddo
    else
       do i=1, numfields
          call header%insert("SPIN"//COOP_STR_OF(i), "0")
       enddo
    endif
    if(present(unit))then
       call header%insert("UNIT", trim(adjustl(unit)))
    else
       call header%insert("UNIT", "Unknown")
    end if
    if(numfields*(numfields+1)/2 .ne. numcls) stop "write_cls: the number of cls must be in the form of n(n+1)/2"
    if(numfields .gt. len_trim(genre)) stop "write_cls: genre must contain all the fields"
    do i=1, numfields
       do j = i, numfields
          call header%insert("TTYPE"//COOP_STR_OF(COOP_MATSYM_INDEX(numfields, i, j)), genre(i:i)//genre(j:j))
          call header%insert("COL"//genre(i:i)//genre(j:j), COOP_STR_OF(COOP_MATSYM_INDEX(numfields, i, j)))
          if(j.ne.i)call header%insert("COL"//genre(j:j)//genre(i:i), COOP_STR_OF(COOP_MATSYM_INDEX(numfields, i, j)))
       enddo
    enddo
    call coop_fits_file_write_binary_table(cls, filename, header)
    call header%free()
#endif
  end subroutine coop_fits_file_write_cls


  subroutine coop_fits_file_create_hdu(this)
    class(coop_fits_file)::this
#if HAS_CFITSIO
    call ftcrhd(this%unit, this%status)
    call ftthdu(this%unit, this%nhdus, this%status)
    call ftghdn(this%unit, this%chdu)
#endif
  end subroutine coop_fits_file_create_hdu

  subroutine coop_fits_file_update_hdu_info(this)
    class(coop_fits_file)::this
#if HAS_CFITSIO
    call ftthdu(this%unit, this%nhdus, this%status)
    call ftghdn(this%unit, this%chdu)
    call ftghdt(this%unit, this%hdutype, this%status)
#endif
  end subroutine coop_fits_file_update_hdu_info

  subroutine coop_fits_file_write_binary_table_d(table, filename, header)
    type(coop_dictionary),optional::header
    COOP_REAL,dimension(:,:)::table
    COOP_UNKNOWN_STRING::filename
    type(coop_fits_file)::fp
    COOP_INT::i
    COOP_INT::status,blocksize,hdutype,tfields,nrows,varidat, column,frow,felem
    COOP_SHORT_STRING::extname
    COOP_SHORT_STRING,dimension(:),allocatable::ttype, tform, tunit
    COOP_FITSIO_CARD::card
#if HAS_CFITSIO
    call system('rm -f '//trim(filename))
    fp%rwmode = 1
    call fp%init(trim(filename))
    if(present(header))then
       call coop_dictionary_lookup(header, "TFIELDS", tfields, size(table,2))
       call coop_dictionary_lookup(header, "NAXIS2", nrows, size(table,1))
       if(tfields .ne. size(table,2))then
          write(*,*) "TFIELDS = ", tfields
          write(*,*) "table has columns: ", size(table,2)
          write(*,*) "fits_file_write_binary_table: matching failed"
          stop
       endif
       if(nrows .ne. size(table,1))then
          write(*,*) "NAXIS2 = ", nrows
          write(*,*) "table has rows: ", size(table,1)
          write(*,*) "fits_file_write_binary_table: matching failed"
          stop
       endif
       call coop_dictionary_lookup(header, "EXTNAME", extname, "COOP_BINARY")
       
    else
       tfields = size(table, 2)
       nrows = size(table, 1)
       extname = "COOP_BINARY"
    endif
    allocate(ttype(tfields), tform(tfields), tunit(tfields))
    if(present(header))then
       do i=1, tfields
          call coop_dictionary_lookup(header, "TTYPE"//COOP_STR_OF(i), ttype(i) , "COLUMN_"//COOP_STR_OF(i))
          call coop_dictionary_lookup(header, "TFORM"//COOP_STR_OF(i), tform(i) , "1D")
          call coop_dictionary_lookup(header, "TUNIT"//COOP_STR_OF(i), tunit(i), " ")
       enddo
    else
       do i=1, tfields
          ttype(i) = "COLUMN_"//COOP_STR_OF(i)
          tform(i) = "1D"
          tunit(i) = " "
       enddo
    endif
    varidat = 0
    call fp%create_binary_table(nrows, tfields, ttype, tform, tunit)
    frow = 1
    felem = 1
    do column = 1, tfields
       select case(trim(tform(column)))
       case("1D", "D")
          call ftpcld(fp%unit, column, frow, felem, nrows, table(:, column), fp%status)
       case("1E", "E")
          call ftpcle(fp%unit, column, frow, felem, nrows, real(table(:, column)), fp%status)
       case("1J", "J")
          call ftpclj(fp%unit, column, frow, felem, nrows, nint(table(:, column)), fp%status)
       case("1I", "I")
          call ftpcli(fp%unit, column, frow, felem, nrows, nint(table(:, column), coop_short_int_length), fp%status)
       case default
          write(*,*) "TFORM"//COOP_STR_OF(column)//"="//trim(tform(column))
          stop "coop_fits_file_write_binary: does not support this TFORM"
       end select
    enddo
    deallocate(ttype, tform, tunit)
    if(present(header))then
       do i = 1, header%n
          if(header%key(i)(1:5).eq. "NAXIS" .or. header%key(i)(1:5).eq. "TTYPE" .or. header%key(i)(1:5).eq. "TFORM" .or. header%key(i)(1:5).eq. "TUNIT")cycle
          select case(trim(header%key(i)))
          case("SIMPLE", "EXTEND", "BITPIX" , "PCOUNT", "GCOUNT", "TFIELDS", "EXTNAME")
             !!do nothing
          case default
             card(1:8) = trim(header%key(i))
             card(9:)= "="//trim(header%val(i))
             call ftprec(fp%unit, card, fp%status)
          end select
       enddo
    endif
    call fp%check_error()
    call fp%close()
#else
    stop "FITSIO is not installed"
#endif
  end subroutine coop_fits_file_write_binary_table_d


  subroutine coop_fits_file_write_binary_table_with_indices_d(indices, table, filename, header)
    type(coop_dictionary),optional::header
    COOP_INT, dimension(:)::indices
    COOP_REAL,dimension(:,:)::table
    COOP_UNKNOWN_STRING::filename
    type(coop_fits_file)::fp
    COOP_INT::i
    COOP_INT::status,blocksize,hdutype,tfields,nrows,varidat, column,frow,felem
    COOP_SHORT_STRING::extname
    COOP_SHORT_STRING,dimension(:),allocatable::ttype, tform, tunit
    COOP_FITSIO_CARD::card
#if HAS_CFITSIO
    if(size(indices) .ne. size(table, 1)) then
       write(*,*) size(indices), size(table,1)
       stop "the size of  pixel indices does not equal to the size of of table"
    endif
    call system('rm -f '//trim(filename))
    fp%rwmode = 1
    call fp%init(trim(filename))
    if(present(header))then
       call coop_dictionary_lookup(header, "TFIELDS", tfields, size(table,2))
       call coop_dictionary_lookup(header, "NAXIS2", nrows, size(table,1))
       if(tfields .ne. size(table,2))then
          write(*,*) "TFIELDS = ", tfields
          write(*,*) "table has columns: ", size(table,2)
          write(*,*) "fits_file_write_binary_table: matching failed"
          stop
       endif
       if(nrows .ne. size(table,1))then
          write(*,*) "NAXIS2 = ", nrows
          write(*,*) "table has rows: ", size(table,1)
          write(*,*) "fits_file_write_binary_table: matching failed"
          stop
       endif
       call coop_dictionary_lookup(header, "EXTNAME", extname, "COOP_BINARY")
       
    else
       tfields = size(table, 2)
       nrows = size(table, 1)
       extname = "COOP_BINARY"
    endif
    allocate(ttype(tfields), tform(tfields), tunit(tfields))
    if(present(header))then
       do i=1, tfields
          call coop_dictionary_lookup(header, "TTYPE"//COOP_STR_OF(i), ttype(i) , "COLUMN_"//COOP_STR_OF(i))
          call coop_dictionary_lookup(header, "TFORM"//COOP_STR_OF(i), tform(i) , "1D")
          call coop_dictionary_lookup(header, "TUNIT"//COOP_STR_OF(i), tunit(i), " ")
       enddo
    else
       do i=1, tfields
          ttype(i) = "COLUMN_"//COOP_STR_OF(i)
          tform(i) = "1D"
          tunit(i) = " "
       enddo
    endif
    varidat = 0
    call fp%create_binary_table(nrows, tfields, ttype, tform, tunit)
    frow = 1
    felem = 1
    call ftpclj(fp%unit, column, frow, felem, nrows, indices, fp%status)
    do column = 1, tfields
       select case(trim(tform(column)))
       case("1D", "D")
          call ftpcld(fp%unit, column, frow, felem, nrows, table(:, column), fp%status)
       case("1E", "E")
          call ftpcle(fp%unit, column, frow, felem, nrows, real(table(:, column)), fp%status)
       case("1J", "J")
          call ftpclj(fp%unit, column, frow, felem, nrows, nint(table(:, column)), fp%status)
       case("1I", "I")
          call ftpcli(fp%unit, column, frow, felem, nrows, nint(table(:, column), coop_short_int_length), fp%status)
       case default
          write(*,*) "TFORM"//COOP_STR_OF(column)//"="//trim(tform(column))
          stop "coop_fits_file_write_binary: does not support this TFORM"
       end select
    enddo
    deallocate(ttype, tform, tunit)
    if(present(header))then
       do i = 1, header%n
          if(header%key(i)(1:5).eq. "NAXIS" .or. header%key(i)(1:5).eq. "TTYPE" .or. header%key(i)(1:5).eq. "TFORM" .or. header%key(i)(1:5).eq. "TUNIT")cycle
          select case(trim(header%key(i)))
          case("SIMPLE", "EXTEND", "BITPIX" , "PCOUNT", "GCOUNT", "TFIELDS", "EXTNAME")
             !!do nothing
          case default
             card(1:8) = trim(header%key(i))
             card(9:)= "="//trim(header%val(i))
             call ftprec(fp%unit, card, fp%status)
          end select
       enddo
    endif
    call fp%check_error()
    call fp%close()
#else
    stop "FITSIO is not installed"
#endif
  end subroutine coop_fits_file_write_binary_table_with_indices_d




  subroutine coop_fits_file_write_binary_table_s(table, filename, header)
    type(coop_dictionary),optional::header
    COOP_SINGLE,dimension(:,:)::table
    COOP_UNKNOWN_STRING::filename
    type(coop_fits_file)::fp
    COOP_INT::i
    COOP_INT::status,blocksize,hdutype,tfields,nrows,varidat, column,frow,felem
    COOP_SHORT_STRING::extname
    COOP_SHORT_STRING,dimension(:),allocatable::ttype, tform, tunit
    COOP_FITSIO_CARD::card
#if HAS_CFITSIO
    call system('rm -f '//trim(filename))
    fp%rwmode = 1
    call fp%init(trim(filename))
    if(present(header))then
       call coop_dictionary_lookup(header, "TFIELDS", tfields, size(table,2))
       call coop_dictionary_lookup(header, "NAXIS2", nrows, size(table,1))
       if(tfields .ne. size(table,2))then
          write(*,*) "TFIELDS = ", tfields
          write(*,*) "table has columns: ", size(table,2)
          write(*,*) "fits_file_write_binary_table: matching failed"
          stop
       endif
       if(nrows .ne. size(table,1))then
          write(*,*) "NAXIS2 = ", nrows
          write(*,*) "table has rows: ", size(table,1)
          write(*,*) "fits_file_write_binary_table: matching failed"
          stop
       endif
       call coop_dictionary_lookup(header, "EXTNAME", extname, "COOP_BINARY")
       
    else
       tfields = size(table, 2)
       nrows = size(table, 1)
       extname = "COOP_BINARY"
    endif
    allocate(ttype(tfields), tform(tfields), tunit(tfields))
    if(present(header))then
       do i=1, tfields
          call coop_dictionary_lookup(header, "TTYPE"//COOP_STR_OF(i), ttype(i) , "COLUMN_"//COOP_STR_OF(i))
          call coop_dictionary_lookup(header, "TFORM"//COOP_STR_OF(i), tform(i) , "1D")
          call coop_dictionary_lookup(header, "TUNIT"//COOP_STR_OF(i), tunit(i), " ")
       enddo
    else
       do i=1, tfields
          ttype(i) = "COLUMN_"//COOP_STR_OF(i)
          tform(i) = "1E"
          tunit(i) = " "
       enddo
    endif
    varidat = 0
    call fp%create_binary_table(nrows, tfields, ttype, tform, tunit)
    frow = 1
    felem = 1
    do column = 1, tfields
       select case(trim(tform(column)))
       case("1D", "D")
          call ftpcld(fp%unit, column, frow, felem, nrows, dble(table(:, column)), fp%status)
       case("1E", "E")
          call ftpcle(fp%unit, column, frow, felem, nrows, table(:, column), fp%status)
       case("1J", "J")
          call ftpclj(fp%unit, column, frow, felem, nrows, nint(table(:, column)), fp%status)
       case("1I", "I")
          call ftpcli(fp%unit, column, frow, felem, nrows, nint(table(:, column), coop_short_int_length), fp%status)
       case default
          write(*,*) "TFORM"//COOP_STR_OF(column)//"="//trim(tform(column))
          stop "coop_fits_file_write_binary: does not support this TFORM"
       end select
    enddo
    deallocate(ttype, tform, tunit)
    if(present(header))then
       do i = 1, header%n
          if(header%key(i)(1:5).eq. "NAXIS" .or. header%key(i)(1:5).eq. "TTYPE" .or. header%key(i)(1:5).eq. "TFORM" .or. header%key(i)(1:5).eq. "TUNIT")cycle
          select case(trim(header%key(i)))
          case("SIMPLE", "EXTEND", "BITPIX" , "PCOUNT", "GCOUNT", "TFIELDS", "EXTNAME")
             !!do nothing
          case default
             card(1:8) = trim(header%key(i))
             card(9:)= "="//trim(header%val(i))
             call ftprec(fp%unit, card, fp%status)
          end select
       enddo
    endif
    call fp%check_error()
    call fp%close()
#else
    stop "FITSIO is not installed"
#endif
  end subroutine coop_fits_file_write_binary_table_s


  subroutine coop_fits_file_write_binary_table_with_indices_s(indices, table, filename, header)
    type(coop_dictionary),optional::header
    COOP_INT, dimension(:)::indices
    COOP_SINGLE,dimension(:,:)::table
    COOP_UNKNOWN_STRING::filename
    type(coop_fits_file)::fp
    COOP_INT::i
    COOP_INT::status,blocksize,hdutype,tfields,nrows,varidat, column,frow,felem
    COOP_SHORT_STRING::extname
    COOP_SHORT_STRING,dimension(:),allocatable::ttype, tform, tunit
    COOP_FITSIO_CARD::card
#if HAS_CFITSIO
    if(size(indices) .ne. size(table, 1)) then
       write(*,*) size(indices), size(table,1)
       stop "the size of  pixel indices does not equal to the size of of table"
    endif
    call system('rm -f '//trim(filename))
    fp%rwmode = 1
    call fp%init(trim(filename))
    if(present(header))then
       call coop_dictionary_lookup(header, "TFIELDS", tfields, size(table,2))
       call coop_dictionary_lookup(header, "NAXIS2", nrows, size(table,1))
       if(tfields .ne. size(table,2))then
          write(*,*) "TFIELDS = ", tfields
          write(*,*) "table has columns: ", size(table,2)
          write(*,*) "fits_file_write_binary_table: matching failed"
          stop
       endif
       if(nrows .ne. size(table,1))then
          write(*,*) "NAXIS2 = ", nrows
          write(*,*) "table has rows: ", size(table,1)
          write(*,*) "fits_file_write_binary_table: matching failed"
          stop
       endif
       call coop_dictionary_lookup(header, "EXTNAME", extname, "COOP_BINARY")
       
    else
       tfields = size(table, 2)
       nrows = size(table, 1)
       extname = "COOP_BINARY"
    endif
    allocate(ttype(tfields), tform(tfields), tunit(tfields))
    if(present(header))then
       do i=1, tfields
          call coop_dictionary_lookup(header, "TTYPE"//COOP_STR_OF(i), ttype(i) , "COLUMN_"//COOP_STR_OF(i))
          call coop_dictionary_lookup(header, "TFORM"//COOP_STR_OF(i), tform(i) , "1D")
          call coop_dictionary_lookup(header, "TUNIT"//COOP_STR_OF(i), tunit(i), " ")
       enddo
    else
       do i=1, tfields
          ttype(i) = "COLUMN_"//COOP_STR_OF(i)
          tform(i) = "1E"
          tunit(i) = " "
       enddo
    endif
    varidat = 0
    call fp%create_binary_table(nrows, tfields, ttype, tform, tunit)
    frow = 1
    felem = 1
    call ftpclj(fp%unit, column, frow, felem, nrows, indices, fp%status)
    do column = 1, tfields
       select case(trim(tform(column)))
       case("1D", "D")
          call ftpcld(fp%unit, column, frow, felem, nrows, dble(table(:, column)), fp%status)
       case("1E", "E")
          call ftpcle(fp%unit, column, frow, felem, nrows, table(:, column), fp%status)
       case("1J", "J")
          call ftpclj(fp%unit, column, frow, felem, nrows, nint(table(:, column)), fp%status)
       case("1I", "I")
          call ftpcli(fp%unit, column, frow, felem, nrows, nint(table(:, column), coop_short_int_length), fp%status)
       case default
          write(*,*) "TFORM"//COOP_STR_OF(column)//"="//trim(tform(column))
          stop "coop_fits_file_write_binary: does not support this TFORM"
       end select
    enddo
    deallocate(ttype, tform, tunit)
    if(present(header))then
       do i = 1, header%n
          if(header%key(i)(1:5).eq. "NAXIS" .or. header%key(i)(1:5).eq. "TTYPE" .or. header%key(i)(1:5).eq. "TFORM" .or. header%key(i)(1:5).eq. "TUNIT")cycle
          select case(trim(header%key(i)))
          case("SIMPLE", "EXTEND", "BITPIX" , "PCOUNT", "GCOUNT", "TFIELDS", "EXTNAME")
             !!do nothing
          case default
             card(1:8) = trim(header%key(i))
             card(9:)= "="//trim(header%val(i))
             call ftprec(fp%unit, card, fp%status)
          end select
       enddo
    endif
    call fp%check_error()
    call fp%close()
#else
    stop "FITSIO is not installed"
#endif
  end subroutine coop_fits_file_write_binary_table_with_indices_s


  subroutine coop_fits_file_create_binary_table(this, nrows, tfields, ttype, tform, tunit)
    class(coop_fits_file)::this
    COOP_INT::nrows, tfields, varidat
    COOP_SHORT_STRING::ttype(tfields), tform(tfields), tunit(tfields), extname
#if HAS_CFITSIO
    varidat = 0
    extname = "COOP_BINARY_TABLE"
    call ftibin(this%unit, nrows, tfields, ttype, tform, tunit, extname, varidat, this%status)
    call this%update_hdu_info()
    call this%load_header()
    call this%check_error()
#endif
  end subroutine coop_fits_file_create_binary_table


  subroutine coop_fits_file_create_ascii_table(this, rowlen, nrows, tfields, ttype, tbcol, tform, tunit)
    class(coop_fits_file)::this
    COOP_INT::rowlen, nrows, tfields
    COOP_INT::tbcol(tfields)
    COOP_SHORT_STRING::ttype(tfields), tform(tfields), tunit(tfields), extname
    extname = "COOP_ASCII_TABLE"
#if HAS_CFITSIO
    call ftitab(this%unit, rowlen, nrows, tfields, ttype, tbcol, tform, tunit, extname, this%status)
    call this%update_hdu_info()
    call this%load_header()
    call this%check_error()
#endif
  end subroutine coop_fits_file_create_ascii_table



  subroutine coop_fits_file_create_image(this, bitpix, naxis, naxes)
    class(coop_fits_file)::this
    COOP_INT::bitpix, naxis
    COOP_INT::naxes(naxis)
#if HAS_CFITSIO
    call ftiimg(this%unit, bitpix, naxis, naxes, this%status)
    call this%update_hdu_info()
    call this%load_header()
    call this%check_error()
#endif
  end subroutine coop_fits_file_create_image

  subroutine coop_fits_file_write_image_2d(image, filename, header)
    type(coop_dictionary),optional::header
    COOP_REAL,dimension(:,:)::image
    COOP_UNKNOWN_STRING::filename
    type(coop_fits_file)::fp
    COOP_INT::i
    COOP_INT::bitpix, naxis, naxes(2), group, fpixel, nelements
    COOP_FITSIO_CARD::card

#if HAS_CFITSIO
    call system('rm -f '//trim(filename))
    if(present(header))then
       call coop_dictionary_lookup(header, "BITPIX", bitpix, -64)
    else
       bitpix = -64
    endif
    call fp%init(filename)
    naxis = 2
    group = 1
    fpixel = 1
    nelements = size(image)
    naxes(1) = size(image,1)
    naxes(2) = size(image,2)
    call fp%create_image(bitpix, naxis, naxes)
    select case(bitpix)
    case(-64)
       call ftpprd(fp%unit, group, fpixel, nelements, image, fp%status)
    case(-32)
       call ftppre(fp%unit, group, fpixel, nelements, real(image), fp%status)
    case default
       write(*,*) "BITPIX = ", bitpix
       write(*,*) "fits_file_write_image: BITPIX must be -32 or -64 for float image"
       stop
    end select
    if(present(header))then
       do i = 1, header%n
          if(header%key(i)(1:5).eq. "NAXIS")cycle
          select case(trim(header%key(i)))
          case("SIMPLE", "NAXIS", "EXTEND", "NAXIS1", "NAXIS2", "NAXIS3","BITPIX" , "PCOUNT", "GCOUNT")
             !!do nothing
          case default
             card(1:8) = trim(header%key(i))
             card(9:)= "="//trim(header%val(i))
             call ftprec(fp%unit, card, fp%status)
          end select
       enddo
    endif
    call fp%check_error()
    call fp%close()
#else
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_write_image_2d

  subroutine coop_fits_file_write_image_3d(image, filename, header)
    type(coop_dictionary),optional::header
    COOP_REAL,dimension(:,:,:)::image
    COOP_UNKNOWN_STRING::filename
    type(coop_fits_file)::fp
    COOP_INT::i
    COOP_FITSIO_CARD::card
    COOP_INT::bitpix, naxis, naxes(3), group, fpixel, nelements
#if HAS_CFITSIO
    call system('rm -f '//trim(filename))
    if(present(header))then
       call coop_dictionary_lookup(header, "BITPIX", bitpix, -64)
    else
       bitpix = -64
    endif
    naxis = 3
    group = 1
    fpixel = 1
    nelements = size(image)
    naxes(1) = size(image,1)
    naxes(2) = size(image,2)
    naxes(3) = size(image,3)
    call fp%init(filename)
    call fp%create_image(bitpix, naxis, naxes)
    select case(bitpix)
    case(-64)
       call ftpprd(fp%unit, group, fpixel, nelements, image, fp%status)
    case(-32)
       call ftppre(fp%unit, group, fpixel, nelements, real(image), fp%status)
    case default
       write(*,*) "BITPIX = ", bitpix
       write(*,*) "fits_file_write_image: BITPIX must be -32 or -64 for float image"
       stop
    end select
    if(present(header))then
       do i = 1, header%n
          if(header%key(i)(1:5).eq. "NAXIS")cycle
          select case(trim(header%key(i)))
          case("SIMPLE", "EXTEND", "BITPIX" , "PCOUNT", "GCOUNT")
             !!do nothing
          case default
             card(1:8) = trim(header%key(i))
             card(9:)= "="//trim(header%val(i))
             call ftprec(fp%unit, card, fp%status)
          end select
       enddo
    endif
    call fp%check_error()
    call fp%close()
#else
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_write_image_3d


  subroutine coop_fits_file_init(this, filename)
    class(coop_fits_file)::this    
    COOP_UNKNOWN_STRING::filename
#if HAS_CFITSIO
    if(this%unit .ne. 0) call this%close()
    this%filename = trim(adjustl(filename))
    this%unit = coop_free_file_unit()
    this%blocksize = 1
    this%status = 0
    call ftinit(this%unit, trim(this%filename), this%blocksize, this%status)
    call this%check_error()
#else
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_init

  subroutine coop_fits_file_close(this)
    class(coop_fits_file)::this
#if HAS_CFITSIO
    call ftclos(this%unit, this%status)
    call this%check_error()
    this%unit = 0
    this%filename = ""
    this%nhdus = 0
    this%chdu = 0
    this%status = 0
    this%rwmode = 0
    this%blocksize= 0
    this%hdutype = 0
    call this%header%free()
#else
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_close

  subroutine coop_fits_file_open(this, filename, ihdu, mode)
    class(coop_fits_file)::this
    COOP_UNKNOWN_STRING, optional::mode
    COOP_UNKNOWN_STRING::filename
    COOP_INT, optional::ihdu
#if HAS_CFITSIO
    if(this%unit .ne. 0) call this%close()
    this%unit = coop_free_file_unit()

    this%filename = trim(adjustl(filename))
    if(present(mode))then
       select case(COOP_UPPER_STR(mode))
       case("RW", "W", "WRITE")
          this%rwmode = 1
       case("R","READ", "READONLY")
          this%rwmode = 0
       case default
          write(*,*) trim(mode)
          write(*,*) "Unknown mode for fits_file_open"
          stop
       end select
    else
       this%rwmode = 0 !!default readonly
    endif
    if(present(ihdu))then
       call ftnopn(this%unit, trim(this%filename)//"["//COOP_STR_OF(ihdu-1)//"]", this%rwmode, this%status)
    else
       call ftopen(this%unit, trim(this%filename), this%rwmode, this%blocksize, this%status)
    endif
    call this%update_hdu_info()
    call this%load_header()
    call this%check_error()
#else
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_open


  subroutine coop_fits_file_open_image(this, filename)
    class(coop_fits_file)::this
    COOP_UNKNOWN_STRING::filename
#if HAS_CFITSIO
    if(this%unit .ne. 0) call this%close()
    this%unit = coop_free_file_unit()

    this%filename = trim(adjustl(filename))
    call ftiopn(this%unit, trim(this%filename), this%rwmode, this%status)
    if(this%status .ne. 0)then
       this%status = 0
       if(this%unit .ne. 0) call this%close()
       this%unit = coop_free_file_unit()
       call ftdopn(this%unit, trim(this%filename), this%rwmode, this%status)
    endif
    if(this%status .ne. 0)then
       this%status = 0
       if(this%unit .ne. 0) call this%close()
       this%unit = coop_free_file_unit()
       call this%open(trim(this%filename))
    endif
    call this%update_hdu_info()
    call this%load_header()
    call this%check_error()
#else
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_open_image


  subroutine coop_fits_file_open_table(this, filename)
    class(coop_fits_file)::this
    COOP_UNKNOWN_STRING::filename
#if HAS_CFITSIO
    if(this%unit .ne. 0) call this%close()
    this%unit = coop_free_file_unit()

    this%filename = trim(adjustl(filename))
    call fttopn(this%unit, trim(this%filename), this%rwmode, this%status)
    if(this%status .ne. 0)then
       this%status = 0
       call this%open(trim(this%filename))
    endif
    call this%check_error()
    call this%update_hdu_info()
    call this%load_header()
    call this%check_error()
#else
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_open_table

  subroutine coop_fits_file_move_to_hdu(this, ihdu)
    class(coop_fits_file)::this
    COOP_INT::ihdu
#if HAS_CFITSIO
    call ftmahd(this%unit, ihdu,  this%hdutype, this%status)
    this%chdu = ihdu
    call this%load_header()
    call this%check_error()
#else
    stop "CFITSIO is not installed"
#endif    
  end subroutine coop_fits_file_move_to_hdu

  subroutine coop_Cls_convert2Pseudo(this, lmin, lmax, kernel)
    class(coop_cls)::this
    COOP_INT::lmin, lmax
    COOP_REAL::kernel(lmin:lmax,lmin:lmax, 4)
    COOP_REAL,dimension(:,:),allocatable::tmp
    if(this%lmin .lt. lmin .or. this%lmax .gt. lmax)then
       write(*,*) "fits_cls2pseudoCls: l range overflow"
    endif
    select case(trim(this%genre))
    case("I", "T")
       this%Cls(this%lmin:this%lmax, 1) = matmul(kernel(this%lmin:this%lmax, this%lmin:this%lmax, coop_pseudoCl_kernel_index_TT),this%Cls(this%lmin:this%lmax, 1))
    case("QU", "EB")
       allocate(tmp(this%lmin:this%lmax, 2))
       this%Cls(this%lmin:this%lmax, coop_EB_index_EB) = matmul(kernel(this%lmin:this%lmax, this%lmin:this%lmax, coop_pseudoCl_kernel_index_EB),this%Cls(this%lmin:this%lmax, coop_EB_index_EB))
       tmp(:,1) = matmul(kernel(this%lmin:this%lmax, this%lmin:this%lmax, coop_pseudoCl_kernel_index_EE_plus_BB),this%Cls(this%lmin:this%lmax, coop_EB_index_EE)+this%Cls(this%lmin:this%lmax, coop_EB_index_BB))
       tmp(:,2) = matmul(kernel(this%lmin:this%lmax, this%lmin:this%lmax, coop_pseudoCl_kernel_index_EE_minus_BB),this%Cls(this%lmin:this%lmax, coop_EB_index_EE)-this%Cls(this%lmin:this%lmax, coop_EB_index_BB))
       this%Cls(this%lmin:this%lmax, coop_EB_index_EE) = (tmp(:,1)+tmp(:,2))/2.d0
       this%Cls(this%lmin:this%lmax, coop_EB_index_BB) = (tmp(:,1)-tmp(:,2))/2.d0
       deallocate(tmp)
    case("IQU", "TQU", "TEB", "IEB")
       allocate(tmp(this%lmin:this%lmax, 2))
       this%Cls(this%lmin:this%lmax, coop_TEB_index_TT) = matmul(kernel(this%lmin:this%lmax, this%lmin:this%lmax, coop_pseudoCl_kernel_index_TT),this%Cls(this%lmin:this%lmax, coop_TEB_index_TT))
       this%Cls(this%lmin:this%lmax, coop_TEB_index_TE) = matmul(kernel(this%lmin:this%lmax, this%lmin:this%lmax, coop_pseudoCl_kernel_index_TE),this%Cls(this%lmin:this%lmax, coop_TEB_index_TE))
       this%Cls(this%lmin:this%lmax, coop_TEB_index_TB) = matmul(kernel(this%lmin:this%lmax, this%lmin:this%lmax, coop_pseudoCl_kernel_index_TB),this%Cls(this%lmin:this%lmax, coop_TEB_index_TB))
       this%Cls(this%lmin:this%lmax, coop_TEB_index_EB) = matmul(kernel(this%lmin:this%lmax, this%lmin:this%lmax, coop_pseudoCl_kernel_index_EB),this%Cls(this%lmin:this%lmax, coop_TEB_index_EB))
       tmp(:,1) = matmul(kernel(this%lmin:this%lmax, this%lmin:this%lmax, coop_pseudoCl_kernel_index_EE_plus_BB),this%Cls(this%lmin:this%lmax, coop_TEB_index_EE)+this%Cls(this%lmin:this%lmax, coop_TEB_index_BB))
       tmp(:,2) = matmul(kernel(this%lmin:this%lmax, this%lmin:this%lmax, coop_pseudoCl_kernel_index_EE_minus_BB),this%Cls(this%lmin:this%lmax, coop_TEB_index_EE)-this%Cls(this%lmin:this%lmax, coop_TEB_index_BB))
       this%Cls(this%lmin:this%lmax, coop_TEB_index_EE) = (tmp(:,1)+tmp(:,2))/2.d0
       this%Cls(this%lmin:this%lmax, coop_TEB_index_BB) = (tmp(:,1)-tmp(:,2))/2.d0
       deallocate(tmp)
    case default
       write(*,*) "genre = "//trim(this%genre)
       stop "This genre is not supported for Cls_convert2Pseudo."
    end select
  end subroutine coop_Cls_convert2Pseudo



end module coop_fitsio_mod


