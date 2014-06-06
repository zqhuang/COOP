module coop_fitswrap_mod
  use coop_wrapper_utils
  implicit none

#include "constants.h"

#define COOP_FITS_IMAGE_DIM 2
  private

  public::coop_fits, coop_fits_image, coop_fits_image_cea

  integer,parameter::sp = kind(1.)
  integer,parameter::dl = kind(1.d0)

  type coop_fits
     COOP_STRING::filename
     type(coop_dictionary):: header
   contains
     procedure::open => coop_fits_open
     procedure::get_header => coop_fits_get_header
     procedure::key_value=>coop_fits_key_value
  end type coop_fits

  type, extends(coop_fits)::coop_fits_image
     COOP_INT::bitpix
     COOP_INT n(COOP_FITS_IMAGE_DIM)
#if COOP_FITS_IMAGE_DIM == 2
     real(sp),dimension(:,:),allocatable::image
#elif COOP_FITS_IMAGE_DIM == 3
     real(sp),dimension(:,:,:),allocatable::image
#endif
     COOP_REAL::transform(COOP_FITS_IMAGE_DIM, COOP_FITS_IMAGE_DIM), center(COOP_FITS_IMAGE_DIM)
   contains
     procedure::get_linear_coordinates => coop_fits_image_get_linear_coordinates
     procedure::get_data => coop_fits_image_get_data
  end type coop_fits_image

  type, extends(coop_fits_image)::coop_fits_image_cea
     contains
       procedure::get_coordinates => coop_fits_image_cea_get_coordinates
  end type coop_fits_image_cea


contains
  
  subroutine coop_fits_open(this, filename)
    class(coop_fits)::this
    COOP_UNKNOWN_STRING::filename
    this%filename = trim(filename)
    call this%get_header()
  end subroutine coop_fits_open

  subroutine coop_fits_get_header(this)
    class(coop_fits)::this
    COOP_LONG_STRING::header
    integer nkeys, i, j, istart, iend
    COOP_REAL delta(COOP_FITS_IMAGE_DIM)
    call coop_fits_print_header(trim(this%filename), header, nkeys)
    istart = 1
    do i=1, nkeys
       j = scan(header(istart:), "=")
       iend = scan(header(istart:), coop_newline)
       j = j+ istart - 1
       iend = iend + istart - 2
       call this%header%insert(header(istart:j-1), header(j+1:iend))
       istart = iend + 2
    enddo
    select type(this)
    class is(coop_fits_image_cea)
       do i=1, COOP_FITS_IMAGE_DIM
          do j=1, COOP_FITS_IMAGE_DIM
             this%transform(i, j) = this%key_value("PC"//trim(coop_num2str(i))//"_"//trim(coop_num2str(j)))
          enddo
          delta(i) = this%key_value("CDELT"//trim(coop_num2str(i))) * coop_SI_degree
          this%center(i) = this%key_value("CRPIX"//trim(coop_num2str(i)))
       enddo
       do i=1, COOP_FITS_IMAGE_DIM
          this%transform(:, i) =  this%transform(:, i) * delta(i)
       enddo
    end select
    select type(this)
    class is(coop_fits_image)
       this%bitpix = nint(this%key_value("BITPIX"))
       if(nint(this%key_value("NAXIS")) .ne. COOP_FITS_IMAGE_DIM)then
          write(*,*) trim(this%filename)//": this is not an "//trim(coop_num2str(COOP_FITS_IMAGE_DIM))//"D image"
          call this%header%print()
          stop
       endif
       do i=1, COOP_FITS_IMAGE_DIM
          this%n(i) = nint(this%key_value("NAXIS"//trim(coop_num2str(i))))
       enddo
       if(any(this%n .eq. 0))then
          write(*,*) trim(this%filename)//": cannot read the dimensions"
          call this%header%print()
          stop
       endif
    end select
  end subroutine coop_fits_get_header

  function coop_fits_key_value(this, key) result(val)
    class(coop_fits)::this
    COOP_UNKNOWN_STRING::key
    COOP_SHORT_STRING::str
    COOP_REAL val
    str = this%header%value(key)
    if(trim(str).ne."")then
       val = coop_str2real(str)
    else
       val = 0.d0
    endif
  end function coop_fits_key_value


  subroutine coop_fits_image_get_linear_coordinates(this, pix, coor)
    class(coop_fits_image)::this
    COOP_INT pix(COOP_FITS_IMAGE_DIM)
    COOP_REAL coor(COOP_FITS_IMAGE_DIM)
    coor = matmul(this%transform, pix - this%center)
  end subroutine coop_fits_image_get_linear_coordinates
  
  subroutine coop_fits_image_cea_get_coordinates(this, pix, coor)
    class(coop_fits_image_cea)::this
    COOP_INT pix(COOP_FITS_IMAGE_DIM)
    COOP_REAL coor(COOP_FITS_IMAGE_DIM)
    call this%get_linear_coordinates(pix, coor)
    coor(2) = asin(coor(2))
  end subroutine coop_fits_image_cea_get_coordinates

  subroutine coop_fits_image_get_data(this)
    class(coop_fits_image)::this
    integer(8) nelements
    integer n, i
#if  COOP_FITS_IMAGE_DIM == 2
    real(dl),dimension(:,:),allocatable::tmp
#elif COOP_FITS_IMAGE_DIM == 3
    real(dl),dimension(:,:,:),allocatable::tmp
#endif
    nelements = this%n(1)
    do i=2, COOP_FITS_IMAGE_DIM
       nelements = nelements*this%n(i)
    enddo
    if(allocated(this%image))then
       if(size(this%image).ne. nelements)then
          deallocate(this%image)
#if COOP_FITS_IMAGE_DIM == 2
          allocate(this%image(this%n(1), this%n(2)))
#elif COOP_FITS_IMAGE_DIM == 3
          allocate(this%image(this%n(1), this%n(2), this%n(3)))
#else
          stop "COOP only support dimension = 2 or 3"
#endif
       endif
    else
#if COOP_FITS_IMAGE_DIM == 2
       allocate(this%image(this%n(1), this%n(2)))
#elif COOP_FITS_IMAGE_DIM == 3
       allocate(this%image(this%n(1), this%n(2), this%n(3)))
#else
       stop "COOP only support dimension = 2 or 3"
#endif
    endif
    select case(this%bitpix)
    case(-32)
       call coop_fits_get_float_data(trim(this%filename), this%image, nelements)
    case(-64)
#if COOP_FITS_IMAGE_DIM == 2
       allocate(tmp(this%n(1), this%n(2)))
#elif COOP_FITS_IMAGE_DIM == 3
       allocate(tmp(this%n(1), this%n(2), this%n(3)))
#else
       stop "COOP only support dimension = 2 or 3"
#endif
       call coop_fits_get_double_data(trim(this%filename), tmp, nelements)
       this%image = real(tmp, sp)
       deallocate(tmp)
    case default
       write(*,*) "Cannot load data for bitpix = "//trim(coop_num2str(this%bitpix))
       stop
    end select

  end subroutine coop_fits_image_get_data


end module coop_fitswrap_mod
