module coop_fitswrap_mod
  use coop_wrapper_utils
  use coop_sphere_mod
  implicit none

#include "constants.h"

  private

  public::coop_fits, coop_fits_image, coop_fits_image_cea

  integer,parameter::sp = kind(1.)
  integer,parameter::dl = kind(1.d0)


  type coop_fits
     COOP_STRING::filename
     type(coop_dictionary):: header
   contains
     procedure::open => coop_fits_open
     procedure::key_value=>coop_fits_key_value
  end type coop_fits

  type, extends(coop_fits)::coop_fits_image
     COOP_INT::bitpix
     COOP_INT::dim
     COOP_INT,dimension(:),allocatable::nside
     COOP_LONG_INT::npix
     real(sp),dimension(:),allocatable::image
     COOP_REAL,allocatable::transform(:, :), center(:)
   contains
     procedure::free => coop_fits_image_free
     procedure::get_linear_coordinates => coop_fits_image_get_linear_coordinates
     procedure::get_data => coop_fits_image_get_data
  end type coop_fits_image



  type, extends(coop_fits_image)::coop_fits_image_cea
     contains
       procedure::pix2ang => coop_fits_image_cea_pix2ang
       procedure::pix2flat => coop_fits_image_cea_pix2flat
  end type coop_fits_image_cea



contains

  subroutine coop_fits_open(this, filename)
    class(coop_fits)::this
    COOP_UNKNOWN_STRING::filename
    if(coop_file_exists(filename))then
       this%filename = trim(filename)
       call coop_fits_get_header(this)
    else
       write(*,*) "The file "//trim(filename)//" does not exist."
    endif
  end subroutine coop_fits_open

  subroutine coop_fits_image_free(this)
    class(coop_fits_image) :: this
    if(allocated(this%image)) deallocate(this%image)
    if(allocated(this%transform))deallocate(this%transform)
    if(allocated(this%center))deallocate(this%center)
    if(allocated(this%nside))deallocate(this%nside)
  end subroutine coop_fits_image_free

  subroutine coop_fits_get_header(this)
    class(coop_fits)::this
    COOP_LONG_STRING::header
    integer nkeys, i, j, istart, iend
    COOP_REAL,dimension(:),allocatable::delta
    call coop_fits_read_header_to_string(trim(this%filename), header, nkeys)
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
    type is (coop_fits)
       return
    class is(coop_fits_image)
       this%bitpix = nint(this%key_value("BITPIX"))
       this%dim = nint(this%key_value("NAXIS"))
       if(this%dim .eq. 0 )then
          call this%header%print()
          write(*,*) "Error: cannot find NAXIS key word in fits file "//trim(this%filename)
          stop
       endif
       call this%free()
       allocate(this%nside(this%dim))
       allocate(this%transform(this%dim, this%dim), this%center(this%dim))
       this%npix = 1
       do i=1, this%dim
          this%nside(i) = nint(this%key_value("NAXIS"//trim(coop_num2str(i))))
          this%npix = this%npix * this%nside(i)
       enddo
       if(any(this%nside .eq. 0))then
          write(*,*) trim(this%filename)//": cannot read the dimensions"
          call this%header%print()
          stop
       endif
       allocate(this%image(0:this%npix-1))
    end select
    !!get transform, center 
    select type(this)
    class is(coop_fits_image_cea)
       if(this%dim .ne. 2) stop "For CEA map the dimension must be 2"
       allocate(delta(this%dim))
       do i=1, this%dim
          do j=1, this%dim
             this%transform(i, j) = this%key_value("PC"//trim(coop_num2str(i))//"_"//trim(coop_num2str(j)))
          enddo
          delta(i) = this%key_value("CDELT"//trim(coop_num2str(i))) * coop_SI_degree
          this%center(i) = this%key_value("CRPIX"//trim(coop_num2str(i)))
       enddo
       do i=1, this%dim
          this%transform(:, i) =  this%transform(:, i) * delta(i)
       enddo
       deallocate(delta)
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
    COOP_INT pix(this%dim)
    COOP_REAL coor(this%dim)
    coor = matmul(this%transform, pix - this%center)
  end subroutine coop_fits_image_get_linear_coordinates


  subroutine coop_fits_image_get_data(this)
    class(coop_fits_image)::this
    real(dl),dimension(:),allocatable::tmp
    select case(this%bitpix)
    case(-32)
       call coop_fits_get_float_data(trim(this%filename), this%image, this%npix)
    case(-64)
       allocate(tmp(0:this%npix-1))
       call coop_fits_get_double_data(trim(this%filename), tmp, this%npix)
       this%image = real(tmp, sp)
       deallocate(tmp)
    case default
       write(*,*) "Cannot load data for bitpix = "//trim(coop_num2str(this%bitpix))
       stop
    end select

  end subroutine coop_fits_image_get_data


!!=======   CEA ============

  subroutine coop_fits_image_cea_pix2ang(this, pix, theta, phi)
    class(coop_fits_image_cea)::this
    COOP_LONG_INT pix
    COOP_REAL theta, phi
    COOP_INT ix, iy
    iy = pix/this%nside(1) 
    ix = pix  - iy*this%nside(1) + 1 - this%center(1)
    iy = iy + 1 - this%center(2)
    phi = this%transform(1,1)* ix + this%transform(1,2)*iy
    theta = coop_pio2 - asin(this%transform(2, 1)*ix + this%transform(2,2)*iy)
  end subroutine coop_fits_image_cea_pix2ang


  subroutine coop_fits_image_cea_pix2flat(this, disc, pix, coor)
    class(coop_fits_image_cea)::this
    COOP_LONG_INT:: pix
    COOP_REAL ::coor(2)
    type(coop_sphere_disc)::disc
    call this%pix2ang(pix, coor(1), coor(2))
    call disc%ang2flat(coor)
  end subroutine coop_fits_image_cea_pix2flat



end module coop_fitswrap_mod
