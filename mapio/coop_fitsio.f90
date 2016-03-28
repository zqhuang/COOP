module coop_fitsio_mod
  use coop_wrapper_utils
  use coop_sphere_mod
  use coop_stacking_mod
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
     procedure::load_int_column => coop_fits_file_load_int_column
     procedure::get_nrows_ncols => coop_fits_file_get_nrows_ncols
     procedure::get_naxes => coop_fits_file_get_naxes
     procedure::get_bitpix => coop_fits_file_get_bitpix
  end type coop_fits_file


contains

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

  subroutine coop_fits_file_load_double_column(this, col, data, bad_value)
    class(coop_fits_file)::this
    COOP_INT::col
    COOP_REAL,dimension(:)::data
    COOP_SINGLE,dimension(:),allocatable::single_data
    COOP_INT,dimension(:),allocatable::int_data
    COOP_SHORT_INT,dimension(:),allocatable::short_int_data
    COOP_REAL, optional::bad_value
    COOP_REAL::nulval
    COOP_SHORT_INT::short_int_nulval
    COOP_INT::ncols, width, repeat, datacode, nrows
    logical anyf
#if HAS_CFITSIO
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
    if(nrows .lt. size(data))then
       write(*,*) trim(this%filename)
       write(*,*) "number of table rows:", nrows
       write(*,*) "Error: want rows ", size(data)
       stop
    endif
    call ftgtcl(this%unit, col, datacode, repeat, width, this%status)
    if(repeat.gt.1)then
       call this%report_error("repeat >1 is not supported in coop_fitsio")
    endif
    select case(datacode)
    case(coop_fitsio_datatype_double)
       call ftgcvd(this%unit, col, 1, 1, size(data), nulval, data, anyf, this%status)
    case(coop_fitsio_datatype_single)
       allocate(single_data(size(data)))
       call ftgcve(this%unit, col, 1, 1, size(data), real(nulval), single_data, anyf, this%status)
       data = single_data
       deallocate(single_data)
    case(coop_fitsio_datatype_int)
       write(*,*) "Warning: integer column in "//trim(this%filename)//" is loaded as double"
       allocate(int_data(size(data)))
       call ftgcvj(this%unit, col, 1, 1, size(data), nint(nulval), int_data, anyf, this%status)
       data = int_data
       deallocate(int_data)
    case(coop_fitsio_datatype_short_int)
       write(*,*) "Warning: short integer column in "//trim(this%filename)//" is loaded as double"
       allocate(short_int_data(size(data)))
       short_int_nulval = nint(nulval)
       call ftgcvi(this%unit, col, 1, 1, size(data), short_int_nulval, short_int_data, anyf, this%status)
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


  subroutine coop_fits_file_load_int_column(this, col, data, bad_value)
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
    logical anyf
#if HAS_CFITSIO
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
    if(nrows .lt. size(data))then
       write(*,*) trim(this%filename)
       write(*,*) "number of table rows:", nrows
       write(*,*) "Error: want rows ", size(data)
       stop
    endif
    call ftgtcl(this%unit, col, datacode, repeat, width, this%status)
    if(repeat.gt.1)then
       call this%report_error("repeat >1 is not supported in coop_fitsio")
    endif
    select case(datacode)
    case(coop_fitsio_datatype_int)
       call ftgcvj(this%unit, col, 1, 1, size(data), nulval, data, anyf, this%status)
    case(coop_fitsio_datatype_single)
       write(*,*) "Warning: float column in "//trim(this%filename)//" is loaded as int"       
       allocate(single_data(size(data)))
       call ftgcve(this%unit, col, 1, 1, size(data), real(nulval), single_data, anyf, this%status)
       data = nint(single_data)
       deallocate(single_data)
    case(coop_fitsio_datatype_double)
       write(*,*) "Warning: double column in "//trim(this%filename)//" is loaded as int"
       allocate(double_data(size(data)))
       call ftgcvd(this%unit, col, 1, 1, size(data), nint(nulval), double_data, anyf, this%status)
       data = nint(double_data)
       deallocate(double_data)
    case(coop_fitsio_datatype_short_int)
       allocate(short_int_data(size(data)))
       short_int_nulval = nint(nulval)
       call ftgcvi(this%unit, col, 1, 1, size(data), short_int_nulval, short_int_data, anyf, this%status)
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
          call ftgpvd(this%unit, 1, 1, size(image), real(bad_value), fimage, anynul, this%status)
       else
          call ftgpvd(this%unit, 1, 1, size(image), 0., fimage, anynul, this%status)
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
          call ftgpvd(this%unit, 1, 1, size(image), real(bad_value), fimage, anynul, this%status)
       else
          call ftgpvd(this%unit, 1, 1, size(image), 0., fimage, anynul, this%status)
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
          call ftgpvd(this%unit, 1, 1, size(image), real(bad_value), fimage, anynul, this%status)
       else
          call ftgpvd(this%unit, 1, 1, size(image), 0., fimage, anynul, this%status)
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
    logical simple, extend
    COOP_INT::i
    COOP_INT::bitpix, naxis, group, fpixel, nelements, naxes(1024)
    COOP_FITSIO_CARD::card
#if HAS_CFITSIO
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
    simple = .true.
    extend = .true.
    group = 1
    fpixel = 1
    call ftphpr(fp%unit, simple, bitpix, naxis, naxes(1:naxis), 0, 1, extend, fp%status)
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

  subroutine coop_fits_file_write_binary_table(table, filename, header)
    type(coop_dictionary),optional::header
    COOP_REAL,dimension(:,:)::table
    COOP_UNKNOWN_STRING::filename
    type(coop_fits_file)::fp
    logical simple, extend
    COOP_INT::i
    COOP_INT::status,readwrite,blocksize,hdutype,tfields,nrows,varidat, column,frow,felem
    COOP_SHORT_STRING::extname
    COOP_SHORT_STRING,dimension(:),allocatable::ttype, tform, tunit
    COOP_SHORT_INT,dimension(:),allocatable::short_int_col
    COOP_FITSIO_CARD::card
#if HAS_CFITSIO
    call fp%init(filename)
    readwrite = 1
    call ftcrhd(fp%unit, fp%status)
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
    call ftphbn(fp%unit, nrows, tfields, ttype, tform, tunit, extname, varidat, fp%status)
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
          if(.not. allocated(short_int_col))allocate(short_int_col(nrows))
          short_int_col = nint(table(:, column))
          call ftpcli(fp%unit, column, frow, felem, nrows, short_int_col, fp%status)
       case default
          write(*,*) "TFORM"//COOP_STR_OF(column)//"="//trim(tform(column))
          stop "coop_fits_file_write_binary: does not support this TFORM"
       end select
    enddo
    COOP_DEALLOC(short_int_col)
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
    stop "CFITSIO is not installed"
#endif
  end subroutine coop_fits_file_write_binary_table


  subroutine coop_fits_file_write_2d_image(image, filename, header)
    type(coop_dictionary),optional::header
    COOP_REAL,dimension(:,:)::image
    COOP_UNKNOWN_STRING::filename
    type(coop_fits_file)::fp
    logical simple, extend
    COOP_INT::i
    COOP_INT::bitpix, naxis, naxes(2), group, fpixel, nelements
    COOP_FITSIO_CARD::card

#if HAS_CFITSIO
    if(present(header))then
       call coop_dictionary_lookup(header, "BITPIX", bitpix, -64)
    else
       bitpix = -64
    endif
    call fp%init(filename)
    simple = .true.
    extend = .true.
    naxis = 2
    group = 1
    fpixel = 1
    nelements = size(image)
    naxes(1) = size(image,1)
    naxes(2) = size(image,2)
    call ftphpr(fp%unit, simple, bitpix, naxis, naxes, 0, 1, extend, fp%status)
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
  end subroutine coop_fits_file_write_2d_image


  subroutine coop_fits_file_write_3d_image(image, filename, header)
    type(coop_dictionary),optional::header
    COOP_REAL,dimension(:,:,:)::image
    COOP_UNKNOWN_STRING::filename
    type(coop_fits_file)::fp
    logical simple, extend
    COOP_INT::i
    COOP_FITSIO_CARD::card
    COOP_INT::bitpix, naxis, naxes(3), group, fpixel, nelements
#if HAS_CFITSIO
    if(present(header))then
       call coop_dictionary_lookup(header, "BITPIX", bitpix, -64)
    else
       bitpix = -64
    endif
    call fp%init(filename)
    simple = .true.
    extend = .true.
    naxis = 2
    group = 1
    fpixel = 1
    nelements = size(image)
    naxes(1) = size(image,1)
    naxes(2) = size(image,2)
    naxes(3) = size(image,3)
    call ftphpr(fp%unit, simple, bitpix, naxis, naxes, 0, 1, extend, fp%status)
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
  end subroutine coop_fits_file_write_3d_image


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
       this%chdu = ihdu
    else
       call ftopen(this%unit, trim(this%filename), this%rwmode, this%blocksize, this%status)
       this%chdu = 1
    endif
    call ftthdu(this%unit, this%nhdus, this%status)
    call ftghdt(this%unit, this%hdutype, this%status)
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
    call ftdopn(this%unit, trim(this%filename), this%rwmode, this%status)
    call ftthdu(this%unit, this%nhdus, this%status)
    call ftghdn(this%unit, this%chdu)
    call ftghdt(this%unit, this%hdutype, this%status)
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
    call ftthdu(this%unit, this%nhdus, this%status)
    call ftghdn(this%unit, this%chdu)
    call ftghdt(this%unit, this%hdutype, this%status)
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



end module coop_fitsio_mod


