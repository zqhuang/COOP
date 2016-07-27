module coop_fitswrap_mod
  use coop_wrapper_utils
  use coop_sphere_mod
  use coop_fitsio_mod
  use coop_healpix_mod
  use coop_stacking_mod
  implicit none

#include "constants.h"

  private

  public::coop_fits, coop_fits_image, coop_fits_image_cea, coop_fits_QU2EB, coop_fits_EB2QU, coop_fits_image_cea_simulate_TEB, coop_fits_image_cea_EB2QU, coop_fits_image_cea_QU2EB, coop_flatsky_maps

  type coop_fits
     COOP_STRING::filename
     type(coop_dictionary):: header
   contains
     procedure::open => coop_fits_open
     procedure::read => coop_fits_open
     procedure::free => coop_fits_free
     procedure::key_value=>coop_fits_key_value
  end type coop_fits

  type, extends(coop_fits)::coop_fits_image
     COOP_INT::bitpix = -64
     COOP_INT::dim = 2
     COOP_INT,dimension(:),allocatable::nside
     COOP_INT::npix = 0
     COOP_REAL,dimension(:),allocatable::image
     COOP_REAL,allocatable::transform(:, :), invtrans(:,:), center(:)
     COOP_INT,allocatable:: icenter(:)
   contains
     procedure::regularize => coop_fits_image_regularize
     procedure::get_linear_coordinates => coop_fits_image_get_linear_coordinates
     procedure::get_data => coop_fits_image_get_data
     procedure::simple_stat => coop_fits_image_simple_stat
  end type coop_fits_image

  type, extends(coop_fits_image)::coop_fits_image_cea
     COOP_REAL smooth_pixsize, pixsize, smooth_dkx, smooth_dky, dkx, dky, dx, dy
     COOP_REAL xmin, xmax, ymin, ymax
     COOP_INT:: smooth_nx, smooth_ny, smooth_npix
     type(coop_sphere_disc)::disc
     COOP_REAL,allocatable::radec_center(:)
     COOP_REAL,dimension(:,:),allocatable:: smooth_image
     COOP_REAL,dimension(:,:),allocatable:: smooth_Q, smooth_U
   contains
     procedure::convert2healpix => coop_fits_image_cea_convert2healpix
     procedure::from_healpix => coop_fits_image_cea_from_healpix
     procedure::get_neighbours => coop_fits_image_cea_get_neighbours
     procedure::write => coop_fits_image_cea_write
     procedure::pix2ixiy => coop_fits_image_cea_pix2ixiy
     procedure::ixiy2pix => coop_fits_image_cea_ixiy2pix
     procedure::pix2ang => coop_fits_image_cea_pix2ang
     procedure::radec2pix => coop_fits_image_cea_radec2pix
     procedure::get_flatmap => coop_fits_image_cea_get_flatmap
     procedure::pix2flat => coop_fits_image_cea_pix2flat
     procedure::cut => coop_fits_image_cea_cut
     procedure::filter => coop_fits_image_cea_filter
     procedure::smooth => coop_fits_image_cea_smooth
     procedure::get_power => coop_fits_image_cea_get_power
     procedure::smooth_flat => coop_fits_image_cea_smooth_flat
     procedure::get_QTUT =>  coop_fits_image_cea_get_QTUT
     procedure::get_QU =>  coop_fits_image_cea_get_QU
     procedure::find_extrema => coop_fits_image_cea_find_extrema
     procedure::stack => coop_fits_image_cea_stack
     procedure::stack2fig => coop_fits_image_cea_stack2fig
     procedure::simulate => coop_fits_image_cea_simulate
     procedure::simulate_flat => coop_fits_image_cea_simulate_flat
     procedure::plot => coop_fits_image_cea_plot
     procedure::EB2QU => coop_fits_image_cea_smooth_EB2QU
     procedure::QU2EB => coop_fits_image_cea_smooth_QU2EB
  end type coop_fits_image_cea

  type coop_flatsky_maps
     COOP_STRING::path=""
     COOP_STRING::filename=""
     COOP_INT::nmaps = 0
     COOP_INT::npix = 0
     COOP_INT,dimension(2)::nside = 0
     logical::has_mask = .false.
     logical::mask_changed = .false.
     COOP_REAL::dx, dy, dkx, dky
     COOP_REAL::total_weight = 0.d0, mask_threshold = 0.d0
     COOP_INT,dimension(:),allocatable::spin
     COOP_SHORT_STRING,dimension(:),allocatable::fields
     COOP_SHORT_STRING,dimension(:),allocatable::units
     COOP_STRING,dimension(:),allocatable::files
     logical,dimension(:),allocatable::map_changed
     COOP_STRING::fmask
     type(coop_fits_image_cea),dimension(:), allocatable::map
     type(coop_fits_image_cea)::mask
     logical,dimension(:),allocatable::unmasked
   contains
     procedure::free => coop_flatsky_maps_free
     procedure::read => coop_flatsky_maps_read
     procedure::read_from_one => coop_flatsky_maps_read_from_one
     procedure::open => coop_flatsky_maps_read
     procedure::write => coop_flatsky_maps_write
     procedure::weighted_sum => coop_flatsky_maps_weighted_sum
     procedure::pixel_values => coop_flatsky_maps_pixel_values
     procedure::is_unmasked => coop_flatsky_maps_is_unmasked
     procedure::is_masked => coop_flatsky_maps_is_masked
     procedure::get_peaks => coop_flatsky_maps_get_peaks
     procedure::get_neighbours => coop_flatsky_maps_get_neighbours
     procedure::get_zeros => coop_flatsky_maps_get_zeros
     procedure::stack_on_peaks => coop_flatsky_maps_stack_on_peaks
     procedure::stack_on_patch => coop_flatsky_maps_stack_on_patch
     procedure::fetch_patch => coop_flatsky_maps_fetch_patch
     procedure::merge => coop_flatsky_maps_merge
     procedure::qu2eb => coop_flatsky_maps_qu2eb
     procedure::eb2qu => coop_flatsky_maps_eb2qu
     procedure::trim_mask => coop_flatsky_maps_trim_mask
     procedure::remove_mono => coop_flatsky_maps_remove_mono
  end type coop_flatsky_maps

contains

  subroutine coop_flatsky_maps_remove_mono(this)
    class(coop_flatsky_maps)::this
    COOP_REAL::mean
    COOP_INT::i
    if(this%has_mask)then
       do i=1, this%nmaps
          mean = this%weighted_sum(this%map(i)%image)/this%total_weight
          this%map(i)%image  = this%map(i)%image - mean
       enddo
    else
       do i=1, this%nmaps
          mean = sum(this%map(i)%image)/this%total_weight
          this%map(i)%image  = this%map(i)%image - mean
       enddo
    endif

   end subroutine coop_flatsky_maps_remove_mono

  subroutine coop_flatsky_maps_trim_mask(this, width)
    class(coop_flatsky_maps)::this
    COOP_REAL::width, rwx, rwy, ry2
    COOP_INT::ix, iy, pix, ibase, iwx, iwy, il, irange
    COOP_INT, dimension(:),allocatable::list
    COOP_REAL, dimension(:),allocatable::weight
    COOP_INT::nlist
    if(width.le.0.d0)return
    if(.not. this%has_mask) return
    iwx = ceiling(width/abs(this%dx))
    iwy = ceiling(width/abs(this%dy))
    rwx = dble(iwx)*1.0000001
    rwy = dble(iwy)*1.0000001
    nlist = 0
    do iy = -iwy, iwy
       irange = nint(iwx*sqrt(1.d0-(iy/rwy)**2))
       nlist= nlist + (2*irange + 1)
    enddo
    allocate(list(nlist), weight(nlist))
    il = 0
    do iy = -iwy, iwy
       irange = nint(iwx*sqrt(1.d0-(iy/rwy)**2))
       ibase = iy*this%nside(1)
       ry2 = (iy/rwy)**2
       do ix = - irange , irange
          il = il + 1
          list(il) = ibase + ix
          weight(il) = exp(-2.d0*((ix/rwx)**2 + ry2))
       enddo
    enddo
    weight = weight/sum(weight)
    
    do iy = 1, this%nside(2)- 1
       ibase = this%nside(1)*iy
       if(iy .lt. iwy .or. iy .ge. this%nside(2) - iwy)then
          this%mask%image(ibase:ibase+this%nside(1)-1) = 0.d0
          cycle
       endif
       this%mask%image(ibase:ibase+iwx-1) = 0.d0
       this%mask%image(ibase+this%nside(1)-iwx:ibase+this%nside(1)-1) = 0.d0
       do ix = iwx, this%nside(1) - iwx - 1
          pix = ibase + ix
          if(this%unmasked(pix))cycle
          this%mask%image(pix+list) = this%mask%image(pix+list) - weight
       enddo
    enddo
    this%unmasked = this%unmasked .and. (this%mask%image.gt.0.9) 
    where(this%unmasked)
       this%mask%image = 1.d0
    elsewhere
       this%mask%image = 0.d0
    end where
    deallocate(list, weight)
    this%mask_changed = .true.
  end subroutine coop_flatsky_maps_trim_mask

  subroutine coop_flatsky_maps_qu2eb(this)
    class(coop_flatsky_maps)::this
    COOP_INT::i
    i = 1
    do while(i.lt. this%nmaps)
       if(this%spin(i).eq.2 .and. this%spin(i+1).eq.2 .and. COOP_UPPER_STR(this%fields(i)).eq."Q" .and. COOP_UPPER_STR(this%fields(i+1)).eq. "U")then
          if(this%has_mask)then
             where(.not. this%unmasked)
                this%map(i)%image = this%map(i)%image * (1.d0-(1.d0-this%mask%image/this%mask_threshold)**2)
                this%map(i+1)%image = this%map(i+1)%image * (1.d0-(1.d0-this%mask%image/this%mask_threshold)**2)
             end where
          endif
          call coop_fits_image_cea_qu2eb(this%map(i), this%map(i+1))
          this%spin(i:i+1) = 0
          this%fields(i) = "E"
          this%fields(i+1) = "B"
          this%map_changed(i)= .true.
          this%map_changed(i+1)= .true.
          i = i + 2
       else
          i = i + 1
       endif
    enddo
  end subroutine coop_flatsky_maps_qu2eb


  subroutine coop_flatsky_maps_eb2qu(this)
    class(coop_flatsky_maps)::this
    COOP_INT::i
    i = 1
    do while(i.lt. this%nmaps)
       if(COOP_UPPER_STR(this%fields(i)).eq."E" .and. COOP_UPPER_STR(this%fields(i+1)).eq."B" .and. this%spin(i).eq.0 .and. this%spin(i+1).eq.0)then
          call coop_fits_image_cea_eb2qu(this%map(i), this%map(i+1))
          this%spin(i:i+1) = 2
          this%fields(i) = "Q"
          this%fields(i+1) = "U"
          this%map_changed(i) = .true.
          this%map_changed(i+1) = .true.
          i = i + 2
       else
          i = i + 1
       endif
    enddo
  end subroutine coop_flatsky_maps_eb2qu

  subroutine coop_flatsky_maps_merge(this, map)
    class(coop_flatsky_maps)::this
    type(coop_flatsky_maps)::map, copy
    COOP_INT::i
    select type(this)
    type is (coop_flatsky_maps)
       if(this%nmaps .le. 0 .or. map%nmaps .le.0) stop "nmaps = 0 error in flatsky_maps_merge"
       if(any(this%nside .ne. map%nside) .or. abs(this%dx/map%dx-1.d0).gt.1.d-5 .or. abs(this%dy/map%dy-1.d0).gt.1.d-5)then
          write(*,*) "map shape error in flatsky_maps_merge"
          stop
       endif
       copy = this
       deallocate(this%map, this%spin, this%fields, this%units, this%files, this%map_changed)
       this%nmaps = copy%nmaps + map%nmaps 
       allocate(this%spin(this%nmaps))
       allocate(this%fields(this%nmaps))
       allocate(this%files(this%nmaps))
       allocate(this%map_changed(this%nmaps))
       allocate(this%units(this%nmaps))
       allocate(this%map(this%nmaps))
       this%spin(1:copy%nmaps) = copy%spin
       this%fields(1:copy%nmaps) = copy%fields
       this%files(1:copy%nmaps) = copy%files
       this%map_changed(1:copy%nmaps) = copy%map_changed
       this%units(1:copy%nmaps) = copy%units
       this%spin(copy%nmaps+1: this%nmaps) = map%spin
       this%fields(copy%nmaps+1: this%nmaps) = map%fields
       this%units(copy%nmaps+1: this%nmaps) = map%units
       this%files(copy%nmaps+1: this%nmaps) = map%files
       this%map_changed(copy%nmaps+1: this%nmaps) = map%map_changed
       if(trim(this%path).ne.trim(map%path))then
          write(*,*) "merging flatsky maps in different directories."
          write(*,*) "creating symbolic link to resolve the problem."
          do i=1, map%nmaps
             if(trim(map%files(i)).ne."NULL")then
                call system("ln -s "//trim(map%path)//trim(map%files(i))//" "//trim(this%path)//trim(map%files(i)))
             endif
          enddo
       endif
       do i=1, copy%nmaps
          this%map(i) = copy%map(i)
       enddo
       do i=1, map%nmaps
          this%map(copy%nmaps+i) = map%map(i)
       enddo
       if(map%has_mask)then
          if(this%has_mask)then
             if(maxval(abs(this%mask%image - map%mask%image)).gt. 1.d-2)then
                this%mask%image = this%mask%image * map%mask%image
                this%mask_changed = .true.
             endif
          else
             this%mask = map%mask
             this%fmask = map%fmask
             this%has_mask = .true.
          endif

       endif
       call copy%free()
    class default
       stop "flatsky_maps_merge does not support extended classes"
    end select

  end subroutine coop_flatsky_maps_merge

  subroutine  coop_flatsky_maps_get_zeros(this, imap, zeros)
    class(coop_flatsky_maps)::this
    COOP_INT::imap
    logical zeros(0:this%npix-1)
    COOP_INT::ix , iy, pix, list(8), ibase, ibase_minus, ibase_plus
    zeros(0:this%nside(1)-1) = .false.
    zeros(this%npix-this%nside(1):this%npix-1) = .false.
    zeros(0:this%npix:this%nside(1)) = .false.
    zeros(this%nside(1)-1:this%npix:this%nside(1)) = .false.
    ibase = 0
    ibase_minus = ibase - this%nside(1)
    ibase_plus = ibase +  this%nside(1)
    do iy = 1, this%nside(2)-2
       ibase_minus = ibase
       ibase = ibase_plus
       ibase_plus = ibase_plus + this%nside(1)
       do ix = 1, this%nside(1) - 2
          list(1) = ibase_minus + ix - 1
          list(2) = ibase_minus + ix
          list(3) = ibase_minus + ix + 1
          list(4) = ibase + ix - 1
          list(5) = ibase + ix + 1
          list(6) = ibase_plus + ix - 1
          list(7) = ibase_plus + ix
          list(8) = ibase_plus + ix + 1
          zeros(ibase + ix) = any(this%map(imap)%image(list) .le. 0.d0) .and. any(this%map(imap)%image(list) .ge. 0.d0)
       enddo
    enddo
  end subroutine coop_flatsky_maps_get_zeros

  function coop_flatsky_maps_pixel_values(this, pix) result(map)
    class(coop_flatsky_maps)::this
    COOP_SINGLE::map(this%nmaps)
    COOP_INT::pix, i
    do i=1, this%nmaps
       map(i) = this%map(i)%image(pix)
    enddo
  end function coop_flatsky_maps_pixel_values

  function coop_flatsky_maps_is_unmasked(this, pix) result(unmasked)
    class(coop_flatsky_maps)::this
    COOP_INT::pix
    logical unmasked
    if(this%has_mask)then
       unmasked = this%unmasked(pix)
    else
       unmasked = .true.
    endif
  end function coop_flatsky_maps_is_unmasked


  function coop_flatsky_maps_is_masked(this, pix) result(masked)
    class(coop_flatsky_maps)::this
    COOP_INT::pix
    logical masked
    if(this%has_mask)then
       masked = .not. this%unmasked(pix)
    else
       masked = .false.
    endif
  end function coop_flatsky_maps_is_masked
  
  function coop_flatsky_maps_weighted_sum(this, array) result(s)
    class(coop_flatsky_maps)::this
    COOP_REAL::s
    COOP_REAL,dimension(:),intent(IN)::array
    if(size(array).ne.this%npix) stop "flatsky_maps_weighted_sum: array size error"
    if(this%has_mask)then
       s = sum(array, mask = this%unmasked)
    else
       s = sum(array)
    end if
  end function coop_flatsky_maps_weighted_sum

  !!write a text file (and maps)
  !!============================
  !! nmaps (integer)
  !! field1 (string)
  !! unit1 (string)
  !! spin1 (integer)
  !! map1 (string)
  !! field2
  !! unit2
  !! spin2
  !! map2 
  !! ...
  !! has_mask (logical)
  !! mask_file (string, present only if has_mask = T)
  !!=============================
  subroutine coop_flatsky_maps_write(this, filename, write_image, write_mask, index_list)
    class(coop_flatsky_maps)::this
    COOP_UNKNOWN_STRING::filename
    COOP_INT, optional::index_list(:)
    type(coop_file)::fp
    logical, optional::write_image, write_mask
    logical wn, wm
    COOP_INT::i, imap
    call fp%open(filename)
    write(fp%unit, "(A)") "#COOP flatsky map file"
    if(present(write_image))then
       wn = write_image
    else
       wn = .true.
    endif
    if(present(write_mask))then
       wm = write_mask
    else
       wm = .false.
    endif
    this%path=trim(coop_file_path_of(filename))
    this%filename=trim(coop_file_name_of(filename))
    if(present(index_list))then
       write(fp%unit, "(I5)") size(index_list)

       do imap = 1, size(index_list)
          i = index_list(imap)
          if(i.gt.this%nmaps) stop "flatsky_map_write: error index overflow"
          write(fp%unit, "(A)") trim(this%fields(i))
          write(fp%unit, "(A)") trim(this%units(i))
          write(fp%unit, "(I5)") this%spin(i)
          if(wn .or. .not. coop_file_exists(trim(this%path)//trim(this%files(i))) .or. this%map_changed(i) )then
             this%files(i) = trim(coop_file_replace_postfix(this%filename, "_MAP"//COOP_STR_OF(imap)//".fits"))
             call this%map(i)%write(trim(this%path)//trim(this%files(i)))
          endif
          write(fp%unit, "(A)") trim(this%files(i))
       enddo

    else
       write(fp%unit, "(I5)") this%nmaps

       do i = 1, this%nmaps
          write(fp%unit, "(A)") trim(this%fields(i))
          write(fp%unit, "(A)") trim(this%units(i))
          write(fp%unit, "(I5)") this%spin(i)
          if(wn .or. .not. coop_file_exists(trim(this%path)//trim(this%files(i))) .or. this%map_changed(i))then
             this%files(i) = trim(coop_file_replace_postfix(this%filename, "_MAP"//COOP_STR_OF(i)//".fits"))
             call this%map(i)%write(trim(this%path)//trim(this%files(i)))
          endif
          write(fp%unit, "(A)") trim(this%files(i))
       enddo
    endif
    if(this%has_mask)then
       write(fp%unit, "(A)") "T"
       if(wm .or. .not.coop_file_exists(trim(this%path)//trim(this%fmask)) .or. this%mask_changed)then
          write(*,*) wm , .not.coop_file_exists(trim(this%path)//trim(this%fmask)), this%mask_changed, trim(this%path)//trim(this%fmask)
          this%fmask = trim(coop_file_replace_postfix(this%filename, "_MASK.fits"))

          call this%mask%write(trim(this%path)//trim(this%fmask))
       endif
       write(fp%unit, "(A)") trim(this%fmask)
    else
       write(fp%unit, "(A)") "F"
    endif
    call fp%close()
  end subroutine coop_flatsky_maps_write

  !!read a text file
  !!============================
  !! nmaps (integer)
  !! field1 (string)
  !! unit1 (string)
  !! spin1 (integer)
  !! map1 (string)
  !! field2
  !! unit2
  !! spin2
  !! map2 
  !! ...
  !! has_mask (logical)
  !! mask_file (string, present only if has_mask = T)
  !!=============================
  subroutine coop_flatsky_maps_read(this, filename)
    class(coop_flatsky_maps)::this
    COOP_UNKNOWN_STRING::filename
    COOP_INT::i
    COOP_REAL,parameter::eps = 1.d-5
    type(coop_file)::fp
    call this%free()
    if(.not. coop_file_exists(filename))then
       write(*,*) trim(filename)//" does not exist"
       stop
    endif
    this%path=coop_file_path_of(filename)
    this%filename=coop_file_name_of(filename)
    call fp%open(filename, "r")
    if(.not.fp%read_int(this%nmaps))call reporterror()
    if(this%nmaps .le. 0 .or. this%nmaps .gt. 20)then
       write(*,*) "flatsky_read: nmaps must be between 1 and 20"
       stop
    endif
    allocate(this%map(this%nmaps))
    allocate(this%spin(this%nmaps))
    allocate(this%fields(this%nmaps))
    allocate(this%files(this%nmaps))
    allocate(this%map_changed(this%nmaps))
    this%map_changed = .false.
    allocate(this%units(this%nmaps))
    i = 0
    do i=1, this%nmaps
       if(.not. fp%read_string(this%fields(i)))call reporterror()
       if(.not. fp%read_string(this%units(i)))call reporterror()
       if(.not. fp%read_int(this%spin(i)))call reporterror()
       if(.not. fp%read_string(this%files(i)))call reporterror()
       call this%map(i)%open(trim(this%path)//trim(this%files(i)))
       if(i.eq.1)then
          this%dx = this%map(1)%dx
          this%dy = this%map(1)%dy
          this%dkx = this%map(1)%dkx
          this%dky = this%map(1)%dky
          this%nside = this%map(1)%nside(1:2)
          this%npix = this%map(1)%npix
       else
          if(abs(this%dx/this%map(i)%dx - 1.d0).gt. eps .or. &
               abs(this%dy / this%map(i)%dy - 1.d0) .gt. eps .or. &
               this%nside(1).ne. this%map(i)%nside(1) .or. this%nside(2).ne.this%map(i)%nside(2))then
             write(*,*) "flatsky_read error: not all maps have the same configuration"
             stop
          endif
       endif
    enddo
    if(.not. fp%read_logical(this%has_mask)) call reporterror()
    if(this%has_mask)then
       if(.not. fp%read_string(this%fmask))call reporterror()
       call this%mask%open(trim(this%path)//trim(this%fmask))
       if(abs(this%dx/this%mask%dx - 1.d0).gt. eps .or. &
            abs(this%dy / this%mask%dy - 1.d0) .gt. eps .or. &
            this%nside(1).ne. this%mask%nside(1) .or. this%nside(2).ne.this%mask%nside(2))then
          write(*,*) "flatsky_read error: not all maps have the same configuration"
          stop
       endif
       this%mask_threshold = coop_flatsky_mask_threshold(this%mask)
       allocate(this%unmasked(this%npix))
       this%unmasked = this%mask%image .gt. this%mask_threshold
       this%total_weight = count(this%unmasked)
    else
       this%total_weight = this%npix
    endif
    call this%remove_mono()
    call fp%close()
  contains
    subroutine reporterror()
      write(*,*) trim(filename)//" does not have the right format"
      stop
    end subroutine reporterror
  end subroutine coop_flatsky_maps_read


!!read from single fits file
  subroutine coop_flatsky_maps_read_from_one(this, filename, mask, nmaps, genre, unit)
    class(coop_flatsky_maps)::this
    COOP_UNKNOWN_STRING::filename
    COOP_UNKNOWN_STRING,optional::mask, genre, unit
    COOP_INT,optional::nmaps
    COOP_INT::i
    COOP_REAL,parameter::eps = 1.d-5
    type(coop_file)::fp
    call this%free()
    if(.not. coop_file_exists(filename))then
       write(*,*) trim(filename)//" does not exist"
       stop
    endif
    this%path=coop_file_path_of(filename)
    this%filename = "NULL"
    if(present(nmaps))then
       if(nmaps .le. 0 .or. nmaps.gt.20) stop "read_from_one: nmaps overflow"
       this%nmaps = nmaps
    else
       this%nmaps = 1
    endif
    allocate(this%map(this%nmaps))
    allocate(this%spin(this%nmaps))
    allocate(this%fields(this%nmaps))
    allocate(this%files(this%nmaps))
    allocate(this%map_changed(this%nmaps))
    this%map_changed = .false.
    allocate(this%units(this%nmaps))
    if(present(unit))then
       this%units = trim(unit)
    else
       this%units = "muK"  !!default
    endif
    if(present(genre))then
       select case(COOP_UPPER_STR(genre))
       case("I", "T", "E", "B", "ZETA")
          this%fields = trim(genre)
          this%spin = 0
       case("Q", "U")
          this%fields = COOP_UPPER_STR(genre)
          this%spin = 2
       case("QU")
          if(this%nmaps .ne. 2) stop "for genre = QU or QU nmaps must be 2"
          this%fields(1) = "Q"
          this%fields(2) = "U"
          this%spin = 2
       case("EB")
          if(this%nmaps .ne. 2) stop "for genre = EB nmaps must be 2"
          this%fields(1) = "E"
          this%fields(2) = "B"
          this%spin = 0
       case("IQU", "TQU")
          if(this%nmaps .ne. 3) stop "for genre = IQU or TQU nmaps must be 3"
          this%fields(1) = "I"
          this%fields(2) = "Q"
          this%fields(3) = "U"
          this%spin = (/ 0, 2, 2 /)
       case("TQTUT")
          if(this%nmaps .ne. 3) stop "for genre = TQTUT  nmaps must be 3"
          this%fields(1) = "T"
          this%fields(2) = "Q_T"
          this%fields(3) = "U_T"
          this%spin = (/ 0, 2, 2 /)
       case("ZQZUZ")
          if(this%nmaps .ne. 3) stop "for genre = ZQZUZ  nmaps must be 3"
          this%fields(1) = "\zeta"
          this%fields(2) = "Q_\zeta"
          this%fields(3) = "U_\zeta"
          this%spin = (/ 0, 2, 2 /)
       case("TEB", "IEB")
          if(this%nmaps .ne. 3) stop "for genre = TEB or IEB nmaps must be 3"
          this%fields(1) = "T"
          this%fields(2) = "E"
          this%fields(3) = "B"
          this%spin = (/ 0, 0, 0 /)
       case default
          write(*,*) "Unknown genre: "//trim(genre)
          stop
       end select
    else  !!try to guess from nmaps
       select case(this%nmaps)
       case(2)
          this%spin = 2
          this%fields(1) = "Q"
          this%fields(2) = "U"
       case(3)
          this%spin = (/ 0 , 2 , 2 /)
          this%fields(1) = "I"
          this%fields(2) = "Q"
          this%fields(2) = "U"
       case default
          this%spin = 0
          this%fields = "I"
       end select
    endif
    this%files(1) = coop_file_name_of(filename)
    call this%map(1)%open(filename)
    this%dx = this%map(1)%dx
    this%dy = this%map(1)%dy
    this%dkx = this%map(1)%dkx
    this%dky = this%map(1)%dky
    this%nside = this%map(1)%nside(1:2)
    this%npix = this%map(1)%npix
    this%has_mask = present(mask)
    if(this%has_mask)then
       this%fmask = coop_file_name_of(mask)
       if(trim(coop_file_path_of(mask)) .ne. trim(this%path))then
          call system("ln -s "//trim(mask)//" "//trim(this%path)//trim(this%fmask))
          write(*,*) "map and mask are not in the same directory; COOP is creating a symbolic link for the fsm format"
       endif

       call this%mask%open(mask)
       if(abs(this%dx/this%mask%dx - 1.d0).gt. eps .or. &
            abs(this%dy / this%mask%dy - 1.d0) .gt. eps .or. &
            this%nside(1).ne. this%mask%nside(1) .or. this%nside(2).ne.this%mask%nside(2))then
          write(*,*) "flatsky_read_from_one error: not all maps have the same configuration"
          stop
       endif
       this%mask_threshold = coop_flatsky_mask_threshold(this%mask)
       allocate(this%unmasked(this%npix))
       this%unmasked = this%mask%image .gt. this%mask_threshold
       this%total_weight = count(this%unmasked)
    else
       this%total_weight = this%npix
    endif
    do i=2, this%nmaps
       this%map(i) = this%map(1)
       this%files(i) = "NULL"
    enddo
  end subroutine coop_flatsky_maps_read_from_one

  

  subroutine coop_flatsky_maps_free(this)
    class(coop_flatsky_maps)::this
    COOP_INT::i
    if(allocated(this%map))then
       do i=1, this%nmaps
          call this%map(i)%free()
       enddo
       deallocate(this%map)
    endif
    COOP_DEALLOC(this%spin)
    COOP_DEALLOC(this%fields)
    COOP_DEALLOC(this%files)
    COOP_DEALLOC(this%map_changed)
    COOP_DEALLOC(this%units)
    COOP_DEALLOC(this%unmasked)
    call this%mask%free()
    this%nmaps = 0
    this%npix = 0
    this%nside = 0
    this%has_mask = .false.
    this%mask_changed = .false.
    this%path=""
    this%filename=""
  end subroutine coop_flatsky_maps_free

  subroutine coop_fits_open(this, filename)
    class(coop_fits)::this
    COOP_UNKNOWN_STRING::filename
    COOP_INT i
    call this%free()
    if(coop_file_exists(filename))then
       this%filename = trim(filename)
       call coop_fits_get_header(this)
    else
       write(*,*) "The file "//trim(filename)//" does not exist."
    endif
    select type(this)
    class is(coop_fits_image)
       call this%get_data()
    end select
  end subroutine coop_fits_open

  subroutine coop_fits_free(this)
    class(coop_fits) :: this
    call this%header%free()
    select type(this)
    class is(coop_fits_image)
       COOP_DEALLOC(this%image)
       COOP_DEALLOC(this%nside)
       COOP_DEALLOC(this%transform)
       COOP_DEALLOC(this%invtrans)       
    end select
    select type(this)
    class is (coop_fits_image_cea)
       COOP_DEALLOC(this%center)
       COOP_DEALLOC(this%icenter)
       COOP_DEALLOC(this%radec_center)    
       COOP_DEALLOC(this%smooth_image)
       COOP_DEALLOC(this%smooth_Q)
       COOP_DEALLOC(this%smooth_U)
    end select
  end subroutine coop_fits_free

  subroutine coop_fits_get_header(this)
    class(coop_fits)::this
    COOP_LONG_STRING::header
    COOP_INT i, j
    COOP_REAL,dimension(:),allocatable::delta, units
    call coop_fits_to_header(this%filename, this%header)
    select type(this)
    type is (coop_fits)
       return
    class is(coop_fits_image)
       this%bitpix = nint(this%key_value("BITPIX"))
       this%dim = nint(this%key_value("NAXIS", 2.d0))
       if(this%dim .eq. 0 )then
          call this%header%print()
          write(*,*) "Error: cannot find NAXIS key word in fits file "//trim(this%filename)
          stop
       endif
       allocate(this%nside(this%dim))
       this%npix = 1
       do i=1, this%dim
          this%nside(i) = nint(this%key_value("NAXIS"//COOP_STR_OF(i)))
          this%npix = this%npix * this%nside(i)
       enddo
       if(any(this%nside .eq. 0))then
          write(*,*) trim(this%filename)//": cannot read the dimensions"
          call this%header%print()
          stop
       endif
       allocate(this%image(0:this%npix-1))
       allocate(this%transform(this%dim, this%dim), this%center(this%dim), this%icenter(this%dim), this%invtrans(this%dim, this%dim))
       this%transform = 0.d0
       this%invtrans = 0.d0       
       do i=1, this%dim
          this%transform(i, i) = 1.d0
          this%invtrans(i, i) = 1.d0          
       enddo
       this%center = 1.d0
       this%icenter = 0
    end select

    !!get transform, center 
    select type(this)
    class is(coop_fits_image_cea)
       if(this%dim .ne. 2) stop "For CEA map the dimension must be 2"
       allocate(this%radec_center(this%dim))
       allocate(units(this%dim))              
       do i=1, this%dim
          if(index(this%header%value("CUNIT"//COOP_STR_OF(i)), "deg").ne.0)then
             units(i) = coop_SI_degree
          elseif(index(this%header%value("CUNIT"//COOP_STR_OF(i)), "arcmin").ne.0)then
             units(i) = coop_SI_arcmin
          else
             write(*,*) "Unknown unit "//trim(this%header%value("CUNIT"//COOP_STR_OF(i)))//"; Using deg as default."
             units(i) = coop_SI_degree
             call this%header%update("CUNIT"//COOP_STR_OF(i), "deg")
          end if
       enddo

       allocate(delta(this%dim))
       do i=1, this%dim
          do j=1, this%dim
             if(i.eq.j)then
                this%transform(i, j) = this%key_value("PC"//COOP_STR_OF(i)//"_"//COOP_STR_OF(j), 1.d0)
             else
                this%transform(i, j) = this%key_value("PC"//COOP_STR_OF(i)//"_"//COOP_STR_OF(j), 0.d0)
             endif
          enddo
          delta(i) = this%key_value("CDELT"//COOP_STR_OF(i)) * units(i)
          if(abs(delta(i)) .lt. 1.d-12)then
             write(*,*) trim(this%filename)//": cannot read delta"
             call this%header%print()
             stop
          endif
          this%center(i) = this%key_value("CRPIX"//COOP_STR_OF(i), 1.d0) 
          this%icenter(i) = nint(this%center(i)) - 1  !!C convention 
          this%radec_center(i) = this%key_value("CRVAL"//COOP_STR_OF(i)) * units(i)
       enddo
       if(abs(abs(delta(1)/delta(2)) - 1.d0) .gt. 1.d-3)then
          print*, delta(1), delta(2)
          stop "pixel size is not the same in x and y directions"
       else
          this%dx = delta(1)
          this%dy = delta(2)
          this%pixsize = sqrt(abs(delta(1)*delta(2)))
       endif
       this%dkx = coop_2pi/(this%dx*this%nside(1))
       this%dky = coop_2pi/(this%dy*this%nside(2))
       do i=1, this%dim
          this%transform(:, i) =  this%transform(:, i) * delta(i)
       enddo
       deallocate(delta, units)
       this%invtrans = this%transform
       call coop_matrix_inverse(this%invtrans)
    class default
       return
    end select
  end subroutine coop_fits_get_header

  function coop_fits_key_value(this, key, default_value) result(val)
    class(coop_fits)::this
    COOP_UNKNOWN_STRING::key
    COOP_SHORT_STRING::str
    COOP_REAL val
    COOP_REAL,optional::default_value
    str = this%header%value(key)
    if(trim(str).ne."")then
       val = coop_str2real(str)
    else
       if(present(default_value))then
          val = default_value
       else
          val = 0.d0
       endif
    endif
  end function coop_fits_key_value


  subroutine coop_fits_image_simple_stat(this,mean,rms, median, fmin, fmax, lower, upper, mask, clsfile, threshold)
    class(coop_fits_image)::this
    COOP_REAL,optional::threshold
    COOP_UNKNOWN_STRING,optional::clsfile
    class(coop_fits_image),optional::mask
    COOP_REAL,optional::mean, rms, lower(3), upper(3), median, fmin, fmax
    COOP_REAL::ave, std, bounds(-3:3), minv, maxv, weight, zerorat, th
    COOP_INT, dimension(:),allocatable::ind_used
    COOP_REAL,dimension(:),allocatable::image_used
    COOP_INT,parameter::n= 80
    COOP_REAL::Cls(n), ells(n)
    COOP_INT::i, nused, j
    type(coop_file)::fp
    if(present(mask))then
       th = coop_flatsky_mask_threshold(mask)
       if(present(threshold))threshold = th
       nused = count(mask%image .gt. th)
       weight = nused
       write(*,*) "mask fraction: "//trim(coop_num2str(weight/mask%npix*100.))//"%"
       ave = sum(this%image, mask=mask%image .gt. th)/weight
       std =  sqrt(sum((this%image-ave)**2, mask=mask%image .gt. th)/weight)
       allocate(ind_used(nused), image_used(nused))
       j = 0
       do i=0, mask%npix-1
          if(mask%image(i) .gt. th)then
             j = j + 1
             ind_used(j) = i
          endif
       enddo
       image_used = this%Image(ind_used)
       call array_get_mult_threshold_double(image_used, nused, COOP_STANDARD_SIGMA_BOUNDS, 7, bounds)
       minv = minval(this%image, mask = mask%image .gt. th)
       maxv = maxval(this%image, mask = mask%image .gt. th)
       zerorat = count(this%image .eq. 0.d0 .and. mask%image.gt.th)/weight
    else
       ave = sum(this%image)/this%npix
       std =  sqrt(sum((this%image-ave)**2/this%npix))
       weight = this%npix
       call array_get_mult_threshold_double(this%image, int(this%npix), COOP_STANDARD_SIGMA_BOUNDS, 7, bounds)
       minv = minval(this%image)
       maxv = maxval(this%image)
       zerorat = count(this%image .eq. 0.d0)/weight
    endif
    if(present(mean))mean =ave
    if(present(rms))rms = std
    write(*,*)
    write(*,*) "=========================================================="
    write(*,*) "Statistics of map :"//trim(this%filename)
    write(*,*) "=========================================================="
    write(*,*) "size: "//COOP_STR_OF(this%nside(1))//" x "//COOP_STR_OF(this%nside(2))
    write(*,*) "mean = ", ave
    write(*,*) "rms = ", std
    write(*,*) "median = ", bounds(0)
    write(*,*) "1sigma lower, upper = ", bounds(-1), bounds(1)
    write(*,*) "2sigma lower, upper = ", bounds(-2), bounds(2)
    write(*,*) "3sigma lower, upper = ", bounds(-3), bounds(3)
    if(present(median))median = bounds(0)
    if(present(lower))lower = bounds(-1:-3:-1)
    if(present(upper))upper = bounds(1:3)
    write(*,*) "min max = ", minv, maxv 
    if(present(fmin))fmin = minv
    if(present(fmax))fmax = maxv
    write(*,*) "zero-value pixels: "//trim(coop_num2str(100.*zerorat,"(F10.3)"))//"%"
    
    if(present(clsfile))then
       if(trim(clsfile).ne."")then
          call fp%open(clsfile)
          select type(this)
             class is(coop_fits_image_cea)
             weight = weight/this%npix
             call coop_set_uniform(n, ells, max(abs(this%dkx), abs(this%dky))*3.d0, min(abs(this%dkx)*this%nside(1), abs(this%dky)*this%nside(2), 3.d4)/10.d0)
             call coop_fits_image_cea_get_power(this, ells, Cls)

             write(fp%unit,"(A)") "# ell    D_l "
             do i=1, n
                write(fp%unit,"(I7, E16.7)") nint(ells(i)), ells(i)**2*Cls(i)/coop_2pi/sqrt(weight) !!taken into account fsky correction
             enddo
          end select
          call fp%close()
       endif
    endif
    write(*,*) "=========================================================="
    write(*,*)
  end subroutine coop_fits_image_simple_stat

  subroutine coop_fits_image_get_linear_coordinates(this, pix, coor)
    class(coop_fits_image)::this
    COOP_INT pix(this%dim)
    COOP_REAL coor(this%dim)
    coor = matmul(this%transform, pix - this%center)
  end subroutine coop_fits_image_get_linear_coordinates


  subroutine coop_fits_image_get_data(this)
    class(coop_fits_image)::this
    type(coop_fits_file)::fp
    call fp%open_image(this%filename)
    call fp%load_image(this%image)
    call fp%close()
  end subroutine coop_fits_image_get_data

  subroutine coop_fits_image_regularize(this, tail)
    class(coop_fits_image) this
    COOP_REAL upper, lower, tail, diff, diff3, bounds(2)
    if(tail .le. 0.d0)return
    call array_get_mult_threshold_double(this%image, this%npix, (/ 1.-tail, tail /), 2, bounds)
    lower = bounds(1)
    upper = bounds(2)
    diff = (upper - lower)/5.d0
    diff3 = diff*3.d0
    where(this%image .lt. lower)
       this%image = lower - diff3*( log(1.d0+log(1.d0+(lower - this%image)/diff)) - 2.d0*log(1.d0+log(1.d0+(lower - this%image)/diff3)) )
    end where
    where(this%image .gt. upper)
       this%image = upper + diff3*(log(1.d0+ log(1.d0+(this%image - upper)/diff))- 2.d0*log(1.d0+ log(1.d0+(this%image - upper)/diff3)))
    end where
  end subroutine coop_fits_image_regularize


!!=======   CEA ============
  subroutine coop_fits_image_cea_radec2pix(this, ra_deg, dec_deg, pix)
    class(coop_fits_image_cea)::this
    COOP_REAL ra_deg, dec_deg, vec(2), ra, dec
    COOP_INT pix
    COOP_INT::ix, iy, ivec(2)
    ra = ra_deg*coop_SI_degree
    dec = dec_deg*coop_SI_degree
    do while(ra .gt. coop_pi)
       ra = ra - coop_2pi
    enddo
    do while(ra .lt. -coop_pi)
       ra = ra + coop_2pi
    enddo
    vec(1) = ra - this%radec_center(1)
    vec(2) = dec - this%radec_center(2)
    vec = matmul(this%invtrans, vec)
    ivec = nint(vec)
    pix = ivec(1) + this%icenter(1) + this%nside(1) * (ivec(2) + this%icenter(2))
    if(pix .ge. this%npix .or. pix .lt. 0)then
       write(*,*) "warning: radec2pix overflow", pix, this%npix
    endif
  end subroutine coop_fits_image_cea_radec2pix

  subroutine coop_fits_image_cea_pix2ixiy(this, pix, ix, iy)
    class(coop_fits_image_cea)::this
    COOP_INT pix
    COOP_INT ix, iy
    iy = pix/this%nside(1) 
    ix = pix  - iy*this%nside(1) - this%icenter(1)
    iy = iy - this%icenter(2)
  end subroutine coop_fits_image_cea_pix2ixiy

  subroutine coop_fits_image_cea_ixiy2pix(this, ix, iy, pix)
    class(coop_fits_image_cea)::this
    COOP_INT pix
    COOP_INT ix, iy
    pix = ix + this%icenter(1) + (iy +  this%icenter(2) ) * this%nside(1) 
  end subroutine coop_fits_image_cea_ixiy2pix
  

  
  subroutine coop_fits_image_cea_pix2ang(this, pix, theta, phi)
    class(coop_fits_image_cea)::this
    COOP_INT pix
    COOP_REAL theta, phi
    COOP_INT ix, iy
    call this%pix2ixiy(pix, ix, iy)
    phi = this%transform(1,1)* ix + this%transform(1,2)*iy + this%radec_center(1)
    theta = coop_pio2 - asin(this%transform(2, 1)*ix + this%transform(2,2)*iy + this%radec_center(2))
  end subroutine coop_fits_image_cea_pix2ang


  subroutine coop_fits_image_cea_pix2flat(this, pix, coor, disc)
    class(coop_fits_image_cea)::this
    COOP_INT:: pix
    COOP_REAL ::coor(2)
    type(coop_sphere_disc), optional::disc
    call this%pix2ang(pix, coor(1), coor(2))
    if(present(disc))then
       call disc%ang2flat(coor)
    else
       call this%disc%ang2flat(coor)
    endif
  end subroutine coop_fits_image_cea_pix2flat

  subroutine coop_fits_image_cea_cut(this, ix, iy, map)
    class(coop_fits_image_cea)::this
    COOP_REAL map(:,:)
    COOP_INT i, j, ix, iy
    do j= 1, size(map, 2) 
       do i= 1, size(map, 1)
          map(i, j) = this%image(ix + i-1 + (iy + j-1) * this%nside(1))
       enddo
    enddo
  end subroutine coop_fits_image_cea_cut

  subroutine coop_fits_image_cea_get_flatmap(this, smooth_pixsize)
    class(coop_fits_image_cea)::this
    COOP_REAL vec(3), vecmean(3), theta, phi, coor(2), coormin(2), coormax(2)
    COOP_REAL smooth_pixsize, weight
    COOP_REAL,dimension(:,:),allocatable::w
    COOP_INT i, j, icoor(2), k
    COOP_INT pix
    this%smooth_pixsize = smooth_pixsize
    vecmean = 0.
    do pix = 0, this%npix - 1
       call this%pix2ang(pix, theta, phi)
       call coop_sphere_ang2vec(theta, phi, vec)
       vecmean = vecmean + vec
    enddo
    vecmean = vecmean/this%npix
    call coop_sphere_vec2ang(vecmean, theta, phi)
    call this%disc%init(theta, phi)
    coormin = 0.
    coormax = 0.
    do i = 0, this%nside(1)-1, this%nside(1)-1
       do j=0, this%nside(2)-1
          pix = i + j*this%nside(1)
          call this%pix2flat(pix, coor)
          coormin = min(coor, coormin)
          coormax = max(coor, coormax)
       enddo
    enddo
    do i = 0, this%nside(1)-1
       do j=0, this%nside(2)-1, this%nside(2)-1
          pix = i + j*this%nside(1)
          call this%pix2flat(pix, coor)
          coormin = min(coor, coormin)
          coormax = max(coor, coormax)
       enddo
    enddo
    this%smooth_nx = nint(max(abs(coormin(1)), abs(coormax(1)))/smooth_pixsize)+2
    this%smooth_ny = nint(max(abs(coormin(2)), abs(coormax(2)))/smooth_pixsize)+2
    COOP_DEALLOC(this%smooth_image)
    allocate(this%smooth_image(-this%smooth_nx:this%smooth_nx, -this%smooth_ny:this%smooth_ny))
    allocate(w(-this%smooth_nx:this%smooth_nx, -this%smooth_ny:this%smooth_ny))
    w = 0.
    this%smooth_image = 0.
    do pix=0, this%npix - 1
       call this%pix2flat(pix, coor)
       coor = coor/this%smooth_pixsize
       icoor = nint(coor)
       do k = -2, 2
          j = icoor(2) + k
          do i=icoor(1)-2 + abs(k), icoor(1) + 2 - abs(k)
             weight = exp(-(coor(1)-i)**2 - (coor(2)-j)**2)
             w(i,j) = w(i,j) + weight
             this%smooth_image(i,j) = this%smooth_image(i, j) + weight * this%image(pix)
          enddo
       enddo
    enddo
    where(w .gt. 0.)
       w = this%smooth_image /w
    end where
    this%smooth_nx = this%smooth_nx - 2
    this%smooth_ny = this%smooth_ny - 2
    this%smooth_dkx = coop_2pi/(this%smooth_pixsize*(this%smooth_nx*2+1))
    this%smooth_dky = coop_2pi/(this%smooth_pixsize*(this%smooth_ny*2+1))
    deallocate(this%smooth_image)
    allocate(this%smooth_image(-this%smooth_nx:this%smooth_nx, -this%smooth_ny:this%smooth_ny))
    this%smooth_image = w(-this%smooth_nx:this%smooth_nx, -this%smooth_ny:this%smooth_ny)
    deallocate(w)
    this%smooth_npix = ( this%smooth_nx * 2 + 1 ) * ( this%smooth_ny * 2 + 1)
    this%xmax = this%smooth_nx * this%smooth_pixsize
    this%ymax = this%smooth_ny * this%smooth_pixsize
    this%xmin = -this%xmax
    this%ymin = -this%ymax
  end subroutine coop_fits_image_cea_get_flatmap

  subroutine coop_fits_image_cea_smooth(this, fwhm, highpass_l1, highpass_l2, lx_cut, ly_cut, lmax, beam)
    class(coop_fits_image_cea)::this
    COOP_REAL::fwhm
    COOP_INT::highpass_l1, highpass_l2, lmax
    COOP_REAL,optional::beam(0:lmax)
    COOP_REAL::hp_omega, sigma2, hp_l2sq
    COOP_INT,optional::lx_cut, ly_cut
    logical::do_lx_cut, do_ly_cut, has_beam
    do_lx_cut = present(lx_cut)
    do_ly_cut = present(ly_cut)
    if(do_lx_cut) do_lx_cut = lx_cut .gt. 0
    if(do_ly_cut) do_ly_cut = ly_cut .gt. 0
    has_beam = present(beam)
    sigma2 = (coop_sigma_by_fwhm*fwhm/coop_sqrt2)**2
    hp_l2sq = highpass_l2**2
    hp_omega = coop_pio2/(highpass_l2 - highpass_l1)
    call this%filter(lmin = highpass_l1, lmax = lmax,  window = smooth_weight)
  contains

    function smooth_weight(kx, ky)
      COOP_REAL kx,ky,k2, k, smooth_weight
      k2 = (kx**2+ky**2)
      k = sqrt(k2)
      smooth_weight = exp(-k2*sigma2)
      if(has_beam) then
         smooth_weight =smooth_weight/max(beam(nint(k)), 1.d-2) !!when beam<0.01 actually no info
      endif
      if(do_lx_cut)then
         if(abs(kx).le.lx_cut)then
            smooth_weight = 0.d0
            return
         endif
         smooth_weight = smooth_weight * cos(coop_pio2*(lx_cut/kx)**4)**2
      endif
      if(do_ly_cut)then
         if(abs(ky).le. ly_cut)then
            smooth_weight = 0.d0
            return
         endif
         smooth_weight = smooth_weight * cos(coop_pio2*(ly_cut/ky)**4)**2
      endif
      if(k2 .lt. hp_l2sq)then
         smooth_weight = smooth_weight * sin((k-highpass_l1)*hp_omega)**2
      endif
    end function smooth_weight
  end subroutine coop_fits_image_cea_smooth

  subroutine coop_fits_image_cea_get_power(this, ells, Cls)
    class(coop_fits_image_cea)::this
    COOP_REAL,dimension(:)::ells
    COOP_REAL,dimension(:)::Cls
    COOP_REAL::weights(size(ells))
    COOP_REAL k2min, k2max, k2, ky2, k
    COOP_INT::n, left, right, mid, i, j
    COOP_COMPLEX,dimension(:,:),allocatable::fk
    n = size(ells)
    if(n.ne.size(Cls)) stop "get_power: size of ells must match size of Cls"
    allocate( fk(0:this%nside(1)/2, 0:this%nside(2)-1))
    call coop_fft_forward(this%nside(1), this%nside(2), this%image, fk)
    k2min = (ells(1))**2*1.000001
    k2max = (ells(n))**2*0.999999
    Cls = 0.d0
    weights = 0.d0
    do j=0, this%nside(2)-1
       ky2 = (min(j, this%nside(2)-j)*this%dky)**2
       if( ky2 .ge. k2max)cycle
       left = 1
       do i = 0, this%nside(1)/2
          k2 = (i*this%dkx)**2 + ky2
          if(k2.ge.k2max) exit
          if(k2 .le. k2min) cycle
          k = sqrt(k2)
          right = n
          do while(right - left .gt. 1)
             mid = (left+right)/2
             if( k .gt. ells(mid))then
                left = mid
             else
                right = mid
             endif
          enddo
          if(left .ge. n .or. left .lt. 1 .or. right-left.ne.1 .or. k.gt. ells(right) .or. k.lt. ells(left)) then
             write(*,*) "error in get_power"
             stop
          endif
          Cls(left) = Cls(left) + (ells(right)-k)/(ells(right)-ells(left))*(abs(fk(i, j))**2)
          weights(left) = weights(left)+ (ells(right)-k)/(ells(right)-ells(left))
          Cls(right) = Cls(right) + (k-ells(left))/(ells(right)-ells(left))*(abs(fk(i, j))**2)
          weights(right) = weights(right)+ (k-ells(left))/(ells(right)-ells(left))

       enddo
    enddo
    deallocate(fk)
    where(weights .gt. 0.d0)
       Cls = Cls/weights
    end where
    Cls = Cls*(this%pixsize**2)/this%npix
  end subroutine coop_fits_image_cea_get_power

  subroutine coop_fits_image_cea_filter(this, lmin, lmax, window)
    class(coop_fits_image_cea)::this
    COOP_INT lmin, lmax, i, j, lx_cut, ly_cut
    external window
    COOP_REAL window
    COOP_REAL k2min, k2max, k2, ky2, ky, kx
    COOP_COMPLEX,dimension(:,:),allocatable::fk
    allocate( fk(0:this%nside(1)/2, 0:this%nside(2)-1))
    call coop_fft_forward(this%nside(1), this%nside(2), this%image, fk)
    k2min = dble(lmin)**2
    k2max = dble(lmax)**2
    do j=0, this%nside(2)-1
       if(j.gt. this%nside(2)-j)then
          ky = (j - this%nside(2))*this%dky
       else
          ky = j*this%dky
       endif
       ky2 = ky**2
       if( ky2 .gt. k2max)then
          fk(:, j) = ( 0.d0, 0.d0 )
          cycle
       endif
       do i = 0, this%nside(1)/2
          kx = i*this%dkx
          k2 = kx**2 + ky2
          if(k2 .lt. k2min .or. k2.gt. k2max)then
             fk(i, j) = ( 0.d0, 0.d0 )
          else
             fk(i, j) = fk(i, j)*window(kx, ky)
          endif
       enddo
    enddo
    call coop_fft_backward(this%nside(1), this%nside(2), fk, this%image)
    deallocate(fk)
  end subroutine coop_fits_image_cea_filter


  subroutine coop_fits_image_cea_smooth_flat(this, lmin, lmax, fk_file)
    class(coop_fits_image_cea)::this
    COOP_UNKNOWN_STRING, optional:: fk_file
    COOP_INT lmin, lmax, i, j, ik
    COOP_REAL k2min, k2max, k2, rn1, rn2, kmin, kmax, omega, j2, kstep, rk, rms, cl, ell
    logical do_qu
    COOP_COMPLEX,dimension(:,:),allocatable::fk, fke, fkb
    type(coop_file)::fkf
    COOP_INT nk
    COOP_REAL,dimension(:),allocatable::karr, Pkarr, warr
    if(allocated(this%smooth_Q) .and. allocated(this%smooth_U))then
       call this%QU2EB()
       do_qu = .true.
    else
       do_qu = .false.
    endif
    allocate( fk(0:this%smooth_nx, 0:this%smooth_ny*2))
    if(do_qu)allocate( fke(0:this%smooth_nx, 0:this%smooth_ny*2), fkb(0:this%smooth_nx, 0:this%smooth_ny*2))
    call coop_fft_forward(this%smooth_nx*2+1, this%smooth_ny*2+1, this%smooth_image, fk)
    if(do_qu)then
       call coop_fft_forward(this%smooth_nx*2+1, this%smooth_ny*2+1, this%smooth_Q, fke)
       call coop_fft_forward(this%smooth_nx*2+1, this%smooth_ny*2+1, this%smooth_U, fkb)
    endif
    rn1 = this%smooth_nx*2+1.d0
    rn2  =this%smooth_ny*2+1.d0
    kmin = (this%smooth_pixsize * lmin/coop_2pi)
    kmax = (this%smooth_pixsize * lmax/coop_2pi)
    omega = coop_pi/(kmax - kmin)
    kmax = min(kmax, coop_sqrt2/2.d0)
    k2min = kmin**2
    k2max = kmax**2
    if(present(fk_file))then
       nk = (kmax-kmin)/this%smooth_pixsize
       kstep = (kmax - kmin)/nk
       allocate(karr(0:nk+1), Pkarr(0:nk+1), warr(0:nk+1))
       call coop_set_uniform(nk+2, karr, kmin-kstep/2.d0, kmax+kstep/2.d0)
       pkarr = 0.d0
       warr = 0.d0
    endif
    do j=0, this%smooth_ny*2
       j2 = (min(j, this%smooth_ny*2+1-j)/rn2)**2
       if( j2 .gt. k2max)then
          fk(0:this%smooth_nx, j) = ( 0.d0, 0.d0 )
          cycle
       endif
       do i = 0, this%smooth_nx
          k2 = j2 + (i/rn1)**2
          if(k2 .lt. k2min .or. k2.gt. k2max)then
             fk(i, j) = ( 0.d0, 0.d0 )
          else
             fk(i, j) = fk(i, j)*sin((sqrt(k2)- kmin)*omega)
             if(do_qu)then
                fke(i,j) = fke(i,j)*sin((sqrt(k2)- kmin)*omega)
                fkb(i,j) = fke(i,j)*sin((sqrt(k2)- kmin)*omega)
             endif
             if(present(fk_file))then
                rk = (sqrt(k2)-kmin)/kstep+0.5d0
                ik = floor(rk)
                rk = rk - ik
                pkarr(ik) = pkarr(ik) + (aimag(fk(i, j))**2 + real(fk(i,j))**2)*(1.-rk)
                warr(ik) = warr(ik) + 1.-rk
                pkarr(ik+1) = pkarr(ik+1) + ( aimag(fk(i, j))**2 + real(fk(i,j))**2 ) * rk
                warr(ik+1) = warr(ik+1) + rk

             endif
          endif
       enddo
    enddo
    if(present(fk_file))then
       where(warr .ne. 0.d0)
          pkarr = pkarr/warr
       end where
       call fkf%open(trim(fk_file))
       rms   = 0
       do i=1, nk
          ell  = sqrt((karr(i)*coop_2pi/this%smooth_pixsize)**2 + 0.25) - 0.5
          cl = pkarr(i)*this%smooth_pixsize**2/this%smooth_npix
          write(fkf%unit, "(2E16.7)")  ell, cl
          rms = rms + cl * (2*ell+1)/coop_4pi* (kstep*coop_2pi/this%smooth_pixsize)
       enddo
       call fkf%close()

    endif
    call coop_fft_backward(this%smooth_nx*2+1, this%smooth_ny*2+1, fk, this%smooth_image)
    if(do_qu)then
       call coop_fft_backward(this%smooth_nx*2+1, this%smooth_ny*2+1, fke, this%smooth_Q)
       call coop_fft_backward(this%smooth_nx*2+1, this%smooth_ny*2+1, fkb, this%smooth_U)
       call this%EB2QU()
       deallocate(fke,fkb)
    endif
    deallocate(fk)
    if(present(fk_file))deallocate(karr, Pkarr, warr)
  end subroutine coop_fits_image_cea_smooth_flat

  subroutine coop_fits_image_cea_find_extrema(this, mask, spot_file, spot_type, radius, nu)
    class(coop_fits_image_cea)::this, mask
    COOP_UNKNOWN_STRING::spot_file, spot_type
    type(coop_file)::cf
    COOP_REAL,optional:: radius, nu
    COOP_REAL angle,rms,mean,weight,cut, nucut
    COOP_INT i, j, irad
    if(present(radius))then
       irad = ceiling(radius/this%smooth_pixsize)
    else
       irad = 1
    endif
    weight = sum(mask%image)
    mean = sum(this%image*mask%image)/weight
    rms = sqrt(sum((this%image-mean)**2*mask%image)/weight)
    if(present(nu))then
       nucut = nu
    else
       nucut = 0.d0
    endif
    call cf%open(trim(spot_file), "w")
    select case(trim(spot_type))
    case("Hot")
       cut = rms*nucut+mean
       do i=-this%smooth_nx+irad, this%smooth_nx - irad
          do j = -this%smooth_ny+irad, this%smooth_ny - irad
             if(  mask%smooth_image(i,j).gt.0.5 .and. this%smooth_image(i, j) .ge. cut)then
                write(cf%unit, "(2I6, E16.7, I6 )") i, j, coop_2pi*coop_random_unit(), min(this%smooth_nx - i, this%smooth_nx + i, this%smooth_ny - j, this%smooth_ny + j)
             endif
          enddo
       enddo
    case("Cold")
       cut = -rms*nucut+mean
       do i=-this%smooth_nx+irad, this%smooth_nx - irad
          do j = -this%smooth_ny+irad, this%smooth_ny - irad
             if(  mask%smooth_image(i,j).gt.0.5 .and. this%smooth_image(i, j) .le. cut)then
                write(cf%unit, "(2I6, E16.7, I6 )") i, j, coop_2pi*coop_random_unit(), min(this%smooth_nx - i, this%smooth_nx + i, this%smooth_ny - j, this%smooth_ny + j)
             endif
          enddo
       enddo
    case("Hot_QTUTOrient")
       if(.not. allocated(this%smooth_Q) .or. .not. allocated(this%smooth_U)) &
            call this%get_QTUT()
       cut = rms*nucut+mean
       do i=-this%smooth_nx+irad, this%smooth_nx - irad
          do j = -this%smooth_ny+irad, this%smooth_ny - irad
             if(  mask%smooth_image(i,j).gt.0.5 .and. this%smooth_image(i, j) .ge. cut)then
                angle = 0.5d0 * COOP_POLAR_ANGLE(this%smooth_Q(i, j), this%smooth_U(i, j))
                write(cf%unit, "(2I6, E16.7, I6 )") i, j, angle, min(this%smooth_nx - i, this%smooth_nx + i, this%smooth_ny - j, this%smooth_ny + j)
             endif
          enddo
       enddo
    case("Cold_QTUTOrient")
       if(.not. allocated(this%smooth_Q) .or. .not. allocated(this%smooth_U)) &
            call this%get_QTUT()
       cut = -rms*nucut+mean
       do i=-this%smooth_nx+irad, this%smooth_nx - irad
          do j = -this%smooth_ny+irad, this%smooth_ny - irad
             if(  mask%smooth_image(i,j).gt.0.5 .and. this%smooth_image(i, j) .le. cut)then
                angle = 0.5d0 * COOP_POLAR_ANGLE(this%smooth_Q(i, j), this%smooth_U(i, j))
                write(cf%unit, "(2I6, E16.7, I6 )") i, j, angle, min(this%smooth_nx - i, this%smooth_nx + i, this%smooth_ny - j, this%smooth_ny + j)
             endif
          enddo
       enddo
    case("Tmax", "Emax", "Bmax")
       cut = mean+rms*nucut
       do i=-this%smooth_nx+irad, this%smooth_nx - irad
          do j = -this%smooth_ny+irad, this%smooth_ny - irad
             if(all(this%smooth_image(i, j) .ge. this%smooth_image(i-1:i+1, j-1:j+1)) .and. mask%smooth_image(i,j).gt.0.5 .and.this%smooth_image(i, j).gt.cut)then
                write(cf%unit, "(2I6, E16.7, I6 )") i, j, coop_2pi*coop_random_unit(), min(this%smooth_nx - i, this%smooth_nx + i, this%smooth_ny - j, this%smooth_ny + j)
             endif
          enddo
       enddo
    case("Tmin", "Emin", "Bmin")
       cut=mean-rms*nucut
       do i=-this%smooth_nx+irad, this%smooth_nx - irad
          do j = -this%smooth_ny+irad, this%smooth_ny - irad
             if(all(this%smooth_image(i, j) .le. this%smooth_image(i-1:i+1, j-1:j+1)) .and. mask%smooth_image(i,j).gt.0.5 .and.this%smooth_image(i, j).lt.cut)then
                write(cf%unit, "(2I6, E16.7, I6 )") i, j, coop_2pi*coop_random_unit(), min(this%smooth_nx - i, this%smooth_nx + i, this%smooth_ny - j, this%smooth_ny + j)
             endif
          enddo
       enddo
    case("Tmax_QTUTOrient")
       cut = mean+rms*nucut
       if(.not. allocated(this%smooth_Q) .or. .not. allocated(this%smooth_U)) &
            call this%get_QTUT()
       do i=-this%smooth_nx+irad, this%smooth_nx - irad
          do j = -this%smooth_ny+irad, this%smooth_ny - irad
             if(all(this%smooth_image(i, j) .ge. this%smooth_image(i-1:i+1, j-1:j+1)) .and. mask%smooth_image(i,j).gt.0.5  .and.this%smooth_image(i, j).gt.cut)then
                angle = 0.5d0 * COOP_POLAR_ANGLE(this%smooth_Q(i, j), this%smooth_U(i, j))
                write(cf%unit, "(2I6, E16.7, I6 )") i, j, angle, min(this%smooth_nx - i, this%smooth_nx + i, this%smooth_ny - j, this%smooth_ny + j)
             endif
          enddo
       enddo
    case("Tmin_QTUTOrient")
       cut = mean-rms*nucut
       if(.not. allocated(this%smooth_Q) .or. .not. allocated(this%smooth_U)) &
            call this%get_QTUT()
       do i=-this%smooth_nx+irad, this%smooth_nx - irad
          do j = -this%smooth_ny+irad, this%smooth_ny - irad
             if(all(this%smooth_image(i, j) .le. this%smooth_image(i-1:i+1, j-1:j+1)) .and. mask%smooth_image(i,j).gt.0.5  .and.this%smooth_image(i, j).lt.cut )then
                angle = 0.5d0 * COOP_POLAR_ANGLE(this%smooth_Q(i, j), this%smooth_U(i, j))
                write(cf%unit, "(2I6, E16.7, I6 )") i, j, angle, min(this%smooth_nx - i, this%smooth_nx + i, this%smooth_ny - j, this%smooth_ny + j)
             endif
          enddo
       enddo
    case("PTmax")
       if(.not. allocated(this%smooth_Q) .or. .not. allocated(this%smooth_U)) &
            call this%get_QTUT()
       do i=-this%smooth_nx+irad, this%smooth_nx - irad
          do j = -this%smooth_ny+irad, this%smooth_ny - irad
             if(all(this%smooth_Q(i, j)**2 + this%smooth_U(i,j)**2 .ge. this%smooth_Q(i-1:i+1, j-1:j+1)**2 + this%smooth_U(i-1:i+1, j-1:j+1)**2 )  .and. mask%smooth_image(i,j).gt.0.5 )then
                angle = 0.5d0 * COOP_POLAR_ANGLE(this%smooth_Q(i, j), this%smooth_U(i, j))
                write(cf%unit, "(2I6, E16.7, I6 )") i, j, angle, min(this%smooth_nx - i, this%smooth_nx + i, this%smooth_ny - j, this%smooth_ny + j)
             endif
          enddo
       enddo
    case("PTmin")
       if(.not. allocated(this%smooth_Q) .or. .not. allocated(this%smooth_U)) &
            call this%get_QTUT()
       do i=-this%smooth_nx+irad, this%smooth_nx - irad
          do j = -this%smooth_ny+irad, this%smooth_ny - irad
             if(all(this%smooth_Q(i, j)**2 + this%smooth_U(i,j)**2 .le. this%smooth_Q(i-1:i+1, j-1:j+1)**2 + this%smooth_U(i-1:i+1, j-1:j+1)**2 ) .and. mask%smooth_image(i,j).gt.0.5)then
                angle = 0.5d0 * COOP_POLAR_ANGLE(this%smooth_Q(i, j), this%smooth_U(i, j))
                write(cf%unit, "(2I6, E16.7, I6 )") i, j, angle, min(this%smooth_nx - i, this%smooth_nx + i, this%smooth_ny - j, this%smooth_ny + j)
             endif
          enddo
       enddo
    case  default
       write(*,*) trim(spot_type)//" :Unknown spot type"
       stop
    end select
    call cf%close()
  end subroutine coop_fits_image_cea_find_extrema

  subroutine coop_fits_image_cea_simulate_flat(this, lmin, lmax, Cls_file)
    class(coop_fits_image_cea)::this
    COOP_UNKNOWN_STRING::Cls_file
    type(coop_file)::fp
    COOP_INT lmin, lmax, i, j, ik, l, ii
    COOP_REAL Cls(3, 3, lmin:lmax), sqrteig(3, lmin:lmax), rot(3,3, lmin:lmax), rin(4)
    COOP_COMPLEX,dimension(:,:,:),allocatable::fk
    COOP_REAL amp, rk
    call coop_random_init()
    if(.not. coop_file_exists(Cls_file))then
       write(*,*) trim(Cls_file)//" does not exist"
       stop
    endif
    if(.not. allocated(this%smooth_image))then
       write(*,*) "call get_flatmap before you do simulate_flat"
       stop
    endif
    if( .not. allocated(this%smooth_Q))allocate(this%smooth_Q(-this%smooth_nx:this%smooth_nx, -this%smooth_ny:this%smooth_ny))
    if( .not. allocated(this%smooth_U))allocate(this%smooth_U(-this%smooth_nx:this%smooth_nx, -this%smooth_ny:this%smooth_ny))
    Cls = 0.d0
    call fp%open(Cls_file, "r")
    do
       read(fp%unit, *, ERR=100, END=100) l, rin
       if(l.lt. lmin) cycle
       if(l.le. lmax)then
          Cls(1, 1, l) = rin(1)!!TT
          Cls(2, 2, l) = rin(2)!!EE
          Cls(3, 3, l) = rin(3)!!BB
          Cls(1, 2, l) = rin(4)!!TE
          Cls(2, 1, l) = rin(4)!!TE
          call coop_matsym_diag(3, 3, Cls(:,:,l), rot(:,:,l))
          do i=1, 3
             sqrteig(i, l) = sqrt(Cls(i, i, l))
          enddo
       endif
       if(l .ge. lmax) exit
    enddo

    amp = sqrt(1.d0/(this%smooth_pixsize**2/this%smooth_npix))
    sqrteig = sqrteig*amp
    do l=lmin, lmax
       sqrteig(:, l) = sqrteig(:, l) * sin(coop_pi*dble(l- lmin)/dble(lmax-lmin))
    enddo
    allocate(fk(0:this%smooth_nx,0:this%smooth_ny*2, 3))
    do i=0, this%smooth_nx
       do j=0, this%smooth_ny*2
          rk  = sqrt((this%smooth_dkx*i)**2 + (this%smooth_dky* min(j, 2*this%smooth_ny - j + 1))**2) + 0.5d0
          ik = floor(rk)
          if(ik.ge.lmin .and. ik .lt. lmax)then
             rk = rk - ik
             do ii=1,3
                fk(i,j,ii) = coop_random_complex_gaussian() * (sqrteig(ii, ik)*(1.d0-rk) + sqrteig(ii, ik+1)*rk)
             enddo
             fk(i, j, :) = matmul((rot(:,:,ik)*(1.d0-rk) + rot(:,:, ik+1)*rk),  fk(i, j, :))
          else
             fk(i, j, :) = 0
          endif
       enddo
    enddo
    call coop_fft_backward(this%smooth_nx*2+1, this%smooth_ny*2+1, fk(:,:,1), this%smooth_image)
    call coop_fft_backward(this%smooth_nx*2+1, this%smooth_ny*2+1, fk(:,:,2), this%smooth_Q)
    call coop_fft_backward(this%smooth_nx*2+1, this%smooth_ny*2+1, fk(:,:,3), this%smooth_U)
    deallocate(fk)
    call fp%close()
    call this%EB2QU()
    return
100  stop "Cl file does not contain valid data"
  end subroutine coop_fits_image_cea_simulate_flat

  subroutine coop_fits_image_cea_stack(this, spot_file, stack_option, nrad, stacked_image, nstack)
    class(coop_fits_image_cea)::this
    COOP_UNKNOWN_STRING::spot_file, stack_option
    type(coop_file)::fp
    COOP_INT i, j, nrad, ii, jj, dis2b, irot, jrot, nstack
    COOP_REAL theta, theta_r, radius, stacked_image(-nrad:nrad,-nrad:nrad)
    COOP_REAL::firot, fjrot, pixv
    logical::accept
    stacked_image = 0.d0
    nstack = 0
    call fp%open(trim(spot_file))
    select case(stack_option)
    case("T")
       do 
          read(fp%unit, *, ERR=100, END=100) i, j, theta, dis2b
          if(dis2b .lt. nrad) cycle
          do ii= -nrad, nrad
             do jj= -nrad, nrad 
                call get_rot(ii, jj, theta, irot, jrot, firot, fjrot, nrad, accept)
                if(accept)then
                   pixv = this%smooth_image(i+ii, j+jj)  
                   stacked_image(irot, jrot) = stacked_image(irot, jrot) + pixv*(1.d0 - firot)*(1.d0-fjrot)
                   stacked_image(irot+1, jrot) = stacked_image(irot+1, jrot) + pixv*firot*(1.d0-fjrot)
                   stacked_image(irot, jrot+1) = stacked_image(irot, jrot+1) + pixv*(1.d0 - firot)*fjrot
                   stacked_image(irot+1, jrot+1) = stacked_image(irot+1, jrot+1) + pixv*firot*fjrot                   
                endif
             enddo
          enddo
          nstack = nstack + 1
       enddo
    case("E")
       do 
          read(fp%unit, *, ERR=100, END=100) i, j, theta, dis2b
          if(dis2b .lt. nrad) cycle
          do ii= -nrad, nrad
             do jj= -nrad, nrad
                call get_rot(ii, jj, theta, irot, jrot, firot, fjrot, nrad, accept)
                if(accept)then
                   pixv = this%smooth_Q(i+ii, j+jj)  
                   stacked_image(irot, jrot) = stacked_image(irot, jrot) + pixv*(1.d0 - firot)*(1.d0-fjrot)
                   stacked_image(irot+1, jrot) = stacked_image(irot+1, jrot) + pixv*firot*(1.d0-fjrot)
                   stacked_image(irot, jrot+1) = stacked_image(irot, jrot+1) + pixv*(1.d0 - firot)*fjrot
                   stacked_image(irot+1, jrot+1) = stacked_image(irot+1, jrot+1) + pixv*firot*fjrot                   
                endif
             enddo
          enddo
          nstack = nstack + 1
       enddo
    case("B")
       do 
          read(fp%unit, *, ERR=100, END=100) i, j, theta, dis2b
          if(dis2b .lt. nrad) cycle
          do ii= -nrad, nrad
             do jj= -nrad, nrad
                call get_rot(ii, jj, theta, irot, jrot, firot, fjrot, nrad, accept)
                if(accept)then
                   pixv = this%smooth_U(i+ii, j+jj)  
                   stacked_image(irot, jrot) = stacked_image(irot, jrot) + pixv*(1.d0 - firot)*(1.d0-fjrot)
                   stacked_image(irot+1, jrot) = stacked_image(irot+1, jrot) + pixv*firot*(1.d0-fjrot)
                   stacked_image(irot, jrot+1) = stacked_image(irot, jrot+1) + pixv*(1.d0 - firot)*fjrot
                   stacked_image(irot+1, jrot+1) = stacked_image(irot+1, jrot+1) + pixv*firot*fjrot                   
                endif
             enddo
          enddo
          nstack = nstack + 1
       enddo
    case("Qr")
       do 
          read(fp%unit, *, ERR=100, END=100) i, j, theta, dis2b
          if(dis2b .lt. nrad) cycle
          do ii= -nrad, nrad
             do jj= -nrad, nrad 
                call get_rot(ii, jj, theta, irot, jrot, firot, fjrot, nrad, accept)
                if(accept)then
                   theta_r = COOP_POLAR_ANGLE(dble(ii), dble(jj))       
                   pixv =  (this%smooth_Q(i+ii, j+jj) * cos(2.*theta_r) + this%smooth_U(i+ii, j+jj)*sin(2.*theta_r))*coop_healpix_QrUrSign                   
                   stacked_image(irot, jrot) = stacked_image(irot, jrot) + pixv*(1.d0 - firot)*(1.d0-fjrot)
                   stacked_image(irot+1, jrot) = stacked_image(irot+1, jrot) + pixv*firot*(1.d0-fjrot)
                   stacked_image(irot, jrot+1) = stacked_image(irot, jrot+1) + pixv*(1.d0 - firot)*fjrot
                   stacked_image(irot+1, jrot+1) = stacked_image(irot+1, jrot+1) + pixv*firot*fjrot                   
                endif
             enddo
          enddo
          nstack = nstack + 1
       enddo
    case("Q")
       do 
          read(fp%unit, *, ERR=100, END=100) i, j, theta, dis2b
          if(dis2b .lt. nrad) cycle          
          do ii= -nrad, nrad
             do jj= -nrad, nrad
                call get_rot(ii, jj, theta, irot, jrot, firot, fjrot, nrad, accept)
                if(accept)then
                   pixv = (this%smooth_Q(i+ii, j+jj) * cos(2.*theta) + this%smooth_U(i+ii, j+jj)*sin(2.*theta)) 
                   stacked_image(irot, jrot) = stacked_image(irot, jrot) + pixv*(1.d0 - firot)*(1.d0-fjrot)
                   stacked_image(irot+1, jrot) = stacked_image(irot+1, jrot) + pixv*firot*(1.d0-fjrot)
                   stacked_image(irot, jrot+1) = stacked_image(irot, jrot+1) + pixv*(1.d0 - firot)*fjrot
                   stacked_image(irot+1, jrot+1) = stacked_image(irot+1, jrot+1) + pixv*firot*fjrot                   
                endif
             enddo
          enddo
          nstack = nstack + 1
       enddo
    case default
       write(*,*) trim(stack_option)//": unknown stack option"
       stop
    end select
100 call fp%close()
    if(nstack .gt. 0) &
         stacked_image = stacked_image / nstack

  contains

    subroutine get_rot(i, j, theta, irot, jrot, firot, fjrot, n, accept)
      COOP_INT::i, j, irot, jrot, n
      COOP_REAL::theta, firot, fjrot, maxn, cost, sint
      logical accept
      cost = cos(theta)
      sint = sin(theta)
      firot = (ii*cost + jj*sint)
      fjrot = (-ii*sint + jj*cost)
      irot = floor(firot)
      jrot = floor(fjrot)
      firot = firot - irot
      fjrot = fjrot - jrot
      maxn = max(abs(irot), abs(irot+1), abs(jrot), abs(jrot+1))
      if(maxn .le. n)then
         accept = .true.
      elseif(maxn.eq.n+1)then  !!maxn = n + 1
         accept = .true.
         if(abs(irot).eq.n+1)then  !!irot = -n-1
            irot = -n
            firot = 0.d0
         elseif(abs(irot+1) .eq. n+1)then
            irot = n-1
            firot = 1.d0
         endif
         if(abs(jrot).eq. n+1)then
            jrot = -n
            fjrot = 0.d0
         elseif(abs(jrot+1).eq.n+1)then
            jrot = n-1
            fjrot = 1.d0
         endif
      else
         accept = .false.
         return
      endif
      
    end subroutine get_rot
    
  end subroutine coop_fits_image_cea_stack


  subroutine coop_fits_image_cea_stack2fig(this, spot_file, stack_option, radius, fig, caption, label, color_table)
    class(coop_fits_image_cea)::this
    COOP_REAL radius
    COOP_UNKNOWN_STRING::spot_file, fig, stack_option
    COOP_UNKNOWN_STRING,optional::caption, label, color_table
    COOP_REAL,dimension(:,:),allocatable::image
    type(coop_asy)::asy
    COOP_INT nstack
    COOP_INT nrad
    if(.not. coop_file_exists(spot_file))then
       write(*,*) "spot file "//trim(spot_file)//" does not exist"
       stop
    endif
    nrad = max(1, floor(radius/this%smooth_pixsize))
    allocate(image(-nrad:nrad, -nrad:nrad))
    call this%stack(spot_file, stack_option, nrad, image, nstack)
    call asy%open(trim(fig))
    if(present(caption))then
       call asy%init(xlabel = "$2\sin{\frac{\theta}{2}} \cos\varphi$", ylabel = "$2\sin{\frac{\theta}{2}} \sin\varphi$", width=7., height=5.5, caption=trim(caption)//"; stacked "//trim(coop_num2str(nstack))//" patches",  nxticks = 5, nyticks = 5)
    else
       call asy%init(xlabel = "$2\sin{\frac{\theta}{2}} \cos\varphi$", ylabel = "$2\sin{\frac{\theta}{2}} \sin\varphi$", width=7., height=5.5, caption="stacked "//trim(coop_num2str(nstack))//" patches", nxticks = 5, nyticks = 5)
    endif
    if(present(color_table))then
       if(present(label))then
          call coop_asy_density(asy, image, -this%smooth_pixsize*nrad, this%smooth_pixsize*nrad, -this%smooth_pixsize*nrad, this%smooth_pixsize*nrad, label = trim(label), color_table=trim(color_table))
       else
          call coop_asy_density(asy, image, -this%smooth_pixsize*nrad, this%smooth_pixsize*nrad, -this%smooth_pixsize*nrad, this%smooth_pixsize*nrad, color_table = trim(color_table))
       endif
    else
       if(present(label))then
          call coop_asy_density(asy, image, -this%smooth_pixsize*nrad, this%smooth_pixsize*nrad, -this%smooth_pixsize*nrad, this%smooth_pixsize*nrad, label = trim(label))
       else
          call coop_asy_density(asy, image, -this%smooth_pixsize*nrad, this%smooth_pixsize*nrad, -this%smooth_pixsize*nrad, this%smooth_pixsize*nrad)
       endif
    endif
    call asy%close()
    deallocate(image)
  end subroutine coop_fits_image_cea_stack2fig

  subroutine coop_fits_image_cea_get_QU(this, Qmap, Umap)
    class(coop_fits_image_cea)::this, Qmap, Umap
    if(this%smooth_nx .ne. Qmap%smooth_nx .or. this%smooth_nx .ne. Umap%smooth_nx .or. this%smooth_ny .ne. Qmap%smooth_ny .or. this%smooth_ny .ne. Umap%smooth_ny .or. abs(this%smooth_pixsize - Qmap%smooth_pixsize).gt. 1.d-6 .or. abs(this%smooth_pixsize - Umap%smooth_pixsize) .gt. 1.d-6)then
       write(*,*) "Cannot merge I, Q, U maps with different configurations."
       stop
    endif
    COOP_DEALLOC(this%smooth_Q)
    COOP_DEALLOC(this%smooth_U)
    allocate( &
         this%smooth_Q(-this%smooth_nx:this%smooth_nx, -this%smooth_ny:this%smooth_ny), &
         this%smooth_U(-this%smooth_nx:this%smooth_nx, -this%smooth_ny:this%smooth_ny) )
    this%smooth_Q = Qmap%smooth_image
    this%smooth_U = Umap%smooth_image
  end subroutine coop_fits_image_cea_get_QU

  subroutine coop_fits_image_cea_get_QTUT(this, qt, ut)
    class(coop_fits_image_cea)::this
    type(coop_fits_image_cea),optional::qt, ut
    COOP_COMPLEX,dimension(:,:,:),allocatable::fk
    COOP_INT::i, j
    COOP_REAL::ky, ky2, kx, kx2, k2, k2min
    if(present(qt) .and. present(ut))then
       allocate(fk(0:this%nside(1)/2,0:this%nside(2)-1, 2))
       call coop_fft_forward(this%nside(1), this%nside(2), this%image, fk(:,:,1))
       k2min = min(this%dkx, this%dky)
       do j=0, this%nside(2)-1
          if(j .gt. this%nside(2)- j)then
             ky = (j-this%nside(2))*this%dky
          else
             ky = j*this%dky
          endif
          ky2 = ky**2
          do i=0, this%nside(1)/2
             kx = this%dkx*i
             kx2 = kx**2
             k2 = kx2+ky2
             if(k2 .lt. k2min)then
                fk(i,j,:) = (0.d0, 0.d0)
             else
                fk(i, j,2) = -2.d0*kx*ky/k2*fk(i,j,1) !!healpix convention
                fk(i, j,1) = (ky2 - kx2)/k2*fk(i,j,1)
             endif
          enddo
       enddo
       call coop_fft_backward(this%nside(1), this%nside(2), fk(:,:,1), qt%image)
       call coop_fft_backward(this%nside(1), this%nside(2), fk(:,:,2), ut%image)
       deallocate(fk)
    else
       if(.not. allocated(this%smooth_image)) stop "you have to call get_flatmap before doing T2QTUT"
       COOP_DEALLOC(this%smooth_Q)
       COOP_DEALLOC(this%smooth_U)
       allocate( &
            this%smooth_Q(-this%smooth_nx:this%smooth_nx, -this%smooth_ny:this%smooth_ny), &
            this%smooth_U(-this%smooth_nx:this%smooth_nx, -this%smooth_ny:this%smooth_ny) )
       this%smooth_Q = this%smooth_image
       this%smooth_U = 0.d0
       call coop_fits_EB2QU(nx = this%smooth_nx*2+1, ny = this%smooth_ny*2+1, Emap = this%smooth_Q, Bmap = this%smooth_U)
    endif
  end subroutine coop_fits_image_cea_get_QTUT


  subroutine coop_fits_EB2QU(nx, ny, Emap, Bmap)
    COOP_INT nx, ny, i, j
    COOP_REAL Emap(nx,ny), Bmap(nx,ny),  kx, ky, k2
    COOP_COMPLEX Qk(0:nx/2, 0:ny-1), Uk(0:nx/2, 0:ny-1), Ek(0:nx/2, 0:ny-1), Bk(0:nx/2, 0:ny-1)
    call coop_fft_forward(nx, ny, Emap, Ek)
    call coop_fft_forward(nx, ny, Bmap, Bk)
    do i=0, nx/2
       do j=0, ny-1
          kx = COOP_REAL_OF(i)/nx
          if(ny - j .lt. j)then
             ky = COOP_REAL_OF(j - ny)/ny
          else
             ky = COOP_REAL_OF(j)/ny
          endif
          k2 = kx**2+ky**2
          if(k2.eq.0.d0)then
             Qk(i, j) = (0.d0, 0.d0)
             Uk(i, j) = (0.d0, 0.d0)
          else
             Qk(i, j) = ((ky**2 - kx**2)*Ek(i,j) - 2.d0*kx*ky*Bk(i,j))/k2
             Uk(i, j) = ((kx**2 - ky**2)*Bk(i,j) - 2.d0*kx*ky*Ek(i,j))/k2 
             !!healpix convention
          endif
       enddo
    enddo
    call coop_fft_backward(nx, ny, Qk, Emap)
    call coop_fft_backward(nx, ny, Uk, Bmap)
  end subroutine coop_fits_EB2QU


  subroutine coop_fits_QU2EB(nx, ny, Qmap, Umap)
    COOP_INT nx, ny, i, j
    COOP_REAL Qmap(nx,ny), Umap(nx,ny), kx, ky, k2
    COOP_COMPLEX Qk(0:nx/2, 0:ny-1), Uk(0:nx/2, 0:ny-1), Ek(0:nx/2, 0:ny-1), Bk(0:nx/2, 0:ny-1)
    call coop_fft_forward(nx, ny, Qmap, Qk)
    call coop_fft_forward(nx, ny, Umap, Uk)
    do i=0, nx/2
       do j=0, ny-1
          kx = COOP_REAL_OF(i)/nx
          if(ny - j .lt. j)then
             ky = COOP_REAL_OF(j - ny)/ny
          else
             ky = COOP_REAL_OF(j)/ny
          endif
          k2 = kx**2 + ky**2
          if(k2.eq.0.d0)then
             Ek(i, j) = 0.d0
             Bk(i, j) = 0.d0
          else
             !!healpix convention
             Ek(i, j) =  ( (ky**2 - kx**2)*Qk(i,j) - 2.d0*kx*ky*Uk(i,j))/k2
             Bk(i, j) =  ( (kx**2 - ky**2)*Uk(i,j) - 2.d0*kx*ky*Qk(i,j))/k2 
          endif
       enddo
    enddo
    call coop_fft_backward(nx, ny, Ek, Qmap)
    call coop_fft_backward(nx, ny, Bk, Umap)
  end subroutine coop_fits_QU2EB

  subroutine coop_fits_image_cea_EB2QU(Emap, Bmap)
    type(coop_fits_image_cea)::Emap, Bmap
    COOP_COMPLEX,dimension(:,:,:),allocatable::fk
    COOP_COMPLEX::tmp
    COOP_INT::i, j
    COOP_REAL::ky, ky2, kx, kx2, k2, k2min
    if(emap%nside(1).ne. bmap%nside(1) .or. emap%nside(2) .ne. bmap%nside(2) .or. abs(emap%dkx/bmap%dkx -1.d0).gt. 1.d-5 .or. abs(emap%dky/bmap%dky-1.d0).gt. 1.d-5)then
       write(*,*) "EB2QU: Emap and Bmap must have the same size"
       stop
    endif
    allocate(fk(0:emap%nside(1)/2,0:emap%nside(2)-1, 2))
    call coop_fft_forward(emap%nside(1), emap%nside(2), emap%image, fk(:,:,1))
    call coop_fft_forward(emap%nside(1), emap%nside(2), bmap%image, fk(:,:,2))
    k2min = min(emap%dkx, emap%dky)
    do j=0, emap%nside(2)-1
       if(j .gt. emap%nside(2)- j)then
             ky = (j-emap%nside(2))*emap%dky
          else
             ky = j*emap%dky
          endif
          ky2 = ky**2
          do i=0, emap%nside(1)/2
             kx = emap%dkx*i
             kx2 = kx**2
             k2 = kx2+ky2
             if(k2 .lt. k2min)then
                fk(i,j,:) = (0.d0, 0.d0)
             else
                tmp = fk(i, j, 2)
                !!healpix convention (U-> -U for IAU convention)
                fk(i, j, 2) = ((kx**2 - ky**2)*fk(i,j,2) - 2.d0*kx*ky*fk(i,j,1))/k2 
                fk(i, j, 1) = ((ky**2 - kx**2)*fk(i,j,1) - 2.d0*kx*ky*tmp)/k2
             endif
          enddo
       enddo
       call coop_fft_backward(emap%nside(1), emap%nside(2), fk(:,:,1), emap%image)
       call coop_fft_backward(emap%nside(1), emap%nside(2), fk(:,:,2), bmap%image)
       deallocate(fk)
  end subroutine coop_fits_image_cea_EB2QU



  subroutine coop_fits_image_cea_QU2EB(Qmap, umap)
    type(coop_fits_image_cea)::Qmap, Umap
    COOP_COMPLEX,dimension(:,:,:),allocatable::fk
    COOP_COMPLEX::tmp
    COOP_INT::i, j
    COOP_REAL::ky, ky2, kx, kx2, k2, k2min
    if(qmap%nside(1).ne. umap%nside(1) .or. qmap%nside(2) .ne. umap%nside(2) .or. abs(qmap%dkx/umap%dkx -1.d0).gt. 1.d-5 .or. abs(qmap%dky/umap%dky-1.d0).gt. 1.d-5)then
       write(*,*) "qu2eb: Qmap and Umap must have the same size"
       stop
    endif
    allocate(fk(0:qmap%nside(1)/2,0:qmap%nside(2)-1, 2))
    call coop_fft_forward(qmap%nside(1), qmap%nside(2), qmap%image, fk(:,:,1))
    call coop_fft_forward(qmap%nside(1), qmap%nside(2), umap%image, fk(:,:,2))
    do j=0, qmap%nside(2)-1
       if(j .gt. qmap%nside(2)- j)then
             ky = (j-qmap%nside(2))*qmap%dky
          else
             ky = j*qmap%dky
          endif
          ky2 = ky**2
          do i=0, qmap%nside(1)/2
             kx = qmap%dkx*i
             kx2 = kx**2
             k2 = kx2+ky2
             if(k2 .eq. 0.d0)then
                fk(i,j,:) = (0.d0, 0.d0)
             else
                tmp = fk(i, j, 2)
                !!healpix convention (U-> -U for IAU convention)
                fk(i, j, 2) = ((kx**2 - ky**2)*fk(i,j,2) - 2.d0*kx*ky*fk(i,j,1))/k2 
                fk(i, j, 1) = ((ky**2 - kx**2)*fk(i,j,1) - 2.d0*kx*ky*tmp)/k2
             endif
          enddo
       enddo
       call coop_fft_backward(qmap%nside(1), qmap%nside(2), fk(:,:,1), qmap%image)
       call coop_fft_backward(qmap%nside(1), qmap%nside(2), fk(:,:,2), umap%image)
       deallocate(fk)
  end subroutine coop_fits_image_cea_QU2EB




  subroutine coop_flatsky_maps_get_peaks(this, sto, max_num)
    class(coop_flatsky_maps)::this
    type(coop_stacking_options)::sto
    COOP_INT,optional::max_num
    COOP_INT::npts
    COOP_REAL::mean, thetaphi(2)
    COOP_INT::i, list(8), nmaps, iskip, ix, iy, ibase, ibase_plus, ibase_minus
    logical,dimension(:),allocatable::zeros1, zeros2
    COOP_INT,parameter::max_npeaks = 1000000
    if(sto%nmaps .ne. this%nmaps) stop "get_peaks: nmaps does not agree"
    if(this%total_weight .lt. 1.d0) stop "get_peaks: no unmasked pixels"
    call sto%free()
    sto%nside = this%nside(1)
    sto%nside2 = this%nside(2)

    if(sto%index_I .ne. 0 .and. (abs(sto%I_lower_nu).lt.coop_stacking_max_threshold .or. abs(sto%I_upper_nu).lt. coop_stacking_max_threshold))then  !!rescale I
       mean = this%weighted_sum(this%map(sto%index_I)%image)/this%total_weight
       sto%sigma_I  = sqrt(this%weighted_sum( (this%map(sto%index_I)%image-mean)**2 )/this%total_weight)
       sto%I_lower = mean + sto%I_lower_nu* sto%sigma_I
       sto%I_upper = mean + sto%I_upper_nu* sto%sigma_I
    endif
    if(sto%index_L .ne. 0 .and. (abs(sto%L_lower_nu).lt.coop_stacking_max_threshold .or. abs(sto%L_upper_nu).lt. coop_stacking_max_threshold) )then  !!rescale L
       mean = this%weighted_sum(this%map(sto%index_L)%image)/this%total_weight
       sto%sigma_L  = sqrt(this%weighted_sum((this%map(sto%index_L)%image-mean)**2)/this%total_weight)
       sto%L_lower = mean+sto%L_lower_nu* sto%sigma_L
       sto%L_upper = mean+sto%L_upper_nu* sto%sigma_L
    endif
    if(sto%index_Q .ne. 0 .and. sto%index_U .ne. 0 .and. (abs(sto%P_lower_nu).lt.coop_stacking_max_threshold .or. abs(sto%P_upper_nu).lt. coop_stacking_max_threshold) )then  !!rescale P2
       sto%sigma_P  = sqrt(this%weighted_sum(this%map(sto%index_Q)%image**2+this%map(sto%index_U)%image**2)/this%total_weight)
       sto%P2_lower = (sto%P_lower_nu * sto%sigma_P)**2
       sto%P2_upper = (sto%P_upper_nu * sto%sigma_P)**2
    endif
    select case(sto%genre)
    case(coop_stacking_genre_saddle, coop_stacking_genre_saddle_Oriented, coop_stacking_genre_col_oriented)
       sto%index_peak = sto%index_I
       if(sto%index_peak .ne. 1 .or. sto%index_L .ne. 4 .or. sto%index_Q .ne. 2 .or. sto%index_U .ne. 3 .or. this%nmaps .lt. 6)then
          stop "get_peaks: wrong configuration for saddle points stacking"
       endif
    case(coop_stacking_genre_Imax, coop_stacking_genre_Imin, coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Imin_Oriented)
       sto%index_peak = sto%index_I
       if(sto%index_peak .le. 0 .or. sto%index_peak .gt. this%nmaps) stop "map index of peak overflow"
    case(coop_stacking_genre_Lmax, coop_stacking_genre_Lmin, coop_stacking_genre_Lmax_Oriented, coop_stacking_genre_Lmin_Oriented)
       sto%index_peak = sto%index_L
       if(sto%index_peak .le. 0 .or. sto%index_peak .gt. this%nmaps) stop "map index of peak overflow"
    case(coop_stacking_genre_random_hot, coop_stacking_genre_random_cold, coop_stacking_genre_random_hot_oriented, coop_stacking_genre_random_cold_oriented)
       if(sto%index_Q .ne. 0 .and. sto%index_U .ne. 0 .and. (abs(sto%P_lower_nu).lt.coop_stacking_max_threshold .or. abs(sto%P_upper_nu).lt. coop_stacking_max_threshold) )then
          sto%index_peak = sto%index_Q
       else
          sto%index_peak = sto%index_I
       endif
    case default
       sto%index_peak = sto%index_Q
    end select


    select case(sto%genre)
    case(coop_stacking_genre_Imax, coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Lmax, coop_stacking_genre_Lmax_Oriented)
       ibase = 0
       ibase_minus = ibase - this%nside(1)
       ibase_plus = ibase +  this%nside(1)
       if(sto%abs_threshold)then
        do iy = 1, this%nside(2)-2
             ibase_minus = ibase
             ibase = ibase_plus
             ibase_plus = ibase_plus + this%nside(1)
             do ix = 1, this%nside(1) - 2
                i  = ibase + ix
                if(this%is_masked(i) .or. sto%reject( this%pixel_values(i) ))cycle
                list(1) = ibase_minus + ix - 1
                list(2) = ibase_minus + ix
                list(3) = ibase_minus + ix + 1
                list(4) = ibase + ix - 1
                list(5) = ibase + ix + 1
                list(6) = ibase_plus + ix - 1
                list(7) = ibase_plus + ix
                list(8) = ibase_plus + ix + 1
                if(this%map(sto%index_peak)%image(i) .gt. 0.d0)then
                   if( all(this%map(sto%index_peak)%image(list) .lt. this%map(sto%index_peak)%image(i)))then
                      call sto%peak_pix%push(i)
                      call this%map(sto%index_peak)%pix2ang(i, thetaphi(1), thetaphi(2))
                      call sto%peak_ang%push(real(thetaphi))
                      call sto%peak_map%push( this%pixel_values(i) )
                   endif
                else
                   if( all(this%map(sto%index_peak)%image(list) .gt. this%map(sto%index_peak)%image(i)))then
                      call sto%peak_pix%push(i)
                      call this%map(sto%index_peak)%pix2ang(i, thetaphi(1), thetaphi(2))
                      call sto%peak_ang%push(real(thetaphi))
                      call sto%peak_map%push( this%pixel_values(i) )
                   endif
                endif
             enddo
          enddo
       else
          do iy = 1, this%nside(2)-2
             ibase_minus = ibase
             ibase = ibase_plus
             ibase_plus = ibase_plus + this%nside(1)
             do ix = 1, this%nside(1) - 2
                i  = ibase + ix
                if(this%is_masked(i) .or. sto%reject( this%pixel_values(i) ))cycle
                list(1) = ibase_minus + ix - 1
                list(2) = ibase_minus + ix
                list(3) = ibase_minus + ix + 1
                list(4) = ibase + ix - 1
                list(5) = ibase + ix + 1
                list(6) = ibase_plus + ix - 1
                list(7) = ibase_plus + ix
                list(8) = ibase_plus + ix + 1
                if( all(this%map(sto%index_peak)%image(list) .lt. this%map(sto%index_peak)%image(i)))then
                   call sto%peak_pix%push(i)
                   call this%map(sto%index_peak)%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push( this%pixel_values(i) )
                endif
             enddo
          enddo
       endif
    case(coop_stacking_genre_Imin, coop_stacking_genre_Imin_Oriented, coop_stacking_genre_Lmin, coop_stacking_genre_Lmin_Oriented)
       ibase = 0
       ibase_minus = ibase - this%nside(1)
       ibase_plus = ibase +  this%nside(1)
       do iy = 1, this%nside(2)-2
          ibase_minus = ibase
          ibase = ibase_plus
          ibase_plus = ibase_plus + this%nside(1)
          do ix = 1, this%nside(1) - 2
             i  = ibase + ix
             if(this%is_masked(i) .or. sto%reject( this%pixel_values(i) ))cycle
             list(1) = ibase_minus + ix - 1
             list(2) = ibase_minus + ix
             list(3) = ibase_minus + ix + 1
             list(4) = ibase + ix - 1
             list(5) = ibase + ix + 1
             list(6) = ibase_plus + ix - 1
             list(7) = ibase_plus + ix
             list(8) = ibase_plus + ix + 1
             if( all(this%map(sto%index_peak)%image(list) .gt. this%map(sto%index_peak)%image(i)))then
                call sto%peak_pix%push(i)
                call this%map(sto%index_peak)%pix2ang(i, thetaphi(1), thetaphi(2))
                call sto%peak_ang%push(real(thetaphi))
                call sto%peak_map%push( this%pixel_values(i) )
             endif
          enddo
       enddo
    case(coop_stacking_genre_random_hot, coop_stacking_genre_random_cold, coop_stacking_genre_random_hot_oriented, coop_stacking_genre_random_cold_oriented)
       npts = 0
       do i=0, this%npix-1
          if(this%is_masked(i) .or. sto%reject( this%pixel_values(i)))cycle
          npts = npts + 1
       enddo
       if(present(max_num))then
          iskip = max(nint(sqrt(dble(npts)/min(max_npeaks, max_num))), 1)
       else
          iskip = max(nint(sqrt(dble(npts)/max_npeaks)), 1)
       endif
       do iy = iskip/2, this%nside(2)-1, iskip
          ibase = iy*this%nside(1)
          do ix=iskip/2, this%nside(1)-1, iskip
             i = ibase + ix
             if(this%is_masked(i) .or. sto%reject( this%pixel_values(i)))cycle
             call sto%peak_pix%push(i)
             call this%map(sto%index_peak)%pix2ang(i, thetaphi(1), thetaphi(2))
             call sto%peak_ang%push(real(thetaphi))
             call sto%peak_map%push( this%pixel_values(i) )
          enddo
       enddo
    case(coop_stacking_genre_Pmax_Oriented)
       ibase = 0
       ibase_minus = ibase - this%nside(1)
       ibase_plus = ibase +  this%nside(1)
       do iy = 1, this%nside(2)-2
          ibase_minus = ibase
          ibase = ibase_plus
          ibase_plus = ibase_plus + this%nside(1)
          do ix = 1, this%nside(1) - 2
             i  = ibase + ix
             if(this%is_masked(i) .or. sto%reject( this%pixel_values(i) ))cycle
             list(1) = ibase_minus + ix - 1
             list(2) = ibase_minus + ix
             list(3) = ibase_minus + ix + 1
             list(4) = ibase + ix - 1
             list(5) = ibase + ix + 1
             list(6) = ibase_plus + ix - 1
             list(7) = ibase_plus + ix
             list(8) = ibase_plus + ix + 1
             if( all(this%map(sto%index_Q)%image(list)**2+this%map(sto%index_U)%image(list)**2 .lt. this%map(sto%index_Q)%image(i)**2 +  this%map(sto%index_U)%image(i)**2 ))then
                call sto%peak_pix%push(i)
                call this%map(sto%index_peak)%pix2ang(i, thetaphi(1), thetaphi(2))
                call sto%peak_ang%push(real(thetaphi))
                call sto%peak_map%push( this%pixel_values(i) )
             endif
          enddo
       enddo
    case(coop_stacking_genre_Pmin_Oriented)
       ibase = 0
       ibase_minus = ibase - this%nside(1)
       ibase_plus = ibase +  this%nside(1)
       do iy = 1, this%nside(2)-2
          ibase_minus = ibase
          ibase = ibase_plus
          ibase_plus = ibase_plus + this%nside(1)
          do ix = 1, this%nside(1) - 2
             i  = ibase + ix
             if(this%is_masked(i) .or. sto%reject( this%pixel_values(i) ))cycle
             list(1) = ibase_minus + ix - 1
             list(2) = ibase_minus + ix
             list(3) = ibase_minus + ix + 1
             list(4) = ibase + ix - 1
             list(5) = ibase + ix + 1
             list(6) = ibase_plus + ix - 1
             list(7) = ibase_plus + ix
             list(8) = ibase_plus + ix + 1
             if( all(this%map(sto%index_Q)%image(list)**2+this%map(sto%index_U)%image(list)**2 .gt. this%map(sto%index_Q)%image(i)**2 +  this%map(sto%index_U)%image(i)**2 ))then
                call sto%peak_pix%push(i)
                call this%map(sto%index_peak)%pix2ang(i, thetaphi(1), thetaphi(2))
                call sto%peak_ang%push(real(thetaphi))
                call sto%peak_map%push( this%pixel_values(i) )
             endif
          enddo
       enddo
    case(coop_stacking_genre_saddle, coop_stacking_genre_saddle_Oriented, coop_stacking_genre_col_oriented)
       allocate(zeros1(0:this%npix-1), zeros2(0:this%npix-1))
       call this%get_zeros(5, zeros1)
       call this%get_zeros(6, zeros2)
       zeros1 = zeros1.and. zeros2
       if(this%has_mask) zeros1 = zeros1 .and. this%unmasked
       do i=0, this%npix-1
          if(zeros1(i) .and. .not. sto%reject(this%pixel_values(i)) .and. this%map(2)%image(i)**2+this%map(3)%image(i)**2 .gt. this%map(4)%image(i)**2 )then
             call sto%peak_pix%push(i)
             call this%map(sto%index_peak)%pix2ang(i, thetaphi(1), thetaphi(2))
             call sto%peak_ang%push(real(thetaphi))
             call sto%peak_map%push( this%pixel_values(i) )
          endif
       enddo
       deallocate(zeros1, zeros2)
    case default
       stop "flatsky_maps_get_peaks: unknown stacking genre"
    end select

  end subroutine coop_flatsky_maps_get_peaks


  subroutine coop_flatsky_maps_stack_on_peaks(this, sto, patch, norm)
    COOP_INT,parameter::n_threads = 8
    class(coop_flatsky_maps)::this
    type(coop_stacking_options)::sto
    type(coop_healpix_patch)::patch
    type(coop_healpix_patch),dimension(n_threads)::p, tmp
    logical,optional::norm
    COOP_INT ithread, i
    if(this%nside(1) .ne. sto%nside .or. this%nside(2) .ne. sto%nside2) stop "flatsky stacking only supports maps with same nside"
    patch%image = 0.d0
    patch%nstack = 0.d0
    patch%nstack_raw = 0
    !!adjust the index of maps
    i = 1
    do while(i.le.patch%nmaps)
       if(patch%tbs%ind(i) .gt. this%nmaps) stop "stack_on_peaks: index overflow"       
       do while(patch%tbs%spin(i) .ne. this%spin(patch%tbs%ind(i)))
          patch%tbs%ind(i)  = patch%tbs%ind(i) + 1
          if(patch%tbs%ind(i) .gt. this%nmaps) stop "stack_on_peaks: cannot find the map with corresponding spin"
       enddo
       if(patch%tbs%spin(i) .ne. 0)then
          if(i.ge.patch%nmaps .or. patch%tbs%ind(i) .ge. this%nmaps) stop "stack_on_peaks: nonzero spin must go in pairs"
          patch%tbs%ind(i+1) = patch%tbs%ind(i) + 1
          i = i+2
       else
          i = i+1
       endif
    enddo
    do ithread=1, n_threads
       p(ithread) = patch
       tmp(ithread) = patch
    enddo
    !$omp parallel do private(i, ithread)
    do ithread = 1, n_threads
       do i=ithread, sto%peak_pix%n, n_threads
          call this%stack_on_patch(center = sto%peak_pix%element(i), angle = sto%rotate_angle(i), patch = p(ithread), tmp_patch = tmp(ithread), factor = sto%norm(i), weight = sto%wnorm(i))
       enddo
    enddo
    !$omp end parallel do

    do ithread = 1, n_threads
       patch%image = patch%image + p(ithread)%image
       patch%nstack = patch%nstack + p(ithread)%nstack
       patch%nstack_raw = patch%nstack_raw + p(ithread)%nstack_raw
       call p(ithread)%free()
       call tmp(ithread)%free()
    enddo
    if(present(norm))then
       if(.not.norm)return
    endif
    if(patch%nstack_raw .ne. 0)then
       do i=1, patch%nmaps
          where(patch%nstack .ne. 0.d0)
             patch%image(:, :, i) = patch%image(:, :, i)/patch%nstack
          end where
       enddo
    else
       write(*,*) "warning: no patches has been found"
       patch%image = 0.d0
    endif
  end subroutine coop_flatsky_maps_stack_on_peaks


  subroutine coop_flatsky_maps_stack_on_patch(this, center, angle, patch, tmp_patch, factor, weight)
    class(coop_flatsky_maps)::this
    COOP_REAL angle
    COOP_INT::center
    type(coop_healpix_patch) patch, tmp_patch
    COOP_REAL::factor, weight
    if(angle .ge. 1.d30)return
    call this%fetch_patch(center = center, angle = angle, patch = tmp_patch)
    if(sum(tmp_patch%nstack*tmp_patch%indisk) .lt. patch%num_indisk_tol)return
    if(angle .gt. coop_4pi)then   !! if angle > 4pi do flip
       call tmp_patch%flipx()
    elseif(angle .lt. -coop_4pi)then
       call tmp_patch%flipy()
    endif
    if(factor .ne. 1.d0)then
       patch%image = patch%image + tmp_patch%image*factor
    else
       patch%image = patch%image + tmp_patch%image
    endif
    if(weight .ne. 1.d0)then
       patch%nstack = patch%nstack + tmp_patch%nstack*weight
    else
       patch%nstack = patch%nstack + tmp_patch%nstack
    endif
    patch%nstack_raw = patch%nstack_raw + tmp_patch%nstack_raw
  end subroutine coop_flatsky_maps_stack_on_patch


  subroutine coop_flatsky_maps_fetch_patch(this, patch, center, angle)
    class(coop_flatsky_maps)::this
    COOP_INT::pix, center
    COOP_REAL angle
    type(coop_healpix_patch) patch
    COOP_INT i, j, k, ix, iy, ix_center, iy_center
    COOP_REAL x, y, r, phi
    COOP_SINGLE qu(2)
    iy_center = center/this%nside(1)
    ix_center = center - iy_center * this%nside(1)
    patch%nstack_raw = 1
    if(all(patch%tbs%spin .eq. 0))then
       do j = -patch%n, patch%n
          do i = -patch%n, patch%n
             x = patch%dr * i
             y = patch%dr * j
             r = sqrt(x**2+y**2)
             phi = COOP_POLAR_ANGLE(x, y) + angle
             ix = nint(r*cos(phi)/this%dx) + ix_center
             iy = nint(r*sin(phi)/this%dy) + iy_center
             if(ix .ge. 0 .and. ix .lt. this%nside(1) .and. iy .ge. 0 .and. iy .lt. this%nside(2))then
                pix = ix  + iy*this%nside(1)
                if(this%is_unmasked(pix))then
                   patch%nstack(i, j) = 1.
                   patch%image(i, j, :) = this%pixel_values(pix)
                else
                   patch%image(i, j, :) = 0.
                   patch%nstack(i, j) = 0.
                endif
             else
                patch%image(i, j, :) = 0.
                patch%nstack(i, j) = 0.
             endif
          enddo
       enddo
    elseif(patch%nmaps .eq. 2 .and. all(patch%tbs%spin .ne. 0))then
       if(all(patch%tbs%local_rotation))then
          do j = -patch%n, patch%n
             do i = -patch%n, patch%n
                x = patch%dr * i
                y = patch%dr * j
                r = sqrt(x**2+y**2)
                phi = COOP_POLAR_ANGLE(x, y) + angle
                ix = nint(r*cos(phi)/this%dx) + ix_center
                iy = nint(r*sin(phi)/this%dy) + iy_center
                if(ix .ge. 0 .and. ix .lt. this%nside(1) .and. iy .ge. 0 .and. iy .lt. this%nside(2))then
                   pix = ix  + iy*this%nside(1)
                   qu(1) = this%map(patch%tbs%ind(1))%image(pix)
                   qu(2) = this%map(patch%tbs%ind(2))%image(pix)
                   call coop_healpix_rotate_qu(qu, phi,patch%tbs%spin(1))
                   qu = qu*coop_healpix_QrUrSign
                   if(this%is_unmasked(pix))then
                      patch%nstack(i, j) = 1.d0
                      patch%image(i, j, 1:2) = qu
                   else 
                      patch%nstack(i, j) = 0.d0
                      patch%image(i, j, 1:2) = 0.d0
                   endif
                else
                   patch%image(i, j, :) = 0.
                   patch%nstack(i, j) = 0.
                endif
             enddo
          enddo
       else
          do j = -patch%n, patch%n
             do i = -patch%n, patch%n
                x = patch%dr * i
                y = patch%dr * j
                r = sqrt(x**2+y**2)
                phi = COOP_POLAR_ANGLE(x, y) + angle
                ix = nint(r*cos(phi)/this%dx) + ix_center
                iy = nint(r*sin(phi)/this%dy) + iy_center
                if(ix .ge. 0 .and. ix .lt. this%nside(1) .and. iy .ge. 0 .and. iy .lt. this%nside(2))then
                   pix = ix  + iy*this%nside(1)
                   qu(1) = this%map(patch%tbs%ind(1))%image(pix)
                   qu(2) = this%map(patch%tbs%ind(2))%image(pix)
                   call coop_healpix_rotate_qu(qu, angle,patch%tbs%spin(1))
                   if(this%is_unmasked(pix))then
                      patch%nstack(i, j) = 1.d0
                      patch%image(i, j, 1:2) = qu
                   else 
                      patch%nstack(i, j) = 0.d0
                      patch%image(i, j, 1:2) = 0.d0
                   endif
                else
                   patch%image(i, j, :) = 0.
                   patch%nstack(i, j) = 0.
                endif
             enddo
          enddo
       endif
    else
       do j = -patch%n, patch%n
          do i = -patch%n, patch%n
             x = patch%dr * i
             y = patch%dr * j
             r = sqrt(x**2+y**2)
             phi = COOP_POLAR_ANGLE(x, y) + angle
             ix = nint(r*cos(phi)/this%dx) + ix_center
             iy = nint(r*sin(phi)/this%dy) + iy_center
             if(ix .ge. 0 .and. ix .lt. this%nside(1) .and. iy .ge. 0 .and. iy .lt. this%nside(2))then
                pix = ix  + iy*this%nside(1)
                if(this%is_masked(pix))then
                   patch%nstack(i,j) = 0.d0
                   patch%image(i, j, :) = 0.d0
                   cycle
                else
                   patch%nstack(i,j) = 1.d0
                endif
             else
                patch%nstack(i,j) = 0.d0
                patch%image(i, j, :) = 0.d0
                cycle
             endif      

             k = 1
             do while(k.le. patch%nmaps)
                if(patch%tbs%spin(k) .eq. 0)then
                   patch%image(i, j, k) = this%map(patch%tbs%ind(k))%image(pix)
                   k = k + 1
                else
                   qu(1) = this%map(patch%tbs%ind(k))%image(pix)
                   qu(2) = this%map(patch%tbs%ind(k+1))%image(pix)
                   if(patch%tbs%local_rotation(k))then
                      call coop_healpix_rotate_qu(qu, phi, patch%tbs%spin(k))
                      qu = qu * coop_healpix_QrUrSign
                   else
                      call coop_healpix_rotate_qu(qu, angle,patch%tbs%spin(k))
                   endif
                   patch%image(i, j, k:k+1) = qu
                   k = k + 2
                end if
             enddo
          enddo
       enddo
    endif
  end subroutine coop_flatsky_maps_fetch_patch

  subroutine coop_fits_image_cea_smooth_EB2QU(this)
    class(coop_fits_image_cea)::this
    if(.not. allocated(this%smooth_Q).or. .not. allocated(this%smooth_U)) stop "EB2QU cannot be done: E, B not allocated yet."
    call coop_fits_EB2QU(this%smooth_nx*2+1, this%smooth_ny*2+1, this%smooth_Q, this%smooth_U)
  end subroutine coop_fits_image_cea_smooth_EB2QU


  subroutine coop_fits_image_cea_smooth_QU2EB(this)
    class(coop_fits_image_cea)::this
    if(.not. allocated(this%smooth_Q).or. .not. allocated(this%smooth_U)) stop "QU2EB cannot be done: Q, U not allocated yet."
    call coop_fits_QU2EB(this%smooth_nx*2+1, this%smooth_ny*2+1, this%smooth_Q, this%smooth_U)
  end subroutine coop_fits_image_cea_smooth_QU2EB




  subroutine coop_fits_image_cea_simulate(this, lmin, lmax, Cls)
    class(coop_fits_image_cea)::this
    COOP_INT lmin, lmax, i, j, ik
    COOP_REAL Cls(lmin:lmax)
    COOP_COMPLEX,dimension(:,:),allocatable::fk
    COOP_REAL amp(lmin:lmax), rk
    amp = sqrt(Cls)/this%pixsize**2
    allocate(fk(0:this%nside(1)/2,0:this%nside(2)-1))
    do i=0, this%nside(1)/2
       do j=0, this%nside(2)/2
          rk  = sqrt((this%dkx*i)**2 + (this%dky*j)**2)
          ik = floor(rk)
          if(ik.ge.lmin .and. ik .lt. lmax)then
             rk = rk - ik
             fk(i,j) = coop_random_complex_gaussian() * (amp(ik)*(1.d0-rk) + amp(ik+1)*rk)
          else
             fk(i, j) = 0
          endif
       enddo
       do j=this%nside(2)/2+1, this%nside(2)-1
          rk  = sqrt((this%dkx*i)**2 + (this%dky*(this%nside(2)-j))**2) 
          ik = floor(rk)
          if(ik.ge.lmin .and. ik .lt. lmax)then
             rk = rk - ik
             fk(i,j) = coop_random_complex_gaussian() * (amp(ik)*(1.d0-rk) + amp(ik+1)*rk)
          else
             fk(i, j) = 0
          endif
       enddo
    enddo

    if(mod(this%nside(2),2) .eq. 0)then
       if(fk(0, this%nside(2)/2) .ne. ( 0.d0, 0.d0 ) )then
          fk(0, this%nside(2)/2) = coop_sqrt2 * real( fk(0, this%nside(2)/2) )
       endif
       if(mod(this%nside(1), 2).eq.0 .and. fk(this%nside(1)/2, this%nside(2)/2) .ne. ( 0.d0, 0.d0 ) ) then
          fk(this%nside(1)/2, this%nside(2)/2) = coop_sqrt2 * real( fk(this%nside(1)/2, this%nside(2)/2) )
       endif
    endif
    if(mod(this%nside(1), 2) .eq. 0)then
       i = this%nside(1)/2
       do j = this%nside(2)/2+1, this%nside(2)-1
          fk(i, j) = conjg(fk(i, this%nside(2)-j))
       enddo
    endif
    call coop_fft_backward(this%nside(1), this%nside(2), fk, this%image)
    deallocate(fk)
  end subroutine coop_fits_image_cea_simulate


  subroutine coop_fits_image_cea_write(this, filename)
    class(coop_fits_image_cea)::this
    COOP_UNKNOWN_STRING,intent(in):: filename
    if(coop_file_exists(filename)) call coop_delete_file(filename)
    call coop_fits_file_write_image(this%image, filename, this%header)
  end subroutine coop_fits_image_cea_write

  subroutine coop_fits_image_cea_plot(this, figname)
    class(coop_fits_image_cea)::this
    COOP_UNKNOWN_STRING::figname
    type(coop_asy)::asy
    if(.not. allocated(this%smooth_image)) stop "call get_flatmap beore plot"
    call asy%open(figname)
    call asy%init(xlabel = "$2\sin{\frac{\theta}{2}}\cos\varphi$", ylabel = "$2\sin{\frac{\theta}{2}}\sin\varphi$", caption = "noise-free simulation", width=7., height=5.5, nxticks = 5, nyticks = 5)
    call coop_asy_density(asy, this%smooth_image, this%xmin, this%xmax, this%ymin, this%ymax, "$I (\mu K)$")
    call asy%close()
    if(allocated(this%smooth_Q))then
       call asy%open(coop_file_add_postfix(figname, "_Q"))
       call asy%init(xlabel = "$2\sin{\frac{\theta}{2}}\cos\varphi$", ylabel = "$2\sin{\frac{\theta}{2}}\sin\varphi$", caption = "noise-free simulation", width=7., height=5.5, nxticks = 5, nyticks = 5)
       call coop_asy_density(asy, this%smooth_Q, this%xmin, this%xmax, this%ymin, this%ymax, "$Q (\mu K)$")
       call asy%close()
    endif

    if(allocated(this%smooth_U))then
       call asy%open(coop_file_add_postfix(figname, "_U"))
       call asy%init(xlabel = "$2\sin{\frac{\theta}{2}}\cos\varphi$", ylabel = "$2\sin{\frac{\theta}{2}}\sin\varphi$", caption = "noise-free simulation", width=7., height=5.5, nxticks = 5, nyticks = 5)
       call coop_asy_density(asy, this%smooth_U, this%xmin, this%xmax, this%ymin, this%ymax, "$U (\mu K)$")
       call asy%close()
    endif

  end subroutine coop_fits_image_cea_plot

  subroutine  coop_fits_image_cea_convert2healpix(this, hp, imap, mask, weights)
    COOP_INT,parameter::n_threads = 8
    class(coop_fits_image_cea)::this
    type(coop_fits_image_cea),optional::weights 
    type(coop_healpix_maps)::hp, mask
    COOP_INT::pix
    COOP_INT:: hpix, imap
    COOP_REAL theta, phi, mean_weights
#ifdef HAS_HEALPIX    
    if(mask%nside .ne. hp%nside .or. mask%ordering .ne. hp%ordering) stop "fits_image_cea_convert2healpix: mask and map must have the same nside and ordering"
    mask%map(:,1) = 0.
    hp%map(:, imap) = 0.
    do pix=1, this%npix-1
       if(present(weights))then
          if(weights%image(pix) .eq. 0.)cycle
       endif
       call this%pix2ang(pix, theta, phi)
       call hp%ang2pix(theta, phi, hpix)
       if(present(weights))then
          mask%map(hpix, 1) = mask%map(hpix, 1) + weights%image(pix)
          hp%map(hpix,imap) = hp%map(hpix, imap) + this%image(pix)*weights%image(pix)
       else
          mask%map(hpix,1) = mask%map(hpix,1) + 1.
          hp%map(hpix,imap) = hp%map(hpix, imap) + this%image(pix)
       endif
    enddo
    mean_weights = 0.5*sqrt(sum(mask%map(:,1)**2)/count(mask%map(:,1).ne.0.))
    where(mask%map(:,1) .ne. mean_weights)
       hp%map(:, imap) = hp%map(:, imap)/mask%map(:,1)
       mask%map(:, 1) = 1.
    elsewhere
       hp%map(:,imap) = 0.
       mask%map(:, 1) = 0.
    end where
#else
    stop "you need to install Healpix"
#endif    
  end subroutine coop_fits_image_cea_convert2healpix



  subroutine  coop_fits_image_cea_from_healpix(this, hp, imap)  !!convert from healpix map
    class(coop_fits_image_cea)::this
    type(coop_healpix_maps)::hp
    COOP_INT::pix, i, j, hpix, imap
    COOP_REAL::theta, phi
#ifdef HAS_HEALPIX    
    !$omp parallel do private(pix, theta, phi, hpix)
    do pix=0, this%npix-1
       call this%pix2ang(pix, theta, phi)
       call hp%ang2pix(theta, phi, hpix)
       this%image(pix) = hp%map(hpix, imap)
    enddo
    !$omp end parallel do
#else
    stop "you need to install Healpix"
#endif    
  end subroutine coop_fits_image_cea_from_healpix


  subroutine coop_fits_image_cea_simulate_TEB(lmin, lmax, Cls, tmap, emap, bmap, qmap, umap)
    type(coop_fits_image_cea)::tmap,emap,bmap
    type(coop_fits_image_cea),optional::qmap,umap
    COOP_INT lmin, lmax, i, j, ik
    COOP_REAL Cls(lmin:lmax,6)
    COOP_REAL::Cls_sqrteig(3,lmin:lmax), Cls_rot(3,3,lmin:lmax)
    COOP_COMPLEX,dimension(:,:,:),allocatable::fk
    COOP_REAL amp(lmin:lmax), ky2, k2max, k2min, k2, kx, ky, kx2
    k2max = lmax**2*0.9999999
    k2min = lmin**2*1.0000001
    call coop_healpix_cls2rot(lmin, lmax, real(cls), cls_sqrteig, cls_rot)
    cls_sqrteig = cls_sqrteig*sqrt(dble(tmap%npix))/tmap%pixsize
    allocate(fk(0:tmap%nside(1)/2,0:tmap%nside(2)-1,3))
    do j=0, tmap%nside(2)-1
       if(j .gt. tmap%nside(2)- j)then
          ky = (j-tmap%nside(2))*tmap%dky
       else
          ky = j*tmap%dky
       endif
       ky2 = ky**2
       if(ky2.ge.k2max)then
          fk(:,j,:) = 0.
          cycle
       endif
       do i=0, tmap%nside(1)/2
          k2 = (tmap%dkx*i)**2 + ky2
          if(k2 .le. k2min .or. k2.ge.k2max)then
             fk(i,j,:)= 0.
             cycle
          endif
          ik = nint(sqrt(k2))
          fk(i,j,:) = matmul(cls_rot(:,:,ik), cls_sqrteig(:,ik)* (/ coop_random_complex_gaussian(), coop_random_complex_gaussian(), coop_random_complex_gaussian() /) ) 
       enddo
    enddo
    if(mod(tmap%nside(2),2) .eq. 0)then
       fk(0, tmap%nside(2)/2,:) = coop_sqrt2 * real( fk(0, tmap%nside(2)/2,:) )
       if(mod(tmap%nside(1), 2).eq.0) &
            fk(tmap%nside(1)/2, tmap%nside(2)/2,:) = coop_sqrt2 * real( fk(tmap%nside(1)/2, tmap%nside(2)/2,:) )
    endif
    if(mod(tmap%nside(1), 2) .eq. 0)then
       i = tmap%nside(1)/2
       do j = tmap%nside(2)/2+1, tmap%nside(2)-1
          fk(i, j,:) = conjg(fk(i, tmap%nside(2)-j,:))
       enddo
    endif
    call coop_fft_backward(tmap%nside(1), tmap%nside(2), fk(:,:,1), tmap%image)
    if(present(Qmap) .and. present(umap))then
       fk(:,:,1) = fk(:,:,2)
       call coop_fft_backward(tmap%nside(1), tmap%nside(2), fk(:,:,1), emap%image)
       fk(:,:,1) = fk(:,:,3)
       call coop_fft_backward(tmap%nside(1), tmap%nside(2), fk(:,:,1), bmap%image)

       fk(:,:,1) = 0.
       do j=0, tmap%nside(2)-1
          if(j .gt. tmap%nside(2)- j)then
             ky = (j-tmap%nside(2))*tmap%dky
          else
             ky = j*tmap%dky
          endif
          ky2 = ky**2
          if(ky2.ge.k2max) cycle
          do i=0, tmap%nside(1)/2
             kx = tmap%dkx*i
             kx2 = kx**2
             k2 = kx2+ky2
             if(k2 .le. k2min .or. k2.ge.k2max)cycle
             fk(i, j,1) = ((ky2 - kx2)*fk(i,j,2) - 2.d0*kx*ky*fk(i,j,3))/k2
             fk(i, j,2) = ((kx2 - ky2)*fk(i,j,3) - 2.d0*kx*ky*fk(i,j,2))/k2 !!healpix convention
          enddo
       enddo

       call coop_fft_backward(tmap%nside(1), tmap%nside(2), fk(:,:,1), qmap%image)
       call coop_fft_backward(tmap%nside(1), tmap%nside(2), fk(:,:,2), umap%image)
    else
       call coop_fft_backward(tmap%nside(1), tmap%nside(2), fk(:,:,2), emap%image)
       call coop_fft_backward(tmap%nside(1), tmap%nside(2), fk(:,:,3), bmap%image)
    
    endif
    deallocate(fk)
  end subroutine coop_fits_image_cea_simulate_TEB


  subroutine coop_fits_image_cea_get_neighbours(this, pix, list, nneigh)
    class(coop_fits_image_cea)::this
    COOP_INT::pix
    COOP_INT::list(8), nneigh
    COOP_INT ix, iy, i
    iy = pix/this%nside(1) 
    ix = pix  - iy*this%nside(1) 
    nneigh = 0
    if(iy .gt. 0)then
       nneigh = nneigh + 1
       list(nneigh) = ix + (iy-1)*this%nside(1)
    endif
    if(iy .lt. this%nside(2)-1)then
       nneigh = nneigh + 1
       list(nneigh) = ix + (iy+1)*this%nside(1)          
    endif

    if(ix .gt. 0)then
       nneigh = nneigh + 1
       list(nneigh) = (ix-1) + iy*this%nside(1)
       if(iy .gt. 0)then
          nneigh = nneigh + 1
          list(nneigh) = (ix-1) + (iy-1)*this%nside(1)
       endif
       if(iy .lt. this%nside(2)-1)then
          nneigh = nneigh + 1
          list(nneigh) = (ix-1) + (iy+1)*this%nside(1)          
       endif
    endif

    if(ix .lt. this%nside(1)-1)then
       nneigh = nneigh + 1
       list(nneigh) = (ix+1) + iy*this%nside(1)
       if(iy .gt. 0)then
          nneigh = nneigh + 1
          list(nneigh) = (ix+1) + (iy-1)*this%nside(1)
       endif
       if(iy .lt. this%nside(2)-1)then
          nneigh = nneigh + 1
          list(nneigh) = (ix+1) + (iy+1)*this%nside(1)          
       endif
    endif
  end subroutine coop_fits_image_cea_get_neighbours


  subroutine coop_flatsky_maps_get_neighbours(this, pix, list, nneigh)
    class(coop_flatsky_maps)::this
    COOP_INT::pix
    COOP_INT::list(8), nneigh
    COOP_INT ix, iy, i
    iy = pix/this%nside(1) 
    ix = pix  - iy*this%nside(1) 
    nneigh = 0
    if(iy .gt. 0)then
       nneigh = nneigh + 1
       list(nneigh) = ix + (iy-1)*this%nside(1)
    endif
    if(iy .lt. this%nside(2)-1)then
       nneigh = nneigh + 1
       list(nneigh) = ix + (iy+1)*this%nside(1)          
    endif

    if(ix .gt. 0)then
       nneigh = nneigh + 1
       list(nneigh) = (ix-1) + iy*this%nside(1)
       if(iy .gt. 0)then
          nneigh = nneigh + 1
          list(nneigh) = (ix-1) + (iy-1)*this%nside(1)
       endif
       if(iy .lt. this%nside(2)-1)then
          nneigh = nneigh + 1
          list(nneigh) = (ix-1) + (iy+1)*this%nside(1)          
       endif
    endif

    if(ix .lt. this%nside(1)-1)then
       nneigh = nneigh + 1
       list(nneigh) = (ix+1) + iy*this%nside(1)
       if(iy .gt. 0)then
          nneigh = nneigh + 1
          list(nneigh) = (ix+1) + (iy-1)*this%nside(1)
       endif
       if(iy .lt. this%nside(2)-1)then
          nneigh = nneigh + 1
          list(nneigh) = (ix+1) + (iy+1)*this%nside(1)          
       endif
    endif
  end subroutine coop_flatsky_maps_get_neighbours


  function coop_flatsky_mask_threshold(mask) result(threshold)
    COOP_INT,parameter::n = 50000
    class(coop_fits_image)::mask
    COOP_REAL::threshold
    COOP_REAL::  mmax, dm, sumc, c(0:n), s
    COOP_INT::   i, j
    mmax = maxval(mask%image)
    dm = mmax/n
    c = 0.d0
    do i = 0, mask%npix-1
       j = nint(mask%image(i)/dm)
       c(j) = c(j) + 1.d0
    enddo
    sumc = sum(c(1:n))*0.05  !!remove 5% of edge with nonzero mask value
    s = 0.d0
    do i = 1, n
       s = s+c(i)
       if(s .ge. sumc)then
          threshold = (i-0.5)*dm
          exit
       endif
    enddo
    threshold = max(threshold, sum(mask%image)/mask%npix*0.05) !!if 5% of the mean is larger take 5% of the mean
  end function coop_flatsky_mask_threshold


end module coop_fitswrap_mod


