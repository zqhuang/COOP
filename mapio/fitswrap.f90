module coop_fitswrap_mod
  use coop_wrapper_utils
  use coop_sphere_mod
  use coop_healpix_mod
  implicit none

#include "constants.h"

  private

  public::coop_fits, coop_fits_image, coop_fits_image_cea, coop_fits_QU2EB, coop_fits_EB2QU

  integer,parameter::sp = kind(1.)
  integer,parameter::dl = kind(1.d0)
  integer,parameter::dlc = kind( (1.d0,1.d0) )

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
     real(dl),dimension(:),allocatable::image
     COOP_REAL,allocatable::transform(:, :), center(:)
     COOP_REAL,allocatable::radec_center(:)
   contains
     procedure::regularize => coop_fits_image_regularize
     procedure::free => coop_fits_image_free
     procedure::get_linear_coordinates => coop_fits_image_get_linear_coordinates
     procedure::get_data => coop_fits_image_get_data
  end type coop_fits_image

  type, extends(coop_fits_image)::coop_fits_image_cea
     COOP_REAL smooth_pixsize, pixsize, smooth_dkx, smooth_dky, dkx, dky
     real(dl) xmin, xmax, ymin, ymax
     COOP_INT:: smooth_nx, smooth_ny, smooth_npix
     type(coop_sphere_disc)::disc
     real(dl),dimension(:,:),allocatable:: smooth_image
     real(dl),dimension(:,:),allocatable:: smooth_Q, smooth_U
   contains
     procedure::convert2healpix => coop_fits_image_cea_convert2healpix
     procedure::write => coop_fits_image_cea_write
     procedure::pix2ang => coop_fits_image_cea_pix2ang
     procedure::get_flatmap => coop_fits_image_cea_get_flatmap
     procedure::pix2flat => coop_fits_image_cea_pix2flat
     procedure::cut => coop_fits_image_cea_cut
     procedure::filter => coop_fits_image_cea_filter
     procedure::smooth_flat => coop_fits_image_cea_smooth_flat
     procedure::get_QTUT =>  coop_fits_image_cea_get_QTUT
     procedure::get_QU =>  coop_fits_image_cea_get_QU
     procedure::find_extrema => coop_fits_image_cea_find_extrema
     procedure::stack => coop_fits_image_cea_stack
     procedure::stack2fig => coop_fits_image_cea_stack2fig
     procedure::simulate => coop_fits_image_cea_simulate
     procedure::simulate_flat => coop_fits_image_cea_simulate_flat
     procedure::plot => coop_fits_image_cea_plot
     procedure::EB2QU => coop_fits_image_cea_EB2QU
     procedure::QU2EB => coop_fits_image_cea_QU2EB
  end type coop_fits_image_cea

  



contains

  subroutine coop_fits_open(this, filename)
    class(coop_fits)::this
    COOP_UNKNOWN_STRING::filename
    COOP_INT i
    if(coop_file_exists(filename))then
       this%filename = trim(filename)
       call coop_convert_to_C_string(this%filename)
       call coop_fits_get_header(this)
    else
       write(*,*) "The file "//trim(filename)//" does not exist."
    endif
    select type(this)
    class is(coop_fits_image)
       call this%get_data()
    end select
  end subroutine coop_fits_open

  subroutine coop_fits_image_free(this)
    class(coop_fits_image) :: this
    if(allocated(this%image)) deallocate(this%image)
    if(allocated(this%transform))deallocate(this%transform)
    if(allocated(this%center))deallocate(this%center)
    if(allocated(this%radec_center))deallocate(this%radec_center)    
    if(allocated(this%nside))deallocate(this%nside)
    select type(this)
    class is(coop_fits_image_cea)
       if(allocated(this%smooth_image))deallocate(this%smooth_image)
       if(allocated(this%smooth_Q))deallocate(this%smooth_Q)
       if(allocated(this%smooth_U))deallocate(this%smooth_U)
    end select
  end subroutine coop_fits_image_free

  subroutine coop_fits_get_header(this)
    class(coop_fits)::this
    COOP_LONG_STRING::header
    COOP_INT nkeys, i, j, istart, iend, ikey
    COOP_REAL,dimension(:),allocatable::delta, units
    call coop_fits_read_header_to_string(this%filename, header, nkeys)
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
       this%dim = nint(this%key_value("NAXIS", 2.d0))
       if(this%dim .eq. 0 )then
          call this%header%print()
          write(*,*) "Error: cannot find NAXIS key word in fits file "//trim(this%filename)
          stop
       endif
       call this%free()
       allocate(this%nside(this%dim))
       allocate(this%transform(this%dim, this%dim), this%center(this%dim), this%radec_center(this%dim))
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
    end select

    !!get transform, center 
    select type(this)
    class is(coop_fits_image_cea)
       if(this%dim .ne. 2) stop "For CEA map the dimension must be 2"
       allocate(units(this%dim))              
       do i=1, this%dim
          if(index(this%header%value("CUNIT"//COOP_STR_OF(i)), "deg").ne.0)then
             units(i) = coop_SI_degree
          elseif(index(this%header%value("CUNIT"//COOP_STR_OF(i)), "arcmin").ne.0)then
             units(i) = coop_SI_arcmin
          else
             write(*,*) "Unknown unit "//trim(this%header%value("CUNIT"//COOP_STR_OF(i)))
             stop
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
          delta(i) = this%key_value("CDELT"//COOP_STR_OF(i)) * coop_SI_degree
          if(abs(delta(i)) .lt. 1.d-12)then
             write(*,*) trim(this%filename)//": cannot read delta"
             call this%header%print()
             stop
          endif
          this%center(i) = this%key_value("CRPIX"//COOP_STR_OF(i), (this%nside(i)+1.d0)/2.d0)
          this%radec_center(i) = this%key_value("CRVAL"//COOP_STR_OF(i)) * coop_SI_degree
       enddo
       if(abs(abs(delta(1)/delta(2)) - 1.d0) .gt. 1.d-3)then
          print*, delta(1), delta(2)
          stop "pixel size is not the same in x and y directions"
       else
          this%pixsize = (abs(delta(1))+abs(delta(2)))/2.d0
       endif
       this%dkx = coop_2pi/this%pixsize/this%nside(1)
       this%dky = coop_2pi/this%pixsize/this%nside(2)
       do i=1, this%dim
          this%transform(:, i) =  this%transform(:, i) * delta(i)
       enddo
       deallocate(delta)
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

  subroutine coop_fits_image_get_linear_coordinates(this, pix, coor)
    class(coop_fits_image)::this
    COOP_INT pix(this%dim)
    COOP_REAL coor(this%dim)
    coor = matmul(this%transform, pix - this%center)
  end subroutine coop_fits_image_get_linear_coordinates


  subroutine coop_fits_image_get_data(this)
    class(coop_fits_image)::this
    real(sp),dimension(:),allocatable::tmp
    select case(this%bitpix)
    case(-64)
       call coop_fits_get_double_data(this%filename, this%image, this%npix)
    case(-32)
       allocate(tmp(0:this%npix-1))
       call coop_fits_get_float_data(this%filename, tmp, this%npix)
       this%image = real(tmp, dl)
       deallocate(tmp)
    case default
       write(*,*) "Cannot load data for bitpix = "//trim(coop_num2str(this%bitpix))
       stop
    end select

  end subroutine coop_fits_image_get_data

  subroutine coop_fits_image_regularize(this, tail)
    class(coop_fits_image) this
    real(dl) upper, lower, tail, diff, diff3
    call array_get_threshold_double(this%image, this%npix, 1.-tail, lower)
    call array_get_threshold_double(this%image, this%npix, tail, upper)
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

  subroutine coop_fits_image_cea_pix2ang(this, pix, theta, phi)
    class(coop_fits_image_cea)::this
    COOP_LONG_INT pix
    COOP_REAL theta, phi
    COOP_INT ix, iy
    iy = pix/this%nside(1) 
    ix = pix  - iy*this%nside(1) + 1 - this%center(1)
    iy = iy + 1 - this%center(2)
    phi = this%transform(1,1)* ix + this%transform(1,2)*iy + this%radec_center(1)
    theta = coop_pio2 - asin(this%transform(2, 1)*ix + this%transform(2,2)*iy + this%radec_center(2))
  end subroutine coop_fits_image_cea_pix2ang


  subroutine coop_fits_image_cea_pix2flat(this, pix, coor, disc)
    class(coop_fits_image_cea)::this
    COOP_LONG_INT:: pix
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
    real(dl) map(:,:)
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
    real(dl),dimension(:,:),allocatable::w
    COOP_INT i, j, icoor(2), k
    COOP_LONG_INT pix
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
    if(allocated(this%smooth_image))deallocate(this%smooth_image)
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


  subroutine coop_fits_image_cea_filter(this, lmin, lmax)
    class(coop_fits_image_cea)::this
    COOP_INT lmin, lmax, i, j
    COOP_REAL k2min, k2max, k2, rn1, rn2, kmin, kmax, omega, j2
    complex(dlc),dimension(:,:),allocatable::fk
    allocate( fk(0:this%nside(1)/2, 0:this%nside(2)-1))
    call coop_fft_forward(this%nside(1), this%nside(2), this%image, fk)
    rn1 = real(this%nside(1), dl)
    rn2  = real(this%nside(2), dl)/rn1
    k2min = (this%pixsize * lmin/coop_2pi*rn1)**2
    k2max = (this%pixsize * lmax/coop_2pi*rn1)**2
    kmin = sqrt(k2min)
    kmax = sqrt(k2max)
    omega = coop_pi/(kmax - kmin)
    do j=0, this%nside(2)-1
       j2 = (min(j, this%nside(2)-j)/ rn2)**2
       if( j2 .gt. k2max)then
          fk(:, j) = ( 0.d0, 0.d0 )
          cycle
       endif
       do i = 0, this%nside(1)/2
          k2 = j2 + real(i)**2
          if(k2 .lt. k2min .or. k2.gt. k2max)then
             fk(i, j) = ( 0.d0, 0.d0 )
          else
             fk(i, j) = fk(i, j)*sin((sqrt(k2)- kmin)*omega)**2
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
    complex(dlc),dimension(:,:),allocatable::fk, fke, fkb
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

  subroutine coop_fits_image_cea_find_extrema(this, mask, spot_file, spot_type, radius, repeat)
    class(coop_fits_image_cea)::this, mask
    COOP_UNKNOWN_STRING::spot_file, spot_type
    type(coop_file)::cf
    COOP_REAL,optional:: radius
    COOP_REAL angle
    COOP_INT i, j, irad, ire
    integer, optional::repeat
    integer::irepeat
    if(present(radius))then
       irad = ceiling(radius/this%smooth_pixsize)
    else
       irad = 1
    endif
    if(present(repeat))then
       irepeat = max(1, repeat)
    else
       irepeat = 1
    endif
    call cf%open(trim(spot_file), "w")
    select case(trim(spot_type))
    case("Tmax", "Emax", "Bmax")
       do i=-this%smooth_nx+irad, this%smooth_nx - irad
          do j = -this%smooth_ny+irad, this%smooth_ny - irad
             if(all(this%smooth_image(i, j) .ge. this%smooth_image(i-1:i+1, j-1:j+1)) .and. mask%smooth_image(i,j).gt.0.5)then
                do ire=1, irepeat
                   write(cf%unit, "(2I6, E16.7, I6 )") i, j, coop_2pi*coop_random_unit(), min(this%smooth_nx - i, this%smooth_nx + i, this%smooth_ny - j, this%smooth_ny + j)
                enddo
             endif
          enddo
       enddo
    case("Tmin", "Emin", "Bmin")
       do i=-this%smooth_nx+irad, this%smooth_nx - irad
          do j = -this%smooth_ny+irad, this%smooth_ny - irad
             if(all(this%smooth_image(i, j) .le. this%smooth_image(i-1:i+1, j-1:j+1)) .and. mask%smooth_image(i,j).gt.0.5)then
                do ire = 1, irepeat
                   write(cf%unit, "(2I6, E16.7, I6 )") i, j, coop_2pi*coop_random_unit(), min(this%smooth_nx - i, this%smooth_nx + i, this%smooth_ny - j, this%smooth_ny + j)
                enddo
             endif
          enddo
       enddo
    case("Tmax_QTUTOrient")
       if(.not. allocated(this%smooth_Q) .or. .not. allocated(this%smooth_U)) &
            call this%get_QTUT()
       do i=-this%smooth_nx+irad, this%smooth_nx - irad
          do j = -this%smooth_ny+irad, this%smooth_ny - irad
             if(all(this%smooth_image(i, j) .ge. this%smooth_image(i-1:i+1, j-1:j+1)) .and. mask%smooth_image(i,j).gt.0.5)then
                angle = 0.5d0 * COOP_POLAR_ANGLE(this%smooth_Q(i, j), this%smooth_U(i, j))
                write(cf%unit, "(2I6, E16.7, I6 )") i, j, angle, min(this%smooth_nx - i, this%smooth_nx + i, this%smooth_ny - j, this%smooth_ny + j)
             endif
          enddo
       enddo
    case("Tmin_QTUTOrient")
       if(.not. allocated(this%smooth_Q) .or. .not. allocated(this%smooth_U)) &
            call this%get_QTUT()
       do i=-this%smooth_nx+irad, this%smooth_nx - irad
          do j = -this%smooth_ny+irad, this%smooth_ny - irad
             if(all(this%smooth_image(i, j) .le. this%smooth_image(i-1:i+1, j-1:j+1)) .and. mask%smooth_image(i,j).gt.0.5 )then
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
       write(*,*) trim(spot_type)

    end select
    call cf%close()
  end subroutine coop_fits_image_cea_find_extrema

  subroutine coop_fits_image_cea_simulate_flat(this, lmin, lmax, Cls_file)
    class(coop_fits_image_cea)::this
    COOP_UNKNOWN_STRING::Cls_file
    type(coop_file)::fp
    COOP_INT lmin, lmax, i, j, ik, l, ii
    COOP_REAL Cls(3, 3, lmin:lmax), sqrteig(3, lmin:lmax), rot(3,3, lmin:lmax), rin(4)
    complex(dlc),dimension(:,:,:),allocatable::fk
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
          call coop_matsymdiag_small(3, Cls(:,:,l), rot(:,:,l))
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
       call asy%init(xlabel = "$2\sin{\frac{\theta}{2}} \cos\varphi$", ylabel = "$2\sin{\frac{\theta}{2}} \sin\varphi$", width=7., height=5.5, caption=trim(caption)//"; stacked "//trim(coop_num2str(nstack))//" patches")
    else
       call asy%init(xlabel = "$2\sin{\frac{\theta}{2}} \cos\varphi$", ylabel = "$2\sin{\frac{\theta}{2}} \sin\varphi$", width=7., height=5.5, caption="stacked "//trim(coop_num2str(nstack))//" patches")
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
    if(allocated(this%smooth_Q))deallocate(this%smooth_Q)
    if(allocated(this%smooth_U))deallocate(this%smooth_U)
    allocate( &
         this%smooth_Q(-this%smooth_nx:this%smooth_nx, -this%smooth_ny:this%smooth_ny), &
         this%smooth_U(-this%smooth_nx:this%smooth_nx, -this%smooth_ny:this%smooth_ny) )
    this%smooth_Q = Qmap%smooth_image
    this%smooth_U = Umap%smooth_image
  end subroutine coop_fits_image_cea_get_QU

  subroutine coop_fits_image_cea_get_QTUT(this)
    class(coop_fits_image_cea)::this
    if(.not. allocated(this%smooth_image)) stop "you have to call get_flatmap before doing T2QTUT"
    if(allocated(this%smooth_Q))deallocate(this%smooth_Q)
    if(allocated(this%smooth_U))deallocate(this%smooth_U)
    allocate( &
         this%smooth_Q(-this%smooth_nx:this%smooth_nx, -this%smooth_ny:this%smooth_ny), &
         this%smooth_U(-this%smooth_nx:this%smooth_nx, -this%smooth_ny:this%smooth_ny) )
    this%smooth_Q = this%smooth_image
    this%smooth_U = 0.d0
    call coop_fits_EB2QU(nx = this%smooth_nx*2+1, ny = this%smooth_ny*2+1, Emap = this%smooth_Q, Bmap = this%smooth_U)
  end subroutine coop_fits_image_cea_get_QTUT

  subroutine coop_fits_EB2QU(nx, ny, Emap, Bmap)
    COOP_INT nx, ny, i, j
    real(dl) Emap(nx,ny), Bmap(nx,ny),  kx, ky, k2
    complex(dlc) Qk(0:nx/2, 0:ny-1), Uk(0:nx/2, 0:ny-1), Ek(0:nx/2, 0:ny-1), Bk(0:nx/2, 0:ny-1)
    call coop_fft_forward(nx, ny, Emap, Ek)
    call coop_fft_forward(nx, ny, Bmap, Bk)
    do i=0, nx/2
       do j=0, ny-1
          kx = real(i, dl)/nx
          if(ny - j .lt. j)then
             ky = real(j - ny, dl)/ny
          else
             ky = real(j, dl)/ny
          endif
          k2 = kx**2+ky**2
          if(k2.eq.0.d0)then
             Qk(i, j) = (0.d0, 0.d0)
             Uk(i, j) = (0.d0, 0.d0)
          else
             Qk(i, j) = -((kx**2 - ky**2)*Ek(i,j) - 2.d0*kx*ky*Bk(i,j))/k2
             Uk(i, j) = -((kx**2 - ky**2)*Bk(i,j) + 2.d0*kx*ky*Ek(i,j))/k2
          endif
       enddo
    enddo
    call coop_fft_backward(nx, ny, Qk, Emap)
    call coop_fft_backward(nx, ny, Uk, Bmap)
  end subroutine coop_fits_EB2QU


  subroutine coop_fits_QU2EB(nx, ny, Qmap, Umap)
    COOP_INT nx, ny, i, j
    real(dl) Qmap(nx,ny), Umap(nx,ny), kx, ky, k2
    complex(dlc) Qk(0:nx/2, 0:ny-1), Uk(0:nx/2, 0:ny-1), Ek(0:nx/2, 0:ny-1), Bk(0:nx/2, 0:ny-1)
    call coop_fft_forward(nx, ny, Qmap, Qk)
    call coop_fft_forward(nx, ny, Umap, Uk)
    do i=0, nx/2
       do j=0, ny-1
          kx = real(i, dl)/nx
          if(ny - j .lt. j)then
             ky = real(j - ny, dl)/ny
          else
             ky = real(j, dl)/ny
          endif
          k2 = kx**2 + ky**2
          if(k2.eq.0.d0)then
             Ek(i, j) = 0.d0
             Bk(i, j) = 0.d0
          else
             Ek(i, j) = - ( (kx**2 - ky**2)*Qk(i,j) + 2.d0*kx*ky*Uk(i,j))/k2
             Bk(i, j) = - ( (kx**2 - ky**2)*Uk(i,j) - 2.d0*kx*ky*Qk(i,j))/k2
          endif
       enddo
    enddo
    call coop_fft_backward(nx, ny, Ek, Qmap)
    call coop_fft_backward(nx, ny, Bk, Umap)
  end subroutine coop_fits_QU2EB


  subroutine coop_fits_image_cea_EB2QU(this)
    class(coop_fits_image_cea)::this
    if(.not. allocated(this%smooth_Q).or. .not. allocated(this%smooth_U)) stop "EB2QU cannot be done: E, B not allocated yet."
    call coop_fits_EB2QU(this%smooth_nx*2+1, this%smooth_ny*2+1, this%smooth_Q, this%smooth_U)
  end subroutine coop_fits_image_cea_EB2QU


  subroutine coop_fits_image_cea_QU2EB(this)
    class(coop_fits_image_cea)::this
    if(.not. allocated(this%smooth_Q).or. .not. allocated(this%smooth_U)) stop "QU2EB cannot be done: Q, U not allocated yet."
    call coop_fits_QU2EB(this%smooth_nx*2+1, this%smooth_ny*2+1, this%smooth_Q, this%smooth_U)
  end subroutine coop_fits_image_cea_QU2EB




  subroutine coop_fits_image_cea_simulate(this, lmin, lmax, Cls)
    class(coop_fits_image_cea)::this
    COOP_INT lmin, lmax, i, j, ik
    COOP_REAL Cls(lmin:lmax)
    complex(dlc),dimension(:,:),allocatable::fk
    COOP_REAL amp(lmin:lmax), rk
    amp = sqrt(Cls/(this%pixsize**2/this%npix))
    allocate(fk(0:this%nside(1)/2,0:this%nside(2)-1))
    do i=0, this%nside(1)/2
       do j=0, this%nside(2)/2
          rk  = sqrt((this%dkx*i)**2 + (this%dky*j)**2) + 0.5d0
          ik = floor(rk)
          if(ik.ge.lmin .and. ik .lt. lmax)then
             rk = rk - ik
             fk(i,j) = coop_random_complex_gaussian() * (amp(ik)*(1.d0-rk) + amp(ik+1)*rk)
          else
             fk(i, j) = 0
          endif
       enddo
       do j=this%nside(2)/2+1, this%nside(2)-1
          rk  = sqrt((this%dkx*i)**2 + (this%dky*(this%nside(2)-j))**2) + 0.5d0
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
    COOP_UNKNOWN_STRING:: filename
    COOP_LONG_STRING::keys, values
    COOP_INT i, j, nkeys
    integer,dimension(:),allocatable::dtypes
    if(coop_file_exists(filename)) call coop_delete_file(filename)
    allocate(dtypes(this%header%n))
    this%filename = trim(filename)
    call coop_convert_to_c_string(this%filename)
    keys = ""
    values = ""
    do i=1, this%header%n
       dtypes(i) = coop_data_type(this%header%val(i))  
       keys = trim(keys)//trim(this%header%key(i))//coop_newline
       values = trim(values)//trim(coop_str_replace(this%header%val(i),"'", ""))//coop_newline
    enddo
    call coop_convert_to_c_string(keys)
    call coop_convert_to_c_string(values)
    nkeys = this%header%n
    call coop_fits_write_image(this%filename, nkeys, keys, values, dtypes, this%image, this%nside(1), this%nside(2))
    deallocate(dtypes)
  end subroutine coop_fits_image_cea_write

  subroutine coop_fits_image_cea_plot(this, figname)
    class(coop_fits_image_cea)::this
    COOP_UNKNOWN_STRING::figname
    type(coop_asy)::asy
    if(.not. allocated(this%smooth_image)) stop "call get_flatmap beore plot"
    call asy%open(figname)
    call asy%init(xlabel = "$2\sin{\frac{\theta}{2}}\cos\varphi$", ylabel = "$2\sin{\frac{\theta}{2}}\sin\varphi$", caption = "noise-free simulation", width=7., height=5.5)
    call coop_asy_density(asy, this%smooth_image, this%xmin, this%xmax, this%ymin, this%ymax, "$I (\mu K)$")
    call asy%close()
    if(allocated(this%smooth_Q))then
       call asy%open(coop_file_add_postfix(figname, "_Q"))
       call asy%init(xlabel = "$2\sin{\frac{\theta}{2}}\cos\varphi$", ylabel = "$2\sin{\frac{\theta}{2}}\sin\varphi$", caption = "noise-free simulation", width=7., height=5.5)
       call coop_asy_density(asy, this%smooth_Q, this%xmin, this%xmax, this%ymin, this%ymax, "$Q (\mu K)$")
       call asy%close()
    endif

    if(allocated(this%smooth_U))then
       call asy%open(coop_file_add_postfix(figname, "_U"))
       call asy%init(xlabel = "$2\sin{\frac{\theta}{2}}\cos\varphi$", ylabel = "$2\sin{\frac{\theta}{2}}\sin\varphi$", caption = "noise-free simulation", width=7., height=5.5)
       call coop_asy_density(asy, this%smooth_U, this%xmin, this%xmax, this%ymin, this%ymax, "$U (\mu K)$")
       call asy%close()
    endif

  end subroutine coop_fits_image_cea_plot

  subroutine  coop_fits_image_cea_convert2healpix(this, hp, imap, mask, lower, upper)  !!convert to healpix map, and mask the bright spots with 5 arcmin circles
    class(coop_fits_image_cea)::this
    type(coop_healpix_maps)::hp, mask
    type(coop_list_integer)::ps
    integer(8)::pix
    COOP_INT:: hpix, imap, i, listpix(0:10000), nlist
    COOP_REAL theta, phi
    COOP_SINGLE,optional::lower, upper
    COOP_SINGLE flower, fupper
#ifdef HAS_HEALPIX    
    if(mask%nside .ne. hp%nside .or. mask%ordering .ne. hp%ordering) stop "fits_image_cea_convert2healpix: mask and map must have the same nside and ordering"
    if(present(lower))then
       flower = lower
    else
       flower = -1.e30
    endif
    if(present(upper))then
       fupper = upper
    else
       fupper = 1.e30
    endif
    mask%map(:,1) = 0.
    hp%map(:, imap) = 0.
    
    call ps%init()
    do pix=0, this%npix-1
       call this%pix2ang(pix, theta, phi)
       call hp%ang2pix(theta, phi, hpix)
       if(this%image(pix).lt.flower .or. this%image(pix).gt.fupper)then
          mask%map(hpix, 1) = -1.
          hp%map(hpix,imap) = 0.
          call ps%push(hpix)
       elseif(mask%map(hpix, 1) .ge. 0.)then
          mask%map(hpix,1) = mask%map(hpix,1) + 1.
          hp%map(hpix,imap) = hp%map(hpix, imap) + this%image(pix)
       endif
    enddo
    where(mask%map(:,1).gt. 0.)
       hp%map(:, imap) = hp%map(:, imap)/mask%map(:,1)
       mask%map(:, 1) = 1.
    end where
    do i = 1, ps%n
       hpix = ps%element(i)
       call hp%query_disc(hpix, coop_SI_arcmin*5.d0, listpix, nlist)
       mask%map(listpix(0:nlist-1), 1) = 0.
       hp%map(listpix(0:nlist-1), :) = 0.
    enddo
    call ps%init()
#else
    stop "you need to install Healpix"
#endif    
  end subroutine coop_fits_image_cea_convert2healpix
  

end module coop_fitswrap_mod
