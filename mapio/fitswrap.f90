module coop_fitswrap_mod
  use coop_wrapper_utils
  use coop_sphere_mod
  implicit none

#include "constants.h"

  private

  public::coop_fits, coop_fits_image, coop_fits_image_cea

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
   contains
     procedure::regularize => coop_fits_image_regularize
     procedure::free => coop_fits_image_free
     procedure::get_linear_coordinates => coop_fits_image_get_linear_coordinates
     procedure::get_data => coop_fits_image_get_data
  end type coop_fits_image



  type, extends(coop_fits_image)::coop_fits_image_cea
     COOP_REAL smooth_pixsize, pixsize, smooth_dkx, smooth_dky
     real(dl) xmin, xmax, ymin, ymax
     COOP_INT:: smooth_nx, smooth_ny, smooth_npix
     type(coop_sphere_disc)::disc
     real(dl),dimension(:,:),allocatable:: smooth_image
     contains
       procedure::pix2ang => coop_fits_image_cea_pix2ang
       procedure::get_flatmap => coop_fits_image_cea_get_flatmap
       procedure::pix2flat => coop_fits_image_cea_pix2flat
       procedure::cut => coop_fits_image_cea_cut
       procedure::filter => coop_fits_image_cea_filter
       procedure::smooth_flat => coop_fits_image_cea_smooth_flat
       procedure::find_extrema => coop_fits_image_cea_find_extrema
       procedure::stack => coop_fits_image_cea_stack
       procedure::stack2fig => coop_fits_image_cea_stack2fig
       procedure::simulate => coop_fits_image_cea_simulate
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
  end subroutine coop_fits_open

  subroutine coop_fits_image_free(this)
    class(coop_fits_image) :: this
    if(allocated(this%image)) deallocate(this%image)
    if(allocated(this%transform))deallocate(this%transform)
    if(allocated(this%center))deallocate(this%center)
    if(allocated(this%nside))deallocate(this%nside)
    select type(this)
    class is(coop_fits_image_cea)
       if(allocated(this%smooth_image))deallocate(this%smooth_image)
    end select
  end subroutine coop_fits_image_free

  subroutine coop_fits_get_header(this)
    class(coop_fits)::this
    COOP_LONG_STRING::header
    integer nkeys, i, j, istart, iend
    COOP_REAL,dimension(:),allocatable::delta
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
       if(abs(abs(delta(1)/delta(2)) - 1.d0) .gt. 1.d-3)then
          print*, delta(1), delta(2)
          stop "pixel size is not the same in x and y directions"
       else
          this%pixsize = (abs(delta(1))+abs(delta(2)))/2.d0
       endif
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
    real(sp),dimension(:),allocatable::tmp
    select case(this%bitpix)
    case(-64)
       call coop_fits_get_float_data(this%filename, this%image, this%npix)
    case(-32)
       allocate(tmp(0:this%npix-1))
       call coop_fits_get_double_data(this%filename, tmp, this%npix)
       this%image = real(tmp, dl)
       deallocate(tmp)
    case default
       write(*,*) "Cannot load data for bitpix = "//trim(coop_num2str(this%bitpix))
       stop
    end select

  end subroutine coop_fits_image_get_data

  subroutine coop_fits_image_regularize(this, tail)
    class(coop_fits_image) this
    real(dl) upper, lower, tail
    call array_get_threshold_double(this%image, this%npix, 1.-tail, lower)
    call array_get_threshold_double(this%image, this%npix, tail, upper)
    where(this%image .lt. lower)
       this%image = lower
    end where
    where(this%image .gt. upper)
       this%image = upper
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
    phi = this%transform(1,1)* ix + this%transform(1,2)*iy
    theta = coop_pio2 + asin(this%transform(2, 1)*ix + this%transform(2,2)*iy)
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
    COOP_REAL k2min, k2max, k2, rn1, rn2, kmin, kmax, omega, j2, kstep, rk
    complex(dlc),dimension(:,:),allocatable::fk
    type(coop_file)::fkf
    COOP_INT nk
    COOP_REAL,dimension(:),allocatable::karr, Pkarr, warr
    allocate( fk(0:this%smooth_nx, 0:this%smooth_ny*2))
    call coop_fft_forward(this%smooth_nx*2+1, this%smooth_ny*2+1, this%smooth_image, fk)
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
             fk(i, j) = fk(i, j)*sin((sqrt(k2)- kmin)*omega)**2
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
    where(warr .ne. 0.d0)
       pkarr = pkarr/warr
    end where
    if(present(fk_file))then
       call fkf%open(trim(fk_file))
       do i=1, nk
          write(fkf%unit, "(2E16.7)")  karr(i)*coop_2pi/this%smooth_pixsize, pkarr(i)*this%smooth_pixsize**2/(this%smooth_nx*2+1.d0)/(this%smooth_ny*2+1.d0)
       enddo
       call fkf%close()
    endif
    call coop_fft_backward(this%smooth_nx*2+1, this%smooth_ny*2+1, fk, this%smooth_image)
    deallocate(fk)
    if(present(fk_file))deallocate(karr, Pkarr, warr)
  end subroutine coop_fits_image_cea_smooth_flat

  subroutine coop_fits_image_cea_find_extrema(this, spotfile, spottype, radius)
    class(coop_fits_image_cea)::this
    COOP_UNKNOWN_STRING::spotfile, spottype
    type(coop_file)::cf
    COOP_REAL,optional:: radius
    COOP_INT i, j, irad
    if(present(radius))then
       irad = ceiling(radius/this%smooth_pixsize)
    else
       irad = 1
    endif
    call cf%open(trim(spotfile), "w")
    select case(trim(spottype))
    case("hot", "Hot", "HOT", "h", "H")
       do i=-this%smooth_nx+irad, this%smooth_nx - irad
          do j = -this%smooth_ny+irad, this%smooth_ny - irad
             if(all(this%smooth_image(i, j) .ge. this%smooth_image(i-1:i+1, j-1:j+1)))then
                write(cf%unit, "(2I6, E16.7, I6 )") i, j, coop_2pi*coop_random_unit(), min(this%smooth_nx - i, this%smooth_nx + i, this%smooth_ny - j, this%smooth_ny + j)
             endif
          enddo
       enddo
    end select
    call cf%close()
  end subroutine coop_fits_image_cea_find_extrema

  subroutine coop_fits_image_cea_simulate(this, lmin, lmax, Cls)
    class(coop_fits_image_cea)::this
    COOP_INT lmin, lmax, i, j, ik
    COOP_REAL Cls(lmin:lmax)
    complex(dlc),dimension(:,:),allocatable::fk
    COOP_REAL amp(lmin:lmax), dkx, dky, rk
    amp = sqrt(Cls*(2.d0*this%smooth_nx+1)*(2.d0*this%smooth_ny+1)/this%smooth_pixsize**2)
    dkx = sqrt((i*coop_2pi/this%smooth_pixsize/(2.d0*this%smooth_nx+1.d0))**2+(j*coop_2pi/this%smooth_pixsize/(2.d0*this%smooth_ny+1.d0))**2)
    allocate(fk(0:this%smooth_nx,0:this%smooth_ny*2))
    do i=0, this%smooth_nx
       do j=0, this%smooth_ny
          rk  = sqrt((this%smooth_dkx*i)**2 + (this%smooth_dky*j)**2) + 0.5d0
          ik = floor(rk)
          if(ik.ge.lmin .and. ik .lt. lmax)then
             rk = rk - ik
             fk(i,j) = coop_random_complex_gaussian() * (amp(ik)*(1.d0-rk) + amp(ik+1)*rk)
          endif
       enddo
       do j=0, this%smooth_ny
          rk  = sqrt((this%smooth_dkx*i)**2 + (this%smooth_dky*(2*this%smooth_ny+1-j))**2) + 0.5d0
          ik = floor(rk)
          if(ik.ge.lmin .and. ik .lt. lmax)then
             rk = rk - ik
             fk(i,j) = coop_random_complex_gaussian() * (amp(ik)*(1.d0-rk) + amp(ik+1)*rk)
          endif
       enddo
    enddo
    deallocate(fk)

  end subroutine coop_fits_image_cea_simulate

  subroutine coop_fits_image_cea_stack(this, spotfile, nrad, stacked_image)
    class(coop_fits_image_cea)::this
    COOP_UNKNOWN_STRING::spotfile
    type(coop_file)::fp
    COOP_INT i, j, nrad, ii, jj, dis2b, irot, jrot
    COOP_REAL theta, radius, stacked_image(-nrad:nrad,-nrad:nrad)
    stacked_image = 0.d0
    call fp%open(trim(spotfile))
    do 
       read(fp%unit, *, ERR=100, END=100) i, j, theta, dis2b
       if(dis2b .lt. nrad) cycle
       do ii= -nrad, nrad
          do jj= -nrad, nrad 
             irot = nint(ii*cos(theta) + jj*sin(theta))
             jrot = nint(-ii*sin(theta) + jj*cos(theta))
             if(abs(irot).le.nrad .and. abs(jrot).le.nrad) &
                  stacked_image(irot, jrot) = stacked_image(irot, jrot) + this%smooth_image(i+ii, j+jj)  
          enddo
       enddo
    enddo
100 call fp%close()
  end subroutine coop_fits_image_cea_stack


  subroutine coop_fits_image_cea_stack2fig(this, spotfile, radius, fig)
    class(coop_fits_image_cea)::this
    COOP_REAL radius
    COOP_UNKNOWN_STRING::spotfile, fig
    COOP_REAL,dimension(:,:),allocatable::image
    type(coop_asy)::asy
    COOP_INT nrad
    nrad = max(1, floor(radius/this%smooth_pixsize))
    allocate(image(-nrad:nrad, -nrad:nrad))
    call this%stack(spotfile, nrad, image)
    call asy%open(trim(fig))
    call asy%init(xlabel = "$2\sin{\frac{\theta}{2}} \cos\varphi$", ylabel = "$2\sin{\frac{\theta}{2}} \sin\varphi$", width=7., height=5.5)
    call coop_asy_density(asy, image, -this%smooth_pixsize*nrad, this%smooth_pixsize*nrad, -this%smooth_pixsize*nrad, this%smooth_pixsize*nrad, label = "$I(\mu K)$")
    call asy%close()
    deallocate(image)
  end subroutine coop_fits_image_cea_stack2fig



end module coop_fitswrap_mod
