module coop_fitswrap_mod
  use coop_wrapper_utils
  use coop_sphere_mod
  use coop_healpix_mod
  implicit none

#include "constants.h"

  private

  public::coop_fits, coop_fits_image, coop_fits_image_cea, coop_fits_QU2EB, coop_fits_EB2QU, coop_fits_image_cea_simulate_TEB


  type coop_fits
     COOP_STRING::filename, fortran_filename
     type(coop_dictionary):: header
   contains
     procedure::open => coop_fits_open
     procedure::free => coop_fits_free
     procedure::key_value=>coop_fits_key_value
  end type coop_fits

  type, extends(coop_fits)::coop_fits_image
     COOP_INT::bitpix
     COOP_INT::dim
     COOP_INT,dimension(:),allocatable::nside
     COOP_LONG_INT::npix
     COOP_REAL,dimension(:),allocatable::image
     COOP_REAL,allocatable::transform(:, :), invtrans(:,:), center(:)     
   contains
     procedure::regularize => coop_fits_image_regularize
     procedure::get_linear_coordinates => coop_fits_image_get_linear_coordinates
     procedure::get_data => coop_fits_image_get_data
     procedure::simple_stat => coop_fits_image_simple_stat
  end type coop_fits_image

  type, extends(coop_fits_image)::coop_fits_image_cea
     COOP_REAL smooth_pixsize, pixsize, smooth_dkx, smooth_dky, dkx, dky
     COOP_REAL xmin, xmax, ymin, ymax
     COOP_INT:: smooth_nx, smooth_ny, smooth_npix
     type(coop_sphere_disc)::disc
     COOP_REAL,allocatable::radec_center(:)     
     COOP_REAL,dimension(:,:),allocatable:: smooth_image
     COOP_REAL,dimension(:,:),allocatable:: smooth_Q, smooth_U
   contains
     procedure::convert2healpix => coop_fits_image_cea_convert2healpix
     procedure::write => coop_fits_image_cea_write
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
     procedure::EB2QU => coop_fits_image_cea_EB2QU
     procedure::QU2EB => coop_fits_image_cea_QU2EB
  end type coop_fits_image_cea


contains

  subroutine coop_fits_open(this, filename)
    class(coop_fits)::this
    COOP_UNKNOWN_STRING::filename
    COOP_INT i
    call this%free()
    if(coop_file_exists(filename))then
       this%fortran_filename = trim(filename)
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

  subroutine coop_fits_free(this)
    class(coop_fits) :: this
    call this%header%free()
    select type(this)
    class is(coop_fits_image)
       if(allocated(this%image)) deallocate(this%image)
       if(allocated(this%nside))deallocate(this%nside)
       if(allocated(this%transform))deallocate(this%transform)
       if(allocated(this%invtrans))deallocate(this%invtrans)       
    end select
    select type(this)
    class is (coop_fits_image_cea)
       if(allocated(this%center))deallocate(this%center)
       if(allocated(this%radec_center))deallocate(this%radec_center)    
       if(allocated(this%smooth_image))deallocate(this%smooth_image)
       if(allocated(this%smooth_Q))deallocate(this%smooth_Q)
       if(allocated(this%smooth_U))deallocate(this%smooth_U)
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
          write(*,*) "Error: cannot find NAXIS key word in fits file "//trim(this%fortran_filename)
          stop
       endif
       allocate(this%nside(this%dim))
       this%npix = 1
       do i=1, this%dim
          this%nside(i) = nint(this%key_value("NAXIS"//COOP_STR_OF(i)))
          this%npix = this%npix * this%nside(i)
       enddo
       if(any(this%nside .eq. 0))then
          write(*,*) trim(this%fortran_filename)//": cannot read the dimensions"
          call this%header%print()
          stop
       endif
       allocate(this%image(0:this%npix-1))
       allocate(this%transform(this%dim, this%dim), this%center(this%dim), this%invtrans(this%dim, this%dim))
       this%transform = 0.d0
       this%invtrans = 0.d0       
       do i=1, this%dim
          this%transform(i, i) = 1.d0
          this%invtrans(i, i) = 1.d0          
       enddo
       this%center = 0.d0
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
             write(*,*) trim(this%fortran_filename)//": cannot read delta"
             call this%header%print()
             stop
          endif
          this%center(i) = this%key_value("CRPIX"//COOP_STR_OF(i), (this%nside(i)+1.d0)/2.d0)
          this%radec_center(i) = this%key_value("CRVAL"//COOP_STR_OF(i)) * units(i)
       enddo
       if(abs(abs(delta(1)/delta(2)) - 1.d0) .gt. 1.d-3)then
          print*, delta(1), delta(2)
          stop "pixel size is not the same in x and y directions"
       else
          this%pixsize = sqrt(abs(delta(1)*delta(2)))
       endif
       this%dkx = coop_2pi/(abs(delta(1))*this%nside(1))
       this%dky = coop_2pi/(abs(delta(2))*this%nside(2))
       do i=1, this%dim
          this%transform(:, i) =  this%transform(:, i) * delta(i)
       enddo
       deallocate(delta)
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


  subroutine coop_fits_image_simple_stat(this,mean,rms, median, fmin, fmax, lower, upper, mask, clsfile)
    class(coop_fits_image)::this
    COOP_UNKNOWN_STRING,optional::clsfile
    class(coop_fits_image),optional::mask
    COOP_REAL,optional::mean, rms, lower(3), upper(3), median, fmin, fmax
    COOP_REAL::ave, std, bounds(-3:3), minv, maxv, weight
    COOP_INT,parameter::n= 80
    COOP_REAL::Cls(n), ells(n)
    COOP_INT::i
    type(coop_file)::fp
    if(present(mask))then
       weight = count(mask%image .gt. 0.d0)
       ave = sum(this%image, mask=mask%image .gt. 0.d0)/weight
       std =  sqrt(sum((this%image-ave)**2, mask=mask%image .gt. 0.d0)/weight)
    else
       ave = sum(this%image)/this%npix
       std =  sqrt(sum((this%image-ave)**2/this%npix))
       weight = this%npix
    endif
    if(present(mean))mean =ave
    if(present(rms))rms = std
    write(*,*)
    write(*,*) "=========================================================="
    write(*,*) "Statistics of map :"//trim(this%fortran_filename)
    write(*,*) "=========================================================="
    write(*,*) "size: "//COOP_STR_OF(this%nside(1))//" x "//COOP_STR_OF(this%nside(2))
    write(*,*) "mean = ", ave
    write(*,*) "rms = ", std
    call array_get_mult_threshold_double(this%image, int(this%npix), COOP_STANDARD_SIGMA_BOUNDS, 7, bounds)
    write(*,*) "median = ", bounds(0)
    write(*,*) "1sigma lower, upper = ", bounds(-1), bounds(1)
    write(*,*) "2sigma lower, upper = ", bounds(-2), bounds(2)
    write(*,*) "3sigma lower, upper = ", bounds(-3), bounds(3)
    if(present(median))median = bounds(0)
    if(present(lower))lower = bounds(-1:-3:-1)
    if(present(upper))upper = bounds(1:3)
    minv = minval(this%image)
    maxv = maxval(this%image)
    write(*,*) "min max = ", minv, maxv 
    if(present(fmin))fmin = minv
    if(present(fmax))fmax = maxv
    write(*,*) "zero-value pixels: "//trim(coop_num2str(100.*count(this%image .eq. 0.d0)/dble(this%npix),"(F10.3)"))//"%"
    
    if(present(clsfile))then
       if(trim(clsfile).ne."")then
          call fp%open(clsfile)
          select type(this)
             class is(coop_fits_image_cea)
             weight = weight/this%npix
             call coop_set_uniform(n, ells, max(this%dkx, this%dky)*3.d0, min(this%dkx*this%nside(1), this%dky*this%nside(2), 3.d4)/10.d0)
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
    COOP_SINGLE,dimension(:),allocatable::tmp
    select case(this%bitpix)
    case(-64)
       call coop_fits_get_double_data(this%filename, this%image, this%npix)
    case(-32)
       allocate(tmp(0:this%npix-1))
       call coop_fits_get_float_data(this%filename, tmp, this%npix)
       this%image = COOP_REAL_OF(tmp)
       deallocate(tmp)
    case default
       write(*,*) "Cannot load data for bitpix = "//trim(coop_num2str(this%bitpix))
       stop
    end select

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
    COOP_LONG_INT pix
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
    pix = ivec(1) + this%center(1)  - 1 + this%nside(1) * (ivec(2) + this%center(2) - 1)
    if(pix .ge. this%npix .or. pix .lt. 0)then
       write(*,*) "warning: radec2pix overflow", pix, this%npix
    endif
  end subroutine coop_fits_image_cea_radec2pix

  
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

  subroutine coop_fits_image_cea_smooth(this, fwhm, highpass_l1, highpass_l2, lowpass_l1, lowpass_l2)
    class(coop_fits_image_cea)::this
    COOP_REAL::fwhm
    COOP_INT::highpass_l1, highpass_l2, lowpass_l1, lowpass_l2
    COOP_REAL::hp_omega, lp_omega, sigma
    sigma = coop_sigma_by_fwhm*fwhm/coop_sqrt2
    hp_omega = coop_pio2/(highpass_l2 - highpass_l1)
    lp_omega = coop_pio2/(lowpass_l2 - lowpass_l1)
    call this%filter(lmin = highpass_l1, lmax = lowpass_l2,  window = smooth_weight)
  contains

    function smooth_weight(k)
      COOP_REAL k, smooth_weight
      smooth_weight = exp(-(k*sigma)**2)
      if(k .lt. highpass_l2)then
         smooth_weight = smooth_weight * sin((k-highpass_l1)*hp_omega)**2
      elseif(k .gt. lowpass_l1)then
         smooth_weight = smooth_weight * sin((lowpass_l2 - k)*lp_omega)**2
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
    Cls = Cls*(this%pixsize**2)**2
  end subroutine coop_fits_image_cea_get_power

  subroutine coop_fits_image_cea_filter(this, lmin, lmax, window)
    class(coop_fits_image_cea)::this
    COOP_INT lmin, lmax, i, j
    external window
    COOP_REAL window
    COOP_REAL k2min, k2max, k2, ky2
    COOP_COMPLEX,dimension(:,:),allocatable::fk
    allocate( fk(0:this%nside(1)/2, 0:this%nside(2)-1))
    call coop_fft_forward(this%nside(1), this%nside(2), this%image, fk)
    k2min = dble(lmin)**2
    k2max = dble(lmax)**2
    do j=0, this%nside(2)-1
       ky2 = (min(j, this%nside(2)-j)*this%dky)**2
       if( ky2 .gt. k2max)then
          fk(:, j) = ( 0.d0, 0.d0 )
          cycle
       endif
       do i = 0, this%nside(1)/2
          k2 = (i*this%dkx)**2 + ky2
          if(k2 .lt. k2min .or. k2.gt. k2max)then
             fk(i, j) = ( 0.d0, 0.d0 )
          else
             fk(i, j) = fk(i, j)*window(sqrt(k2))
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

  subroutine  coop_fits_image_cea_convert2healpix(this, hp, imap, mask, hits, lower, upper)  !!convert to healpix map, and mask the bright spots with 5 arcmin circles
    class(coop_fits_image_cea)::this
    type(coop_fits_image_cea),optional::hits !!if this presents, cut off points with hits less than 30% of the mean hits
    type(coop_healpix_maps)::hp, mask
    type(coop_list_integer)::ps
    integer(8)::pix
    COOP_INT:: hpix, imap, i, listpix(0:10000), nlist
    COOP_REAL theta, phi, meanhits
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
       if(present(hits))then
          if(hits%image(pix) .le. 0.)cycle
       endif
       call hp%ang2pix(theta, phi, hpix)
       if(this%image(pix).lt.flower .or. this%image(pix).gt.fupper)then
          mask%map(hpix, 1) = -1.
          hp%map(hpix,imap) = 0.
          call ps%push(hpix)
       elseif(mask%map(hpix, 1) .ge. 0.)then
          if(present(hits))then
             mask%map(hpix, 1) = mask%map(hpix,1)+hits%image(pix)
             hp%map(hpix,imap) = hp%map(hpix, imap) + this%image(pix)*hits%image(pix)
          else
             mask%map(hpix,1) = mask%map(hpix,1) + 1.
             hp%map(hpix,imap) = hp%map(hpix, imap) + this%image(pix)
          endif
       endif
    enddo
    meanhits = sum(mask%map(:,1), mask= mask%map(:,1).gt.0.)/count(mask%map(:,1).gt.0.)*0.25  !!discard pixels with small obs time
    where(mask%map(:,1).gt. meanhits)
       hp%map(:, imap) = hp%map(:, imap)/mask%map(:,1)
       mask%map(:, 1) = 1.
    elsewhere
       hp%map(:,imap) = 0.
       mask%map(:, 1) = 0.
    end where
    do i = 1, ps%n
       hpix = ps%element(i)
       call hp%query_disc(hpix, coop_SI_arcmin*5.d0, listpix, nlist)
       mask%map(listpix(0:nlist-1), 1) = 0.
       hp%map(listpix(0:nlist-1), imap) = 0.
    enddo
    call ps%init()
#else
    stop "you need to install Healpix"
#endif    
  end subroutine coop_fits_image_cea_convert2healpix



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
    cls_sqrteig = cls_sqrteig/tmap%pixsize**2
    allocate(fk(0:tmap%nside(1)/2,0:tmap%nside(2)-1,3))
    do j=0, tmap%nside(2)-1
       ky2 = (tmap%dky*min(j, tmap%nside(2)-j))**2
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
    call coop_fft_backward(tmap%nside(1), tmap%nside(2), fk(:,:,2), emap%image)
    call coop_fft_backward(tmap%nside(1), tmap%nside(2), fk(:,:,3), bmap%image)
    if(present(Qmap) .and. present(umap))then
       do j=0, tmap%nside(2)-1
          if(j .lt. tmap%nside(2)- j)then
             ky = (j-tmap%nside(2))*tmap%dky
          else
             ky = j*tmap%dky
          endif
          ky2 = ky**2
          if(ky2.ge.k2max) cycle
          do i=0, tmap%nside(1)/2
             kx = i*tmap%dkx
             kx2 = kx**2
             k2 = kx2+ky2
             if(k2 .le. k2min .or. k2.ge.k2max)cycle
             fk(i, j,1) = -((kx2 - ky2)*fk(i,j,2) - 2.d0*kx*ky*fk(i,j,3))/k2
             fk(i, j,2) = -((kx2 - ky2)*fk(i,j,3) + 2.d0*kx*ky*fk(i,j,2))/k2
          enddo
       enddo
       call coop_fft_backward(tmap%nside(1), tmap%nside(2), fk(:,:,1), qmap%image)
       call coop_fft_backward(tmap%nside(1), tmap%nside(2), fk(:,:,2), umap%image)
    endif
    deallocate(fk)
  end subroutine coop_fits_image_cea_simulate_TEB





end module coop_fitswrap_mod
