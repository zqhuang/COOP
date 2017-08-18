module coop_gstack_mod
  use coop_wrapper_utils
  use coop_sphere_mod
  use coop_fitsio_mod
  use coop_healpix_mod
  use coop_stacking_mod
  use coop_fitswrap_mod
  implicit none

#include "constants.h"


  !!general type of fits maps  
  type coop_general_stack_points
     type(coop_dictionary)::header
     COOP_STRING::fmt
     COOP_STRING::map, map_pt, map_nu, map_e, map_orient, map_sym,  point_type, orient, symmetry, output, mask, prefix
     COOP_REAL::fwhm, fwhm_pt, fwhm_nu, fwhm_e, fwhm_orient, fwhm_sym
     COOP_REAL::nu_min, nu_max, e_min, e_max
     COOP_REAL,dimension(:,:),allocatable::points
     logical::check_nu = .false.
     logical::check_e = .false.
     logical::check_orient = .false.
     logical::check_sym = .false.
     logical::xup = .false.
     logical::yup = .false.     
     COOP_REAL::randrot = coop_2pi
   contains
     procedure::init => coop_general_stack_points_init     
     procedure::free => coop_general_stack_points_free
     procedure::read => coop_general_stack_points_read
     procedure::open => coop_general_stack_points_read     
     procedure::write => coop_general_stack_points_write
  end type coop_general_stack_points

contains

    
  subroutine coop_general_stack_points_init(this, inifile)
    class(coop_general_stack_points)::this
    COOP_UNKNOWN_STRING::inifile
    type(coop_list_realarr)::pts
    COOP_INT::i, idot
    COOP_INT::nneigh, list(8)
    type(coop_healpix_maps)::hmap_pt, hmap_nu, hmap_e, hmap_orient, hmap_sym, hmask, zeros1, zeros2
    type(coop_fits_image_cea)::fmap_pt, fmap_nu, fmap_e, fmap_orient, fmap_sym, fmask
    COOP_REAL:: rms, imin, imax, e2min, e2max, e2, mean, summask
    COOP_REAL::arr(0:5)
    call this%free()
    call coop_load_dictionary(inifile, this%header)
    call coop_dictionary_lookup(this%header, "format", this%fmt, "HEALPIX")
    call coop_dictionary_lookup(this%header, "map", this%map)
    call coop_dictionary_lookup(this%header, "output", this%output)
    call coop_dictionary_lookup(this%header, "mask", this%mask)        
    idot = scan(this%map, ".", .true.)
    if(idot .eq. 0) call coop_return_error("map must be a fits file")
    if(this%map(idot+1:idot+4) .ne. "fits") call coop_return_error("map must be a fits file")
    if(.not. coop_file_exists( trim(this%map) )) call coop_return_error("map "//trim(this%map)//" does not exist.")
    this%prefix = this%map(1:idot-1)
    call coop_dictionary_lookup(this%header, "point_type", this%point_type, "MAX")
    call coop_dictionary_lookup(this%header, "orient", this%orient, "RANDOM")
    call coop_dictionary_lookup(this%header, "symmetry", this%symmetry, "SYMMETRIC")
    !!read the fwhm's
    call coop_dictionary_lookup(this%header, "fwhm", this%fwhm, 0.d0)
    this%fwhm = max(this%fwhm, 0.d0)
    call coop_dictionary_lookup(this%header, "fwhm_pt", this%fwhm_pt, this%fwhm)
    this%fwhm_pt = max(this%fwhm, this%fwhm_pt)
    call coop_dictionary_lookup(this%header, "fwhm_nu", this%fwhm_nu, this%fwhm)
    this%fwhm_nu = max(this%fwhm, this%fwhm_nu)
    call coop_dictionary_lookup(this%header, "fwhm_e", this%fwhm_e, this%fwhm)
    this%fwhm_e = max(this%fwhm_e, this%fwhm)
    call coop_dictionary_lookup(this%header, "fwhm_orient", this%fwhm_orient, this%fwhm)
    this%fwhm_orient = max(this%fwhm_orient, this%fwhm)    
    call coop_dictionary_lookup(this%header, "fwhm_sym", this%fwhm_sym, this%fwhm)
    this%fwhm_sym = max(this%fwhm_sym, this%fwhm)

    call coop_dictionary_lookup(this%header, "map_pt", this%map_pt, "")
    call coop_dictionary_lookup(this%header, "map_nu", this%map_nu,"")
    call coop_dictionary_lookup(this%header, "map_e", this%map_e,"")
    call coop_dictionary_lookup(this%header, "map_orient", this%map_orient,"")
    call coop_dictionary_lookup(this%header, "map_sym", this%map_sym,"")
    
    if(trim(this%map_pt).eq."")then
       select case(trim(this%point_type))
       case("MAX", "MIN", "ALL")
          this%map_pt = trim(this%prefix)//"_AMPLITUDE_fwhm"//COOP_FILESTR_OF(this%fwhm_pt)//"a.fits"
       case("SADDLE")
          this%map_pt = trim(this%prefix)//"_SADDLE_fwhm"//COOP_FILESTR_OF(this%fwhm_pt)//"a.fits"        
       case default
          call coop_return_error("point type "//trim(this%point_type)//" is not supported")
       end select
    endif

    



    call coop_dictionary_lookup(this%header, "nu_min", this%nu_min, -1.d30)
    call coop_dictionary_lookup(this%header, "nu_max", this%nu_max, 1.d30)
    if(this%nu_min .ge. this%nu_max)then
       call coop_return_error("setting nu_min >= nu_max is silly...")
    endif
    call coop_dictionary_lookup(this%header, "e_min", this%e_min, -1.d30)
    call coop_dictionary_lookup(this%header, "e_max", this%e_max, 1.d30)
    if(this%e_min .ge. this%e_max)then
       call coop_return_error("setting e_min >= e_max is silly...")
    endif
    
    this%check_nu = this%nu_min .gt. -9999.d0 .or. this%nu_max .lt. 9999.d0
    if(this%check_nu .and. trim(this%map_nu).eq."") this%map_nu = trim(this%prefix)//"_AMPLITUDE_fwhm"//COOP_FILESTR_OF(this%fwhm_nu)//"a.fits"

    this%check_e =  this%e_min .gt. -9999.d0 .or. this%e_max .lt. 9999.d0
    if(this%check_e .and. trim(this%map_e) .eq. "") this%map_e = trim(this%prefix)//"_ECCENTRICITY_fwhm"//COOP_FILESTR_OF(this%fwhm_e)//"a.fits"
    
    this%check_orient = .true.    
    select case(trim(this%orient))
    case("RANDOM")
       this%randrot = coop_2pi
       this%check_orient = .false.
    case("ORIGINAL")
       this%randrot = 0.d0
       this%check_orient = .false.
    case("HESSIAN")
       this%randrot = 0.d0
       if(trim(this%map_orient).eq."")this%map_orient = trim(this%prefix)//"_ORIENT_HESSIAN_fwhm"//COOP_FILESTR_OF(this%fwhm_orient)//"a.fits"       
    case("QU")
       this%randrot = 0.d0
       if(trim(this%map_orient).eq."")this%map_orient = trim(this%prefix)//"_ORIENT_QU_fwhm"//COOP_FILESTR_OF(this%fwhm_orient)//"a.fits"
    case default
       call coop_return_error("orient = "//trim(this%orient)//" is not supported")
    end select

    this%check_sym = .true.
    select case(trim(this%symmetry))
    case("SYMMETRIC")
       this%check_sym = .false.
       this%xup = .false.
       this%yup = .false.
    case("X_UP")
       this%xup = .true.
       this%yup = .false.
    case("Y_UP")
       this%xup = .false.
       this%yup = .true.      
    case("XY_UP")
       this%xup = .true.
       this%yup = .true.       
    case default
       call coop_return_error("symmetry = "//trim(this%symmetry)//" is not supported")
    end select
    if(this%check_sym .and. trim(this%map_sym).eq."") this%map_sym = trim(this%prefix)//"_PARITY_fwhm"//COOP_FILESTR_OF(this%fwhm_sym)//"a.fits"
    if(this%check_sym .and. this%fwhm_sym .lt. this%fwhm_pt)then
       write(*,*) "Warning: fwhm for parity selection < fwhm for point-type selection?"
    endif


    call pts%init()
    
    select case(trim(this%fmt))
    case("HEALPIX")
       !!load original map
       if(trim(this%point_type).eq."SADDLE")then
          call hmap_pt%read(this%map, nmaps_wanted=2, nmaps_to_read = 1)
       else
          call hmap_pt%read(this%map, nmaps_wanted=1)
       endif
       hmap_pt%spin = 0       
       call hmap_pt%map2alm(index_list=(/ 1 /) )
       call hmap_pt%convert2nested()
       !!load mask
       select case(trim(this%mask))
       case("", "NONE")
          call hmask%init(nside = hmap_pt%nside, nmaps = 1, genre="MASK", nested = .true.)
          hmask%map(:,1) = 1.
       case("AUTOMATIC")
          if(coop_file_exists(trim(this%prefix)//"_AUTOMASK.fits"))then
             call hmask%read(trim(this%prefix)//"_AUTOMASK.fits")
             call hmask%convert2nested()
          else
             call hmask%init(nside = hmap_pt%nside, nmaps = 1, genre="MASK", nested = .true.)
             do i=0, hmask%npix - 1
                call neighbours_nest(hmap_pt%nside, i, list, nneigh)
                if(all(hmap_pt%map(list(1:nneigh),1).eq.hmap_pt%map(i,1)))then
                   hmask%map(i,1)=0.
                else
                   hmask%map(i,1) = 1.
                endif
             enddo
             call hmask%write(trim(this%prefix)//"_AUTOMASK.fits")
          endif
       case default
          call hmask%read(this%mask)
          call hmask%convert2nested()
       end select
       !!load amplitude selection map
       if(this%check_nu)then
          write(*,*) "doing amplitude selection"
          if(coop_file_exists(this%map_nu))then
             call hmap_nu%read(this%map_nu)
             call hmap_nu%convert2nested()
          else             
             call hmap_nu%init(nside = hmap_pt%nside, nmaps = 3, spin = (/ 0, 1, 1 /), lmax = hmap_pt%lmax)
             hmap_nu%alm(:, :, 1) = hmap_pt%alm(:, :, 1)
             hmap_nu%alm(:, :, 2) = hmap_pt%alm(:, :, 1)
             hmap_nu%alm(:, :, 3) = 0.
             call hmap_nu%filter_alm(fwhm = sqrt(this%fwhm_nu**2 - this%fwhm**2)*coop_SI_arcmin, index_list = (/ 1 /) )
             call hmap_nu%filter_alm(fwhm = sqrt(this%fwhm_nu**2 - this%fwhm**2)*coop_SI_arcmin, lpower=1.d0, index_list = (/ 2 /) )             
             call hmap_nu%alm2map()
             call hmap_nu%convert2nested()
             call hmap_nu%write(this%map_nu)
          endif
          summask = sum(dble(hmask%map(:,1)))
          mean = sum(dble(hmap_nu%map(:,1)*hmask%map(:,1)))/summask
          rms = sqrt(sum(dble((hmap_nu%map(:,1)-mean)**2*hmask%map(:,1)))/summask)
          imin = rms*this%nu_min + mean
          imax = rms*this%nu_max + mean
       endif
       !!load eccentricity selection map
       if(this%check_e)then
          write(*,*) "doing eccentricity selection"          
          if(coop_file_exists(this%map_e))then
             call hmap_e%read(this%map_e)
             call hmap_e%convert2nested()             
          else
             call hmap_e%init(nside = hmap_pt%nside, nmaps = 3, spin = (/ 0, 2, 2 /), lmax = hmap_pt%lmax)
             hmap_e%alm(:, :, 1) = hmap_pt%alm(:, :, 1)             
             call hmap_e%filter_alm(fwhm = sqrt(this%fwhm_e**2 - this%fwhm**2)*coop_SI_arcmin, lpower=2.d0, index_list= (/ 1 /) )
             hmap_e%alm(:, :, 2) = hmap_e%alm(:, :, 1)
             hmap_e%alm(:, :, 3) = 0.
             call hmap_e%alm2map()
             call hmap_e%convert2nested()             
             call hmap_e%write(this%map_e)
          endif
       endif
       !!load orientation selection map
       if(this%check_orient)then
          write(*,*) "getting orientation"
          if( coop_file_exists(this%map_orient))then
             call hmap_orient%read(this%map_orient)
             call hmap_orient%convert2nested()
          else
             call hmap_orient%init(nside = hmap_pt%nside, nmaps = 2, spin = (/ 2, 2 /), lmax = hmap_pt%lmax)
             hmap_orient%alm(:,:,1) = hmap_pt%alm(:, :, 1)
             if(trim(this%orient).eq."HESSIAN")then
                call hmap_orient%filter_alm(fwhm = sqrt(this%fwhm_orient**2 - this%fwhm**2)*coop_SI_arcmin, lpower=2.d0, index_list= (/ 1/) )
             endif
             hmap_orient%alm(:, :, 2) = 0.             
             call hmap_orient%alm2map()
             call hmap_orient%convert2nested()
             call hmap_orient%write(this%map_orient)
          endif
       endif
       !!load parity selection map
       if(this%check_sym)then
          write(*,*) "doing parity selection"          
          if(coop_file_exists(this%map_sym))then
             call hmap_sym%read(this%map_sym)
             call hmap_sym%convert2nested()
          else
             call hmap_sym%init(nside = hmap_pt%nside, nmaps = 2, spin = (/ 1, 1 /), lmax=hmap_pt%lmax )
             hmap_sym%alm(:, :, 1) = hmap_pt%alm(:, :, 1)
             call hmap_sym%filter_alm(fwhm=sqrt(this%fwhm_sym**2 - this%fwhm**2)*coop_SI_arcmin, lpower=1.d0, index_list = (/ 1 /) )
             hmap_sym%alm(:, :, 2) = 0.
             call hmap_sym%alm2map()
             call hmap_sym%convert2nested()             
             call hmap_sym%write(this%map_sym)
          endif
       endif
       write(*,*) "scanning the map "//trim(this%map_pt)
       !!finally let's also smooth the point-type selection map, if necessary
       if(trim(this%point_type).eq."SADDLE")then
          call hmap_pt%filter_alm(fwhm = sqrt(this%fwhm_pt**2 - this%fwhm**2)*coop_SI_arcmin, lpower = 1.d0, index_list = (/ 1 /) )
          hmap_pt%spin = 1
          hmap_pt%alm(:,:,2)=0.
          call hmap_pt%alm2map()
          call hmap_pt%write(trim(this%map_pt))                 
       elseif(this%fwhm_pt .gt. this%fwhm)then
          call hmap_pt%filter_alm(fwhm = sqrt(this%fwhm_pt**2 - this%fwhm**2)*coop_SI_arcmin, index_list = (/ 1 /) )
          call hmap_pt%alm2map()
          call hmap_pt%write(trim(this%map_pt))                 
       endif
       call hmap_pt%convert2nested()


       !!make the selection
       select case(trim(this%point_type))
       case("MAX")
          do i=0, hmap_pt%npix-1
             if(hmask%map(i,1).ge. 0.5)then
                call neighbours_nest(hmap_pt%nside, i, list, nneigh)
                if(all(hmap_pt%map(list(1:nneigh), 1).lt.hmap_pt%map(i, 1)))then
                   call check_point()
                endif
             endif
          enddo
       case("MIN")
          do i=0, hmap_pt%npix-1
             if(hmask%map(i,1).ge. 0.5)then
                call neighbours_nest(hmap_pt%nside, i, list, nneigh)
                if(all(hmap_pt%map(list(1:nneigh), 1).gt.hmap_pt%map(i, 1)))then
                   call check_point()
                endif
             endif
          enddo
       case("SADDLE")
          call hmap_pt%zeros(1, zeros1, hmask)
          call hmap_pt%zeros(2, zeros2, hmask)
          do i=0, hmap_pt%npix-1
             if(zeros1%map(i, 1) .gt. 0.5  .and. zeros2%map(i, 1) .gt. 0.5  .and. hmask%map(i,1).gt.0.5)then
                call check_point()
             endif
          enddo
          call zeros1%free()
          call zeros2%free()
       case("ALL")
          do i=0, hmap_pt%npix-1
             if(hmask%map(i,1).gt.0.5)then
                call check_point()
             endif
          enddo
       end select
       !!free the maps
       call hmap_pt%free()
       call hmap_nu%free()
       call hmap_e%free()
       call hmap_orient%free()
       call hmap_sym%free()
    case("RA-DEC")
       stop "RA-DEC option has not been implemented yet."
    case default
       call coop_return_error("Format "//trim(this%fmt)//" is not supported.")
    end select
    write(*,*) "found "//COOP_STR_OF(pts%n)//" stacking points"
    allocate(this%points(6, pts%n))
    do i=1, pts%n
       this%points(:, i) = pts%element(i)
    enddo
    call pts%free()

  contains
    subroutine check_point()
      if(this%check_nu)then
         if(hmap_nu%map(i, 1) .lt. imin .or. hmap_nu%map(i, 1) .gt. imax)  return
      endif
      if(this%check_e)then
         e2 = (hmap_e%map(i,2)**2+hmap_e%map(i, 3)**2)/hmap_e%map(i,1)**2
         if( e2 .lt. e2min .or. e2 .gt. e2max) return
      endif
      call hmap_pt%pix2ang(i, arr(1), arr(2))
      if(this%check_orient)then
         arr(3) = COOP_POLAR_ANGLE(dble(hmap_orient%map(i,1)), dble(hmap_orient%map(i, 2)))/2.d0 + coop_rand01()*coop_pi
      else
         call random_number(arr(3))
         arr(3) = this%randrot * arr(3)
      endif
      if(this%check_sym)then
         if(hmap_sym%map(i, 2)*cos(arr(3))-hmap_sym%map(i,1)*sin(arr(3)) .gt. 0.d0 .and. this%xup)then
            arr(4) = -1.d0
         else
            arr(4) = 1.d0
         endif
         if(hmap_sym%map(i, 1)*cos(arr(3))+hmap_sym%map(i,2)*sin(arr(3)) .lt. 0.d0 .and. this%yup)then
            arr(5)= -1.d0
         else
            arr(5) = 1.d0
         endif
      else
         arr(4:5) = 1.d0
      endif
      arr(0) = dble(i)
      call pts%push( real(arr))                 
    end subroutine check_point
    
  end subroutine coop_general_stack_points_init


  
  subroutine coop_general_stack_points_free(this)
    class(coop_general_stack_points)::this
    call this%header%free()
    COOP_DEALLOC(this%points)
  end subroutine coop_general_stack_points_free

  subroutine coop_general_stack_points_read(this)
    class(coop_general_stack_points)::this
  end subroutine coop_general_stack_points_read

  subroutine coop_general_stack_points_write(this, filename)
    class(coop_general_stack_points)::this
    COOP_UNKNOWN_STRING,optional::filename
    if(present(filename))then
       call coop_fits_file_write_image_2d(this%points, filename, this%header)
    else
       call coop_fits_file_write_image_2d(this%points, trim(this%output), this%header)
    endif
  end subroutine coop_general_stack_points_write
  
end module coop_gstack_mod


