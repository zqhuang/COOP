module coop_stacking_mod
  use coop_wrapper_firstorder
#ifdef HAS_HEALPIX
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  use udgrade_nr
#endif
  implicit none
#include "constants.h"

  COOP_REAL,parameter::coop_stacking_max_threshold = 10.d0
  COOP_INT,parameter::coop_stacking_genre_Null = 0
  COOP_INT,parameter::coop_stacking_genre_Imax = 1
  COOP_INT,parameter::coop_stacking_genre_Imin = 2
  COOP_INT,parameter::coop_stacking_genre_Imax_Oriented = 3
  COOP_INT,parameter::coop_stacking_genre_Imin_Oriented = 4
  COOP_INT,parameter::coop_stacking_genre_Lmax = 5
  COOP_INT,parameter::coop_stacking_genre_Lmin = 6
  COOP_INT,parameter::coop_stacking_genre_Lmax_Oriented = 7
  COOP_INT,parameter::coop_stacking_genre_Lmin_Oriented = 8
  COOP_INT,parameter::coop_stacking_genre_Pmax_Oriented = 9
  COOP_INT,parameter::coop_stacking_genre_Pmin_Oriented = 10


  type coop_stacking_options
     COOP_INT::genre = coop_stacking_genre_Null
     COOP_INT::nmaps = 0
     COOP_INT::nside = 0
     COOP_INT::index_I = 0
     COOP_INT::index_Q = 0
     COOP_INT::index_U = 0
     COOP_INT::index_L = 0
     COOP_STRING::caption
     COOP_SINGLE::I_lower = -1.e20
     COOP_SINGLE::I_upper = 1.e20
     COOP_SINGLE::L_lower = -1.e20
     COOP_SINGLE::L_upper = 1.e20
     COOP_SINGLE::P2_lower = -1.e20
     COOP_SINGLE::P2_upper = 1.e20
     COOP_SINGLE::I_lower_nu = -1.e20
     COOP_SINGLE::I_upper_nu = 1.e20
     COOP_SINGLE::L_lower_nu = -1.e20
     COOP_SINGLE::L_upper_nu = 1.e20
     COOP_SINGLE::P_lower_nu = -1.e20
     COOP_SINGLE::P_upper_nu = 1.e20
     COOP_SINGLE::P2byI2_lower = -1.e20
     COOP_SINGLE::P2byI2_upper = 1.e20
     COOP_SINGLE::P2byL2_lower = -1.e20
     COOP_SINGLE::P2byL2_upper = 1.e20
     COOP_INT::threshold_option = 0
     logical::addpi = .true.  !!randomly add pi on polarization directions
     logical::nested = .false.
     type(coop_list_integer)::peak_pix
     type(coop_list_realarr)::peak_ang
     type(coop_list_realarr)::peak_map
   contains
     procedure::export => coop_stacking_options_export
     procedure::import => coop_stacking_options_import     
     procedure::init => coop_stacking_options_init
     procedure::free => coop_stacking_options_free
     procedure::reject=>coop_stacking_options_reject
     procedure::pix => coop_stacking_options_pix
     procedure::export_pix => coop_stacking_options_export_pix
     procedure::export_ang => coop_stacking_options_export_ang     
     procedure::rotate_angle => coop_stacking_options_rotate_angle
     procedure::peak_r => coop_stacking_options_peak_r
     procedure::peak_e => coop_stacking_options_peak_e
     procedure::peak_get_angle_r_e => coop_stacking_options_peak_get_angle_r_e
  end type coop_stacking_options


  type coop_to_be_stacked
     COOP_INT::nmaps = 0
     COOP_INT,dimension(:),allocatable::ind
     COOP_INT,dimension(:),allocatable::spin
     COOP_STRING,dimension(:),allocatable::label
     logical, dimension(:),allocatable::headless_vector !!headless vector
     logical, dimension(:),allocatable::local_rotation !!local rotation
     COOP_REAL, dimension(:),allocatable::zmin !!lowerbound
     COOP_REAL, dimension(:),allocatable::zmax !!upperbound
   contains     
     procedure::init => coop_to_be_stacked_init
     procedure::free => coop_to_be_stacked_free
  end type coop_to_be_stacked


contains

  
  subroutine coop_stacking_options_export(this, filename)
    class(coop_stacking_options)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file)fp
    COOP_INT i
    call fp%open(filename, "u")
    write(fp%unit) this%genre, this%nmaps, this%nside, this%index_I, this%index_Q, this%index_U, this%index_L
    write(fp%unit) this%I_lower, this%I_upper, this%L_lower, this%L_upper, this%P2_lower, this%P2_upper, this%I_lower_nu, this%I_upper_nu, this%L_lower_nu, this%L_upper_nu, this%P_lower_nu, this%P_upper_nu, this%P2byI2_lower, this%P2byI2_upper, this%P2byL2_lower, this%P2byL2_upper
    write(fp%unit) this%genre
    write(fp%unit) this%threshold_option, this%addpi, this%nested
    do i=1, this%peak_pix%n
       write(fp%unit) this%peak_pix%element(i), this%peak_ang%element(i), this%peak_map%element(i)
    enddo
    call fp%close()
  end subroutine coop_stacking_options_export

  subroutine coop_stacking_options_import(this, filename)
    class(coop_stacking_options)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file) fp
    COOP_INT pix
    COOP_SINGLE::thetaphi(2)
    COOP_SINGLE,dimension(:),allocatable::map
    if(.not. coop_file_exists(filename))then
       write(*,*) "stack option file "//trim(adjustl(filename))//" does not exist"
       stop
    endif
    call this%free()
    call fp%open(filename, "ur")
    read(fp%unit) this%genre, this%nmaps, this%nside, this%index_I, this%index_Q, this%index_U, this%index_L
    read(fp%unit) this%I_lower, this%I_upper, this%L_lower, this%L_upper, this%P2_lower, this%P2_upper, this%I_lower_nu, this%I_upper_nu, this%L_lower_nu, this%L_upper_nu, this%P_lower_nu, this%P_upper_nu, this%P2byI2_lower, this%P2byI2_upper, this%P2byL2_lower, this%P2byL2_upper
    read(fp%unit) this%genre
    read(fp%unit) this%threshold_option, this%addpi, this%nested
    allocate(map(this%nmaps))
    do i=1, this%peak_pix%n
       read(fp%unit) pix, thetaphi, map
       call this%peak_pix%push(pix)
       call this%peak_ang%push(thetaphi)
       call this%peak_map%push(map)
    enddo
    call fp%close()
    deallocate(map)
  end subroutine coop_stacking_options_import

  subroutine coop_stacking_options_free(this)
    class(coop_stacking_options)::this
    call this%peak_pix%init()
    call this%peak_ang%init()
    call this%peak_map%init()
  end subroutine coop_stacking_options_free

  subroutine coop_to_be_stacked_free(this, nmaps)
    class(coop_to_be_stacked)::this
    COOP_INT,optional::nmaps
    if(allocated(this%ind))deallocate(this%ind)
    if(allocated(this%spin))deallocate(this%spin)
    if(allocated(this%label))deallocate(this%label)
    if(allocated(this%zmin))deallocate(this%zmin)        
    if(allocated(this%zmax))deallocate(this%zmax)    
    if(allocated(this%headless_vector)) deallocate(this%headless_vector)
    if(allocated(this%local_rotation)) deallocate(this%local_rotation)
    if(present(nmaps))then
       this%nmaps = nmaps
       allocate(this%ind(this%nmaps), this%spin(this%nmaps), this%label(this%nmaps), this%local_rotation(this%nmaps), this%headless_vector(this%nmaps), this%zmin(this%nmaps), this%zmax(this%nmaps))       
    else
       this%nmaps = 0
    endif
    this%zmin = 1.1e31
    this%zmax = -1.1e31
    this%headless_vector = .false.
    this%local_rotation = .false.
    this%label = ""
    this%spin = 0
    this%ind = 1
  end subroutine coop_to_be_stacked_free

  subroutine coop_to_be_stacked_init(this, str)
    class(coop_to_be_stacked)::this
    COOP_UNKNOWN_STRING::str
    type(coop_list_string)::l, subl
    COOP_STRING::line
    COOP_INT i, j
    select case(trim(adjustl(str)))
    case("I", "T", "E", "B")
       call this%free(1)
       this%label = "$"//trim(adjustl(str))//"(\mu K)$"
    case("zeta", "Z", "\zeta")
       call this%free(1)
       this%label = "$10^5\zeta$"
    case("QU")
       call this%free(2)
       this%ind = (/ 1, 2 /)       
       this%spin = 2
       this%label(1) =  "$Q(\mu K)$"
       this%label(2) =  "$U(\mu K)$"
       this%headless_vector = .true.
    case("QrUr")
       call this%free(2)
       this%ind = (/ 1, 2 /)
       this%spin = 2
       this%label(1) =  "$Q_r(\mu K)$"
       this%label(2) =  "$U_r(\mu K)$"
       this%headless_vector = (/ .true., .false. /)
       this%local_rotation = .true.
    case("QTUT")
       call this%free(2)
       this%ind = (/ 1, 2 /)
       this%spin = 2
       this%label(1) =  "$Q_T(\mu K)$"
       this%label(2) =  "$U_T(\mu K)$"
       this%headless_vector = .true.
    case("QLTULT")
       call this%free(2)
       this%ind = (/ 1, 2 /)
       this%spin = 2
       this%label(1) =  "$Q_{\nabla^2 T}(\mu K/\mathrm{rad}^2)$"
       this%label(2) =  "$U_{\nabla^2 T}(\mu K/\mathrm{rad}^2)$"
       this%headless_vector = .true.
    case default  !! label1@index1@spin1@headless_vector1@local_roataion1@zmin1@zmax1::label2:index2:spin2@headless_vector2@local_roataion2@zmin1@zmax1
       call coop_string_to_list(str, l, "::")
       call this%free(l%n)
       do i = 1, l%n
          call coop_string_to_list(l%element(i), subl, "@")
          this%label(i) = trim(subl%element(1))
          if(subl%n .ge. 2)then
             call subl%get_element(2, line)
             read(line, *) this%ind(i)
             if(subl%n .ge. 3)then
                call subl%get_element(3, line)
                read(line, *) this%spin(i)
                if(subl%n .ge. 4)then
                   this%headless_vector(i) = (trim(subl%element(4)) .eq. "T" .and. this%spin(i) .eq. 2)
                   if(subl%n .ge. 5)then
                      this%local_rotation(i) = (trim(subl%element(5)) .eq. "T" .and. this%spin(i) .eq. 2)
                      if(subl%n .ge. 6)then
                         call subl%get_element(6, line)
                         read(line, *) this%zmin(i)
                         if(subl%n .ge. 7)then
                            call subl%get_element(7, line)
                            read(line, *) this%zmax(i)
                         endif
                      endif
                   endif
                endif                
             endif
          endif
       enddo
       call l%init()
       call subl%init()
    end select
    !!check spin pairs
    i = 1
    do 
       if(this%spin(i) .ne. 0)then
          if( i .ge. this%nmaps) stop "to_be_stacked_init: nonzero spin must go in pairs"
          if( this%spin(i+1) .ne. this%spin(i) .or. this%ind(i+1).ne. this%ind(i)+1)stop "to_be_stacked_init: nonzero spin must go in pairs"
          i = i + 2
       else
          i = i + 1
       endif
    enddo
  end subroutine coop_to_be_stacked_init

  subroutine coop_stacking_options_init(this, domax, peak_name, Orient_name, nmaps)
    class(coop_stacking_options)::this
    COOP_UNKNOWN_STRING::peak_name, orient_name
    COOP_SHORT_STRING::p
    logical domax
    COOP_INT::nmaps
    call this%free()
    p = trim(adjustl(peak_name))
    if(trim(adjustl(orient_name)).eq. "NULL" .or. trim(adjustl(orient_name)) .eq. "random")then
       if(domax)then
          if(len_trim(p).gt. 9)then
             if(p(1:9).eq."$\nabla^2")then
                this%genre = coop_stacking_genre_Lmax
             else
                this%genre = coop_stacking_genre_Imax
             endif
          else
             this%genre = coop_stacking_genre_Imax
          endif
          this%caption = trim(p)//" maxima"
       else
          if(len_trim(p) .gt. 9)then
             if(p(1:9) .eq. "$\nabla^2")then
                this%genre = coop_stacking_genre_Lmin
             else
                this%genre = coop_stacking_genre_Imin
             endif
          else
             this%genre = coop_stacking_genre_Imin
          endif
          this%caption = trim(p)//" minima"          
       endif
    else
       if(domax)then
          if(p(1:2) .eq. "$P")then
             this%genre = coop_stacking_genre_Pmax_Oriented             
          else
             if(len_trim(p).gt. 9)then
                if(p(1:9).eq."$\nabla^2")then
                   this%genre = coop_stacking_genre_Lmax_Oriented
                else
                   this%genre = coop_stacking_genre_Imax_Oriented               
                endif
             else
                this%genre = coop_stacking_genre_Imax_Oriented    
             endif
          endif
          this%caption = trim(p)//" maxima, "//trim(adjustl(Orient_name))//" oriented"          
       else
          if(p(1:2).eq."$P")then
             this%genre = coop_stacking_genre_Pmin_Oriented             
          else
             if(len_trim(p) .gt. 9)then
                if(p(1:9) .eq. "$\nabla^2")then
                   this%genre = coop_stacking_genre_Lmin_Oriented
                else
                   this%genre = coop_stacking_genre_Imin_Oriented
                endif
             else
                this%genre = coop_stacking_genre_Imin_Oriented
             endif
          endif
          this%caption = trim(adjustl(peak_name))//" minima, "//trim(adjustl(Orient_name))//" oriented"                    
       endif
    endif
    this%nmaps = nmaps
    !!now set default indices
    select case(nmaps)
    case(1)
       select case(this%genre)
       case(coop_stacking_genre_Imax, coop_stacking_genre_Imin)
          this%index_I = 1
          this%threshold_option = 4          
       case(coop_stacking_genre_Lmax, coop_stacking_genre_Lmin)
          this%index_L = 1
          this%threshold_option = 1
       case default
          stop "nmaps does not match stacking option"
       end select
    case(2)
       select case(this%genre)
       case(coop_stacking_genre_Pmax_Oriented, coop_stacking_genre_Pmin_Oriented)
          this%index_Q = 1
          this%index_U = 2
          this%threshold_option = 2          
       case(coop_stacking_genre_Imax, coop_stacking_genre_Imin, coop_stacking_genre_Lmax, coop_stacking_genre_Lmin)
          this%index_I = 1
          this%index_L = 2
          this%threshold_option = 5          
       case default
          stop "nmaps does not match stacking option"          
       end select
    case(3)
       select case(this%genre)
       case(coop_stacking_genre_Pmax_Oriented, coop_stacking_genre_Pmin_Oriented, coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Imin_Oriented, coop_stacking_genre_Imax, coop_stacking_genre_Imin)
          this%index_I = 1
          this%index_Q = 2
          this%index_U = 3
          this%threshold_option = 6          
       case(coop_stacking_genre_Lmax, coop_stacking_genre_Lmin, coop_stacking_genre_Lmax_Oriented, coop_stacking_genre_Lmin_Oriented)  
          this%index_L = 1
          this%index_Q = 2
          this%index_U = 3
          this%threshold_option = 3
       end select
    case(4:)
       this%index_I = 1
       this%index_Q = 2
       this%index_U = 3
       this%index_L = 4
       this%threshold_option = 7       
    end select
  end subroutine coop_stacking_options_init


  function coop_stacking_options_reject(this, map) result(rej)
    class(coop_stacking_options)::this
    logical rej
    COOP_SINGLE::map(:)
    COOP_REAL::p2
    select case(this%threshold_option)
    case(7)
       P2 = map(this%index_Q)**2+map(this%index_U)**2 
       rej = (map(this%index_I) > this%I_upper .or. map(this%index_I) < this%I_lower .or. P2> this%P2_upper .or. P2 < this%P2_lower .or.  map(this%index_L) > this%L_upper .or. map(this%index_L) < this%L_lower .or. P2 > this%P2byI2_upper*map(this%index_I)**2 .or. P2 < this%P2byI2_lower * map(this%index_I)**2 .or. P2 > this%P2byL2_upper*map(this%index_L)**2 .or. P2 < this%P2byL2_lower*map(this%index_L)**2)
    case(6)
       P2 = map(this%index_Q)**2+map(this%index_U)**2        
       rej = (map(this%index_I) > this%I_upper .or. map(this%index_I) < this%I_lower .or. P2 > this%P2_upper .or. P2 < this%P2_lower .or. P2 > this%P2byI2_upper*map(this%index_I)**2 .or. P2 < this%P2byI2_lower * map(this%index_I)**2)
    case(5)
       rej = ( map(this%index_I) > this%I_upper .or. map(this%index_I) < this%I_lower  .or.  map(this%index_L) > this%L_upper .or. map(this%index_L) < this%L_lower)
    case(4)
       rej = (map(this%index_I) > this%I_upper .or. map(this%index_I) < this%I_lower)
    case(3)
       P2 = map(this%index_Q)**2+map(this%index_U)**2               
       rej = ( P2 > this%P2_upper .or. P2 < this%P2_lower .or.  map(this%index_L) > this%L_upper .or. map(this%index_L) < this%L_lower .or. P2 > this%P2byL2_upper*map(this%index_L)**2 .or. P2 < this%P2byL2_lower*map(this%index_L)**2)
    case(2)
       P2 = map(this%index_Q)**2+map(this%index_U)**2                      
       rej = (P2 > this%P2_upper .or. P2 < this%P2_lower)
    case(1)
       rej = ( map(this%index_L) > this%L_upper .or.  map(this%index_L) < this%L_lower)
    case default
       rej  = .false.
    end select
  end function coop_stacking_options_reject

  function coop_stacking_options_pix(this, nside, i) result(pix)
    class(coop_stacking_options)::this
    COOP_INT nside, i, pix
    COOP_SINGLE::thetaphi(2)
#ifdef HAS_HEALPIX    
    if(this%nside .eq. nside)then
       pix = this%peak_pix%element(i)
       return
    endif
    call this%peak_ang%get_element(i, thetaphi)
    if(this%nested)then
       call ang2pix_nest(nside, dble(thetaphi(1)), dble(thetaphi(2)), pix)
    else
       call ang2pix_ring(nside, dble(thetaphi(1)), dble(thetaphi(2)), pix)
    endif
#endif    
  end function coop_stacking_options_pix


  subroutine coop_stacking_options_export_pix(this, filename)
    class(coop_stacking_options)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file)::fp
    COOP_INT :: i
    call fp%open(filename, "w")
    do i = 1, this%peak_pix%n
       write(fp%unit, "(I8)") this%peak_pix%element(i)
    enddo
    call fp%close()
  end subroutine coop_stacking_options_export_pix

  subroutine coop_stacking_options_export_ang(this, filename)
    class(coop_stacking_options)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file)::fp
    COOP_INT :: i
    logical dor
    call fp%open(filename, "w")
    if(this%addpi)then
       do i = 1, this%peak_ang%n
          write(fp%unit, "(3E14.5)") this%peak_ang%element(i),  this%rotate_angle(i) + coop_rand01()*coop_pi
       enddo       
    else
       do i = 1, this%peak_ang%n
          write(fp%unit, "(3E14.5)") this%peak_ang%element(i),  this%rotate_angle(i)
       enddo
    endif
    call fp%close()
  end subroutine coop_stacking_options_export_ang
  

  function coop_stacking_options_rotate_angle(this, i) result(angle)
    class(coop_stacking_options)::this
    COOP_REAL angle
    COOP_SINGLE::map(this%nmaps)
    COOP_INT i
    if(i.gt. this%peak_pix%n) stop "rotate_angle: pix overflow"
    select case(this%genre)
    case(coop_stacking_genre_Imax, coop_stacking_genre_Imin, coop_stacking_genre_Lmax, coop_stacking_genre_Lmin, coop_stacking_genre_Null)
       call random_number(angle)
       angle = angle*coop_2pi
    case(coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Imin_Oriented, coop_stacking_genre_Lmax_Oriented, coop_stacking_genre_Lmin_Oriented, coop_stacking_genre_Pmax_Oriented, coop_stacking_genre_Pmin_Oriented)
       call this%peak_map%get_element(i, map)
       angle = COOP_POLAR_ANGLE(map(this%index_Q), map(this%index_U))/2.d0
    case default
       stop "rotate_angle: unknown genre"
    end select
  end function coop_stacking_options_rotate_angle

  function coop_stacking_options_peak_r(this, i) result(r)
    COOP_REAL, parameter::max_radius = coop_pio2
    class(coop_stacking_options)::this
    COOP_INT i
    COOP_REAL r
    COOP_SINGLE::map(this%nmaps)
    call this%peak_map%get_element(i, map)
    r = min(sqrt(abs(map(this%index_I))/max(abs(map(this%index_L)),1.d-20)), max_radius)
  end function coop_stacking_options_peak_r

  function coop_stacking_options_peak_e(this, i) result(e)
    COOP_REAL, parameter::max_radius = coop_pio2
    class(coop_stacking_options)::this
    COOP_INT i
    COOP_REAL e
    COOP_SINGLE::map(this%nmaps)
    call this%peak_map%get_element(i, map)
    e = min(sqrt(map(this%index_Q)**2 + map(this%index_U)**2)/max(abs(map(this%index_L)),1.d-20), 1.d0)
  end function coop_stacking_options_peak_e

  subroutine coop_stacking_options_peak_get_angle_r_e(this, i, angle, r, e)
    COOP_REAL, parameter::max_radius = coop_pio2
    class(coop_stacking_options)::this
    COOP_REAL angle,r,e
    COOP_SINGLE::map(this%nmaps)
    COOP_INT i
    if(i.gt. this%peak_pix%n) stop "rotate_angle: pix overflow"
    call this%peak_map%get_element(i, map)    
    select case(this%genre)
    case(coop_stacking_genre_Imax, coop_stacking_genre_Imin)
       call random_number(angle)
       angle = angle*coop_2pi
       e = 0.d0
       if(this%index_I.gt.0 .and. this%index_L.gt.0)then
          r = min(sqrt(abs(map(this%index_I))/max(abs(map(this%index_L)),1.d-20)), max_radius)
       else
          r = 0.d0
       endif
    case(coop_stacking_genre_Lmax, coop_stacking_genre_Lmin)
       call random_number(angle)
       angle = angle*coop_2pi
       e = 0.d0
       r = 0.d0
    case(coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Imin_Oriented)
       angle = COOP_POLAR_ANGLE(map(this%index_Q), map(this%index_U))
       if(this%index_L .gt. 0)then
          e = min(sqrt(map(this%index_Q)**2 + map(this%index_U)**2)/max(abs(map(this%index_L)),1.d-20), 0.99d0)
       else
          e = min(sqrt(map(this%index_Q)**2 + map(this%index_U)**2)/max(abs(map(this%index_I)),1.d-20), 0.99d0)
       endif
       if(this%index_I.gt.0 .and. this%index_L.gt.0)then
          r = min(sqrt(abs(map(this%index_I))/max(abs(map(this%index_L)),1.d-20)), max_radius)
       else
          r = 0.d0
       endif
    case(coop_stacking_genre_Lmax_Oriented, coop_stacking_genre_Lmin_Oriented, coop_stacking_genre_Pmax_Oriented, coop_stacking_genre_Pmin_Oriented)
       angle = COOP_POLAR_ANGLE(map(this%index_Q), map(this%index_U))
       e = 0.d0
       r = 0.d0
    case default
       stop "rotate_angle: unknown genre"
    end select
    r = min(sqrt(abs(map(this%index_I))/max(abs(map(this%index_L)),1.d-20)), max_radius)
  end subroutine coop_stacking_options_peak_get_angle_r_e
  
end module coop_stacking_mod
