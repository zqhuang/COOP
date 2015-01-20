module coop_stacking_mod
  use coop_wrapper_firstorder
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
     logical::nested = .false.
     type(coop_list_integer)::peak_pix
     type(coop_list_realarr)::peak_ang
     type(coop_list_realarr)::peak_map
   contains
     procedure::init => coop_stacking_options_init
     procedure::reject=>coop_stacking_options_reject
     procedure::export_pix => coop_stacking_options_export_pix
     procedure::export_ang => coop_stacking_options_export_ang     
     procedure::rotate_angle => coop_stacking_options_rotate_angle
     procedure::peak_r => coop_stacking_options_peak_r
     procedure::peak_e => coop_stacking_options_peak_e
     procedure::peak_get_angle_r_e => coop_stacking_options_peak_get_angle_r_e
  end type coop_stacking_options


contains

  subroutine coop_stacking_options_init(this, domax, peak_name, Orient_name, nmaps)
    class(coop_stacking_options)::this
    COOP_UNKNOWN_STRING::peak_name, orient_name
    COOP_SHORT_STRING::p
    logical domax
    COOP_INT::nmaps
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
       case(coop_stacking_genre_Lmax, coop_stacking_genre_Lmin)
          this%index_L = 1
       case default
          stop "nmaps does not match stacking option"
       end select
    case(2)
       select case(this%genre)
       case(coop_stacking_genre_Pmax_Oriented, coop_stacking_genre_Pmin_Oriented)
          this%index_Q = 1
          this%index_U = 2
       case(coop_stacking_genre_Imax, coop_stacking_genre_Imin, coop_stacking_genre_Lmax, coop_stacking_genre_Lmin)
          this%index_I = 1
          this%index_L = 2
       case default
          stop "nmaps does not match stacking option"          
       end select
    case(3)
       select case(this%genre)
       case(coop_stacking_genre_Pmax_Oriented, coop_stacking_genre_Pmin_Oriented, coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Imin_Oriented, coop_stacking_genre_Imax, coop_stacking_genre_Imin)
          this%index_I = 1
          this%index_Q = 2
          this%index_U = 3
       case(coop_stacking_genre_Lmax, coop_stacking_genre_Lmin, coop_stacking_genre_Lmax_Oriented, coop_stacking_genre_Lmin_Oriented)  
          this%index_L = 1
          this%index_Q = 2
          this%index_U = 3
       end select
    case(4:)
       this%index_I = 1
       this%index_Q = 2
       this%index_U = 3
       this%index_L = 4
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

  subroutine coop_stacking_options_export_ang(this, filename, addpi)
    class(coop_stacking_options)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file)::fp
    COOP_INT :: i
    logical, optional::addpi
    logical dor
    call fp%open(filename, "w")
    if(present(addpi))then
       dor = addpi
    else
       dor = .false.
    endif
    if(dor)then
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
