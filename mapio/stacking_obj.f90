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
  COOP_INT,parameter::coop_stacking_genre_Pmax_Oriented = 5
  COOP_INT,parameter::coop_stacking_genre_Pmin_Oriented = 6

  type coop_stacking_options
     COOP_INT::genre = coop_stacking_genre_Null
     COOP_INT::nmaps = 0
     COOP_INT::index_I = 0
     COOP_INT::index_Q = 0
     COOP_INT::index_U = 0
     COOP_INT::index_L = 0
     COOP_STRING::caption
     COOP_SINGLE::I_lower = -1.e31
     COOP_SINGLE::I_upper = 1.e31
     COOP_SINGLE::L_lower = -1.e31
     COOP_SINGLE::L_upper = 1.e31
     COOP_SINGLE::P2_lower = 0.
     COOP_SINGLE::P2_upper = 1.e31
     COOP_SINGLE::P2byI2_lower = -1.e31
     COOP_SINGLE::P2byI2_upper = 1.e31
     COOP_SINGLE::P2byL2_lower = -1.e31
     COOP_SINGLE::P2byL2_upper = 1.e31
   contains
     procedure::init => coop_stacking_options_init
     procedure::accept=>coop_stacking_options_accept
  end type coop_stacking_options


contains

  subroutine coop_stacking_options_init(this, domax, peak_name, Orient_name)
    class(coop_stacking_options)::this
    COOP_UNKNOWN_STRING::peak_name, orient_name
    COOP_SHORT_STRING::p
    logical domax
    p = trim(adjustl(peak_name))
    if(trim(adjustl(orient_name)).eq."" .or. trim(adjustl(orient_name)) .eq. "random")then
       if(domax)then
          this%genre = coop_stacking_genre_Imax
          this%caption = trim(p)//" maxima"
       else
          this%genre = coop_stacking_genre_Imin
          this%caption = trim(p)//" mainima"          
       endif
    else
       if(domax)then
          if(p(1:2).eq."$P")then
             this%genre = coop_stacking_genre_Pmax_Oriented             
          else
             this%genre = coop_stacking_genre_Imax_Oriented
          endif
          this%caption = trim(p)//" maxima, "//trim(adjustl(Orient_name))//" oriented"          
       else
          if(p(1:2).eq."$P")then
             this%genre = coop_stacking_genre_Pmin_Oriented             
          else
             this%genre = coop_stacking_genre_Imin_Oriented
          endif
          this%caption = trim(adjustl(peak_name))//" minima, "//trim(adjustl(Orient_name))//" oriented"                    
       endif
    endif
       
  end subroutine coop_stacking_options_init


  function coop_stacking_options_accept(this, I, P2, L) result(acc)
    class(coop_stacking_options)::this
    logical acc
    COOP_REAL, optional::I, P2, L
    if(present(I))then
       if(present(P2))then
          if(present(L))then
             acc = (I .lt. this%I_upper .and. I .gt. this%I_lower .and. P2 .lt. this%P2_upper .and. P2 .gt. this%P2_lower .and.  L .lt. this%L_upper .and. L .gt. this%L_lower .and. P2 .lt. this%P2byI2_upper*I**2 .and. P2 .gt. this%P2byI2_lower * I**2 .and. P2 .lt. this%P2byL2_upper*L**2 .and. P2 .gt. this%P2byL2_lower*L**2)
          else
             acc = (I .lt. this%I_upper .and. I .gt. this%I_lower .and. P2 .lt. this%P2_upper .and. P2 .gt. this%P2_lower .and. P2 .lt. this%P2byI2_upper*I**2 .and. P2 .gt. this%P2byI2_lower * I**2)
          endif
       else
          if(present(L))then
             acc = (I .lt. this%I_upper .and. I .gt. this%I_lower  .and.  L .lt. this%L_upper .and. L .gt. this%L_lower)
          else
             acc = (I .lt. this%I_upper .and. I .gt. this%I_lower)
          endif
       endif
    else
       if(present(P2))then
          if(present(L))then
             acc = ( P2 .lt. this%P2_upper .and. P2 .gt. this%P2_lower .and.  L .lt. this%L_upper .and. L .gt. this%L_lower .and. P2 .lt. this%P2byL2_upper*L**2 .and. P2 .gt. this%P2byL2_lower*L**2)
          else
             acc = (P2 .lt. this%P2_upper .and. P2 .gt. this%P2_lower)
          endif
       else
          if(present(L))then
             acc = (L .lt. this%L_upper .and. L .gt. this%L_lower)
          else
             acc = .true.
          endif
       endif
    endif
  end function coop_stacking_options_accept
  
end module coop_stacking_mod
