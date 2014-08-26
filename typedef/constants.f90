module coop_constants_mod
  implicit none
#include "constants.h"
!! universal

#include "constants_datatypes.h"
#include "constants_char.h"
#include "constants_math.h"
#include "constants_physics.h"
#include "constants_settings.h"
  
  COOP_STRING, parameter:: coop_version = "0.0"

  COOP_INT::coop_global_feedback = 1  !!0 no feedback, 1 key info,  2 verbose, 3 debug 

contains

  subroutine coop_feedback(info, feedback_level, action)
    COOP_UNKNOWN_STRING info
    COOP_INT, optional:: feedback_level
    COOP_UNKNOWN_STRING,optional::action
    character ans
    if(present(feedback_level))then
       if(coop_global_feedback .ge. feedback_level)then
          write(*,*) trim(info)
       endif
    else
       if(coop_global_feedback .ge. 0)then
          write(*,*) trim(info)
       endif
    endif
    if(present(action))then
       select case(trim(action))
       case("stop", "STOP")
          stop
       case("pause", "PAUSE")
          write(*,*) "continue? (Y/N)"
          read(*,*) ans
          if(ans .ne. "Y" .and. ans .ne. "y") stop
       end select
    endif
  end subroutine coop_feedback


end module coop_constants_mod
