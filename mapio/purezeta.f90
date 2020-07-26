program test
#if HAS_HEALPIX  
  use coop_wrapper_utils
  use coop_wrapper_firstorder  
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  implicit none
#include "constants.h"
  COOP_STRING::inpfile, outfile, weight_option, action, maskfile
  COOP_REAL::fwhm_arcmin
  logical want_fluc
  COOP_INT::ninp, hp_l1, hp_l2
  type(coop_healpix_maps)::map, mask
  call coop_MPI_init()
  call coop_random_init()
  if(iargc().le.0)then
     write(*,"(A)") "Syntax:"
     write(*,"(A)") "./TE2Z -inp INPUT_MAP -fwhm FWHM_ARCMIN -out OUTPUT_FILE -fluc [F|T] -action [T2Z|TE2Z|E2Z] -weight [vis|latevis|earlyvis|recomb_slice|reion_slice] [-mask MASKFILE -hpl1 L1 -hpl2 L2]"
     stop
  endif
  call coop_get_command_line_argument(key = "inp", arg = inpfile)
  call coop_get_command_line_argument(key = "out", arg = outfile)
  call coop_get_command_line_argument(key = "mask", arg = maskfile, default="")
  call coop_get_command_line_argument(key = "weight", arg = weight_option, default = "vis")
  call coop_get_command_line_argument(key = "action", arg = action, default = "T2Z")    
  call coop_get_command_line_argument(key = "fwhm", arg = fwhm_arcmin)  
  call coop_get_command_line_argument(key = "fluc", arg = want_fluc, default = .true.)
  call coop_get_command_line_argument(key = "hpl1", arg = hp_l1, default = 0)
  call coop_get_command_line_argument(key = "hpl2", arg = hp_l2, default = hp_l1)
  if(hp_l1 .gt. 0) write(*,*) "Doing highpass filtering: ", hp_l1, hp_l2
  call map%read(inpfile, nmaps_wanted = 7)
  if(trim(maskfile).ne."")then
     call mask%read(maskfile, nmaps_wanted = 1)
     select case(trim(action))
     case("T2Z")
        call map%t2zeta(fwhm_arcmin = fwhm_arcmin, want_unconstrained = want_fluc, weight_option = trim(weight_option), mask = mask, hp_l1 = hp_l1, hp_l2 = hp_l2)
     case("TE2Z")
        call map%te2zeta(fwhm_arcmin = fwhm_arcmin, want_unconstrained = want_fluc, weight_option = trim(weight_option), hp_l1 = hp_l1, hp_l2 = hp_l2, mask = mask)
     case("E2Z")
        call map%e2zeta(fwhm_arcmin = fwhm_arcmin, want_unconstrained = want_fluc, weight_option = trim(weight_option), hp_l1 = hp_l1, hp_l2 = hp_l2, mask = mask)
     case default
        write(*,*) "action = "//trim(action)
        write(*,*) "action = T2Z/TE2Z/E2Z"
        stop "unknown action"
     end select
     
  else
     select case(trim(action))
     case("T2Z")
        call map%t2zeta(fwhm_arcmin = fwhm_arcmin, want_unconstrained = want_fluc, weight_option = trim(weight_option))
     case("TE2Z")
        call map%te2zeta(fwhm_arcmin = fwhm_arcmin, want_unconstrained = want_fluc, weight_option = trim(weight_option), hp_l1 = hp_l1, hp_l2 = hp_l2)
     case("E2Z")
        call map%e2zeta(fwhm_arcmin = fwhm_arcmin, want_unconstrained = want_fluc, weight_option = trim(weight_option), hp_l1 = hp_l1, hp_l2 = hp_l2)
     case default
        write(*,*) "action = "//trim(action)
        write(*,*) "action = T2Z/TE2Z/E2Z"
        stop "unknown action"
     end select
  endif
  if(want_fluc)then
     call map%write(outfile)
  else
     call map%write(outfile, index_list = (/ 1 /) )
  endif
  call coop_MPI_Finalize()
#else
  stop "You need healpix to compile this."
#endif  
end program test
