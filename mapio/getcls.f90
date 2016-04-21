program ClFromMap
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools

  implicit none
#include "constants.h"

  COOP_STRING::fmap1, fmap2, clfile, fmask
  type(coop_healpix_maps)::map1, map2, mask
  type(coop_file)::fp
  logical::has_mask, has_map2, want_i
  COOP_INT::l, lmax, lmax_mask, lmin, nmaps1, nmaps2
  COOP_INT::nmaps
  COOP_REAL,dimension(:,:),allocatable::Cl_Pseudo, Cl
  COOP_REAL,dimension(:,:,:),allocatable::kernel
  COOP_REAL,dimension(:),allocatable::Cl_mask

  if(iargc().lt.2)then
     write(*,*) "Syntax:"
     write(*,*) "./MAP2CLS -map1 ... [-nmaps1 ...] [-map2 ... -nmaps2 ...] [-mask ...] -lmax ... [-lmax_mask ..] [-lmin ...]  -out ..."
     stop
  endif
  call coop_get_command_line_argument(key = "map1", arg = fmap1)
  call coop_get_command_line_argument(key = "nmaps1", arg = nmaps1, default = 0)
  call coop_get_command_line_argument(key = "nmaps2", arg = nmaps2, default = 0)
  call coop_get_command_line_argument(key = "lmax", arg = lmax)
  call coop_get_command_line_argument(key = "lmin", arg = lmin, default = 2)
  if(lmin .gt. lmax) stop "lmin > lmax?"
  call coop_get_command_line_argument(key = "map2", arg = fmap2, default = "")
  has_map2 = trim(fmap2).ne.""
  call coop_get_command_line_argument(key = "mask", arg = fmask, default = "")
  has_mask = trim(fmask) .ne. ""
  if(has_mask)then
     call coop_get_command_line_argument(key = "lmax_mask", arg = lmax_mask, default = 200)
  endif
  call coop_get_command_line_argument(key = "out", arg = clfile)
  if(nmaps1 .gt. 0)then
     call map1%read(fmap1, nmaps_wanted = nmaps1)
  else
     call map1%read(fmap1)
  endif
  if(has_map2)then
     if(nmaps2 .gt. 0)then
        call map2%read(fmap2, nmaps_wanted = nmaps2)
     else
        call map2%read(fmap2)
     endif
  endif
  call fp%open(clfile)


  if(has_mask)then
     select case(map1%nmaps)
     case(1) !!T map
        if(map1%spin(1).ne.0) stop "nmaps = 1 but not an I map"
     case(2)
        want_i = .false.
        if(map1%spin(1).ne.2 .or. map1%spin(2).ne.2) stop "nmaps = 2 but not a QU map"
     case(3)
        want_i = .true.
        if(map1%spin(1).ne.0 .or. map1%spin(2).ne.2 .or. map1%spin(3).ne.2) stop "nmaps = 3 but not an IQU map"
     case default
        stop "GetCl only works with I, QU, IQU maps"
     end select
     call mask%read(fmask, nmaps_wanted = 1)
     allocate(cl_mask(0:lmax_mask))
     call mask%map2alm(lmax = lmax_mask)
     cl_mask = mask%cl(0:lmax_mask, 1)
     if(map1%nmaps .eq. 1)then
        allocate(kernel(lmin:lmax, lmin:lmax, 1), Cl_pseudo(lmin:lmax, 1), Cl(lmin:lmax, 1))
        call coop_pseudoCl_get_kernel(lmax_mask, Cl_mask, lmin, lmax, kernel(:,:,1))
     else
        allocate(kernel(lmin:lmax, lmin:lmax, 4), Cl_pseudo(lmin:lmax, 6), Cl(lmin:lmax, 6))
        call coop_pseudoCl_get_kernel_pol(lmax_mask, Cl_mask, lmin, lmax, kernel, want_i)
     endif


     call map1%apply_mask(mask)
     call map1%map2alm(lmax = lmax)
     if(has_map2)then
        call map2%apply_mask(mask)
        call map2%map2alm(lmax = lmax)
        call map1%get_cls(map2)
     else
        call map1%get_cls()
     endif
     select case(map1%nmaps)
     case(1) !!T map
        cl_pseudo(lmin:lmax,1) = map1%cl(lmin:lmax, 1)
        call coop_pseudoCl2Cl(lmin = lmin, lmax = lmax, Cl_pseudo = Cl_pseudo(lmin:lmax,1), kernel = kernel(lmin:lmax, lmin:lmax, 1),  Cl = Cl(lmin:lmax, 1))
        do l = lmin, lmax
           write(fp%unit, "(I6, E16.7)") l, l*(l+1.d0)/coop_2pi*Cl(l, 1)
        enddo
     case(2)
        cl_pseudo(lmin:lmax, coop_TEB_index_EE) = map1%cl( lmin:lmax, COOP_MATSYM_INDEX(2, 1, 1) )
        cl_pseudo(lmin:lmax, coop_TEB_index_BB) = map1%cl( lmin:lmax,  COOP_MATSYM_INDEX(2, 2, 2) )
        cl_pseudo(lmin:lmax, coop_TEB_index_EB) = map1%cl( lmin:lmax,  COOP_MATSYM_INDEX(2, 1, 2) )
        call coop_pseudoCl2Cl_pol(lmin = lmin, lmax = lmax, Cl_pseudo = Cl_pseudo, kernel = kernel,  Cl = Cl, want_i = .false.)
        write(fp%unit, "(A8, 3A16)") "# ell ", " EE ",  " BB ", " EB "
        do l = lmin, lmax
           write(fp%unit, "(I8, 3E16.7)") l, l*(l+1.d0)/coop_2pi*Cl(l, coop_TEB_index_EE), l*(l+1.d0)/coop_2pi*Cl(l, coop_TEB_index_BB), l*(l+1.d0)/coop_2pi*Cl(l, coop_TEB_index_EB)
        enddo
     case(3)
        cl_pseudo(lmin:lmax, 1:6) = map1%cl(lmin:lmax, 1:6)
        call coop_pseudoCl2Cl_pol(lmin = lmin, lmax = lmax, Cl_pseudo = Cl_pseudo, kernel = kernel,  Cl = Cl)
        write(fp%unit, "(A8, 6A16)") "# ell ", " TT ", " EE ",  " BB ", " TE ", " EB ", " TB "
        do l = lmin, lmax
           write(fp%unit, "(I8, 6E16.7)") l, l*(l+1.d0)/coop_2pi*Cl(l, :)
        enddo

     end select
  else
     call map1%map2alm(lmax = lmax)
     if(has_map2)then
        call map2%map2alm(lmax = lmax)
        call map1%get_cls(map2)
     else
        call map1%get_cls()
     endif
     do l = lmin, lmax
        write(fp%unit, "(I6, "//COOP_STR_OF(map1%nmaps*(map1%nmaps+1)/2)//"E16.7)") l, map1%Cl(l, :)*(l*(l+1.)/coop_2pi)
     enddo
  endif
  call fp%close()

end program ClFromMap
