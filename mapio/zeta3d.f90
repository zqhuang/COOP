module coop_zeta3d
  use coop_wrapper_firstorder
  use coop_healpix_mod
  implicit none
#include "constants.h"


contains

  subroutine coop_3d_generate_cmb(cosmology, fnl, lmax, nside, prefix)
    class(coop_cosmology_firstorder)::cosmology
    external fnl
    COOP_UNKNOWN_STRING::prefix
    COOP_REAL,dimension(:,:, :),allocatable::trans
    COOP_INT::l, lmax, nside, m1, m2, m3
    logical::need_compute
    type(coop_file)fp
    type(coop_zeta_shell),dimension(:),allocatable::shells
    call cosmology%compute_source(0)
    call coop_set_default_zeta_r(cosmology, cosmology%source(0))
    allocate(trans(cosmology%source(0)%nsrc, coop_zeta_nr, 2:lmax))
    if(coop_file_exists(trim(prefix)//"_trans.dat"))then
       call fp%open(trim(prefix)//"_trans.dat", "u")
       read(fp%unit) m1, m2, m3
       if(m1 .eq. cosmology%source(0)%nsrc .and. m2.eq.coop_zeta_nr .and. m3.eq.lmax)then
          read(fp%unit) trans
          call fp%close()
          need_compute = .false.
       else
          need_compute = .true.
       endif
    else
       need_compute = .true.
    endif
    if(need_compute)then
       !$omp parallel do
       do l=2, lmax
          call coop_get_zeta_trans_l(cosmology%source(0),  l, coop_zeta_nr, coop_zeta_r, trans(:,:,l))
       enddo
       !$omp end parallel do
       call fp%open(trim(prefix)//"_trans.dat", "u")
       write(fp%unit)  cosmology%source(0)%nsrc, coop_zeta_nr, lmax
       write(fp%unit) trans
       call fp%close()
    endif
    allocate(shells(coop_zeta_nr))
    deallocate(trans)
  end subroutine coop_3d_generate_cmb
  



end module coop_zeta3d
