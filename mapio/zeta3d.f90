Module coop_zeta3d_mod
  use coop_wrapper_firstorder
  use coop_healpix_mod
  implicit none
#include "constants.h"

  private

  COOP_REAL,dimension(:,:, :),allocatable::coop_zeta3d_trans
  type(coop_zeta_shell),dimension(:),allocatable::coop_zeta3d_shells

  public::coop_zeta3d_generate_cmb, coop_zeta3d_cleanup

contains

  subroutine coop_zeta3d_cleanup()
    COOP_INT::i
    COOP_DEALLOC(coop_zeta3d_trans)
    if(allocated(coop_zeta3d_shells))then
       do i=1, size(coop_zeta3d_shells)
          call coop_zeta3d_shells(i)%free()
       enddo
       deallocate(coop_zeta3d_shells)
    endif
  end subroutine coop_zeta3d_cleanup

  subroutine coop_zeta3d_generate_cmb(cosmology, fnl, lmax, nside, prefix, mapprefix, fwhm_arcmin)
    class(coop_cosmology_firstorder)::cosmology
    external fnl
    COOP_UNKNOWN_STRING::prefix
    COOP_UNKNOWN_STRING::mapprefix
    COOP_REAL, optional::fwhm_arcmin
    COOP_INT::l, lmax, nside, m1, m2, m3, i
    logical::need_compute
    type(coop_file)fp
    type(coop_healpix_maps)::hm, zeta
    COOP_REAL, dimension(:), allocatable::vis
    COOP_REAL::fwhm, rat, tau
    if(present(fwhm_arcmin))then
       fwhm = fwhm_arcmin 
    else
       fwhm = 60.d0  !!default 1 degree resolution
    endif
    call hm%init(nside = nside, nmaps = 2, genre = "TE", lmax = lmax)
#if DO_ZETA_TRANS
    call zeta%init(nside = nside, nmaps = 2, genre = "ZETA", lmax = lmax)
#else
    call zeta%init(nside = nside, nmaps = 1, genre = "ZETA", lmax = lmax)
#endif
    allocate(vis(0:lmax))
    hm%alm = 0.
    zeta%alm = 0.
    if(allocated(coop_zeta3d_shells) .and. allocated(coop_zeta3d_trans))then  !!prepare the Gaussian field
       do i=1, coop_zeta_nr
          call coop_zeta3d_shells(i)%map_project(fnl, lmax, hm%alm(0:lmax, 0:lmax, 1), coop_zeta3d_trans(i, 1, 0:lmax)*coop_zeta_dr(i))
          call coop_zeta3d_shells(i)%map_project(fnl, lmax, hm%alm(0:lmax, 0:lmax, 2), coop_zeta3d_trans(i, 2, 0:lmax)*coop_zeta_dr(i))
          tau = cosmology%tau0 - coop_zeta_r(i)
          if(tau .gt. 0.d0)then
             vis = cosmology%visofa(cosmology%aoftau(tau))
             call coop_zeta3d_shells(i)%map_project(fnl, lmax, zeta%alm(0:lmax, 0:lmax, 1), vis*coop_zeta_dr(i))
          endif
       enddo
    else
       call cosmology%compute_source(0)
       call coop_set_default_zeta_r(cosmology, cosmology%source(0))
       allocate(coop_zeta3d_trans(coop_zeta_nr, cosmology%source(0)%nsrc, 0:lmax))
       coop_zeta3d_trans = 0.d0
       if(coop_file_exists(trim(prefix)//"_trans_"//COOP_STR_OF(lmax)//".dat"))then
          call coop_feedback("Reading transfer data file.")
          call fp%open(trim(prefix)//"_trans_"//COOP_STR_OF(lmax)//".dat", "u")
          read(fp%unit) m1, m2, m3
          if(m1 .eq. cosmology%source(0)%nsrc .and. m2.eq.coop_zeta_nr .and. m3.eq.lmax)then
             read(fp%unit) coop_zeta3d_trans
             call fp%close()
             need_compute = .false.
          else
             need_compute = .true.
          endif
       else
          need_compute = .true.
       endif
       if(need_compute)then
          call coop_feedback("Computing transfer function.")
          !$omp parallel do
          do l=2, lmax
             call coop_get_zeta_trans_l(cosmology%source(0),  l, coop_zeta_nr, coop_zeta_r, coop_zeta3d_trans(:,:,l))
             coop_zeta3d_trans(:,2,l) = coop_zeta3d_trans(:,2,l)*sqrt((l+2.)*(l+1.)*l*(l-1.))
          enddo
          !$omp end parallel do
          call fp%open(trim(prefix)//"_trans_"//COOP_STR_OF(lmax)//".dat", "u")
          write(fp%unit)  cosmology%source(0)%nsrc, coop_zeta_nr, lmax
          write(fp%unit) coop_zeta3d_trans
          call fp%close()
       endif
       allocate(coop_zeta3d_shells(coop_zeta_nr))
       do i=1, coop_zeta_nr
          coop_zeta3d_shells(i)%r = coop_zeta_r(i)
          coop_zeta3d_shells(i)%dr = coop_zeta_dr(i)
          rat = (coop_zeta3d_shells(i)%r/cosmology%distlss/0.8)**8
          if(rat .lt. 1.d0)then
             call coop_zeta3d_shells(i)%set_lmax(min(ceiling(30*(1.d0-rat)+lmax*rat), lmax))
          else
             call coop_zeta3d_shells(i)%set_lmax(lmax)
          endif
       enddo
       if(coop_file_exists(trim(prefix)//"_3D_"//COOP_STR_OF(lmax)//".dat"))then
          call coop_feedback("Reading 3D file.")
          call fp%open(trim(prefix)//"_3D_"//COOP_STR_OF(lmax)//".dat", "u")
          do i=1, coop_zeta_nr
             read(fp%unit) coop_zeta3d_shells(i)%alm_real
             call coop_zeta3d_shells(i)%map_project(fnl, lmax, hm%alm(0:lmax, 0:lmax, 1), coop_zeta3d_trans(i, coop_index_source_T, 0:lmax)*coop_zeta_dr(i))
             call coop_zeta3d_shells(i)%map_project(fnl, lmax, hm%alm(0:lmax, 0:lmax, 2), coop_zeta3d_trans(i, coop_index_source_E, 0:lmax)*coop_zeta_dr(i))
#if DO_ZETA_TRANS
             call coop_zeta3d_shells(i)%map_project(fnl, lmax, zeta%alm(0:lmax, 0:lmax, 2), coop_zeta3d_trans(i, coop_index_source_zeta, 0:lmax)*coop_zeta_dr(i))
#endif
             tau = cosmology%tau0 - coop_zeta_r(i)
             if(tau .gt. 0.d0)then
                vis = cosmology%visofa(cosmology%aoftau(tau))
                vis(0) = 0
                vis(1) = 0
                call coop_zeta3d_shells(i)%map_project(fnl, lmax, zeta%alm(0:lmax, 0:lmax, 1), vis*coop_zeta_dr(i))
             endif
          enddo
          call fp%close()
       else
          call coop_feedback("Generating 3D maps now")
          call coop_generate_3Dzeta(cosmology, coop_zeta_nr, coop_zeta3d_shells)
          call coop_feedback("Projecting into a CMB map")
          call fp%open(trim(prefix)//"_3D_"//COOP_STR_OF(lmax)//".dat", "u")
          do i=1, coop_zeta_nr
             write(fp%unit) coop_zeta3d_shells(i)%alm_real
             call coop_zeta3d_shells(i)%map_project(fnl, lmax, hm%alm(0:lmax, 0:lmax, 1), coop_zeta3d_trans(i, coop_index_source_T, 0:lmax)*coop_zeta_dr(i))
             call coop_zeta3d_shells(i)%map_project(fnl, lmax, hm%alm(0:lmax, 0:lmax, 2), coop_zeta3d_trans(i, coop_index_source_E, 0:lmax)*coop_zeta_dr(i))
             tau = cosmology%tau0 - coop_zeta_r(i)
             if(tau .gt. 0.d0)then
                vis = cosmology%visofa(cosmology%aoftau(tau))
                call coop_zeta3d_shells(i)%map_project(fnl, lmax, zeta%alm(0:lmax, 0:lmax, 1), vis*coop_zeta_dr(i))
             endif
#if DO_ZETA_TRANS
             call coop_zeta3d_shells(i)%map_project(fnl, lmax, zeta%alm(0:lmax, 0:lmax, 2), coop_zeta3d_trans(i, coop_index_source_zeta, 0:lmax)*coop_zeta_dr(i))
#endif
          enddo
          call fp%close()
       endif
    endif
    do l = 2,  lmax
       hm%alm(0:l, l, :) =  hm%alm(0:l, l, :)*coop_Gaussian_filter(fwhm_arcmin  =fwhm, l = l)
       zeta%alm(0:l, l, :) =  zeta%alm(0:l, l, :)*coop_Gaussian_filter(fwhm_arcmin  =fwhm, l = l)
    enddo
    hm%alm = hm%alm * (1.e6 * COOP_DEFAULT_TCMB) !!switch to muK
    zeta%alm = zeta%alm * (1.e5)  !!swith to unit 1.e-5
    call hm%get_Cls()
    call fp%open(trim(mapprefix)//"_"//COOP_STR_OF(lmax)//"_Cls.txt")
    do l = 2, hm%lmax
       write(fp%unit, "(I5, 10E16.7)") l, hm%Cl(l, :)*(l+1.)*l/coop_2pi
    enddo
    call fp%close()
    call hm%alm2map()
    call hm%write(trim(mapprefix)//"_"//COOP_STR_OF(lmax)//"_TE.fits")
    call hm%free()
    call zeta%alm2map()
    call zeta%write(trim(mapprefix)//"_"//COOP_STR_OF(lmax)//"_zeta.fits")
    call zeta%free()
    deallocate(vis)
  end subroutine coop_zeta3d_generate_cmb
  



end module coop_zeta3d_mod
