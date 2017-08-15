Module coop_zeta3d_mod
  use coop_wrapper_firstorder
  use coop_healpix_mod
  use coop_fitsio_mod
  implicit none
#include "constants.h"

  private

  COOP_REAL,dimension(:,:, :),allocatable::coop_zeta3d_trans
  type(coop_zeta_shell),dimension(:),allocatable::coop_zeta3d_shells

  public::coop_zeta3d_generate_cmb, coop_zeta3d_cleanup, coop_zeta3d_shells, coop_zeta3d_trans, coop_load_zeta3d_trans, coop_dump_zeta3d_trans

contains


  subroutine coop_load_zeta3d_trans(filename, lmax_loaded)
    COOP_UNKNOWN_STRING::filename
    COOP_INT::lmax_loaded, nr, nsrc, lmax
    COOP_REAL,dimension(:,:,:),allocatable::tmp
    type(coop_fits_file)::fits
    type(coop_file)::fp
    if(trim(coop_file_postfix_of(filename)).eq."fits")then  !!write a fits file
       call fits%open_image(filename)
       call coop_dictionary_lookup(fits%header, "LMAX", lmax)
       call coop_dictionary_lookup(fits%header, "NRBINS", nr)
       call coop_dictionary_lookup(fits%header, "NSRC", nsrc)

       if( size(coop_zeta3d_trans, 1) .ne. nsrc .or. size(coop_zeta3d_trans, 2) .ne. nr)then
          write(*,*) trim(filename)
          write(*,*) size(coop_zeta3d_trans, 1),nsrc
          write(*,*) size(coop_zeta3d_trans, 2) .ne. nr
          stop "load_zeta3d_trans: size mismatch"
       endif
       if(lmax+1 .le.  size(coop_zeta3d_trans,3))then
          call fits%load_image_3d(coop_zeta3d_trans(:,:,0:lmax))
          lmax_loaded = lmax
          call fits%close()
          return
       else
          allocate(tmp(nr, nsrc, 0:lmax))
          call fits%load_image_3d(tmp)
          coop_zeta3d_trans(:,:,:) = tmp(:,:,0:size(coop_zeta3d_trans,3)-1)
          deallocate(tmp)
          call fits%close()
          lmax_loaded = size(coop_zeta3d_trans,3)-1
          return
       endif


    else
       call fp%open(filename, "ur")
       read(fp%unit)  nsrc, nr, lmax
       if( size(coop_zeta3d_trans, 1) .ne. nsrc .or. size(coop_zeta3d_trans, 2) .ne. nr)then
          write(*,*) size(coop_zeta3d_trans, 1), nsrc
          write(*,*) size(coop_zeta3d_trans, 2), nr
          write(*,*) trim(filename)
          stop "load_zeta3d_trans: size mismatch"
       endif
       if(lmax+1 .lt.  size(coop_zeta3d_trans,3))then
          write(*,*) "reading up to lmax = ", lmax
          read(fp%unit) coop_zeta3d_trans(:,:,0:lmax)
          lmax_loaded = lmax
          call fp%close()
          return
       else
          allocate(tmp(nr, nsrc, 0:lmax))
          read(fp%unit) tmp
          coop_zeta3d_trans(:,:,:) = tmp(:,:,0:size(coop_zeta3d_trans,3)-1)
          deallocate(tmp)
          call fp%close()
          lmax_loaded = size(coop_zeta3d_trans,3)-1
          return
       endif
    endif
  end subroutine coop_load_zeta3d_trans


  subroutine coop_dump_zeta3d_trans(filename)
    COOP_UNKNOWN_STRING::filename
    type(coop_file)::fp
    type(coop_dictionary)::header
    if(trim(coop_file_postfix_of(filename)).eq."fits")then  !!write a fits file
       call header%insert("NRBINS", COOP_STR_OF(size(coop_zeta3d_trans,1)))
       call header%insert("NSRC", COOP_STR_OF(size(coop_zeta3d_trans, 2)))
       call header%insert("LMAX", COOP_STR_OF(size(coop_zeta3d_trans,3)-1))
       call coop_fits_file_write_image_3d(coop_zeta3d_trans, filename, header)
       call header%free()
    else
       call fp%open(filename, "u")
       write(fp%unit)  size(coop_zeta3d_trans, 1), size(coop_zeta3d_trans, 2), size(coop_zeta3d_trans,3)-1
       write(fp%unit) coop_zeta3d_trans
       call fp%close()
    endif
  end subroutine coop_dump_zeta3d_trans


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

  subroutine coop_zeta3d_generate_cmb(hm, zeta, cosmology, fnl, lmax, nside, transfile, prefix3d, prefixmap, fwhm_arcmin, writefile)
    class(coop_cosmology_firstorder)::cosmology
    external fnl
    logical writefile
    COOP_UNKNOWN_STRING::transfile, prefix3d
    COOP_UNKNOWN_STRING::prefixmap
    COOP_REAL, optional::fwhm_arcmin
    COOP_INT::l, lmax, nside, m1, m2, m3, i, nlc, isrc, lmax_loaded
    COOP_INT, dimension(:),allocatable::computed_ells
    COOP_REAL, dimension(:),allocatable::computed_trans,computed_trans2,  rls
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
    call hm%init(nside = nside, nmaps = 2, genre = "TE", lmax = lmax, unit="muK")
    call zeta%init(nside = nside, nmaps = 1, genre = "ZETA", lmax = lmax, unit="E-5")

    allocate(vis(0:lmax))
    hm%alm = 0.
    zeta%alm = 0.
    if(allocated(coop_zeta3d_shells) .and. allocated(coop_zeta3d_trans))then  !!prepare the Gaussian field
       do i=1, coop_zeta_nr
          tau = cosmology%tau0 - coop_zeta_r(i)
          if(tau .gt. 0.d0)then
             vis = cosmology%visofa(cosmology%aoftau(tau))
          else
             vis = 0.
          endif
          call coop_zeta3d_shells(i)%map_project(fnl, lmax, hm%alm(0:lmax, 0:lmax, 1), coop_zeta3d_trans(i, coop_index_source_T, 0:lmax)*coop_zeta_dr(i), alm_total2 = hm%alm(0:lmax, 0:lmax, 2), weight2 =  coop_zeta3d_trans(i, coop_index_source_E, 0:lmax)*coop_zeta_dr(i), alm_total3 =   zeta%alm(0:lmax, 0:lmax, 1), weight3 = vis*coop_zeta_dr(i))
       enddo
    else
       call cosmology%compute_source(0)
       call coop_set_default_zeta_r(cosmology, cosmology%source(0))
       allocate(coop_zeta3d_trans(coop_zeta_nr, cosmology%source(0)%nsrc, 0:lmax))
       coop_zeta3d_trans = 0.d0
       if(coop_file_exists(transfile))then
          call coop_feedback("Reading transfer data file.")
          call coop_load_zeta3d_trans(transfile, lmax_loaded)
       else
          lmax_loaded = 1
       endif
       if(lmax_loaded .lt. lmax)then
          call coop_feedback("Computing transfer function.")
          nlc = coop_nl_range(lmin = lmax_loaded+1, lmax = lmax)
          allocate(computed_ells(nlc), computed_trans(nlc), computed_trans2(nlc), rls(nlc))
          call coop_set_ells(computed_ells, lmax_loaded+1, lmax)
          rls = computed_ells
          !$omp parallel do 
          do i = 1, nlc
             call coop_get_zeta_trans_l(cosmology%source(0),  computed_ells(i), coop_zeta_nr, coop_zeta_r, coop_zeta3d_trans(:,:,computed_ells(i)))
             coop_zeta3d_trans(:,coop_index_source_E, computed_ells(i)) = coop_zeta3d_trans(:,coop_index_source_E,computed_ells(i))*sqrt((computed_ells(i)+2.)*(computed_ells(i)+1.)*computed_ells(i)*(computed_ells(i)-1.))
          enddo
          !$omp end parallel do
          do isrc = 1, cosmology%source(0)%nsrc
             do i = 1, coop_zeta_nr
                computed_trans = coop_zeta3d_trans(i, isrc, computed_ells)
                call coop_spline(nlc, rls, computed_trans, computed_trans2)
                do l = lmax_loaded+1, lmax
                   call coop_splint(nlc, rls, computed_trans, computed_trans2, dble(l), coop_zeta3d_trans(i, isrc, l))
                enddo
             enddo
          enddo
          call coop_dump_zeta3d_trans(transfile)
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
       if(coop_file_exists(trim(prefix3d)//"_3D_"//COOP_STR_OF(lmax)//".dat"))then
          call coop_feedback("Reading 3D file.")
          call fp%open(trim(prefix3d)//"_3D_"//COOP_STR_OF(lmax)//".dat", "u")
          do i=1, coop_zeta_nr
             read(fp%unit) coop_zeta3d_shells(i)%alm_real
             tau = cosmology%tau0 - coop_zeta_r(i)
             if(tau .gt. 0.d0)then
                vis = cosmology%visofa(cosmology%aoftau(tau))
             else
                vis = 0.
             endif
             call coop_zeta3d_shells(i)%map_project(fnl, lmax, hm%alm(0:lmax, 0:lmax, 1), coop_zeta3d_trans(i, coop_index_source_T, 0:lmax)*coop_zeta_dr(i), alm_total2 = hm%alm(0:lmax, 0:lmax, 2), weight2 =  coop_zeta3d_trans(i, coop_index_source_E, 0:lmax)*coop_zeta_dr(i), alm_total3 =   zeta%alm(0:lmax, 0:lmax, 1), weight3 = vis*coop_zeta_dr(i))
          enddo
          call fp%close()
       else
          call coop_feedback("Generating 3D maps now")
          call coop_generate_3Dzeta(cosmology, coop_zeta_nr, coop_zeta3d_shells)
          call coop_feedback("Projecting into a CMB map")
          if(writefile)call fp%open(trim(prefix3d)//"_3D_"//COOP_STR_OF(lmax)//".dat", "u")
          do i=1, coop_zeta_nr
             if(writefile)write(fp%unit) coop_zeta3d_shells(i)%alm_real
             tau = cosmology%tau0 - coop_zeta_r(i)
             if(tau .gt. 0.d0)then
                vis = cosmology%visofa(cosmology%aoftau(tau))
             else
                vis = 0.
             endif
             call coop_zeta3d_shells(i)%map_project(fnl, lmax, hm%alm(0:lmax, 0:lmax, 1), coop_zeta3d_trans(i, coop_index_source_T, 0:lmax)*coop_zeta_dr(i), alm_total2 = hm%alm(0:lmax, 0:lmax, 2), weight2 =  coop_zeta3d_trans(i, coop_index_source_E, 0:lmax)*coop_zeta_dr(i), alm_total3 =   zeta%alm(0:lmax, 0:lmax, 1), weight3 = vis*coop_zeta_dr(i))
          enddo
          if(writefile)call fp%close()
       endif
    endif
    do l = 2,  lmax
       hm%alm(l,0:l, :) =  hm%alm(l,0:l,:)*coop_Gaussian_filter(fwhm_arcmin  =fwhm, l = l)
       zeta%alm(l, 0:l, :) =  zeta%alm(l,0:l,:)*coop_Gaussian_filter(fwhm_arcmin  =fwhm, l = l)
    enddo
    hm%alm = hm%alm * (1.e6 * COOP_DEFAULT_TCMB) !!switch to muK
    zeta%alm = zeta%alm * (1.e5)  !!swith to unit 1.e-5
    call hm%get_Cls()
    call fp%open(trim(prefixmap)//"_"//COOP_STR_OF(lmax)//"_Cls.txt")
    do l = 2, hm%lmax
       write(fp%unit, "(I5, 10E16.7)") l, hm%Cl(l, :)*(l+1.)*l/coop_2pi
    enddo
    call fp%close()
    call hm%alm2map()
    call zeta%alm2map()
    if(writefile)then
       call hm%write(trim(prefixmap)//"_"//COOP_STR_OF(lmax)//"_TE.fits")
       call zeta%write(trim(prefixmap)//"_"//COOP_STR_OF(lmax)//"_zeta.fits")
    endif
    deallocate(vis)
  end subroutine coop_zeta3d_generate_cmb

end module coop_zeta3d_mod
