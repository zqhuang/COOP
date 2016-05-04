program map
  use coop_healpix_mod
  use coop_wrapper_firstorder
  use coop_wrapper_utils
#ifdef HAS_HEALPIX  
  use healpix_types
  use alm_tools
  use pix_tools
  use head_fits
  use fitstools
#endif  
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX  
  COOP_INT::fail_times = 0
  COOP_INT,parameter::nmax = 8
  character(len=80), dimension(64) :: header
  COOP_STRING,dimension(nmax)::fin
  COOP_STRING::fout, fbeam, field
  COOP_INT,dimension(:),allocatable::nmaps_in, ordering_in, nside_in
  COOP_INT nin, i, j, l, il, l1, l2, pix
  real,dimension(:,:),allocatable:: map_tmp
  COOP_LONG_INT npixtot
  type(coop_list_integer),dimension(nmax)::indices_wanted
  COOP_STRING::inline
  COOP_INT num_maps_wanted, npix, k
  logical::inline_mode = .false.
  type(coop_healpix_maps) hgm, hgm2, mask
  COOP_REAL fwhm, scal, threshold, eig_p, eig_n, pol_amp
  COOP_SINGLE::sigmas(100), maxabs, fsky
  type(coop_file)::fp

  nin = 1
  inline_mode =  (iargc() .gt. 0)
  if(.not. inline_mode)then
     write(*,*) "options are: SPLIT; PARTIAL2FULL; SMOOTH; GSMOOTH; DOBEAM; MULTIPLY;I2TQTUT;I2TQUL;I2TQULDD;I2TQUL6D;IQU2TEB;T2ZETA; IQU2ZETA; QU2ZETA; SCALE;INFO;ADD;SUBTRACT;MAKEMASK; SHUFFLE; HIGHPASS; LOWPASS; LOG; EXP; LOGIQU; EXPIQU; GAUSSIANIZE"
  endif
  do while(nin .le. nmax)
     if(inline_mode)then
        fin(nin) = trim(coop_inputArgs(nin))
     else
        write(*,*) "Enter input file and press Enter (or just press Enter key to finish):"
        read(*,'(A)') fin(nin)
     endif
     select case(trim(fin(nin)))
     case("INFO")
        do i=1, nin
           npixtot = getsize_fits(trim(fin(i)), nmaps = hgm%nmaps, nside = hgm%nside, ordering = hgm%ordering)
           write(*,"(A)")"**********************************"
           write(*,"(A)") trim(fin(i))//":"
           write(*, "(A,I5)") "nmaps = "//COOP_STR_OF(hgm%nmaps)
           write(*, "(A,I5)") "nside = "//COOP_STR_OF(hgm%nside)
           if(hgm%ordering .eq. COOP_RING)then
              write(*, "(A)") "ordering: ring"
           elseif(hgm%ordering .eq. COOP_NESTED)then
              write(*, "(A)") "ordering: nested"
           else
              write(*, "(A)") "ordering: unknown"
           endif
        enddo
     case("MAKEMASK")
        if(inline_mode)then
           inline = coop_inputArgs(nin+1)
           read(inline, *) threshold
           fout  = trim(coop_inputArgs(nin+2))
        else           
           if(nin.ge.3) stop "MAKEMASK option can only be applied to 1 map"
           call hgm%read(fin(1), nmaps_wanted = 1)
           write(*,*) "map min:", minval(hgm%map(:,1))
           write(*,*) "map max:", maxval(hgm%map(:,1))
           write(*,*) "Enter the threshold:"
           read(*,*) threshold
           write(*,*) "Enter the output file name:"
           read(*,'(A)') fout
        endif
        !$omp parallel do
        do i=0, hgm%npix-1
           if(hgm%map(i, 1).gt. threshold)then
              hgm%map(i, 1) = 0.
           else
              hgm%map(i, 1) = 1.
           endif
        enddo
        !$omp end parallel do
        call hgm%write(trim(fout))
        call hgm%free()
        goto 500
     case("SPLIT")
        nin = nin -1
        do i=1, nin
           call coop_healpix_split(trim(fin(i)))
        enddo
        print*, "maps are all split"
        goto 500
     case("PARTIAL2FULL")
        nin = nin -1
        do i=1, nin
           call hgm%read_simple(fin(i))
           call hgm%write(fin(i))
        enddo
        print*, "maps are all converted to fullsky"
        goto 500
     case("DOBEAM")
        if(inline_mode)then
           fbeam = trim(coop_inputArgs(nin+1))
        else
           write(*,*) "Enter the beam file:"
           read(*, "(A)") fbeam
        endif
        fbeam = trim(adjustl(fbeam))
        if( .not. coop_file_exists(fbeam))then
           write(*,"(A)") "beam file "//trim(fbeam)//" does not exist"
           if(inline_mode)stop
           do while(.not. coop_file_exists(fbeam))
              write(*,*) "Enter the beam file:"
              read(*, "(A)") fbeam
              fbeam = adjustl(fbeam)
           enddo
        endif
        do i=1, nin-1
           call hgm%read(trim(fin(i)))        
           call fp%open(trim(fbeam), 'r')
           call hgm%map2alm()
           do l=2, hgm%lmax
              read(fp%unit, *) il, scal
              if(il.ne.l) stop "beam file broken"
              hgm%alm(l, :, :) = hgm%alm(l, :, :)*scal
           enddo
           call fp%close()
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_beam")))
        enddo
        print*, "beam applied"
        goto 500
     case("HIGHPASS")
        if(inline_mode)then
           inline = coop_inputArgs(nin+1)
           read(inline, *) l1
           inline = coop_inputArgs(nin+2)
           read(inline, *) l2
        else
           write(*,*) "Enter l1, l2"
           read(*,*) l1, l2
        endif
        nin = nin -1
        do i=1, nin
           call hgm%read(trim(fin(i)))
           call hgm%map2alm()
           hgm%alm(0:l1, :,:) = 0.
           do l=l1+1, l2-1
              hgm%alm(l, :, :) = hgm%alm(l, :, :)*sin(dble(l-l1)/(l2-l1)*coop_pio2)
           enddo
           call hgm%alm2map()
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_hp_"//COOP_STR_OF(l1)//"_"//COOP_STR_OF(l2))))

        enddo

        print*, "maps are all highpassed"
        goto 500
     case("LOWPASS")
        if(inline_mode)then
           inline = coop_inputArgs(nin+1)
           read(inline, *) l1
           inline = coop_inputArgs(nin+2)
           read(inline, *) l2
        else
           write(*,*) "Enter l1, l2"
           read(*,*) l1, l2
        endif
        nin = nin -1
        do i=1, nin
           call hgm%read(trim(fin(i)))
           call hgm%map2alm()
           hgm%alm(l2:hgm%lmax, :,:) = 0.
           do l=l1+1, l2-1
              hgm%alm(l, :, :) = hgm%alm(l, :, :)*sin(dble(l2-l)/(l2-l1)*coop_pio2)
           enddo
           call hgm%alm2map()
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_lp_"//COOP_STR_OF(l1)//"_"//COOP_STR_OF(l2))))

        enddo

        print*, "maps are all lowpassed"
        goto 500
     case("SMOOTH")
        if(inline_mode)then
           inline = coop_inputArgs(nin+1)
           read(inline, *) fwhm
        else
           write(*,*) "Enter the fwhm in arcmin:"
           read(*,*) fwhm
        endif
        fwhm = fwhm*coop_SI_arcmin
        nin = nin -1
        do i=1, nin
           call coop_healpix_smooth_mapfile(trim(fin(i)), fwhm)
        enddo
        print*, "maps are all smoothed"
        goto 500
     case("GSMOOTH")
        if(inline_mode)then
           inline = coop_inputArgs(nin+1)
           read(inline, *) fwhm
           inline = coop_inputArgs(nin+2)
           read(inline, *) l1
           inline = coop_inputArgs(nin+3)
           read(inline, *) l2
        else
           write(*,*) "Enter the fwhm in arcmin:"
           read(*,*) fwhm
           write(*, *) "Eneter highpass l1"
           read(*, *) l1
           write(*, *) "Eneter highpass l2"
           read(*, *) l2
        endif
        fwhm = fwhm*coop_SI_arcmin
        nin = nin -1
        do i=1, nin
           call coop_healpix_gsmooth_mapfile(trim(fin(i)), fwhm = fwhm, highpass_l1 = l1, highpass_l2 = l2)
        enddo
        print*, "maps are all gsmoothed"
        goto 500
     case("SCALE")
        if(inline_mode)then
           inline = coop_inputArgs(nin+1)
           read(inline, *) scal
        else           
           write(*,*) "Enter the scale"
           read(*,*) scal
        endif
        nin  = nin - 1
        do i=1, nin
           call hgm%read(trim(fin(i)))
           hgm%map = hgm%map * scal
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_scal")))
        enddo
        call hgm%free
        goto 500
!!$     case("GAUSSIANIZE")
!!$        if(inline_mode)then
!!$           inline = coop_inputArgs(nin+1)
!!$        else
!!$           write(*,*) "Enter the mask (use NONE for no mask)"
!!$           read(*, "(A)") inline
!!$        endif
!!$        nin = nin - 1        
!!$        if(trim(adjustl(inline)).eq. "NONE" .or. trim(adjustl(inline)).eq. "")then  !!no mask
!!$           do i = 1, nin
!!$              call hgm%read(trim(fin(i)))
!!$              if( (hgm%nmaps .eq. 1 .and. hgm%spin(1) .eq. 0) .or. (hgm%nmaps .eq. 3 .and. hgm%spin(1) .eq. 0 .and. hgm%iq .eq. 2) )then
!!$                 write(*,*) "GAUSSIANIZING "//trim(fin(i))
!!$                 call hgm%map2alm()
!!$                 call hgm%simulate()
!!$                 call hgm%write(trim(coop_file_add_postfix(fin(i), "_gauss_sim")))
!!$              else
!!$                 write(*,*) "GAUSSIANIZE not supported for "//trim(fin(i))
!!$              endif
!!$           end do
!!$        else
!!$           call mask%read(trim(adjustl(inline)))
!!$           fsky = sum(dble(mask%map(:, 1)))/mask%npix           
!!$           do i = 1, nin
!!$              call hgm%read(trim(fin(i)))
!!$              if( (hgm%nmaps .eq. 1 .and. hgm%spin(1) .eq. 0) .or. (hgm%nmaps .eq. 3 .and. hgm%spin(1) .eq. 0 .and. hgm%iq .eq. 2) )then
!!$                 write(*,*) "GAUSSIANIZING "//trim(fin(i))                 
!!$                 call hgm%apply_mask(mask = mask)
!!$                 call hgm%map2alm()
!!$                 hgm%cl(0:1, :) = 0.
!!$                 hgm%cl = hgm%cl/fsky
!!$                 call hgm%simulate()
!!$                 call hgm%write(trim(coop_file_add_postfix(fin(i), "_gauss_sim")))
!!$              else
!!$                 write(*,*) "GAUSSIANIZE not supported for "//trim(fin(i))
!!$              endif
!!$           end do
!!$        endif
!!$        call hgm%free()
!!$        goto 500
     case("LOG")
        nin  = nin - 1
        do i=1, nin
           call hgm%read(trim(fin(i)))
           do j = 1, hgm%nmaps
              scal = sqrt(sum(hgm%map(:, j)**2)/hgm%npix)*1.e-12
              hgm%map(:, j) = log(max(hgm%map(:, j), 0.)+scal)
           enddo
           call hgm%set_fields("LNI")
           call hgm%set_units("1")
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_ln")))
        enddo
        call hgm%free
        goto 500
     case("EXP")
        nin  = nin - 1
        do i=1, nin
           call hgm%read(trim(fin(i)))
           hgm%map = exp(hgm%map)
           call hgm%set_fields("I")
           call hgm%set_units("muK")
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_exp")))
        enddo
        call hgm%free
        goto 500
     case("LOGIQU")
        nin  = nin - 1
        do i=1, nin
           call hgm%read(trim(fin(i)), nmaps_wanted = 3)
           scal = sqrt(sum(hgm%map(:, 1)**2)/hgm%npix)*1.d-12
           fail_times = 0
           !$omp parallel do private(eig_p, eig_n, pol_amp, pix) reduction(+:fail_times)
           do pix = 0, hgm%npix-1
              if(hgm%map(pix,1) .le. 0.d0)then
                 fail_times = fail_times + 1
                 hgm%map(pix,1) = scal
              endif
              pol_amp = sqrt(hgm%map(pix, 2)**2 + hgm%map(pix, 3)**2)
              if(pol_amp .le. 0.d0)then
                 hgm%map(pix,1) = log(hgm%map(pix,1))
                 hgm%map(pix,2:3) = 0.
              elseif(pol_amp .ge. hgm%map(pix,1)*0.99)then
                 hgm%map(pix,1) = log(hgm%map(pix,1))
                 hgm%map(pix,2:3) = 0.
                 fail_times = fail_times + 1
              else
                 eig_p = log(hgm%map(pix, 1) + pol_amp)
                 eig_n = log(hgm%map(pix, 1) - pol_amp)
                 hgm%map(pix, 1) = 0.5d0 * (eig_p + eig_n)
                 hgm%map(pix, 2:3) = hgm%map(pix, 2:3) * (0.5d0 * (eig_p - eig_n)/pol_amp)
              endif
           enddo
           !$omp end parallel do
           call hgm%set_units("1")
           call hgm%set_field(1, "LNI")
           call hgm%set_field(2, "LNQ")
           call hgm%set_field(3, "LNU")           
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_ln")))           
           write(*,*) trim(fin(i))//" number of negative pixels (cannot take log): "//COOP_STR_OF(fail_times)
        enddo
        call hgm%free
        goto 500                
     case("EXPIQU")
        nin  = nin - 1
        do i=1, nin
           call hgm%read(trim(fin(i)), nmaps_wanted = 3)
           !$omp parallel do private(eig_p, eig_n, pol_amp, pix)
           do pix = 0, hgm%npix-1
              pol_amp = sqrt(hgm%map(pix, 2)**2 + hgm%map(pix, 3)**2)
              eig_p = exp(hgm%map(pix, 1) + pol_amp)
              eig_n = exp(hgm%map(pix, 1) - pol_amp)
              hgm%map(pix, 1) = 0.5d0 * (eig_p + eig_n)
              hgm%map(pix, 2:3) = hgm%map(pix, 2:3) * (0.5d0 * (eig_p - eig_n)/pol_amp)
           enddo
           !$omp end parallel do
           call hgm%set_units("muK")
           call hgm%set_field(1, "I")
           call hgm%set_field(2, "Q")
           call hgm%set_field(3, "U")           
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_exp")))   
        enddo
        call hgm%free
        goto 500                
     case("MAX")
        if(inline_mode)then
           inline = coop_inputArgs(nin+1)
           read(inline, *) maxabs
        else           
           write(*,*) "Enter the maximum pixel value"
           read(*,*) maxabs
        endif
        nin  = nin - 1
        do i=1, nin
           call hgm%read(trim(fin(i)))
           hgm%map = max(min(hgm%map, maxabs), -maxabs)
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_maxcut"))) 
        enddo
        call hgm%free
        goto 500        
     case("MULTIPLY")
        if(inline_mode)then
           fout = trim(coop_inputArgs(nin+1))
        else
           fout = ""
           do while(trim(fout).eq."")
              write(*,*) "Enter the output file name: "
              read(*,'(A)') fout
           enddo           
        endif
        nin  = nin - 1
        if(nin.lt.2) stop "nmaps<2, cannot multiply"
        call hgm%read(trim(fin(1)))
        call hgm2%read(trim(fin(2)))
        if(hgm2%nside .ne. hgm%nside) stop "map with different resolution cannot be multiplied"
        if(hgm2%nmaps .eq. 1)then
           do i=1, hgm%nmaps
              hgm%map(:, i) = hgm%map(:, i) * hgm2%map(:,1)
           enddo
           call hgm%write(trim(fout))
        elseif(hgm%nmaps.eq.1)then
           do i=1, hgm2%nmaps
              hgm2%map(:, i) = hgm2%map(:, i) * hgm%map(:,1)
           enddo
           call hgm2%write(trim(fout))
        elseif(hgm%nmaps.eq.hgm2%nmaps)then
           do i=1, hgm%nmaps
              hgm%map(:, i) = hgm%map(:, i) * hgm2%map(:,i)
           enddo
           call hgm%write(trim(fout))
        else
           stop "different number of maps can not be matched"
        endif
        call hgm%free()
        call hgm2%free()
        goto 500
     case("ADD")
        if(inline_mode)then
           fout = trim(coop_inputArgs(nin+1))           
        else
           fout = ""
           do while(trim(fout).eq."")
              write(*,*) "Enter the output file name: "
              read(*,'(A)') fout
           enddo
        endif
        nin  = nin - 1
        if(nin.lt.2) stop "nmaps<2, cannot add"
        call hgm%read(trim(fin(1)))
        call hgm2%read(trim(fin(2)))
        if(hgm2%nside .ne. hgm%nside) stop "map with different resolution cannot be added"
        if(hgm%nmaps.eq.hgm2%nmaps)then
           do i=1, hgm%nmaps
              hgm%map(:, i) = hgm%map(:, i) + hgm2%map(:,i)
           enddo
           call hgm%write(trim(fout))
        else
           stop "different number of maps can not be matched"
        endif
        call hgm%free()
        call hgm2%free()
        goto 500
     case("SUBTRACT")
        if(inline_mode)then
           fout = trim(coop_inputArgs(nin+1))           
        else
           fout = ""
           do while(trim(fout).eq."")
              write(*,*) "Enter the output file name: "
              read(*,'(A)') fout
           enddo
        endif
        nin  = nin - 1
        if(nin.lt.2) stop "nmaps<2, cannot subtract"        
        call hgm%read(trim(fin(1)))
        call hgm2%read(trim(fin(2)))
        if(hgm2%nside .ne. hgm%nside) stop "map with different resolution cannot be subtracted"
        if(hgm%nmaps.eq.hgm2%nmaps)then
           do i=1, hgm%nmaps
              hgm%map(:, i) = hgm%map(:, i) - hgm2%map(:,i)
           enddo
           call hgm%write(trim(fout))
        else
           stop "different number of maps can not be matched"
        endif
        call hgm%free()
        call hgm2%free()
        goto 500
     case("I2TQUL")
        nin = nin -1
        do i=1, nin
           call hgm%read(trim(fin(i)), nmaps_wanted = 4, nmaps_to_read = 1 )
           call hgm%get_QUL()
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_conv_TQUL")))
        enddo
        print*, "maps are all converted to TQUL"
        goto 500
     case("I2TQULDD")
        nin = nin -1
        do i=1, nin
           call hgm%read(trim(fin(i)), nmaps_wanted = 6, nmaps_to_read = 1 )
           call hgm%get_QULDD()
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_conv_TQULDD")))
        enddo
        print*, "maps are all converted to TQULDD"
        goto 500
     case("I2TQUL4D")
        if(inline_mode)then
           inline = coop_inputArgs(nin+1)
           read(inline, *) fwhm
        else           
           write(*,*) "Enter FWHM in arcmin:"
           read(*,*) fwhm
        end if
        fwhm = coop_SI_arcmin*fwhm
        nin = nin -1
        do i=1, nin
           call hgm%read(trim(fin(i)), nmaps_wanted = 8, nmaps_to_read = 1 )
           call hgm%get_QUL4D(fwhm, tophat = .true.)
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_conv_TQUL4D")))
        enddo
        print*, "maps are all converted to TQUL4D"
        goto 500                        
     case("I2TQUL6D")
        nin = nin -1
        do i=1, nin
           call hgm%read(trim(fin(i)), nmaps_wanted = 10, nmaps_to_read = 1 )
           call hgm%get_QUL6D()
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_conv_TQUL6D")))
        enddo
        print*, "maps are all converted to TQUL6D"
        goto 500                
     case("I2TQTUT")
        nin = nin -1
        do i=1, nin
           call hgm%read(trim(fin(i)), nmaps_wanted = 3, nmaps_to_read = 1 )
           call hgm%get_QU()
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_conv_TQTUT")))
        enddo
        print*, "maps are all converted to TQTUT"
        goto 500
     case("IQU2TEB")
        nin = nin -1
        do i=1, nin
           call hgm%read(trim(fin(i)))
           call hgm%qu2EB()
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_conv_TEB")))
        enddo
        print*, "maps are all converted to TEB"
        goto 500
     case("IQU2ZETA")
        if(inline_mode)then
           inline = coop_inputArgs(nin+1)
           read(inline, *) fwhm
        else
           write(*,*) "Enter the fwhm in arcmin:"
           read(*,*) fwhm
        endif
        if(trim(coop_inputArgs(nin+2)).eq."SINGLE_SLICE")then
           write(*,*) "Single slice zeta not supproted now"
           stop
        else
           write(*,*) "Doing weighted zeta"
        endif
        nin = nin -1
        do i=1, nin
           call hgm%read(trim(fin(i)), nmaps_wanted = 3)
           call hgm%te2zeta(fwhm_arcmin = fwhm, want_unconstrained = .true.)
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_conv_ZETA")), index_list = (/ 1, 2, 3 /) )
        enddo
        print*, "maps are all converted to zeta"
        goto 500
     case("QU2ZETA")
        if(inline_mode)then
           inline = coop_inputArgs(nin+1)
           read(inline, *) fwhm
        else
           write(*,*) "Enter the fwhm in arcmin:"
           read(*,*) fwhm
        endif
        if(trim(coop_inputArgs(nin+2)).eq."SINGLE_SLICE")then
           write(*,*) "Single slice zeta not supported"
        else
           write(*,*) "Doing weighted zeta"
        endif
        
        nin = nin -1
        do i=1, nin
           call hgm%read(trim(fin(i)), nmaps_wanted = 3, nmaps_to_read = 2)
           call hgm%E2zeta(fwhm_arcmin = fwhm, want_unconstrained = .true.)
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_conv_ZETA")), index_list = (/ 1, 2, 3/) )
        enddo
        print*, "maps are all converted to zeta"
        goto 500
     case("T2ZETA")
        if(inline_mode)then
           inline = coop_inputArgs(nin+1)
           read(inline, *) fwhm
        else
           write(*,*) "Enter the fwhm in arcmin:"
           read(*,*) fwhm
        endif        
        if(trim(coop_inputArgs(nin+2)).eq."SINGLE_SLICE")then
           write(*,*) "Single slice zeta not supported now"
        else
           write(*,*) "Doing weighted zeta"
        endif
        nin = nin -1
        do i=1, nin
           call hgm%read(trim(fin(i)), nmaps_wanted = 3, nmaps_to_read = 1)
           call hgm%t2zeta(fwhm_arcmin = fwhm, want_unconstrained = .true.)
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_conv_ZETA")), index_list = (/ 1, 2, 3/) )
        enddo
        print*, "maps are all converted to zeta"
        goto 500                        
     case("", "SHUFFLE")
        nin = nin - 1
        exit
     end select
     if(coop_file_exists(fin(nin)))then
        nin = nin+1
     else
        write(*,*) "The file "//trim(fin(nin))//" does not exist!"
        if(inline_mode) stop
     endif
  enddo
  if(inline_mode)stop "inline mode does not support shuffle"
10 write(*,*) "Enter the output file name: "
  read(*,'(A)') fout
  if(trim(fout).eq."")goto 10
  allocate(nmaps_in(nin), ordering_in(nin), nside_in(nin))
  do i=1,nin
     npixtot = getsize_fits(trim(fin(i)), nmaps = nmaps_in(i), nside = nside_in(i), ordering =ordering_in(i))
  enddo
  if(any(nside_in .ne. nside_in(1)))then
     print*, "Error: nside of all the input maps are not the same."
     goto 500
  endif
  num_maps_wanted = 0
  npix = nside2npix(nside_in(1))
  do i=1,nin
     write(*,*) trim(coop_num2str(nmaps_in(i)))//" maps are found in file "//trim(fin(i))//", enter the indices of the maps that you want to keep (numbers between 1 to "//trim(coop_num2str(nmaps_in(i)))//" seperated by spaces):"
     read(*, "(A)") inline
     write(*,*) trim(inline)
     call coop_string_to_list(inline, indices_wanted(i))
     num_maps_wanted = num_maps_wanted + indices_wanted(i)%n
  enddo
  write(*,*) "totally "//trim(coop_num2str(num_maps_wanted))//" maps"
  call hgm%init(nside = nside_in(1), nmaps = num_maps_wanted, genre = "I")
  j = 1
  do i=1, nin
     call hgm2%read(trim(fin(i)), nmaps_wanted =  indices_wanted(i)%element(indices_wanted(i)%n))
     call hgm2%convert2ring()
     do k = 1, indices_wanted(i)%n
        hgm%map(:, j) = hgm2%map(:, indices_wanted(i)%element(k))
        j = j + 1
     enddo
  enddo
  write(*,*) "enter the genre of each map"  
  do i=1, hgm%nmaps
     write(*,*) "map "//COOP_STR_OF(i)//":"
     read(*,*) field
     call hgm%set_field(i, field)
  enddo
  write(*,*) "enter the unit of each map"  
  do i=1, hgm%nmaps
     write(*,*) "map "//COOP_STR_OF(i)//":"
     read(*,*) field
     call hgm%set_unit(i, field)
  enddo
  call coop_delete_file(trim(fout))
  call hgm%write(trim(fout))
  goto 500
450 stop "end of file"  
500 continue
#else
  stop "you need to install healpix"
#endif
end program map
