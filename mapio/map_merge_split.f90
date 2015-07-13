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
  integer,parameter::nmax = 8
  character(len=80), dimension(64) :: header
  COOP_STRING,dimension(nmax)::fin
  COOP_STRING::fout, fbeam
  integer,dimension(:),allocatable::nmaps_in, ordering_in, nside_in
  integer nin, i, j, l, il, l1, l2
  real,dimension(:,:),allocatable:: map_in  
  real,dimension(:,:),allocatable:: map_out, map_tmp
  integer(8) npixtot
  type(coop_list_integer),dimension(nmax)::indices_wanted
  COOP_STRING::inline
  integer num_maps_wanted, npix, k
  logical::inline_mode = .false.
  type(coop_healpix_maps) hgm, hgm2
  COOP_REAL fwhm, scal, threshold
  COOP_SINGLE::sigmas(100)
  type(coop_file)::fp

  nin = 1
  inline_mode =  (iargc() .gt. 0)
  if(.not. inline_mode)then
     write(*,*) "options are: SPLIT; SMOOTH; DOBEAM; MULTIPLY;I2TQTUT;I2TQUL;IQU2TEB;T2ZETA; IQU2ZETA; QU2ZETA; SCALE;INFO;ADD;SUBTRACT;MAKEMASK; SHUFFLE; HIGHPASS; LOWPASS"
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
           call hgm%write(trim(fin(i)))
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
        elseif(hgm2%nmaps.eq.1)then
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
           hgm2 = hgm
           call hgm2%map2alm(index_list = (/ 1 /) )
           hgm2%alm(0:1, :, 1) = 0.
           do l = 2, hgm2%lmax
              hgm2%alm(l, :, 1) = hgm2%alm(l, :, 1)*(l*(l+1.d0))
           enddo
           call hgm2%alm2map(index_list = (/ 1 /) )
           call hgm2%get_QU()
           hgm%map(:, 2:3) = hgm2%map(:,2:3)
           hgm%map(:, 4) = hgm2%map(:, 1)
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_converted_to_TQUL")))
        enddo
        print*, "maps are all converted to TQUL"
        goto 500
     case("I2TQTUT")
        nin = nin -1
        do i=1, nin
           call hgm%read(trim(fin(i)), nmaps_wanted = 3, nmaps_to_read = 1 )
           call hgm%get_QU()
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_converted_to_TQTUT")))
        enddo
        print*, "maps are all converted to TQTUT"
        goto 500
     case("IQU2TEB")
        nin = nin -1
        do i=1, nin
           call hgm%read(trim(fin(i)))
           call hgm%qu2EB()
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_converted_to_TEB")))
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
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_converted_to_ZETA")), index_list = (/ 1, 2, 3 /) )
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
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_converted_to_ZETA")), index_list = (/ 1, 2, 3/) )
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
           call hgm%write(trim(coop_file_add_postfix(fin(i), "_converted_to_ZETA")), index_list = (/ 1, 2, 3/) )
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
     write(*,*) trim(coop_num2str(nmaps_in(i)))//" maps are found in file "//trim(fin(i))//", enter the indices of the maps that you want to keep (numbers between 1 to "//trim(coop_num2str(nmaps_in(i)))//" seperated with spaces):"
     read(*, "(A)") inline
     write(*,*) trim(inline)
     call coop_string_to_list(inline, indices_wanted(i))
     num_maps_wanted = num_maps_wanted + indices_wanted(i)%n
  enddo
  write(*,*) "totally "//trim(coop_num2str(num_maps_wanted))//" maps"
  allocate(map_out(0:npix-1, num_maps_wanted))
  j = 1
  do i=1, nin
     allocate(map_tmp(0:npix-1, nmaps_in(i)))
     call input_map(trim(fin(i)), map_tmp, npix, nmaps_in(i), fmissval = 0.)
     if(ordering_in(i) .eq. COOP_NESTED) call convert_nest2ring(nside_in(i), map_tmp)
     do k = 1, indices_wanted(i)%n
        map_out(:, j) = map_tmp(:, indices_wanted(i)%element(k))
        j = j + 1
     enddo
     deallocate(map_tmp)
  enddo
  call write_minimal_header(header,dtype = 'MAP', nside=nside_in(1), order=COOP_RING, creator='Zhiqi Huang') !!ring order
  call coop_delete_file(trim(fout))
  call output_map(map_out, header, trim(fout))
  goto 500
450 stop "end of file"  
500 continue
#else
  stop "you need to install healpix"
#endif
end program map
