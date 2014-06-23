program map
  use coop_healpix_mod
  use coop_wrapper_utils
  use healpix_types
  use alm_tools
  use pix_tools
  use head_fits
  use fitstools
  implicit none
#include "constants.h"
  
  integer,parameter::nmax = 8
  character(len=80), dimension(64) :: header
  character(LEN=1024),dimension(nmax)::fin
  character(LEN=1024)::fout
  integer,dimension(:),allocatable::nmaps_in, ordering_in, nside_in
  integer nin, i, j
  real,dimension(:,:),allocatable:: map_in  
  real,dimension(:,:),allocatable:: map_out, map_tmp
  integer(8) npixtot
  type(coop_list_integer),dimension(nmax)::indices_wanted
  character(LEN=1024)::inline
  integer num_maps_wanted, npix, k
  type(coop_healpix_maps) hgm, hgm2
  real*8 fwhm, scal
  write(*,*) "options are: SPLIT; SMOOTH; MULTIPLY;IQU2TQTUT;IQU2TEB;SCALE"
  nin = 1
  do while(nin .le. nmax)
     write(*,*) "Enter input file and press Enter (or just press Enter key to finish):"
     read(*,'(A)') fin(nin)
     select case(trim(fin(nin)))
     case("SPLIT")
        nin = nin -1
        do i=1, nin
           call coop_healpix_split(trim(fin(i)))
        enddo
        print*, "maps are all split"
        goto 500
     case("SMOOTH")
        write(*,*) "Enter the fwhm in arcmin:"
        read(*,*) fwhm
        fwhm = fwhm*coop_SI_arcmin
        nin = nin -1
        do i=1, nin
           call coop_healpix_smooth_mapfile(trim(fin(i)), fwhm)
        enddo
        print*, "maps are all smoothed"
        goto 500
     case("SCALE")
        write(*,*) "Enter the scale"
        read(*,*) scal
        nin  = nin - 1
        do i=1, nin
           call hgm%read(trim(fin(i)))
           hgm%map = hgm%map * scal
           call hgm%write(trim(fin(i)))
        enddo
        call hgm%free
        goto 500
     case("MULTIPLY")
        nin  = nin - 1
        if(nin.lt.2) stop "nmaps<2, cannot multiply"
        fout = ""
        do while(trim(fout).eq."")
           write(*,*) "Enter the output file name: "
           read(*,'(A)') fout
        enddo
        call hgm%read(trim(fin(1)), nmaps_wanted = 1)
        do i=2, nin
           call hgm2%read(trim(fin(i)), nmaps_wanted = 1)
           if(hgm2%nside .ne. hgm%nside) stop "map with different resolution cannot be multiplied"
           
           hgm%map(:, 1) = hgm%map(:, 1) * hgm2%map(:,1)
        enddo
        call hgm%write(trim(fout))
        call hgm%free()
        call hgm2%free()
        goto 500
     case("IQU2TQTUT")
        nin = nin -1
        do i=1, nin
           call hgm%read(trim(fin(i)))
           call hgm%iqu2TQTUT()
           call hgm%write(trim(coop_file_add_postfix(fin(i), "converted_to_TQTUT")))
        enddo
        print*, "maps are all converted to TQTUT"
        goto 500
     case("IQU2TEB")
        nin = nin -1
        do i=1, nin
           call hgm%read(trim(fin(i)))
           call hgm%iqu2TEB()
           call hgm%write(trim(coop_file_add_postfix(fin(i), "converted_to_TEB")))
        enddo
        print*, "maps are all converted to TEB"
        goto 500
     case("")
        nin = nin - 1
        exit
     end select
     if(coop_file_exists(fin(nin)))then
        nin = nin+1
     else
        write(*,*) "The file "//trim(fin(nin))//" does not exist!"
     endif
  enddo
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
500 continue
end program map
