program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, meanmap, V0map, rmsmap, mask
  COOP_REAL::r_deg, numin, numax
  COOP_REAL,dimension(:),allocatable::nu, V0
  COOP_REAL::global_mean, global_rms, summ, rat, Amp
  COOP_INT::nside, nnu, i
  COOP_REAL:: nuc, Gauss_dAdnu, dAdnu
  COOP_STRING::map_file, prefix, mask_file
  type(coop_dynamic_array_integer)::inds
  type(coop_file)::fp
  COOP_INT::nlist
  COOP_INT,allocatable::listpix(:)
  logical::gf
  call coop_MPI_init()
  if(iargc().le.0)then
     write(*,"(A)") "----------------------------------------------------------"     
     write(*,"(A)") "Syntax:"
     write(*,"(A)") "./MM0 -map MAP_FILE -prefix PREFIX"
     write(*,"(A)") "----------------------------------------------------------"
     write(*,"(A)") "other optional inputs are:"
     write(*,"(A)") "-mask -MASK_FILE"     
     write(*,"(A)") "-numin NU_MIN"
     write(*,"(A)") "-numax NU_MAX"
     write(*,"(A)") "-nnu N_NU"
     write(*,"(A)") "-radius RAIDUS_IN_DEGREE [0.]"
     write(*,"(A)") "-nside NSIDE [2]"
     write(*,"(A)") "-gf [T|F]"
     stop
  endif
  call coop_random_init()
  call coop_get_command_line_argument(key = "map", arg = map_file)
  call coop_get_command_line_argument(key = "mask", arg = mask_file, default = "NONE")  
  call coop_get_command_line_argument(key = "numin", arg = numin, default = -2.1d0)
  call coop_get_command_line_argument(key = "numax", arg = numax, default = 2.1d0)
  call coop_get_command_line_argument(key = "nnu", arg = nnu, default = 22)
  call coop_get_command_line_argument(key = "radius", arg = r_deg, default = 0.d0)
  call coop_get_command_line_argument(key = "nside", arg = nside, default = 2)
  call coop_get_command_line_argument(key = "prefix", arg = prefix)
  call coop_get_command_line_argument(key = "gf", arg = gf, default = .false.)  
  
  call map%read(map_file)
  allocate(nu(nnu), V0(nnu))
  call coop_set_uniform(nnu, nu, numin, numax)

  if(nside .gt. 0)then
     call meanmap%init(nside = nside, nmaps = 1, genre = "I", nested = .true.)
     call rmsmap%init(nside = nside, nmaps = 1, genre = "I", nested = .true.)
     call V0map%init(nside = nside, nmaps = nnu, genre = "I", nested = .true.)




     call map%scan_local_minkowski0(1, nu, meanmap, rmsmap, V0map, r_deg, do_gaussian_fit = gf)  
!!$     call meanmap%write(trim(adjustl(prefix))//"_mean.fits")
!!$     call rmsmap%write(trim(adjustl(prefix))//"_rms.fits")
!!$     call V0map%write(trim(adjustl(prefix))//"_V0.fits")

     if(trim(mask_file) .ne. "NONE")then
        call mask%read(mask_file, nested = .true.)
        if(mask%nside .ne. V0map%nside)stop "mask has the wrong nside"
        call V0map%apply_mask(mask)     
        do i=1, nnu
           V0(i) = sum(V0map%map(:, i))/sum(mask%map(:, 1))
        enddo
     else
        if(mask%nside .ne. V0map%nside)then
           do i=1, nnu
              V0(i) = sum(V0map%map(:, i))/V0map%npix
           enddo
        endif
     endif
  else
     if(trim(mask_file) .ne. "NONE")then
        call mask%read(mask_file, nested = .true.)
        if(mask%nside .ne. map%nside)then
           write(*,*) "For nside = 0 you need to make sure mask and map have the same nside"
           stop
        endif
        call map%apply_mask(mask)
        summ = sum(dble(mask%map(:, 1)))
        if(gf)then
           call inds%get_indices( mask%map(:, 1) .gt. 0.5,  start_index = 0)
           call coop_fit_gaussian(dble(map%map(inds%i, 1)),  max(inds%n/200, 20), global_mean, global_rms, Amp)
        else
           global_mean = sum(dble(map%map(:, 1)))/summ
           global_rms = sqrt(sum(dble((map%map(:, 1) - global_mean)**2*mask%map(:,1)))/summ)
        endif
        where (mask%map(:, 1) .lt. 0.5)
           map%map(:, 1) = global_mean - global_rms * 100. !!just put out of the box
        end where
     else
        summ = map%npix
        if(gf)then
           call coop_fit_gaussian(dble(map%map(:,1)), max(map%npix/200, 20), global_mean, global_rms, Amp)
        else
           global_mean = sum(dble(map%map(:,1)))/summ
           global_rms = sqrt(sum(dble((map%map(:,1)-global_mean)**2))/summ)
        endif
     endif
     do i=1, nnu
        V0(i) = count(map%map(:,1) .gt. global_mean + global_rms * nu(i))/summ
     enddo
  endif
  call fp%open(trim(adjustl(prefix))//"_V0.txt", "w")
  do i=1, nnu-1
     nuc = (nu(i)+nu(i+1))/2.d0
     gauss_dAdnu = exp(-nuc**2/2.d0)/sqrt(coop_2pi)
     dAdnu = -(V0(i+1)-V0(i))/(nu(i+1)-nu(i))
     rat = dAdnu/gauss_dAdnu
     if(dAdnu .le. 0.)then
        dAdnu = gauss_dAdnu
     endif
     write(fp%unit, "(3E16.7)") nuc, dAdnu*log(dAdnu/gauss_dAdnu), rat
  enddo
  call fp%close()
  call coop_MPI_finalize()
end program test
