program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, meanmap, V1map, rmsmap, mask
  COOP_REAL::r_deg, numin, numax
  COOP_REAL,dimension(:),allocatable::nu, V1
  COOP_REAL::global_mean, global_rms, summ
  COOP_INT::nside, nnu, i
  COOP_REAL:: nuc, Gauss_dAdnu, dAdnu
  COOP_STRING::map_file, prefix, mask_file
  type(coop_file)::fp
  call coop_MPI_init()
  if(iargc().le.0)then
     write(*,"(A)") "----------------------------------------------------------"     
     write(*,"(A)") "Syntax:"
     write(*,"(A)") "./MM1 -map MAP_FILE -prefix PREFIX"
     write(*,"(A)") "----------------------------------------------------------"
     write(*,"(A)") "other optional inputs are:"
     write(*,"(A)") "-mask -MASK_FILE"     
     write(*,"(A)") "-numin NU_MIN"
     write(*,"(A)") "-numax NU_MAX"
     write(*,"(A)") "-nnu N_NU"
     write(*,"(A)") "-radius RAIDUS_IN_DEGREE [0.]"
     write(*,"(A)") "-nside NSIDE [2]"     
     stop
  endif
  call coop_random_init()
  call coop_get_command_line_argument(key = "map", arg = map_file)
  call coop_get_command_line_argument(key = "mask", arg = mask_file, default = "NONE")  
  call coop_get_command_line_argument(key = "numin", arg = numin, default = 0.d0)
  call coop_get_command_line_argument(key = "numax", arg = numax, default = 2.d0)
  call coop_get_command_line_argument(key = "nnu", arg = nnu, default = 20)
  call coop_get_command_line_argument(key = "radius", arg = r_deg, default = 0.d0)
  call coop_get_command_line_argument(key = "nside", arg = nside, default = 2)
  call coop_get_command_line_argument(key = "prefix", arg = prefix)
  
  call map%read(map_file)
  allocate(nu(nnu), V1(nnu))
  call coop_set_uniform(nnu, nu, numin, numax)

  if(nside .gt. 0)then
     call meanmap%init(nside = nside, nmaps = 1, genre = "I", nested = .true.)
     call rmsmap%init(nside = nside, nmaps = 1, genre = "I", nested = .true.)
     call V1map%init(nside = nside, nmaps = nnu, genre = "I", nested = .true.)




     call map%scan_local_minkowski0(1, nu, meanmap, rmsmap, V1map, r_deg)  
     call meanmap%write(trim(adjustl(prefix))//"_mean.fits")
     call rmsmap%write(trim(adjustl(prefix))//"_rms.fits")
     call V1map%write(trim(adjustl(prefix))//"_V1.fits")

     if(trim(mask_file) .ne. "NONE")then
        call mask%read(mask_file, nested = .true.)
        if(mask%nside .ne. V1map%nside)stop "mask has the wrong nside"
        call V1map%apply_mask(mask)     
        do i=1, nnu
           V1(i) = sum(V1map%map(:, i))/sum(mask%map(:, 1))
        enddo
     else
        if(mask%nside .ne. V1map%nside)then
           do i=1, nnu
              V1(i) = sum(V1map%map(:, i))/V1map%npix
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
        global_mean = sum(dble(map%map(:, 1)))/summ
        global_rms = sqrt(sum(dble((map%map(:, 1) - global_mean)**2*mask%map(:,1)))/summ)
        where (mask%map(:, 1) .lt. 0.5)
           map%map(:, 1) = global_mean - global_rms * 100. !!just put out of the box
        end where
     else
        summ = map%npix
        global_mean = sum(dble(map%map(:,1)))/summ
        global_rms = sqrt(sum(dble((map%map(:,1)-global_mean)**2))/summ)
     endif
     do i=1, nnu
        V1(i) = count(map%map(:,1) .gt. global_mean + global_rms * nu(i))/summ
     enddo
  endif
  call fp%open(trim(adjustl(prefix))//".txt", "w")
  do i=1, nnu-1
     nuc = (nu(i)+nu(i+1))/2.d0
     gauss_dAdnu = exp(-nuc**2/2.d0)/sqrt(coop_2pi)
     dAdnu = -(V1(i+1)-V1(i))/(nu(i+1)-nu(i))     
     write(fp%unit, *) nuc, dAdnu*log(dAdnu/gauss_dAdnu)
  enddo
  call fp%close()
  call coop_MPI_finalize()
end program test
