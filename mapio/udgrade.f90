program udg
  use coop_wrapper_utils
  use coop_healpix_mod
  use udgrade_nr
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  implicit none
#include "constants.h"
  character(LEN=1024)::f_map_in, f_map_out, line
  integer nside_in, nside_out, npix_in, npix_out
  real::fmissval = 0.
  logical::pessimistic = .false.
  real,dimension(:,:),allocatable::map_in, map_out
  integer ordering, nmaps
  integer(8) npixtot
  character(LEN=80)::header(64)

  write(*,*) "usage:"
  write(*,*) "./UDG map_in map_out nside_out pessimistic fmissval"
  if(iargc().lt.5) stop
  f_map_in = coop_InputArgs(1)
  f_map_out = coop_InputArgs(2) 
  if(trim(f_map_in).eq."") stop
  line = coop_InputArgs(3) 
  read(line, *) nside_out
  line = coop_InputArgs(4)
  if(trim(line).ne."")then
     read(line, *) pessimistic
  endif
  line = coop_InputArgs(5)
  if(trim(line).ne."")then
     read(line, *) fmissval
  endif
     
  npixtot = getsize_fits(trim(f_map_in), nmaps = nmaps, nside = nside_in, ordering = ordering)
  npix_in = nside2npix(nside_in)
  npix_out = nside2npix(nside_out)
  allocate(map_in(0:npix_in-1, nmaps))
  allocate(map_out(0:npix_out-1, nmaps))
  call input_map(trim(f_map_in), map_in, npix_in, nmaps, fmissval = 0.)
  if(ordering .eq. COOP_NESTED)then
     call udgrade_nest(map_in, nside_in, map_out, nside_out, fmissval, pessimistic)
  else
     call udgrade_ring(map_in, nside_in, map_out, nside_out, fmissval, pessimistic)
  endif
  call write_minimal_header(header,dtype = 'MAP', nside=nside_out, order = ordering, creator='Zhiqi Huang')
  call coop_delete_file(trim(f_map_out))
  call output_map(map_out, header, trim(f_map_out))
end program udg
