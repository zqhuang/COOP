program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  COOP_INT, parameter::nbins = 50
  COOP_REAL, parameter::maxnu = 5.d0
  type(coop_flatsky_maps)::map
  COOP_STRING::filename, output, xlabel, ylabel, caption
  COOP_INT::sig, i, loc
  type(coop_asy)::fig
  COOP_REAL::rms, mean, dx, x(-nbins:nbins), counts(-nbins:nbins)
  if(iargc() .lt. 2)then
     write(*,*) "./F1D -map MAPFILE -out OUTPUTFILE [-sig MAP_INDEX -xlabel ... -ylabel ...]"
     stop
  endif
  call coop_get_command_line_argument(key = 'map', arg = filename)
  call coop_get_command_line_argument(key = 'out', arg = output)
  call coop_get_command_line_argument(key = 'ylabel', arg = ylabel, default= '$P(x)$')
  call coop_get_command_line_argument(key = 'xlabel', arg = xlabel, default = '$x$')
  call coop_get_command_line_argument(key = 'caption', arg = caption, default  = '')
  call coop_get_command_line_argument(key = 'sig', arg = sig, default = 1)

  call map%read(filename)
  if(sig .lt. 1 .or. sig .gt. map%nmaps)then
     write(*,*) "-sig out of range"
     stop
  endif
  call fig%open(output)
  call fig%init(xlabel = trim(xlabel), ylabel = trim(ylabel), caption = trim(caption))

  mean = map%weighted_sum(map%map(sig)%image)/map%total_weight
  map%map(sig)%image =   map%map(sig)%image - mean
  rms =  sqrt(map%weighted_sum(map%map(sig)%image**2)/map%total_weight)
  map%map(sig)%image =  map%map(sig)%image/rms
  dx = maxnu/nbins
  counts = 0.d0
  x = dx* (/ (i, i = -nbins, nbins) /)
  do i=0, map%npix-1
     if(map%unmasked(i))then
        loc = nint(map%map(sig)%image(i)/dx)
        if(loc .ge. -nbins .and. loc.le. nbins)then
           counts(loc) = counts(loc)+ 1.d0
        endif
     endif
  enddo
  counts = counts/dx/map%total_weight
  call fig%plot(x = x, y = counts, color = "darkred", linewidth = 2., linetype = "solid")
  counts = exp(-x**2/2.d0)/sqrt(coop_2pi)
  call fig%plot(x = x, y = counts, color = "skyblue", linewidth = 2., linetype = "dotted")
  call fig%close()  
end program test
