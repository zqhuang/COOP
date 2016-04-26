program PlotF
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_asy)::fig
  COOP_INT,parameter::nsamples = 720
  COOP_STRING::xlabel, ylabel, filename, fcolor, bcolor, btype, output, caption, xvar, yvar, legend
  COOP_LONG_STRING::line
  COOP_SINGLE::bwidth, xsize, ysize,  xmin,xmax, ymin,ymax, xlegend, ylegend
  COOP_INT::i, ixvar, iyvar, nkeys, ikey, ix, iy, ncontours, ic, itheta, npt, lcols
  COOP_REAL,dimension(:),allocatable::covline
  COOP_REAL::cov(2,2), mean(2), eigs(2)
  COOP_REAL::x(nsamples), y(nsamples), r, rho, theta(nsamples), sqrteigs(2), vec(2)
  COOP_REAL::cls(3) = (/ 0.683d0, 0.954d0, 0.997d0 /)
  type(coop_file)::fp
  type(coop_list_string)::ls
  logical::has_legend 
  if(iargc() < 2)then
     write(*,*) "Syntax:"
     write(*,*) "./PLOTF -out OUTPUT -xvar X -yvar Y -cov1 COVFILE1 -color1 COLOR1  [-linetype1 LINETYPE1 -linewidth1 LINEWIDTH1 -legend1 LEGEND1 -fillcolor1 FILLCOLOR1 -width FIG_WIDTH_INCH -height FIG_HEIGHT_INCH -xlabel XLABEL -ylabel YLABE -caption CAPTION -ncontours NUM_CONTOURS(default 2) -cov2 ... -color2 ...-xlegend LEGEND_X_RELATIVE_POSITION[0-1] -ylegend LEGEND_Y_RELATIVE_POSITION[0-1] -legend_cols LEGEND_COLUMNS -xmin XMIN -xmax XMAX -ymin YMIN -ymax YMAX] "
     write(*,*) "Examples of colors:  black,  skyblue, red, RGB:255:100:100,  GRAY:120"
     write(*,*) "Examples of linetype: solid, dashed, dotted"

     stop
  endif
  call coop_set_uniform(nsamples, theta, coop_pi/nsamples, coop_2pi*(nsamples-0.5)/nsamples)
  call coop_get_command_line_argument(key = 'out', arg = output)
  call coop_get_command_line_argument(key = 'xvar', arg = xvar)
  call coop_get_command_line_argument(key = 'yvar', arg = yvar)
  call coop_get_command_line_argument(key = 'ncontours', arg = ncontours, default = 2)
  call coop_get_command_line_argument(key = 'caption', arg = caption, default='NULL')
  call coop_get_command_line_argument(key = 'xlabel', arg = xlabel, default='NULL')
  call coop_get_command_line_argument(key = 'ylabel', arg = ylabel, default='NULL')
  call coop_get_command_line_argument(key = 'width', arg = xsize, default= 8.)
  call coop_get_command_line_argument(key = 'height', arg = ysize, default= 6.)
  call coop_get_command_line_argument(key = 'xmin', arg = xmin, default= 1.1e31)
  call coop_get_command_line_argument(key = 'xmax', arg = xmax, default= -1.1e31)
  call coop_get_command_line_argument(key = 'ymin', arg = ymin, default= 1.1e31)
  call coop_get_command_line_argument(key = 'ymax', arg = ymax, default= -1.1e31)
  call coop_get_command_line_argument(key = 'xlegend', arg = xlegend, default= 0.5)
  call coop_get_command_line_argument(key = 'ylegend', arg = ylegend, default= 0.95)
  call coop_get_command_line_argument(key = 'legend_cols', arg = lcols, default = 1)

  call fig%open(output)
  call fig%init(xlabel = xlabel, ylabel = ylabel, caption = caption, width= xsize, height = ysize)
  has_legend = .false.
  do i = 1, 5
     call coop_get_command_line_argument(key = 'cov'//COOP_STR_OF(i), arg = filename, default = '')
     if(trim(filename) == '')cycle
     call coop_get_command_line_argument(key = 'fillcolor'//COOP_STR_OF(i), arg = fcolor, default = "invisible")
     call coop_get_command_line_argument(key = 'color'//COOP_STR_OF(i), arg = bcolor)
     call coop_get_command_line_argument(key = 'legend'//COOP_STR_OF(i), arg = legend, default = '')
     call coop_get_command_line_argument(key = 'linetype'//COOP_STR_OF(i), arg = btype, default="solid")
     call coop_get_command_line_argument(key = 'linewidth'//COOP_STR_OF(i), arg = bwidth, default=1.)
     call fp%open(filename)
     read(fp%unit, "(A)")line
     nkeys = (len_trim(line)-1)/16 + 1
     ixvar = 0
     iyvar = 0
     do ikey = 1, nkeys
        if(trim(line((ikey-1)*16+1:ikey*16)) == trim(xvar)) ixvar = ikey
        if(trim(line((ikey-1)*16+1:ikey*16)) == trim(yvar)) iyvar = ikey
     enddo
     if(ixvar .eq. 0)then
        write(*,*) "the variable"//trim(xvar)//" is not found in file "//trim(filename)
        stop
     endif
     if(iyvar .eq. 0)then
        write(*,*) "the variable "//trim(yvar)//" is not found in file "//trim(filename)
        stop
     endif
     allocate(covline(nkeys))
     read(fp%unit, *) covline
     mean(1) = covline(ixvar)
     mean(2) = covline(iyvar)
     do ikey   = 1, nkeys
        read(fp%unit, *) covline
        if(ikey .eq. ixvar)then
           cov(1,1) = covline(ixvar)*1.000001
           cov(1,2) = covline(iyvar)
        endif
        if(ikey .eq. iyvar)then
           cov(2,1) = covline(ixvar)
           cov(2,2) = covline(iyvar)*1.0000001
        endif
     enddo
     deallocate(covline)
     call fp%close()
     call coop_matsym_sqrt_small(2, cov)
     sqrteigs = sqrt(eigs)
     do ic = 1, ncontours
        r = sqrt(-2.d0*log(coop_IncompleteGamma(0.5d0, dble(ic)**2/2.d0)/coop_sqrtpi))
        do itheta  = 1, nsamples
           vec(1) = r*cos(theta(itheta))
           vec(2) = r*sin(theta(itheta))
           vec = matmul(cov, vec) + mean
           x(itheta) = vec(1) 
           y(itheta) = vec(2)
        enddo
        if(ic .eq. 1 .and. legend .ne.'')then
           call fig%contour(x = x, y = y, colorfill = fcolor, linecolor = bcolor, linetype = btype, linewidth = bwidth, legend=legend)
           has_legend = .true.
        else
           call fig%contour(x = x, y = y, colorfill = fcolor, linecolor = bcolor, linetype = btype, linewidth = bwidth)
        endif
     enddo
  enddo
  if(has_legend)call fig%legend(xratio = xlegend, yratio = ylegend, cols = lcols)
  call fig%close()
end program PlotF
