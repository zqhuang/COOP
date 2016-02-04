program TestNpeak
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_asy)::fig
  COOP_STRING::file, output, caption, xlabel, ylabel, legend, interp, attach
  character(len = 5) tmp
  COOP_SHORT_STRING::color, linetype
  COOP_SINGLE::linewidth, legendx, legendy
  COOP_INT::i, legendcols
  COOP_STRING::xcol, ycol
  logical::has_legend
  if(iargc().lt.2)then
     write(*,*) "============================================================================"
     write(*,*) "Syntax:"
     write(*,*) "./MPLOT -xlabel XLABEL -ylabel YLABEL -file1 FILENAME1 [-color1 COLOR1 -linetype1 LINETYPE1 -linewidth1 LINEWIDTH1 -legend1 LENGEND1 -interp1 INTERPOLATE METHOD1 -xcol1 COLUMN_FOR_X1 -ycol1 COLUMN_FOR_Y1] -out OUTPUT_FILENAME"
     write(*,*) "============================================================================"
     write(*,*) "Other optional inputs are:"
     write(*,*) "-xlog F/T"
     write(*,*) "-ylog F/T"
     write(*,*) "-xmin XMIN"
     write(*,*) "-xmax XMAX"
     write(*,*) "-ymin YMIN"
     write(*,*) "-ymax YMAX"
     write(*,*) "-caption CAPTION"
     write(*,*) "-width FIGURE_WIDTH_IN_INCHES"
     write(*,*) "-height FIGURE_HEIGHT_IN_INCHES"
     write(*,*) "-clip  CLIP_THE_FIGURE (only when xmin, xmax, ymin, ymax are given)"
     write(*,*) "-legendx LEGEND_X_POSITION (from 0 to 1)"
     write(*,*) "-legendy LEGEND_Y_POSITION (from 0 to 1)"
     write(*,*) "-legendcols LEGEND_COLUMNS (1 or 2)"
     write(*,*) "-file2 FILENAME2 [-color2 ... -linetype2 ... -linewidth2 ...-legend2 ...-interp2 ...]"
     write(*,*) "-file3 FILENAME3 [-color3 ... -linetype3 ... -linewidth3 ...-legend3 ... -interp3 ...]"
     write(*,*) "...up to 12 files"
     write(*,*) "-attach1 ATTACH_FILENAME1"
     write(*,*) "-attach2 ATTACH_FILENAME2"
     write(*,*) "-attach3 ATTACH_FILENAME3"
     write(*,*) "...up to 12 files"
     write(*,*) "============================================================================"
     write(*,*) "xlabel and ylabel are latex expressions, for example: -xlabel '$x$' -ylabel 'a function of x' "
     write(*,*) "color can be red, green, ... or RGB:255:0:100 or GRAY:134"
     write(*,*) "linetype can be solid, dashed, dotted, dashdotted, longdashed, longdashdotted"
     write(*,*) "linewidth can be any real number of order 1."
     write(*,*) "interp can be LinearLinear, LinearLog, LogLinear, LogLog"
     write(*,*) "You can use gnuplot style for ycol, such as -ycol '$2 * $3 + $4'"
     write(*,*) "============================================================================"
     write(*,*) "Usage for -attach1, -attach2..."
     write(*,*) "Sometimes you want to plot columns in two files. "
     write(*,*) "Example: you have a file X with 4 columns and file Y with 3 columns. With the -attach argument you can combine -file1 with -attache1 to (virtually) have a combined file with 7 columns."
     write(*,*) "If you want to plot the 2nd column of X against the 1st column of Y (which is the 5th column in the combined file), you do:"
     write(*,*) "MPLOT -file1 X -attach1 Y -xcol 2 -ycol 5"
     write(*,*) "If you want to plot the sum of the 2nd column of X and the 3rd column of Y (that is the 7th column in the combined file) against the 1st column of X. You do:"
     write(*,*) "MPLOT -file1 X -attach1 Y -xcol1 1 -ycol1 '$2 + $7'"
     write(*,*) "============================================================================"
     stop
  endif
  call coop_get_command_line_argument(key ="out" , arg = output, default="coop_figure.txt")
  if(trim(output).eq."")then
     stop "argument -out must be a file name"
  endif
  call fig%open(output)
  call coop_get_command_line_argument(key = "width", arg = fig%width, default =  coop_asy_default_width)
  call coop_get_command_line_argument(key ="height", arg =  fig%height, default =  coop_asy_default_height)
  write(fig%unit, "(2F12.2)") fig%width, fig%height
  call coop_get_command_line_argument(key ="caption" , arg = caption, default= "" )
  if(trim(adjustl(caption)) .eq."")then
     write(fig%unit, "(A)") "NULL"
  else
     write(fig%unit, "(A)") trim(adjustl(caption))
  endif
  call coop_get_command_line_argument(key ="xlabel" , arg = xlabel , default= "")
  call coop_get_command_line_argument(key ="ylabel" , arg = ylabel , default= "")
  if(trim(adjustl(xlabel)).eq."")then
     write(fig%unit, "(A)") "NULL"
  else
     write(fig%unit, "(A)") trim(adjustl(xlabel))
  endif
  if(trim(adjustl(ylabel)).eq."")then
     write(fig%unit, "(A)") "NULL"
  else
     write(fig%unit, "(A)") trim(adjustl(ylabel))
  endif

  call coop_get_command_line_argument(key ="xlog" , arg = fig%xlog, default= .false. )
  call coop_get_command_line_argument(key ="ylog" , arg = fig%ylog, default= .false.)
  if(fig%xlog)then
     tmp(1:2) = "1 "
  else
     tmp(1:2) = "0 "
  endif
  if(fig%ylog)then
     tmp(3:4) = "1 "
  else
     tmp(3:4) = "0 "
  endif
  tmp(5:5) = "0"
  write(fig%unit, "(A)") tmp
  call coop_get_command_line_argument(key ="clip" , arg =  fig%clip, default=.false.)
  call coop_get_command_line_argument(key = "xmin", arg=fig%xmin, default = 1.1e30)
  fig%adjust_xmin = (abs(fig%xmin) .gt. 1.e30)
  call coop_get_command_line_argument(key = "xmax", arg=fig%xmax, default = -1.1e30)
  fig%adjust_xmax = (abs(fig%xmax) .gt. 1.e30)

  call coop_get_command_line_argument(key = "ymin", arg=fig%ymin, default = 1.1e30)
  fig%adjust_ymin = (abs(fig%ymin) .gt. 1.e30)
  call coop_get_command_line_argument(key = "ymax", arg=fig%ymax, default = -1.1e30)
  fig%adjust_ymax = (abs(fig%ymax) .gt. 1.e30)
  if(fig%adjust_xmin .or. fig%adjust_xmax .or. fig%adjust_ymin .or. fig%adjust_ymax) fig%clip = .false.

  if(fig%clip)then
     write(fig%unit, "(A)") "1"
  else
     write(fig%unit, "(A)") "0"
  endif

  write(fig%unit,"(2G14.5)") fig%xmin, fig%xmax
  write(fig%unit,"(2G14.5)") fig%ymin, fig%ymax
  write(fig%unit,"(I5)") 0   !!0 means any number of blocks
  fig%color(1) = "black"
  fig%color(2) = "red"
  fig%color(3) = "blue"
  fig%color(4) = "cyan"
  fig%color(5) = "violet"
  fig%color(6) = "orange"
  fig%color(7) = "skyblue"
  fig%color(8) = "gray"
  fig%color(9) = "brown"
  fig%color(10) = "green"
  fig%color(11) = "pink"
  fig%color(12) = "yellow"
  fig%linewidth(1:3) = 2.
  fig%linewidth(4:6) = 1.5
  fig%linewidth(7:12) = 1.
  fig%linetype(1) = "solid"
  fig%linetype(2) = "dotted"
  fig%linetype(3) = "dashed"
  fig%linetype(4) = "dashdotted"
  fig%linetype(5) = "longdashed"
  fig%linetype(6) = "longdashdotted"
  fig%linetype(7:12) = "solid"
  has_legend = .false.
  do i = 1, 12
     call coop_get_command_line_argument(key="file"//COOP_STR_OF(i), arg=file, default="")
     call coop_get_command_line_argument(key="attach"//COOP_STR_OF(i), arg=attach, default="")
     if(trim(file).eq."")cycle
     call coop_get_command_line_argument(key="color"//COOP_STR_OF(i), arg=color, default = fig%color(i))
     call coop_get_command_line_argument(key="linetype"//COOP_STR_OF(i), arg=linetype, default=fig%linetype(i))
     call coop_get_command_line_argument(key="linewidth"//COOP_STR_OF(i), arg=linewidth, default=fig%linewidth(i))
     call coop_get_command_line_argument(key="legend"//COOP_STR_OF(i), arg=legend, default="")
     call coop_get_command_line_argument(key="interp"//COOP_STR_OF(i), arg=interp, default="")
     call coop_get_command_line_argument(key="xcol"//COOP_STR_OF(i), arg=xcol, default="")
     call coop_get_command_line_argument(key="ycol"//COOP_STR_OF(i), arg=ycol, default="")
     if(trim(attach).ne."")then
        if(trim(legend).eq."")then
           call fig%plot_file(filename = file, interpolate = interp, xcol= xcol, ycol = ycol, color=color, linetype=linetype, linewidth= linewidth, filename2=attach)
        else
           call fig%plot_file(filename = file, interpolate = interp, xcol = xcol, ycol = ycol, color=color, linetype = linetype, linewidth = linewidth, legend = legend, filename2 = attach)
           has_legend = .true.
        endif
     else
        if(trim(legend).eq."")then
           call fig%plot_file(filename = file, interpolate = interp, xcol= xcol, ycol = ycol, color=color, linetype=linetype, linewidth= linewidth)
        else
           call fig%plot_file(filename = file, interpolate = interp, xcol = xcol, ycol = ycol, color=color, linetype = linetype, linewidth = linewidth, legend = legend)
           has_legend = .true.
        endif
     endif
  enddo
  if(has_legend)then
     call coop_get_command_line_argument(key="legendx", arg= legendx, default = 0.5)
     call coop_get_command_line_argument(key="legendy", arg=legendy, default = 0.92)
     call coop_get_command_line_argument(key="legendcols", arg=legendcols, default = 1)
     call fig%legend(legendx, legendy, legendcols)
  endif
  call fig%close()
  call system("fasy.sh "//trim(output))
end program TestNpeak
