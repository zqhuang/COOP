module coop_asy_mod
  use coop_wrapper_typedef
  use coop_interpolation_mod
  use coop_file_mod
  use coop_list_mod
  implicit none
  private

  public::coop_asy, coop_asy_path, coop_asy_error_bar, coop_asy_interpolate_curve, coop_asy_gray_color, coop_asy_rgb_color, coop_asy_label, coop_asy_legend, coop_asy_dot, coop_asy_line, coop_asy_labels, coop_asy_dots, coop_asy_lines, coop_asy_contour, coop_asy_curve, coop_asy_density,  coop_asy_topaxis, coop_asy_rightaxis, coop_asy_clip, coop_asy_plot_function, coop_asy_plot_likelihood, coop_asy_curve_from_file, coop_asy_path_from_array, coop_asy_histogram, coop_asy_band


#include "constants.h"

  COOP_INT, parameter::sp = coop_single_real_length
  COOP_INT, parameter::dl = coop_real_length
  COOP_SINGLE::coop_asy_default_width = 4.8
  COOP_SINGLE::coop_asy_default_height = 3.9
  integer, parameter::coop_asy_num_line_types = 12

  type, extends(coop_file) :: coop_asy
     COOP_SINGLE  xmin, xmax, ymin, ymax, width, height
     logical::xlog, ylog
     COOP_SHORT_STRING, dimension(coop_asy_num_line_types)::color
     COOP_SHORT_STRING, dimension(coop_asy_num_line_types)::linetype
     COOP_SINGLE , dimension(coop_asy_num_line_types)::linewidth
   contains
     procedure::init => coop_asy_init
     procedure::write_coor => coop_asy_write_coor
     procedure::write_limits => coop_asy_write_limits
     procedure::curve => coop_asy_curve_d
     procedure::interpolate_curve => coop_asy_interpolate_curve_d
     procedure::plot => coop_asy_curve_d
     procedure::legend => coop_asy_legend_relative
     procedure::add_legend => coop_asy_add_legend
     procedure::line => coop_asy_line_d
     procedure::lines => coop_asy_lines_d
     procedure::label => coop_asy_label_relative
     procedure::dot => coop_asy_dot_d
     procedure::dots => coop_asy_dots_d
     procedure::contour => coop_asy_contour_d
     procedure::band => coop_asy_band_d
     procedure::density => coop_asy_density_d
     procedure::xrel => coop_asy_xrel
     procedure::yrel => coop_asy_yrel
     procedure::expand => coop_asy_expand
     procedure::arrow => coop_asy_arrow_d
     procedure::arrows => coop_asy_arrows_d
  end type coop_asy



  COOP_INT , parameter::coop_asy_path_max_nclosed = 4096

  interface coop_asy_arrow
     module procedure coop_asy_arrow_d, coop_asy_arrow_s
  end interface coop_asy_arrow

  interface coop_asy_arrows
     module procedure coop_asy_arrows_d, coop_asy_arrows_s
  end interface coop_asy_arrows
  
  
  interface coop_asy_band
     module procedure coop_asy_band_s, coop_asy_band_d
  end interface coop_asy_band

  interface coop_asy_error_bar
     module procedure coop_asy_error_bar_s, coop_asy_error_bar_d
  end interface coop_asy_error_bar

  interface coop_asy_interpolate_curve
     module procedure coop_asy_interpolate_curve_d, coop_asy_interpolate_curve_s
  end interface coop_asy_interpolate_curve

  interface coop_asy_gray_color
     module procedure coop_asy_gray_color_s, coop_asy_gray_color_d
  end interface coop_asy_gray_color

  interface coop_asy_rgb_color
     module procedure coop_asy_rgb_color_s, coop_asy_rgb_color_d
  end interface coop_asy_rgb_color

  interface coop_asy_label
     module procedure coop_asy_label_s, coop_asy_label_d
  end interface coop_asy_label

  interface coop_asy_legend
     module procedure coop_asy_legend_s, coop_asy_legend_d, coop_asy_legend_default, coop_asy_legend_location
  end interface coop_asy_legend

  interface coop_asy_dot
     module procedure coop_asy_dot_s, coop_asy_dot_d
  end interface coop_asy_dot

  interface coop_asy_clip
     module procedure coop_asy_clip_s, coop_asy_clip_d
  end interface coop_asy_clip


  interface coop_asy_line
     module procedure coop_asy_line_s, coop_asy_line_d
  end interface coop_asy_line

  interface coop_asy_labels
     module procedure coop_asy_labels_s, coop_asy_labels_d
  end interface coop_asy_labels

  interface coop_asy_dots
     module procedure coop_asy_dots_s, coop_asy_dots_d
  end interface coop_asy_dots

  interface coop_asy_lines
     module procedure coop_asy_lines_s, coop_asy_lines_d
  end interface coop_asy_lines

  interface coop_asy_contour
     module procedure coop_asy_contour_s, coop_asy_contour_d, coop_asy_contour_mult_s, coop_asy_contour_mult_d,  coop_asy_contour_arr_d, coop_asy_contour_arr_s, coop_asy_contour_path
  end interface coop_asy_contour

  interface coop_asy_curve
     module procedure coop_asy_curve_s, coop_asy_curve_d
  end interface coop_asy_curve

  interface coop_asy_plot_likelihood
     module procedure coop_asy_plot_likelihood_s, coop_asy_plot_likelihood_d
  end interface coop_asy_plot_likelihood


  interface coop_asy_density
     module procedure coop_asy_density_s, coop_asy_density_d, coop_asy_irregular_density_s, coop_asy_irregular_density_d
  end interface coop_asy_density

  interface coop_asy_path_append
     module procedure coop_asy_path_append_s,  coop_asy_path_append_d
  end interface coop_asy_path_append

  interface coop_asy_topaxis
     module procedure coop_asy_topaxis_s, coop_asy_topaxis_d
  end interface coop_asy_topaxis

  interface coop_asy_rightaxis
     module procedure coop_asy_rightaxis_s, coop_asy_rightaxis_d
  end interface coop_asy_rightaxis

  type coop_asy_path
     COOP_INT  nclosed
     COOP_INT ,dimension(:),allocatable:: length
     type(coop_list_realarr) l
  end type coop_asy_path



contains



  subroutine coop_asy_init(fp,  xmin, xmax, ymin, ymax, width, height, caption, xlabel, ylabel, xlog, ylog, zlog, doclip, nblocks)
    class(coop_asy) fp
    COOP_SINGLE ,optional:: width, height
    COOP_UNKNOWN_STRING,optional:: caption, xlabel, ylabel
    COOP_INT ,optional::nblocks
    logical,optional::xlog, ylog, zlog, doclip
    COOP_SINGLE ,optional:: xmin, xmax, ymin, ymax
    character(len = 5) tmp
    if(present(width))then
       fp%width = width
    else
       fp%width = coop_asy_default_width
    endif
    if(present(height))then
       fp%height = height
    else
       fp%height = coop_asy_default_height
    endif
    write(fp%unit, "(2F12.2)") fp%width, fp%height
    if(present(caption))then
       if(caption.eq."")then
          write(fp%unit, "(A)") "NULL"
       else
          write(fp%unit, "(A)") trim(caption)
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(xlabel))then
       if(trim(xlabel).eq."")then
          write(fp%unit, "(A)") "NULL"
       else
          write(fp%unit, "(A)") trim(xlabel)
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(ylabel))then
       if(trim(ylabel).eq."")then
          write(fp%unit, "(A)") "NULL"
       else
          write(fp%unit, "(A)") trim(ylabel)
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(xlog))then
       fp%xlog  =xlog
       if(xlog)then
          tmp(1:2) = "1 "
       else
          tmp(1:2) = "0 "
       endif
    else
       fp%xlog = .false.
       tmp(1:2) = "0 "
    endif
    if(present(ylog))then
       fp%ylog = ylog
       if(ylog)then
          tmp(3:4) = "1 "
       else
          tmp(3:4) = "0 "
       endif
    else
       fp%ylog = .false.
       tmp(3:4) = "0 "
    endif
    if(present(zlog))then
       if(zlog)then
          tmp(5:5) = "1"
       else
          tmp(5:5) = "0"
       endif
    else
       tmp(5:5) = "0"
    endif
    write(fp%unit, "(A)") tmp
    if(present(doclip))then
       if(doclip)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"
    endif
    if(present(xmin))then
       fp%xmin = xmin
    else
       fp%xmin = 1.1e30
    endif
    if(present(xmax))then
       fp%xmax = xmax
    else
       fp%xmax = -1.1e30
    endif
    if(present(ymin))then
       fp%ymin = ymin
    else
       fp%ymin = 1.1e30
    endif
    if(present(ymax))then
       fp%ymax = ymax
    else
       fp%ymax = -1.1e30
    endif
    write(fp%unit,"(2G14.5)") fp%xmin, fp%xmax
    write(fp%unit,"(2G14.5)") fp%ymin, fp%ymax
    if(present(nblocks))then
       write(fp%unit,"(I5)") nblocks    !!plot only n blocks
    else
       write(fp%unit,"(I5)") 0   !!0 means any number of blocks
    endif
    fp%color(1) = "black"
    fp%color(2) = "red"
    fp%color(3) = "blue"
    fp%color(4) = "green"
    fp%color(5) = "violet"
    fp%color(6) = "cyan"
    fp%color(7) = "skyblue"
    fp%color(8) = "orange"
    fp%color(9) = "brown"
    fp%color(10) = "gray"
    fp%color(11) = "pink"
    fp%color(12) = "yellow"
    fp%linewidth(1:3) = 1.2
    fp%linewidth(4:6) = 0.9
    fp%linewidth(7:12) = 0.6
    fp%linetype(1) = "solid"
    fp%linetype(2) = "dotted"
    fp%linetype(3) = "dashed"
    fp%linetype(4) = "dashdotted"
    fp%linetype(5) = "longdashed"
    fp%linetype(6) = "longdashdotted"
    fp%linetype(7:12) = "solid"
  end subroutine coop_asy_init

  subroutine coop_asy_write_limits(fp, xmin, xmax, ymin, ymax)
    class(coop_asy) fp
    COOP_SINGLE  xmin, xmax, ymin, ymax
    write(fp%unit, "(2G14.5)") xmin, xmax
    write(fp%unit, "(2G14.5)") ymin, ymax
    fp%xmin = min(xmin, fp%xmin)
    fp%xmax = min(xmax, fp%xmax)
    fp%ymin = min(ymin, fp%ymin)
    fp%ymax = min(ymax, fp%ymax)
  end subroutine coop_asy_write_limits

  subroutine coop_asy_write_coor(fp, x, y, x2, y2)
    class(coop_asy) fp
    COOP_SINGLE  x, y
    COOP_SINGLE ,optional:: x2, y2
    if(present(x2) .and. present(y2))then
       write(fp%unit, "(4G14.5)") x, y, x2, y2
       fp%xmin = min(x2, fp%xmin)
       fp%xmax = max(x2, fp%xmax)
       fp%ymin = min(y2, fp%ymin)
       fp%ymax = max(y2, fp%ymax)
    else
       write(fp%unit, "(2G14.5)") x, y
    endif
    fp%xmin = min(x, fp%xmin)
    fp%xmax = max(x, fp%xmax)
    fp%ymin = min(y, fp%ymin)
    fp%ymax = max(y, fp%ymax)
  end subroutine coop_asy_write_coor

  subroutine coop_asy_dot_d(fp, x, y, color, symbol)
    class(coop_asy) fp
    COOP_INT  n
    COOP_REAL  x,y
    COOP_UNKNOWN_STRING,optional:: color
    COOP_UNKNOWN_STRING,optional::symbol
    write(fp%unit, "(A)") "DOTS"
    write(fp%unit, "(A)") "1"
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(symbol))then
       write(fp%unit, "(A)") trim(symbol)
    else
       write(fp%unit, "(A)") "dot"
    endif
    call fp%write_coor(real(x, sp), real(y, sp))
  end subroutine coop_asy_dot_d


  subroutine coop_asy_dot_s(fp, x, y, color, symbol)
    class(coop_asy) fp
    COOP_INT  n
    COOP_SINGLE  x,y
    COOP_UNKNOWN_STRING,optional:: color
    COOP_UNKNOWN_STRING,optional::symbol
    write(fp%unit, "(A)") "DOTS"
    write(fp%unit, "(A)") "1"
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(symbol))then
       write(fp%unit, "(A)") trim(symbol)
    else
       write(fp%unit, "(A)") "dot"
    endif
    call fp%write_coor(x, y)
  end subroutine coop_asy_dot_s



  subroutine coop_asy_dots_d(fp, x, y, color, symbol)
    class(coop_asy) fp
    COOP_INT  n,i
    COOP_REAL ,dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional:: color
    COOP_UNKNOWN_STRING,optional::symbol
    n = coop_getdim( "coop_asy_dot_block", size(x), size(y))
    write(fp%unit, "(A)") "DOTS"
    write(fp%unit, "(I8)") n
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(symbol))then
       write(fp%unit, "(A)") trim(symbol)
    else
       write(fp%unit, "(A)") "dot"
    endif
    do i=1,n
       call fp%write_coor(real(x(i),sp), real(y(i),sp))
    enddo
  end subroutine coop_asy_dots_d

  subroutine coop_asy_dots_s(fp, x, y, color, symbol)
    class(coop_asy) fp
    COOP_INT  n,i
    COOP_SINGLE ,dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional:: color
    COOP_UNKNOWN_STRING,optional::symbol
    n = coop_getdim( "coop_asy_dot_block", size(x), size(y))
    write(fp%unit, "(A)") "DOTS"
    write(fp%unit, "(I8)") n
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(symbol))then
       write(fp%unit, "(A)") trim(symbol)
    else
       write(fp%unit, "(A)") "dot"
    endif
    do i=1,n
       call fp%write_coor(x(i),y(i))
    enddo
  end subroutine coop_asy_dots_s


  subroutine coop_asy_line_d(fp, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy) fp
    COOP_REAL ::xstart, ystart, xend, yend
    COOP_UNKNOWN_STRING, optional::color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING lineproperty
    write(fp%unit, "(A)") "LINES"
    write(fp%unit, "(A)") "1"
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    call fp%write_coor(real(xstart, sp), real(ystart, sp), real(xend, sp), real(yend, sp))
  end subroutine coop_asy_line_d


  subroutine coop_asy_line_s(fp, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy) fp
    real::xstart, ystart, xend, yend
    COOP_UNKNOWN_STRING, optional::color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING lineproperty
    write(fp%unit, "(A)") "LINES"
    write(fp%unit, "(A)") "1"
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    call fp%write_coor(xstart, ystart, xend, yend)
  end subroutine coop_asy_line_s
  
  subroutine coop_asy_lines_d(fp, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy) fp
    COOP_REAL ,dimension(:),intent(IN)::xstart, ystart, xend, yend
    COOP_UNKNOWN_STRING, optional::color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING lineproperty
    COOP_INT  n, i
    write(fp%unit, "(A)") "LINES"
    n = coop_getdim("coop_asy_lines", size(xstart), size(ystart), size(xend), size(yend))
    write(fp%unit, "(I8)") n
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    do i = 1, n
       call fp%write_coor(real(xstart(i),sp), real(ystart(i),sp), real(xend(i), sp), real(yend(i), sp))
    enddo
  end subroutine coop_asy_lines_d


  subroutine coop_asy_lines_s(fp, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy) fp
    COOP_SINGLE ,dimension(:),intent(IN)::xstart, ystart, xend, yend
    COOP_UNKNOWN_STRING, optional::color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING lineproperty
    COOP_INT  n, i
    write(fp%unit, "(A)") "LINES"
    n = coop_getdim("coop_asy_lines", size(xstart), size(ystart), size(xend), size(yend))
    write(fp%unit, "(I8)") n
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    do i = 1, n
       call fp%write_coor( xstart(i), ystart(i), xend(i), yend(i))
    enddo
  end subroutine coop_asy_lines_s

  subroutine coop_asy_add_legend(fp,  legend,  color, linetype, linewidth)
    class(coop_asy)::fp
    COOP_UNKNOWN_STRING::legend
    COOP_UNKNOWN_STRING, optional :: color, linetype
    COOP_SINGLE, optional::linewidth
    COOP_STRING lineproperty
    write(fp%unit, "(A)") "LEGEND"
    write(fp%unit, "(A)") "VIRTUAL"
    if(trim(adjustl(legend)) .ne. "")then
       write(fp%unit, "(A)") trim(adjustl(legend))
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    else
       lineproperty = trim(lineproperty)//"_8"
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    fp%ymax = fp%ymax + (fp%ymax-fp%ymin)*0.01
  end subroutine coop_asy_add_legend

  subroutine coop_asy_interpolate_curve_d(fp, xraw, yraw, interpolate, color, linetype, linewidth, legend)
    COOP_INT ,parameter::n = 256
    COOP_REAL  x(n), y(n),  minx, maxx, dx
    COOP_INT  w(n)
    COOP_REAL ,dimension(:),allocatable::xx, yy, yy2
    class(coop_asy) fp
    COOP_REAL ,dimension(:),intent(IN)::xraw, yraw
    COOP_UNKNOWN_STRING interpolate
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_UNKNOWN_STRING, optional::legend
    COOP_INT  i, j, m, nn, npt
    COOP_STRING lineproperty
    logical do_raw
    m = Coop_getdim("coop_asy_interpolate_curve", size(xraw), size(yraw))
    minx = minval(xraw)
    maxx = maxval(xraw)
    y = 0.
    w = 0
    if(minx .ge. maxx)then
       npt = 1
       x(1) = minx
       select case(trim(interpolate))
       case("LinearLog", "LogLog")
          y(1) = exp(sum(log(yraw))/m)
       case default
          y(1) = sum(yraw)/m
       end select
       goto 100
    endif
    npt = n
    do_raw = .false.
    select case(trim(interpolate))
    case("LinearLinear")
       call coop_set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((xraw(i) - minx)/dx)+1
          y(j) = y(j) + yraw(i)
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
    case("LinearLog")
       call coop_set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((xraw(i) - minx)/dx)+1
          y(j) = y(j) + log(yraw(i))
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       y = exp(y)
    case("LogLinear")
       minx = log(minx)
       maxx = log(maxx)
       call coop_set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((log(xraw(i)) - minx)/dx)+1
          y(j) = y(j) + yraw(i)
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       x = exp(x)
    case("LogLog")
       minx = log(minx)
       maxx = log(maxx)
       call coop_set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((log(xraw(i)) - minx)/dx)+1
          y(j) = y(j) + log(yraw(i))
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       x = exp(x)
       y = exp(y)
    case default
       do_raw = .true.
       npt = m
    end select

100 write(fp%unit, "(A)") "CURVE"
    write(fp%unit, "(I8)") npt
    if(present(legend))then
       if(trim(legend).ne."")then
          write(fp%unit, "(A)") trim(legend)
       else
          write(fp%unit, "(A)") "NULL"
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    write(fp%unit, "(A)") "0"
    if(do_raw)then
       do i=1,npt
          call fp%write_coor( real(xraw(i), sp), real(yraw(i), sp) )
       enddo
    else
       do i=1,npt
          call fp%write_coor(real(x(i), sp), real(y(i), sp))
       enddo
    endif
  end subroutine coop_asy_interpolate_curve_d

  subroutine coop_asy_plot_likelihood_d(fp, xraw, yraw, color, linetype, linewidth, legend, left_tail, right_tail)
    class(coop_asy) fp
    COOP_REAL ,dimension(:),intent(IN)::xraw, yraw
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_UNKNOWN_STRING, optional::legend
    COOP_INT  i, j, m, nn, npt, imax
    COOP_STRING lineproperty
    logical, optional::left_tail, right_tail
    logical do_left_tail, do_right_tail
    COOP_REAL  ymax
    m = Coop_getdim("coop_asy_fit_likliehood", size(xraw), size(yraw))
    imax = coop_maxloc(COOP_REAL_OF(yraw))
    ymax = yraw(imax)
    npt = m 
    if(imax.gt.1)  npt = npt + 1
    if(imax .lt. m) npt = npt + 1
    do_left_tail = .false.
    do_right_tail = .false.
    if(present(left_tail))then
       if(left_tail .and. yraw(1) .gt. ymax/100.d0 .and. yraw(1).lt. ymax/10.d0)then
          do_left_tail = .true.
          npt = npt + 1
       endif
    endif
    if(present(right_tail))then
       if(right_tail .and. yraw(m) .gt. ymax/100.d0 .and. yraw(m).lt. ymax/10.d0)then
          do_right_tail = .true.
          npt = npt+1
       endif
    endif
100 write(fp%unit, "(A)") "CURVE"
    write(fp%unit, "(I8)") npt
    if(present(legend))then
       if(trim(legend).ne."")then
          write(fp%unit, "(A)") trim(legend)
       else
          write(fp%unit, "(A)") "NULL"
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    write(fp%unit, "(A)") "0"
    if(do_left_tail)then
       call fp%write_coor(real(xraw(1)*2.-xraw(2), sp), 0._sp)
    endif
    do i=1,imax-1
       call fp%write_coor( real(xraw(i), sp), real(yraw(i), sp) )
    enddo
    if(imax.gt.1) &
         call fp%write_coor( real((xraw(i)+xraw(i-1))/2., sp), real(yraw(i)*0.75+yraw(i-1)*0.25, sp) )
    call fp%write_coor( real(xraw(imax), sp), real(yraw(imax), sp) )
    if(imax .lt. m) &
         call fp%write_coor( real((xraw(i)+xraw(i+1))/2., sp), real(yraw(i)*0.75+yraw(i+1)*0.25, sp) )
    do i=imax+1, m
       call fp%write_coor( real(xraw(i), sp), real(yraw(i), sp) )
    enddo
    if(do_right_tail)then
       call fp%write_coor(real(xraw(m)*2.-xraw(m-1), sp), 0._sp)
    endif
  end subroutine coop_asy_plot_likelihood_d

  subroutine coop_asy_interpolate_curve_s(fp, xraw, yraw, interpolate, color, linetype, linewidth, legend)
    COOP_INT ,parameter::n = 256
    COOP_REAL  x(n), y(n),  minx, maxx, dx
    COOP_INT  w(n)
    COOP_REAL ,dimension(:),allocatable::xx, yy, yy2
    class(coop_asy) fp
    COOP_SINGLE ,dimension(:),intent(IN)::xraw, yraw
    COOP_UNKNOWN_STRING interpolate
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_UNKNOWN_STRING, optional::legend
    COOP_INT  i, j, m, nn, npt
    COOP_STRING lineproperty
    logical do_raw
    m = Coop_getdim("coop_asy_interpolate_curve", size(xraw), size(yraw))
    minx = minval(xraw)
    maxx = maxval(xraw)
    y = 0.
    w = 0
    if(minx .ge. maxx)then
       npt = 1
       x(1) = minx
       select case(trim(interpolate))
       case("LinearLog", "LogLog")
          y(1) = exp(sum(log(yraw))/m)
       case default
          y(1) = sum(yraw)/m
       end select
       goto 100
    endif
    npt = n
    do_raw = .false.
    select case(trim(interpolate))
    case("LinearLinear")
       call coop_set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((xraw(i) - minx)/dx)+1
          y(j) = y(j) + yraw(i)
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
    case("LinearLog")
       call coop_set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((xraw(i) - minx)/dx)+1
          y(j) = y(j) + log(yraw(i))
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m,xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       y = exp(y)
    case("LogLinear")
       minx = log(minx)
       maxx = log(maxx)
       call coop_set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((log(xraw(i)) - minx)/dx)+1
          y(j) = y(j) + yraw(i)
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       x = exp(x)
    case("LogLog")
       minx = log(minx)
       maxx = log(maxx)
       call coop_set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((log(xraw(i)) - minx)/dx)+1
          y(j) = y(j) + log(yraw(i))
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       x = exp(x)
       y = exp(y)
    case default
       do_raw = .true.
       npt = m
    end select

100 write(fp%unit, "(A)") "CURVE"
    write(fp%unit, "(I8)") npt
    if(present(legend))then
       if(trim(legend).ne."")then
          write(fp%unit, "(A)") trim(legend)
       else
          write(fp%unit, "(A)") "NULL"
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    write(fp%unit, "(A)") "0"
    if(do_raw)then
       do i=1,npt
          call fp%write_coor( xraw(i), yraw(i) )
       enddo
    else
       do i=1,npt
          call fp%write_coor( real(x(i), sp), real(y(i), sp) )
       enddo
    endif
  end subroutine coop_asy_interpolate_curve_s



  subroutine coop_asy_plot_likelihood_s(fp, xraw, yraw, color, linetype, linewidth, legend, left_tail, right_tail)
    class(coop_asy) fp
    COOP_SINGLE ,dimension(:),intent(IN)::xraw, yraw
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_UNKNOWN_STRING, optional::legend
    COOP_INT  i, j, m, nn, npt, imax
    COOP_STRING lineproperty
    logical  do_left_tail, do_right_tail
    logical, optional::left_tail, right_tail
    real ymax
    m = Coop_getdim("coop_asy_plot_likelihood", size(xraw), size(yraw))
    imax = coop_maxloc(COOP_REAL_OF(yraw))
    ymax = yraw(imax)
    npt = m
    if(imax .gt. 1) npt = npt + 1
    if(imax .lt. m) npt = npt + 1
    do_left_tail = .false.
    do_right_tail = .false.
    if(present(left_tail))then
       if(left_tail .and. yraw(1) .gt. ymax/100.d0 .and. yraw(1).lt. ymax/10.d0)then
          do_left_tail = .true.
          npt = npt + 1
       endif
    endif
    if(present(right_tail))then
       if(right_tail .and. yraw(m) .gt. ymax/100.d0 .and. yraw(m).lt. ymax/10.d0)then
          do_right_tail = .true.
          npt = npt+1
       endif
    endif
100 write(fp%unit, "(A)") "CURVE"
    write(fp%unit, "(I8)") npt
    if(present(legend))then
       if(trim(legend).ne."")then
          write(fp%unit, "(A)") trim(legend)
       else
          write(fp%unit, "(A)") "NULL"
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    write(fp%unit, "(A)") "0"
    if(do_left_tail) call fp%write_coor(xraw(1)*2.-xraw(2), 0.)
    do i=1,imax-1
       call fp%write_coor( xraw(i), yraw(i) )
    enddo
    if(imax.gt.1) &
         call fp%write_coor( (xraw(i)+xraw(i-1))/2.,  yraw(i)*0.75+yraw(i-1)*0.25 )
    call fp%write_coor(xraw(imax), yraw(imax))
    if(imax .lt. m) &
         call fp%write_coor( (xraw(i)+xraw(i+1))/2., yraw(i)*0.75+yraw(i+1)*0.25 )
    do i=imax+1, m
       call fp%write_coor( xraw(i), yraw(i) )
    enddo
    if(do_right_tail) call fp%write_coor(xraw(m)*2.-xraw(m-1), 0.)
  end subroutine coop_asy_plot_likelihood_s



  subroutine coop_asy_curve_from_file(fp, fname, interpolate, xcol, ycol, color, linetype, linewidth, legend)
    COOP_INT ,parameter::n = 256
    COOP_INT , optional::xcol, ycol
    COOP_REAL  x(n), y(n),  minx, maxx, dx
    COOP_INT  w(n)
    COOP_REAL ,dimension(:),allocatable::xx, yy, yy2
    class(coop_asy) fp
    type(coop_file) fp2
    COOP_REAL ,dimension(:),allocatable::xraw, yraw, line
    COOP_UNKNOWN_STRING interpolate, fname
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_UNKNOWN_STRING, optional::legend
    COOP_INT  i, j, m, nn, npt, ncols, xl, yl
    COOP_STRING lineproperty
    logical do_raw
    if(.not. coop_file_exists(fname))then
       call coop_return_error("coop_asy_curve_from_file", "The data file "//trim(fname)//" does not exist", "stop")
    endif
    m = coop_file_numlines(fname)
    ncols = coop_file_numColumns(fname)
    if(ncols .eq. 0)then
       call coop_return_error("coop_asy_curve_from_file", "The data file "//trim(fname)//" is empty", "stop")
    endif
    if(present(xcol))then
       if(xcol .gt. ncols) stop "coop_asy_curve_from_file: xcol > ncol"
       xl = xcol
    else
       x = min(ncols-1, 1)
    endif
    if(present(ycol))then
       if(ycol .gt. ncols) stop "coop_asy_curve_from_file: ycol > ncol"
       yl = ycol
    else
       yl = min(ncols, 2)
    endif
    allocate(line(0:ncols))
    allocate(xraw(m), yraw(m))
    call fp2%open(trim(fname), "r")
    do i = 1, m
       if(fp2%read_real_array(line(1:ncols)))then
          line(0) = i
          xraw(i) = line(xl)
          yraw(i) = line(yl)
       else
          call coop_return_error("coop_asy_curve_from_file", trim(fname)//" has a wrong format", "stop")
       endif
    enddo
    call fp2%close()
    minx = minval(xraw)
    maxx = maxval(xraw)
    y = 0.
    w = 0
    if(minx .ge. maxx)then
       npt = 1
       x(1) = minx
       select case(trim(interpolate))
       case("LinearLog", "LogLog")
          y(1) = exp(sum(log(yraw))/m)
       case default
          y(1) = sum(yraw)/m
       end select
       goto 100
    endif
    npt = n
    do_raw = .false.
    select case(trim(interpolate))
    case("LinearLinear")
       call coop_set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((xraw(i) - minx)/dx)+1
          y(j) = y(j) + yraw(i)
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
    case("LinearLog")
       call coop_set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((xraw(i) - minx)/dx)+1
          y(j) = y(j) + log(yraw(i))
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       y = exp(y)
    case("LogLinear")
       minx = log(minx)
       maxx = log(maxx)
       call coop_set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((log(xraw(i)) - minx)/dx)+1
          y(j) = y(j) + yraw(i)
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       x = exp(x)
    case("LogLog")
       minx = log(minx)
       maxx = log(maxx)
       call coop_set_uniform(n, x, minx, maxx)
       dx = (maxx - minx)/(n-1.)
       do i=1, m
          j = nint((log(xraw(i)) - minx)/dx)+1
          y(j) = y(j) + log(yraw(i))
          w(j) = w(j) + 1
       enddo
       where(w.ne.0)
          y = y/w
       end where
       if(any(w.eq.0))then
          m = count(w.ne.0)
          allocate(xx(m), yy(m), yy2(m))
          j = 1
          do i=1, n
             if(w(i).ne.0)then
                xx(j) = x(i)
                yy(j) = y(i)
                j = j + 1
             endif
          enddo
          call coop_spline(m, xx, yy, yy2)
          do i=1, n
             if(w(i).eq.0)then
                call coop_splint(m, xx, yy, yy2, x(i), y(i))
             endif
          enddo
          deallocate(xx, yy, yy2)
       endif
       x = exp(x)
       y = exp(y)
    case default
       do_raw = .true.
       npt = m
    end select

100 write(fp%unit, "(A)") "CURVE"
    write(fp%unit, "(I8)") npt
    if(present(legend))then
       if(trim(legend).ne."")then
          write(fp%unit, "(A)") trim(legend)
       else
          write(fp%unit, "(A)") "NULL"
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    write(fp%unit, "(A)") "0"
    if(do_raw)then
       do i=1,npt
          call fp%write_coor( real(xraw(i), sp), real(yraw(i), sp))
       enddo
    else
       do i=1,npt
          call fp%write_coor(real(x(i), sp), real(y(i), sp))
       enddo
    endif
    deallocate(xraw, yraw, line)
  end subroutine coop_asy_curve_from_file


  subroutine coop_asy_curve_d(fp, x, y, smooth, color, linetype, linewidth, legend)
    class(coop_asy) fp
    logical,optional::smooth
    COOP_REAL ,dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_UNKNOWN_STRING, optional::legend
    COOP_INT  i,n
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_curve", size(x), size(y))
    write(fp%unit, "(A)") "CURVE"
    write(fp%unit, "(I8)") n
    if(present(legend))then
       if(trim(legend).ne."")then
          write(fp%unit, "(A)") trim(legend)
       else
          write(fp%unit, "(A)") "NULL"
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    do i=1,n
       call fp%write_coor(real(x(i), sp), real(y(i), sp))
    enddo
  end subroutine coop_asy_curve_d

  subroutine coop_asy_curve_s(fp, x, y, smooth, color, linetype, linewidth, legend)
    class(coop_asy) fp
    logical,optional::smooth
    COOP_SINGLE ,dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_INT  i,n
    COOP_UNKNOWN_STRING,optional::legend
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_curve", size(x), size(y))
    write(fp%unit, "(A)") "CURVE"
    write(fp%unit, "(I8)") n
    if(present(legend))then
       if(trim(legend).ne."")then
          write(fp%unit, "(A)") trim(legend)
       else
          write(fp%unit, "(A)") "NULL"
       endif
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    do i=1,n
       call fp%write_coor( x(i), y(i) )
    enddo
  end subroutine coop_asy_curve_s

  subroutine coop_asy_labels_d(fp, labels, x, y, color, alignment)
    COOP_STRING, dimension(:),intent(IN)::labels
    COOP_REAL ,dimension(:),intent(IN)::x, y
    class(coop_asy) fp
    COOP_UNKNOWN_STRING,optional::color
    COOP_UNKNOWN_STRING,optional::alignment
    COOP_INT  n, i
    if(present(alignment))then
       select case(trim(adjustl(alignment)))
       case("left", "LEFT", "Left", "l", "L")
          write(fp%unit, "(A)") "LEFTLABELS"
       case("right", "RIGHT", "Right", "r", "R")
          write(fp%unit, "(A)") "RIGHTLABELS"
       case default
          write(fp%unit, "(A)") "LABELS"          
       end select
    else
       write(fp%unit, "(A)") "LABELS"
    endif    
    n = coop_getdim("coop_asy_labels", size(x), size(y), size(labels))
    write(fp%unit,"(I8)") n
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    do i=1,n
       call fp%write_coor(real(x(i), sp), real(y(i), sp))
       if(trim(labels(i)).eq."")then
          write(fp%unit, "(A)") "NULLL"
       else
          write(fp%unit, "(A)") trim(labels(i))
       endif
    enddo
  end subroutine coop_asy_labels_d

  subroutine coop_asy_labels_s(fp, labels, x, y, color, alignment)
    COOP_STRING, dimension(:),intent(IN)::labels
    COOP_SINGLE ,dimension(:),intent(IN)::x, y
    class(coop_asy) fp
    COOP_UNKNOWN_STRING,optional::color
    COOP_UNKNOWN_STRING,optional::alignment
    COOP_INT  n, i
    if(present(alignment))then
       select case(trim(alignment))
       case("left","LEFT","Left","l","L")
          write(fp%unit, "(A)") "LEFTLABELS"
       case("right","RIGHT","Right","r","R")
          write(fp%unit, "(A)") "RIGHTLABELS"
       case default
          write(fp%unit, "(A)") "LABELS"          
       end select
    else
       write(fp%unit, "(A)") "LABELS"
    endif
    n = coop_getdim("coop_asy_labels", size(x), size(y), size(labels))
    write(fp%unit,"(I8)") n
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    do i=1,n
       call fp%write_coor( x(i), y(i))
       if(trim(labels(i)).eq."")then
          write(fp%unit, "(A)") "NULLL"
       else
          write(fp%unit, "(A)") trim(labels(i))
       endif
    enddo
  end subroutine coop_asy_labels_s


  subroutine coop_asy_label_d(fp, label, x, y, color, alignment)
    class(coop_asy) fp
    COOP_UNKNOWN_STRING label
    COOP_UNKNOWN_STRING,optional::color
    COOP_UNKNOWN_STRING,optional::alignment
    COOP_REAL  x, y
    if(present(alignment))then
       select case(trim(alignment))
       case("left","LEFT","Left","l","L")
          write(fp%unit, "(A)") "LEFTLABELS"
       case("right","RIGHT","Right","r","R")
          write(fp%unit, "(A)") "RIGHTLABELS"
       case default
          write(fp%unit, "(A)") "LABELS"          
       end select
    else
       write(fp%unit, "(A)") "LABELS"
    endif
    write(fp%unit, "(A)") "1"
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    call fp%write_coor( real(x, sp), real(y, sp) )
    if(trim(label).eq."")then
       write(fp%unit, "(A)") "NULL"
    else
       write(fp%unit, "(A)") trim(label)
    endif
  end subroutine coop_asy_label_d

  subroutine coop_asy_label_relative(fp, label, xratio, yratio, color, alignment)
    class(coop_asy) fp
    COOP_UNKNOWN_STRING label
    COOP_UNKNOWN_STRING,optional::color
    COOP_UNKNOWN_STRING,optional::alignment
    COOP_SINGLE  xratio, yratio
    if(present(alignment))then
       select case(trim(alignment))
       case("left","LEFT","Left","l","L")
          write(fp%unit, "(A)") "LEFTLABELS"
       case("right","RIGHT","Right","r","R")
          write(fp%unit, "(A)") "RIGHTLABELS"
       case default
          write(fp%unit, "(A)") "LABELS"          
       end select
    else
       write(fp%unit, "(A)") "LABELS"
    endif
    write(fp%unit, "(A)") "1"
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    call fp%write_coor(fp%xrel(xratio), fp%yrel(yratio))
    if(trim(label).eq."")then
       write(fp%unit, "(A)") "NULL"
    else
       write(fp%unit, "(A)") trim(label)
    endif
  end subroutine coop_asy_label_relative

  
  subroutine coop_asy_label_s(fp, label, x, y, color, alignment)
    class(coop_asy) fp
    COOP_UNKNOWN_STRING label
    COOP_UNKNOWN_STRING,optional::color
    COOP_UNKNOWN_STRING,optional::alignment
    COOP_SINGLE  x, y
    if(present(alignment))then
       select case(trim(alignment))
       case("left","LEFT","Left","l","L")
          write(fp%unit, "(A)") "LEFTLABELS"
       case("right","RIGHT","Right","r","R")
          write(fp%unit, "(A)") "RIGHTLABELS"
       case default
          write(fp%unit, "(A)") "LABELS"          
       end select
    else
       write(fp%unit, "(A)") "LABELS"
    endif    
    write(fp%unit, "(A)") "1"
    if(present(color))then
       write(fp%unit, "(A)") trim(color)
    else
       write(fp%unit, "(A)") "black"
    endif
    call fp%write_coor( x, y)
    if(trim(label).eq."")then
       write(fp%unit, "(A)") "NULL"
    else
       write(fp%unit, "(A)") trim(label)
    endif
  end subroutine coop_asy_label_s



  subroutine coop_asy_contour_d(fp, x, y, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) fp
    logical,optional::smooth
    COOP_REAL ,dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_INT  i,n
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_contour", size(x), size(y))
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "1" !!type 1 contour
    if(present(colorfill))then
       write(fp%unit, "(A)") trim(colorfill)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(linecolor).or. present(linetype).or.present(linewidth))then
       if(present(linecolor))then
          lineproperty=trim(linecolor)
       else
          lineproperty = "black"
       endif
       if(present(linetype))then
          lineproperty = trim(lineproperty)//"_"//trim(linetype)
       else
          lineproperty = trim(lineproperty)//"_solid"
       endif
       if(present(linewidth))then
          lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
       endif
       write(fp%unit, "(A)") trim(lineproperty)
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(fp%unit, "(A)") "1"  !!1 path 
    write(fp%unit, "(I8)") n
    do i=1,n
       call fp%write_coor(real(x(i), sp), real(y(i), sp))
    enddo
  end subroutine coop_asy_contour_d

  subroutine coop_asy_contour_s(fp, x, y, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) fp
    logical,optional::smooth
    COOP_SINGLE ,dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_INT  i,n
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_contour", size(x), size(y))
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "1"   !!type I contour
    if(present(colorfill))then
       write(fp%unit, "(A)") trim(colorfill)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(linecolor).or. present(linetype).or.present(linewidth))then
       if(present(linecolor))then
          lineproperty=trim(linecolor)
       else
          lineproperty = "black"
       endif
       if(present(linetype))then
          lineproperty = trim(lineproperty)//"_"//trim(linetype)
       else
          lineproperty = trim(lineproperty)//"_solid"
       endif
       if(present(linewidth))then
          lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
       endif
       write(fp%unit, "(A)") trim(lineproperty)
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(fp%unit, "(A)") "1"  !!1 path 
    write(fp%unit, "(I8)") n
    do i=1,n
       call fp%write_coor( x(i), y(i) )
    enddo
  end subroutine coop_asy_contour_s


  subroutine coop_asy_contour_path(fp, path, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) fp
    logical,optional::smooth
    type(coop_asy_path) path
    COOP_SINGLE  xy(2)
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_INT  i, ipath, pl
    COOP_STRING lineproperty
    if(path%nclosed.eq.0 .or. (.not. path%l%isinit()) .or. path%l%n .eq. 0)return
    if(path%l%dim .ne.2) stop "coop_asy_contour_path: wrong dimension of list in path"
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "1"   !!type I contour
    if(present(colorfill))then
       write(fp%unit, "(A)") trim(colorfill)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(linecolor).or. present(linetype).or.present(linewidth))then
       if(present(linecolor))then
          lineproperty=trim(linecolor)
       else
          lineproperty = "black"
       endif
       if(present(linetype))then
          lineproperty = trim(lineproperty)//"_"//trim(linetype)
       else
          lineproperty = trim(lineproperty)//"_solid"
       endif
       if(present(linewidth))then
          lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
       endif
       write(fp%unit, "(A)") trim(lineproperty)
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(fp%unit, "(I8)") path%nclosed
    i = 1
    do ipath = 1, path%nclosed
       write(fp%unit, "(I8)") path%length(ipath)
       do pl = 1, path%length(ipath)
          call path%l%get_element(i, xy)
          call fp%write_coor(xy(1), xy(2))
          i = i + 1
       enddo
    enddo
  end subroutine coop_asy_contour_path

  subroutine coop_asy_contour_mult_d(fp, x, y, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) fp
    logical,optional::smooth
    COOP_REAL ,dimension(:,:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_INT  i,n, m,ipath
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_contour", size(x,1), size(y,1))
    m = coop_getdim("coop_asy_contour", size(x,2), size(y,2))
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "1"   !!type I contour
    if(present(colorfill))then
       write(fp%unit, "(A)") trim(colorfill)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(linecolor).or. present(linetype).or.present(linewidth))then
       if(present(linecolor))then
          lineproperty=trim(linecolor)
       else
          lineproperty = "black"
       endif
       if(present(linetype))then
          lineproperty = trim(lineproperty)//"_"//trim(linetype)
       else
          lineproperty = trim(lineproperty)//"_solid"
       endif
       if(present(linewidth))then
          lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
       endif
       write(fp%unit, "(A)") trim(lineproperty)
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(fp%unit, "(I8)") m
    do ipath = 1, m
       write(fp%unit, "(I8)") n
       do i=1,n
          call fp%write_coor(real(x(i, ipath), sp), real(y(i, ipath), sp))
       enddo
    enddo
  end subroutine coop_asy_contour_mult_d

  subroutine coop_asy_contour_mult_s(fp, x, y, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) fp
    logical,optional::smooth
    COOP_SINGLE ,dimension(:,:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_INT  i,n, m,ipath
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_contour", size(x,1), size(y,1))
    m = coop_getdim("coop_asy_contour", size(x,2), size(y,2))
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "1"   !!type I contour
    if(present(colorfill))then
       write(fp%unit, "(A)") trim(colorfill)
    else
       write(fp%unit, "(A)") "black"
    endif
    if(present(linecolor).or. present(linetype).or.present(linewidth))then
       if(present(linecolor))then
          lineproperty=trim(linecolor)
       else
          lineproperty = "black"
       endif
       if(present(linetype))then
          lineproperty = trim(lineproperty)//"_"//trim(linetype)
       else
          lineproperty = trim(lineproperty)//"_solid"
       endif
       if(present(linewidth))then
          lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
       endif
       write(fp%unit, "(A)") trim(lineproperty)
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(fp%unit, "(I8)") m
    do ipath = 1, m
       write(fp%unit, "(I8)") n
       do i=1,n
          call fp%write_coor( x(i, ipath), y(i, ipath) )
       enddo
    enddo
  end subroutine coop_asy_contour_mult_s


  subroutine coop_asy_contour_arr_d(fp, z, xmin, xmax, ymin, ymax, cvals, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) fp
    COOP_REAL ,dimension(:,:)::z
    COOP_REAL  xmin, xmax, ymin, ymax
    COOP_REAL  cvals(:)
    COOP_INT  nx, ny, i, nc
    COOP_STRING,dimension(:)::colorfill
    COOP_STRING,dimension(:)::linecolor, linetype
    COOP_SINGLE ,dimension(:),optional::linewidth
    logical smooth
    nx = size(z, 1)
    ny = size(z, 2)
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "2" !!type 2 contour
    call fp%write_limits( real(xmin, sp), real(xmax, sp), real(ymin, sp), real(ymax, sp))
    if(present(linewidth))then
       nc = coop_getdim("coop_asy_contour_arr_d", size(cvals), size(colorfill), size(linecolor), size(linetype), size(linewidth))
    else
       nc = coop_getdim("coop_asy_contour_arr_d", size(cvals), size(colorfill), size(linecolor), size(linetype))
    endif
    write(fp%unit, "(I8)") nc
    write(fp%unit, "("//trim(coop_num2str(nc))//"G14.5)") cvals
    do i = 1, nc
       write(fp%unit, "(A)") trim(colorfill(i))
       if(present(linewidth))then
          write(fp%unit, "(A)") trim(coop_asy_linestyle(linecolor(i), linetype(i), linewidth(i)))
       else
          write(fp%unit, "(A)") trim(coop_asy_linestyle(linecolor(i), linetype(i)))
       endif
    enddo
    if(smooth)then
       write(fp%unit, "(A)") "1"
    else
       write(fp%unit, "(A)") "0"
    endif
    write(fp%unit, "(2I8)") nx, ny
    do i=1,nx
       write(fp%unit, "("//trim(coop_num2str(ny))//"G14.5)") real(z(i,:))
    enddo
  end subroutine coop_asy_contour_arr_d
  


  subroutine coop_asy_contour_arr_s(fp, z, xmin, xmax, ymin, ymax, cvals, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) fp
    COOP_SINGLE ,dimension(:,:)::z
    COOP_SINGLE  xmin, xmax, ymin, ymax
    COOP_SINGLE  cvals(:)
    COOP_INT  nx, ny, i, nc
    COOP_STRING,dimension(:)::colorfill
    COOP_STRING,dimension(:)::linecolor, linetype
    COOP_SINGLE ,dimension(:),optional::linewidth
    logical smooth
    nx = size(z, 1)
    ny = size(z, 2)
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "2" !!type 2 contour
    call fp%write_limits(xmin, xmax, ymin, ymax)
    if(present(linewidth))then
       nc = coop_getdim("coop_asy_contour_arr_s", size(cvals), size(colorfill), size(linecolor), size(linetype), size(linewidth))
    else
       nc = coop_getdim("coop_asy_contour_arr_s", size(cvals), size(colorfill), size(linecolor), size(linetype))
    endif
    write(fp%unit, "(I8)") nc
    write(fp%unit, "("//trim(coop_num2str(nc))//"G14.5)") cvals
    do i = 1, nc
       write(fp%unit, "(A)") trim(colorfill(i))
       if(present(linewidth))then
          write(fp%unit, "(A)") trim(coop_asy_linestyle(linecolor(i), linetype(i), linewidth(i)))
       else
          write(fp%unit, "(A)") trim(coop_asy_linestyle(linecolor(i), linetype(i)))
       endif
    enddo
    if(smooth)then
       write(fp%unit, "(A)") "1"
    else
       write(fp%unit, "(A)") "0"
    endif
    write(fp%unit, "(2I8)") nx, ny
    do i=1,nx
       write(fp%unit, "("//trim(coop_num2str(ny))//"G14.5)") z(i,:)
    enddo
  end subroutine coop_asy_contour_arr_s




  subroutine coop_asy_density_d(fp, z, xmin, xmax, ymin, ymax, label, zmin, zmax, color_table)
    class(coop_asy) fp
    COOP_REAL ,dimension(:,:)::z
    COOP_REAL  xmin, xmax, ymin, ymax
    COOP_REAL ,optional::zmin, zmax
    COOP_UNKNOWN_STRING, optional::label
    COOP_REAL  minz, maxz
    COOP_INT  nx, ny, i
    COOP_UNKNOWN_STRING, optional::color_table
    nx = size(z, 1)
    ny = size(z, 2)
    write(fp%unit, "(A)") "DENSITY"
    if(present(color_table))then
       if(trim(color_table).ne."")then
          write(fp%unit, "(A)") trim(color_table)
       else
          write(fp%unit, "(A)") "Rainbow"
       endif
    else
       write(fp%unit, "(A)") "Rainbow"
    endif

    if(present(label))then
       if(trim(label).ne."")then
          write(fp%unit,"(A)") trim(label)
       else
          write(fp%unit,"(A)") "NULL"
       endif
    else
       write(fp%unit,"(A)") "NULL"
    endif
    call fp%write_limits(real(xmin, sp), real(xmax, sp), real(ymin, sp), real(ymax, sp))
    if(present(zmin))then
       minz = zmin
    else
       call coop_array_get_threshold(z, 0.99d0, minz)
    endif
    if(present(zmax))then
       maxz = zmax
    else
       call coop_array_get_threshold(z, 0.01d0, maxz)
    endif

    write(fp%unit, "(2G14.5)") minz, maxz
    write(fp%unit, "(A)") "0"
    write(fp%unit, "(2I8)") nx, ny
    do i=1,nx
       write(fp%unit, "("//trim(coop_num2str(ny))//"G14.5)") real(z(i,:))
    enddo
  end subroutine coop_asy_density_d

  subroutine coop_asy_density_s(fp, z, xmin, xmax, ymin, ymax, label, zmin, zmax, color_table)
    class(coop_asy) fp
    COOP_SINGLE ,dimension(:,:)::z
    COOP_SINGLE  xmin, xmax, ymin, ymax
    COOP_SINGLE ,optional::zmin, zmax
    COOP_SINGLE  minz, maxz
    COOP_INT  nx, ny, i
    COOP_UNKNOWN_STRING, optional::label, color_table
    nx = size(z, 1)
    ny = size(z, 2)
    write(fp%unit, "(A)") "DENSITY"
    if(present(color_table))then
       if(trim(color_table).ne."")then
          write(fp%unit,"(A)") trim(color_table)
       else
          write(fp%unit,"(A)") "Rainbow"
       endif
    else
       write(fp%unit,"(A)") "Rainbow"
    endif
    if(present(label))then
       if(trim(label).ne."")then
          write(fp%unit,"(A)") trim(label)
       else
          write(fp%unit,"(A)") "NULL"
       endif
    else
       write(fp%unit,"(A)") "NULL"
    endif
    call fp%write_limits(xmin, xmax, ymin, ymax)
    if(present(zmin))then
       minz = zmin
    else
       call coop_array_get_threshold(z, 0.99, minz)
    endif
    if(present(zmax))then
       maxz = zmax
    else
       call coop_array_get_threshold(z, 0.01, maxz)
    endif
    write(fp%unit, "(2G14.5)") minz, maxz
    write(fp%unit, "(A)") "0"
    write(fp%unit, "(2I8)") nx, ny
    do i=1,nx
       write(fp%unit, "("//trim(coop_num2str(ny))//"G14.5)") z(i,:)
    enddo
  end subroutine coop_asy_density_s


 subroutine coop_asy_irregular_density_d(fp, x, y, z, label, xmin, xmax, ymin, ymax, zmin, zmax, color_table)
    class(coop_asy) fp
    COOP_REAL ,dimension(:)::x, y, z
    COOP_UNKNOWN_STRING, optional::label, color_table
    COOP_REAL ,optional::xmin, xmax, ymin, ymax, zmin, zmax
    COOP_INT  i, n
    COOP_REAL  minx, miny, minz, maxx, maxy, maxz
    n = coop_getdim("coop_asy_irregular_density", size(x), size(y), size(z))
    if(present(xmin))then
       minx = xmin
    else
       call coop_array_get_threshold(x, 0.99d0, minx)
    endif
    if(present(xmax))then
       maxx = xmax
    else
       call coop_array_get_threshold(x, 0.01d0, maxx)
    endif
    if(present(ymin))then
       miny = ymin
    else
       call coop_array_get_threshold(y, 0.99d0, miny)
    endif
    if(present(ymax))then
       maxy = ymax
    else
       call coop_array_get_threshold(y, 0.01d0, maxy)
    endif
    if(present(zmin))then
       minz = zmin
    else
       call coop_array_get_threshold(z, 0.99d0, minz)
    endif
    if(present(zmax))then
       maxz = zmax
    else
       call coop_array_get_threshold(z, 0.01d0, maxz)
    endif
    write(fp%unit,"(A)") "DENSITY"
    if(present(color_table))then
       if(trim(color_table).ne."")then
          write(fp%unit,"(A)") trim(color_table)
       else
          write(fp%unit,"(A)") "Rainbow"
       endif
    else
       write(fp%unit,"(A)") "Rainbow"
    endif
    if(present(label))then
       if(trim(label).ne."")then
          write(fp%unit,"(A)") trim(label)
       else
          write(fp%unit,"(A)") "NULL"
       endif
    else
       write(fp%unit,"(A)") "NULL"
    endif
    call fp%write_limits(real(minx, sp), real(maxx, sp), real(miny, sp), real(maxy, sp))
    write(fp%unit, "(2G14.5)") minz, maxz
    write(fp%unit,"(A)") "1"
    write(fp%unit,"(I10)") n
    do i=1,n
       write(fp%unit, "(3G14.5)") real(x(i)), real(y(i)), real(z(i))
    enddo
  end subroutine coop_asy_irregular_density_d

  subroutine coop_asy_irregular_density_s(fp, x, y, z, label, xmin, xmax, ymin, ymax, zmin, zmax, color_table)
    class(coop_asy) fp
    COOP_SINGLE ,dimension(:)::x, y, z
    COOP_UNKNOWN_STRING, optional::label, color_table
    COOP_SINGLE ,optional::xmin, xmax, ymin, ymax, zmin, zmax
    COOP_INT  i, n
    COOP_SINGLE  minx, miny, minz, maxx, maxy, maxz
    n = coop_getdim("coop_asy_irregular_density", size(x), size(y), size(z))
    if(present(xmin))then
       minx = xmin
    else
       call coop_array_get_threshold(x, 0.99, minx)
    endif
    if(present(xmax))then
       maxx = xmax
    else
       call coop_array_get_threshold(x, 0.01, maxx)
    endif
    if(present(ymin))then
       miny = ymin
    else
       call coop_array_get_threshold(y,  0.99, miny)
    endif
    if(present(ymax))then
       maxy = ymax
    else
       call coop_array_get_threshold(y, 0.01, maxy)
    endif
    if(present(zmin))then
       minz = zmin
    else
       call coop_array_get_threshold(z, 0.99, minz)
    endif
    if(present(zmax))then
       maxz = zmax
    else
       call coop_array_get_threshold(z, 0.01, maxz)
    endif
    write(fp%unit,"(A)") "DENSITY"

    if(present(color_table))then
       if(trim(color_table).ne."")then
          write(fp%unit,"(A)") trim(color_table)
       else
          write(fp%unit,"(A)") "Rainbow"
       endif
    else
       write(fp%unit,"(A)") "Rainbow"
    endif

    if(present(label))then
       if(trim(label).ne."")then
          write(fp%unit,"(A)") trim(label)
       else
          write(fp%unit,"(A)") "NULL"
       endif
    else
       write(fp%unit,"(A)") "NULL"
    endif
    call fp%write_limits(minx, maxx, miny, maxy) !real(minx, sp), real(maxx, sp), real(miny, sp), real(maxy, sp))
    write(fp%unit, "(2G14.5)") minz, maxz
    write(fp%unit,"(A)") "1"
    write(fp%unit,"(I10)") n
    do i=1,n
       write(fp%unit, "(3G14.5)") x(i), y(i), z(i)
    enddo
  end subroutine coop_asy_irregular_density_s

  subroutine coop_asy_clip_d(fp, x, y, smooth)
    class(coop_asy) fp
    logical,optional::smooth
    COOP_REAL ,dimension(:),intent(IN)::x,y
    COOP_INT  i,n
    n = coop_getdim("coop_asy_clip", size(x), size(y))
    write(fp%unit, "(A)") "CLIP"
    write(fp%unit, "(I8)") n
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    do i=1,n
       call fp%write_coor(real(x(i), sp), real(y(i), sp))
    enddo

  end subroutine coop_asy_clip_d

  subroutine coop_asy_clip_s(fp, x, y, smooth)
    class(coop_asy) fp
    logical,optional::smooth
    COOP_SINGLE ,dimension(:),intent(IN)::x,y
    COOP_INT  i,n
    n = coop_getdim("coop_asy_clip", size(x), size(y))
    write(fp%unit, "(A)") "CLIP"
    write(fp%unit, "(I8)") n
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    do i=1,n
       call fp%write_coor( x(i), y(i))
    enddo
  end subroutine coop_asy_clip_s

  function coop_asy_linestyle(color, linetype, linewidth)
    COOP_UNKNOWN_STRING color
    COOP_UNKNOWN_STRING,optional::linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING coop_asy_linestyle
    if(trim(color).eq."" .or. trim(color).eq."NULL")then
       coop_asy_linestyle = "NULL"
       return
    endif
    if(present(linetype))then
       coop_asy_linestyle = trim(color)//"_"//trim(linetype)
    else
       coop_asy_linestyle = trim(color)//"_solid"
    endif
    if(present(linewidth))then
       coop_asy_linestyle = trim(coop_asy_linestyle)//"_"//trim(coop_num2str(linewidth))
    endif
  end function coop_asy_linestyle

  subroutine coop_asy_path_append_s(path, x, y)
    COOP_SINGLE  x, y
    type(coop_asy_path) path
    call path%l%push( (/ x, y /) )
  end subroutine coop_asy_path_append_s

  subroutine coop_asy_path_append_d(path, x, y)
    COOP_REAL  x, y
    type(coop_asy_path) path
    call coop_asy_path_append_s(path, real(x), real(y))
  end subroutine coop_asy_path_append_d

  subroutine coop_asy_path_close(path)
    type(coop_asy_path) path
    if(allocated(path%length))then
       path%nclosed = path%nclosed + 1
       if(path%nclosed .gt. coop_asy_path_max_nclosed) stop "coop_asy_path_close: too many closed paths"
       path%length(path%nclosed) = path%l%n - sum(path%length(1:path%nclosed-1))
    else
       allocate(path%length(coop_asy_path_max_nclosed))
       path%nclosed = 1
       path%length(1) = path%l%n
    endif
  end subroutine coop_asy_path_close

  subroutine coop_asy_path_clear(path)
    type(coop_asy_path) path
    if(allocated(path%length))deallocate(path%length)
    call path%l%init()
  end subroutine coop_asy_path_clear

  !!generate a path that encloses the region of f(x, y) > threshold
  subroutine coop_asy_path_2dcl(path, f, xmin, xmax, ymin, ymax, threshold, resolution)
    COOP_INT ,parameter::default_resolution = 64
    COOP_INT ,optional:: resolution
    type(coop_asy_path) path
    external f
    COOP_SINGLE  f
    COOP_SINGLE  xmin, xmax, ymin, ymax, threshold, dx, dy, dxby2, dyby2
    COOP_INT  n
    logical,dimension(:,:),allocatable::above
    COOP_INT ,dimension(:,:,:),allocatable::lines
    COOP_INT  i, j, imin, jmin, last
    call coop_asy_path_clear(path) !!clear the path
    if(present(resolution))then
       n = resolution
    else
       n = default_resolution
    endif
    dx = (xmax-xmin)/n
    dxby2 = dx/2.
    dy = (ymax-ymin)/n
    dyby2 = dy/2.
    allocate(above(0:n+1, 0:n+1),lines(2, -1:n, -1:n))
    above(0,:) = .false.
    above(:,0) = .false.
    above(n+1,:) = .false.
    above(:, n+1) = .false.
    do j = 1, n
       do i=1,n
          above(i, j) = (f(xmin+dx*(i-0.5), ymin+dy*(j-0.5)) .gt. threshold)
       enddo
    enddo
    lines = 0
    do j=0, n
       do i=0, n
          if(above(i+1, j+1) .and. .not. above(i+1, j))then
             lines(1, i, j) = 1
          elseif(.not. above(i+1, j+1) .and. above(i+1, j))then
             lines(1, i, j)= -1
          endif
          if(above(i, j+1) .and. .not. above(i+1, j+1))then
             lines(2, i, j) =  1
          elseif(.not. above(i, j+1) .and. above(i+1, j+1))then
             lines(2, i, j) = -1
          endif
       enddo
    enddo
    jmin = 0
    imin = 0
100 continue
    do 
       if(lines(1, imin, jmin).eq.0)then
          imin = imin + 1
          if(imin .gt. n)then
             jmin = jmin+1
             if(jmin .gt. n)then
                deallocate(above, lines)
                return
             endif
             imin = 0
          endif
       else
          exit
       endif
    enddo
    j = jmin
    if(lines(1, imin, jmin) .eq. 1)then
       lines(1, imin, jmin) = 0
       i = imin
       call insert_point()
       i =  i + 1
       call insert_point()
       last = -1
    else
       lines(1, imin, jmin) = 0
       imin = imin + 1
       i = imin
       call insert_point()
       i = i - 1
       call insert_point()
       last = 1
    endif

    do 
       select case(last)
       case(1)
          if(lines(2, i, j-1) .eq. -1)then
             lines(2, i, j-1) = 0
             j = j - 1
             last = 2
          elseif(lines(1, i-1, j) .eq. -1)then
             lines(1, i-1, j) = 0
             i = i - 1             
          else
             lines(2, i, j) = 0
             j = j + 1
             last = -2
          endif
       case(-1)
          if(lines(2, i, j) .eq. 1)then
             lines(2, i, j) = 0
             j = j + 1
             last = -2
          elseif(lines(1, i, j) .eq. 1)then
             lines(1, i, j) = 0
             i = i + 1
          else
             lines(2, i, j-1) = 0
             j = j -1
             last =2
          endif
       case(2)
          if(lines(1, i, j) .eq. 1)then
             lines(1, i, j) = 0
             i = i + 1
             last  = -1
          elseif(lines(2, i, j-1).eq.-1)then
             lines(2, i, j-1) = 0
             j = j - 1
          else
             lines(1, i-1, j) = 0
             i = i - 1
             last = 1
          endif
       case default !!actually must be -2
          if(lines(1, i-1, j).eq. -1)then
             lines(1, i-1, j) = 0
             i = i-1
             last = 1
          elseif(lines(2, i, j) .eq. 1)then
             lines(2, i, j) = 0
             j = j + 1
          else 
             lines(1, i, j) = 0
             i = i + 1
             last = -1
          endif
       end select
       if(i.eq.imin .and. j.eq.jmin)then !!this happens in finite steps
          call coop_asy_path_close(path)
          goto 100 
       else
          call insert_point()
       endif
    enddo

  contains

    subroutine insert_point()
      COOP_SINGLE  x, y
      COOP_SINGLE  fp(2), s
      x =  xmin+i*dx
      y = ymin+j*dy
      if(i.eq.0 .or. i.eq.n)then
         fp(1) = 0 
      else
         fp(1) = (f(x+dxby2, y) - f(x-dxby2, y))
      endif
      if(j.eq.0 .or. j.eq.n)then
         fp(2) = 0
      else
         fp(2) =  (f(x, y+dyby2) - f(x, y-dyby2))
      endif
      s= sum(fp**2)
      if(s.gt.0)then
         s = (threshold - f(x, y))/s
         x = min(xmax, max(x + max(min(fp(1)*s, 0.5), -0.5)*dx, xmin))
         y = min(ymax, max(y + max(min(fp(2)*s, 0.5), -0.5)*dy, ymin))
      endif
      call coop_asy_path_append(path, x, y)
    end subroutine insert_point

  end subroutine coop_asy_path_2dcl


  subroutine coop_asy_path_from_array_center(path, f, xmin, xmax, ymin, ymax, threshold)
    type(coop_asy_path) path
    COOP_SINGLE  f(:,:)
    COOP_SINGLE  xmin, xmax, ymin, ymax, threshold, dx, dy
    COOP_INT  nx, ny, n
    nx = size(f, 1)
    ny = size(f, 2)
    n = min(max(nx, ny,16)*2, 256) !!higher resolution is useless since we do bilinear interpolation
    dx = (xmax-xmin)/nx
    dy = (ymax- ymin)/ny
    call coop_asy_path_2dcl(path, finterp, xmin, xmax, ymin, ymax, threshold, n)
  contains
    
    function  finterp(x, y)
      COOP_SINGLE  x, y, finterp
      COOP_INT  i, j, ip1, jp1
      COOP_SINGLE  ri, rj
      ri = (x-xmin)/dx + 0.5
      rj = (y-ymin)/dy + 0.5
      i = floor(ri)
      j = floor(rj)
      ri = ri - i
      rj = rj - j
      ip1 = i + 1
      jp1 = j + 1
      i = min(max(i,1), nx)
      ip1 = min(max(ip1,1), nx)
      j = min(max(j, 1), ny)
      jp1 = min(max(jp1, 1), ny)
      finterp = (f(i, j) * (1.-ri) + f(ip1, j)*ri)*(1.-rj)+(f(i, jp1)*(1.-ri) + f(ip1, jp1)*ri)*rj
    end function finterp
    
  end subroutine coop_asy_path_from_array_center



  subroutine coop_asy_path_from_array(path, f, xmin, xmax, ymin, ymax, threshold)
    type(coop_asy_path) path
    COOP_SINGLE  f(:,:)
    COOP_SINGLE  xmin, xmax, ymin, ymax, threshold, dx, dy
    COOP_INT  nx, ny, n
    nx = size(f, 1)
    ny = size(f, 2)
    n = min(max(nx, ny,16)*2, 256) !!higher resolution is useless since we do bilinear interpolation
    dx = (xmax-xmin)/(nx-1)
    dy = (ymax- ymin)/(ny-1)
    call coop_asy_path_2dcl(path, finterp, xmin, xmax, ymin, ymax, threshold, n)
  contains
    
    function  finterp(x, y)
      COOP_SINGLE  x, y, finterp
      COOP_INT  i, j
      COOP_SINGLE  ri, rj
      ri = (x-xmin)/dx + 1.
      rj = (y-ymin)/dy + 1.
      i = min(max(floor(ri), 1), nx-1)
      j = min(max(floor(rj), 1), ny-1)
      ri = ri - i
      rj = rj - j
      finterp = (f(i, j)*(1.-ri) + f(i+1, j)*ri)*(1.-rj)+(f(i, j+1)*(1.-ri) + f(i+1, j+1)*ri)*rj
    end function finterp
    
  end subroutine coop_asy_path_from_array


  subroutine coop_asy_legend_default(fp)
    class(coop_asy) fp
    call coop_asy_legend_location(fp, "N")
  end subroutine coop_asy_legend_default

  subroutine coop_asy_legend_location(fp, loc, cols)
    class(coop_asy) fp
    COOP_UNKNOWN_STRING :: loc
    COOP_INT , optional::cols
    write(fp%unit, "(A)") "LEGEND"
    if(trim(loc).ne."")then
       write(fp%unit, "(A)") trim(loc)
    else
       write(fp%unit, "(A)") "E"
    end if
    if(present(cols))then
       write(fp%unit, "(I8)")  cols
    else
       write(fp%unit, "(I8)")  1
    endif
  end subroutine coop_asy_legend_location

  subroutine coop_asy_legend_d(fp, x, y, cols)
    class(coop_asy) fp
    COOP_REAL x, y
    COOP_INT ,optional::cols
    write(fp%unit, "(A)") "LEGEND"
    write(fp%unit, "(A)") "NULL"
    write(fp%unit, "(2G15.4)") x, y
    if(present(cols))then
       write(fp%unit, "(I8)")  cols
    else
       write(fp%unit, "(I8)")  1
    endif
  end subroutine coop_asy_legend_d

  subroutine coop_asy_legend_s(fp, x, y, cols)
    class(coop_asy) fp
    COOP_SINGLE  x, y
    COOP_INT ,optional::cols
    write(fp%unit, "(A)") "LEGEND"
    write(fp%unit, "(A)") "NULL"
    write(fp%unit, "(2G15.4)") x, y
    if(present(cols))then
       write(fp%unit, "(I8)")  cols
    else
       write(fp%unit, "(I8)")  1
    endif
  end subroutine coop_asy_legend_s

  subroutine coop_asy_legend_relative(fp, xratio, yratio, cols)
    class(coop_asy) fp
    COOP_SINGLE  xratio, yratio, x, y
    COOP_INT ,optional::cols
    write(fp%unit, "(A)") "LEGEND"
    write(fp%unit, "(A)") "NULL"
    write(fp%unit, "(2G15.4)") fp%xrel(xratio), fp%yrel(yratio)
    if(present(cols))then
       write(fp%unit, "(I8)")  cols
    else
       write(fp%unit, "(I8)")  1
    endif
  end subroutine coop_asy_legend_relative
  

  subroutine coop_asy_topaxis_s(fp, xmin, xmax, islog,label)
    COOP_SINGLE  xmin, xmax
    logical islog
    class(coop_asy) fp
    COOP_UNKNOWN_STRING label
    write(fp%unit, "(A)")"EXTRA_AXIS"
    write(fp%unit, "(A)") "top"
    if(trim(label).ne."")then
       write(fp%unit, "(A)") trim(label)
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(islog)then
       write(fp%unit, "(A)") "1"
    else
       write(fp%unit, "(A)") "0"
    endif
    write(fp%unit, "(2G14.5)") xmin, xmax
  end subroutine coop_asy_topaxis_s

  subroutine coop_asy_topaxis_d(fp, xmin, xmax, islog,label)
    COOP_REAL  xmin, xmax
    logical islog
    class(coop_asy) fp
    COOP_UNKNOWN_STRING label
    call coop_asy_topaxis_s(fp, real(xmin), real(xmax), islog,label)
  end subroutine coop_asy_topaxis_d


  subroutine coop_asy_rightaxis_s(fp, ymin, ymax, islog, label)
    COOP_SINGLE  ymin, ymax
    logical islog
    class(coop_asy) fp
    COOP_UNKNOWN_STRING label
    write(fp%unit, "(A)")"EXTRA_AXIS"
    write(fp%unit, "(A)") "right"
    if(trim(label).ne."")then
       write(fp%unit, "(A)") trim(label)
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(islog)then
       write(fp%unit, "(A)") "1"
    else
       write(fp%unit, "(A)") "0"
    endif
    write(fp%unit, "(2G14.5)") ymin, ymax
  end subroutine coop_asy_rightaxis_s

  subroutine coop_asy_rightaxis_d(fp, ymin, ymax, islog, label)
    COOP_REAL  ymin, ymax
    logical islog
    COOP_UNKNOWN_STRING label
    class(coop_asy) fp
    call coop_asy_rightaxis_s(fp, real(ymin), real(ymax), islog,label)
  end subroutine coop_asy_rightaxis_d

  subroutine coop_asy_error_bar_d(fp, x, y, dy_minus, dy_plus, dx_minus, dx_plus, color)
    class(coop_asy) fp
    COOP_REAL  x, y
    COOP_REAL ,optional::dy_minus, dy_plus, dx_minus, dx_plus
    COOP_UNKNOWN_STRING, optional::color
    if(present(color))then
       call coop_asy_label(fp, "$\bullet$", x, y, color)
       if(present(dy_minus))then
          call coop_asy_line(fp, x, y, x, y-dy_minus, linewidth = 1., color=color)
          call coop_asy_label(fp, "-", x, y-dy_minus, color=color)
       endif
       if(present(dy_plus))then
          call coop_asy_line(fp, x, y, x, y+dy_plus, linewidth = 1., color= color)
          call coop_asy_label(fp, "-", x, y+dy_plus, color=color)
       endif
       if(present(dx_minus))then
          call coop_asy_line(fp, x, y, x-dx_minus, y, linewidth = 1., color=color)
          call coop_asy_label(fp, "{\tiny $|$}", x-dx_minus, y, color=color)
       endif
       if(present(dx_plus))then
          call coop_asy_line(fp, x, y, x+dx_plus, y, linewidth = 1., color=color)
          call coop_asy_label(fp, "{\tiny $|$}", x+dx_plus, y, color=color)
       endif
    else
       call coop_asy_label(fp, "$\bullet$", x, y)
       if(present(dy_minus))then
          call coop_asy_line(fp, x, y, x, y-dy_minus, linewidth = 1.)
          call coop_asy_label(fp, "-", x, y-dy_minus)
       endif
       if(present(dy_plus))then
          call coop_asy_line(fp, x, y, x, y+dy_plus, linewidth = 1.)
          call coop_asy_label(fp, "-", x, y+dy_plus)
       endif
       if(present(dx_minus))then
          call coop_asy_line(fp, x, y, x-dx_minus, y, linewidth = 1.)
          call coop_asy_label(fp, "{\tiny $|$}", x-dx_minus, y)
       endif
       if(present(dx_plus))then
          call coop_asy_line(fp, x, y, x+dx_plus, y, linewidth = 1.)
          call coop_asy_label(fp, "{\tiny $|$}", x+dx_plus, y)
       endif
    endif
  end subroutine coop_asy_error_bar_d

  subroutine coop_asy_error_bar_s(fp, x, y, dy_minus, dy_plus, dx_minus, dx_plus, color)
    class(coop_asy) fp
    COOP_SINGLE  x, y
    COOP_SINGLE ,optional::dy_minus, dy_plus, dx_minus, dx_plus
    COOP_UNKNOWN_STRING, optional::color
    if(present(color))then
       call coop_asy_label(fp, "$\bullet$", x, y, color)
       if(present(dy_minus))then
          call coop_asy_line(fp, x, y, x, y-dy_minus, linewidth = 1., color=color)
          call coop_asy_label(fp, "-", x, y-dy_minus, color=color)
       endif
       if(present(dy_plus))then
          call coop_asy_line(fp, x, y, x, y+dy_plus, linewidth = 1., color= color)
          call coop_asy_label(fp, "-", x, y+dy_plus, color=color)
       endif
       if(present(dx_minus))then
          call coop_asy_line(fp, x, y, x-dx_minus, y, linewidth = 1., color=color)
          call coop_asy_label(fp, "{\tiny $|$}", x-dx_minus, y, color=color)
       endif
       if(present(dx_plus))then
          call coop_asy_line(fp, x, y, x+dx_plus, y, linewidth = 1., color=color)
          call coop_asy_label(fp, "{\tiny $|$}", x+dx_plus, y, color=color)
       endif
    else
       call coop_asy_label(fp, "$\bullet$", x, y)
       if(present(dy_minus))then
          call coop_asy_line(fp, x, y, x, y-dy_minus, linewidth = 1.)
          call coop_asy_label(fp, "-", x, y-dy_minus)
       endif
       if(present(dy_plus))then
          call coop_asy_line(fp, x, y, x, y+dy_plus, linewidth = 1.)
          call coop_asy_label(fp, "-", x, y+dy_plus)
       endif
       if(present(dx_minus))then
          call coop_asy_line(fp, x, y, x-dx_minus, y, linewidth = 1.)
          call coop_asy_label(fp, "{\tiny $|$}", x-dx_minus, y)
       endif
       if(present(dx_plus))then
          call coop_asy_line(fp, x, y, x+dx_plus, y, linewidth = 1.)
          call coop_asy_label(fp, "{\tiny $|$}", x+dx_plus, y)
       endif
    endif
  end subroutine coop_asy_error_bar_s


  function coop_asy_rgb_color_s(r, g, b) result(color)
    COOP_SINGLE  r, g, b
    COOP_SHORT_STRING color
    color = "RGB:"//trim(coop_num2str(nint(min(1., max(0., r))*255.*8.)/8.))//":"//trim(coop_num2str(nint(min(1., max(0., g))*255.*8.)/8.))//":"//trim(coop_num2str(nint(min(1., max(0., b))*255.*8.)/8.))
  end function coop_asy_rgb_color_s

  function coop_asy_rgb_color_d(r,g,b) result(color)
    COOP_REAL  r, g, b
    COOP_SHORT_STRING color
    color = coop_asy_rgb_color_s(real(r), real(g), real(b))
  end function coop_asy_rgb_color_d


  function coop_asy_gray_color_s(gray) result(color)
    COOP_SINGLE  gray
    COOP_SHORT_STRING color
    color = "GRAY:"//trim(coop_num2str(nint(min(1., max(0., gray))*255.*8.)/8.))
  end function coop_asy_gray_color_s


  function coop_asy_gray_color_d(gray) result(color)
    COOP_REAL  gray
    COOP_SHORT_STRING color
    color = coop_asy_gray_color_s(real(gray))
  end function coop_asy_gray_color_d

  subroutine coop_asy_histogram(x, nbins, filename)
    COOP_UNKNOWN_STRING::filename
    COOP_REAL x(:), dx
    COOP_INT nbins, iloc, i
    COOP_REAL xmin, xmax, c(nbins), xb(nbins)
    type(coop_asy)::fig
    call coop_array_get_threshold(x, 0.997d0, xmin)
    call coop_array_get_threshold(x, 0.003d0, xmax)
    dx = (xmax-xmin)/nbins
    if(dx .eq. 0.d0) stop "data with single point cannot be histogrammed"
    call coop_set_uniform(nbins, xb, xmin + dx/2.d0, xmax-dx/2.d0)
    c = 0.
    call fig%open(trim(filename))
    call fig%init(xlabel = "x", ylabel = "P")
    do i=1, size(x)
       iloc = nint((x(i) - xmin)/dx+0.5)
       if(iloc.gt.0 .and. iloc .le.nbins) c(iloc) = c(iloc)+1
    enddo
    call coop_asy_curve(fig, xb, c)
    call fig%close()
  end subroutine coop_asy_histogram

  subroutine coop_asy_plot_function(func, filename, xlabel, ylabel)
    type(coop_function)func
    COOP_UNKNOWN_STRING::filename
    COOP_UNKNOWN_STRING,optional::xlabel, ylabel
    type(coop_asy)::asy
    COOP_INT, parameter::n = 512
    COOP_REAL:: x(n), y(n)
    integer i
    call asy%open(trim(filename))
    if(present(xlabel))then
       if(present(ylabel))then
          call asy%init(xlabel = xlabel, ylabel = ylabel, xlog = func%xlog, ylog = func%ylog)
       else
          call asy%init(xlabel = xlabel, ylabel = "$f(x)$", xlog = func%xlog, ylog = func%ylog)
       endif
    else
       if(present(ylabel))then
          call asy%init(xlabel = "$x$", ylabel = ylabel, xlog = func%xlog, ylog = func%ylog)
       else
          call asy%init(xlabel = "$x$", ylabel = "$f(x)$", xlog = func%xlog, ylog = func%ylog)
       endif
    endif
    call coop_set_uniform(n, x,func%xmin, func%xmax)
    if(func%xlog) x = exp(x)
    do i=1, n
       y(i) = func%eval(x(i))
    enddo
    call coop_asy_curve(asy, x, y)
    call asy%close()
  end subroutine coop_asy_plot_function



  subroutine coop_asy_band_d(fp, x, ylower, yupper, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) fp
    logical,optional::smooth
    COOP_REAL ,dimension(:),intent(IN)::x,ylower, yupper
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_INT  i,n
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_contour", size(x), size(ylower), size(yupper))
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "1" !!type 1 contour
    if(present(colorfill))then
       write(fp%unit, "(A)") trim(colorfill)
    else
       write(fp%unit, "(A)") "lightgray"
    endif
    if(present(linecolor).or. present(linetype).or.present(linewidth))then
       if(present(linecolor))then
          lineproperty=trim(linecolor)
       else
          lineproperty = "black"
       endif
       if(present(linetype))then
          lineproperty = trim(lineproperty)//"_"//trim(linetype)
       else
          lineproperty = trim(lineproperty)//"_solid"
       endif
       if(present(linewidth))then
          lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
       endif
       write(fp%unit, "(A)") trim(lineproperty)
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(fp%unit, "(A)") "1"  !!1 path 
    write(fp%unit, "(I8)") n*2
    if(x(2).gt. x(1))then
       do i = 1, n
          call fp%write_coor(real(x(i), sp), real(ylower(i), sp))
       enddo
       do i = n, 1, -1
          call fp%write_coor(real(x(i), sp), real(yupper(i), sp))
       enddo
    else
       do i = n, 1, -1
          call fp%write_coor(real(x(i), sp), real(ylower(i), sp))
       enddo
       do i = 1, n
          call fp%write_coor(real(x(i), sp), real(yupper(i), sp))
       enddo
    endif
  end subroutine coop_asy_band_d

  subroutine coop_asy_band_s(fp, x, ylower, yupper, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) fp
    logical,optional::smooth
    COOP_SINGLE ,dimension(:),intent(IN)::x,ylower, yupper
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_INT  i,n
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_contour", size(x), size(ylower), size(yupper))
    write(fp%unit, "(A)") "CONTOUR"
    write(fp%unit, "(A)") "1" !!type 1 contour
    if(present(colorfill))then
       write(fp%unit, "(A)") trim(colorfill)
    else
       write(fp%unit, "(A)") "lightgray"
    endif
    if(present(linecolor).or. present(linetype).or.present(linewidth))then
       if(present(linecolor))then
          lineproperty=trim(linecolor)
       else
          lineproperty = "black"
       endif
       if(present(linetype))then
          lineproperty = trim(lineproperty)//"_"//trim(linetype)
       else
          lineproperty = trim(lineproperty)//"_solid"
       endif
       if(present(linewidth))then
          lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
       endif
       write(fp%unit, "(A)") trim(lineproperty)
    else
       write(fp%unit, "(A)") "NULL"
    endif
    if(present(smooth))then
       if(smooth)then
          write(fp%unit, "(A)") "1"
       else
          write(fp%unit, "(A)") "0"
       endif
    else
       write(fp%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(fp%unit, "(A)") "1"  !!1 path 
    write(fp%unit, "(I8)") n*2
    if(x(2).gt. x(1))then
       do i = 1, n
          call fp%write_coor(real(x(i), sp), real(ylower(i), sp))
       enddo
       do i = n, 1, -1
          call fp%write_coor(real(x(i), sp), real(yupper(i), sp))
       enddo
    else
       do i = n, 1, -1
          call fp%write_coor(real(x(i), sp), real(ylower(i), sp))
       enddo
       do i = 1, n
          call fp%write_coor(real(x(i), sp), real(yupper(i), sp))
       enddo
    endif
  end subroutine coop_asy_band_s

  function coop_asy_xrel(fp, xratio) result(xrel)
    class(coop_asy)::fp
    COOP_SINGLE  xratio
    COOP_SINGLE  xrel
    if(fp%xlog)then
       xrel = fp%xmax**xratio*fp%xmin**(1.d0-xratio)
    else
       xrel = xratio*fp%xmax + (1.d0-xratio)*fp%xmin
    endif
  end function coop_asy_xrel

  function coop_asy_yrel(fp, yratio) result(yrel)
    class(coop_asy)::fp
    COOP_SINGLE  yratio
    COOP_SINGLE  yrel
    if(fp%ylog)then
       yrel = fp%ymax**yratio*fp%ymin**(1.d0-yratio)
    else
       yrel = yratio*fp%ymax + (1.d0-yratio)*fp%ymin
    endif
  end function coop_asy_yrel
  
  subroutine coop_asy_expand(fp, xl, xr, yl, yr)
    class(coop_asy) fp
    COOP_SINGLE xl, xr, yl, yr, dx, dy
    write(fp%unit, "(A)") "EXPAND"
    write(fp%unit, "(4G14.5)")  xl, xr, yl, yr
    dx = fp%xmax - fp%xmin
    dy = fp%ymax - fp%ymin
    fp%xmin = fp%xmin - dx*xl
    fp%xmax = fp%xmax + dx*xr
    fp%ymin = fp%ymin - dy*yl
    fp%ymax = fp%ymax + dy*yr
  end subroutine coop_asy_expand

  subroutine coop_asy_arrow_d(fp, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy)::fp
    COOP_REAL xstart, ystart, xend, yend
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING lineproperty
    write(fp%unit, "(A)") "ARROWS"
    write(fp%unit, "(A)") "1"
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    call fp%write_coor(real(xstart, sp), real(ystart, sp), real(xend, sp), real(yend, sp))    
  end subroutine coop_asy_arrow_d

  subroutine coop_asy_arrow_s(fp, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy)::fp
    COOP_SINGLE xstart, ystart, xend, yend
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING lineproperty
    write(fp%unit, "(A)") "ARROWS"
    write(fp%unit, "(A)") "1"
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    call fp%write_coor(xstart, ystart, xend, yend)
  end subroutine coop_asy_arrow_s

  subroutine coop_asy_arrows_d(fp, n, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy)::fp
    COOP_INT n, i
    COOP_REAL xstart(n), ystart(n), xend(n), yend(n)
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING lineproperty
    if(n.gt. 9999)then
       write(*, '(A, I10)') "n = ", n
       stop "cannot plot so many arrows"
    endif
    write(fp%unit, "(A)") "ARROWS"
    write(fp%unit, "(I5)") n
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    do i=1, n
       call fp%write_coor(real(xstart(i), sp), real(ystart(i), sp), real(xend(i), sp), real(yend(i), sp))
    enddo
  end subroutine coop_asy_arrows_d
  
  subroutine coop_asy_arrows_s(fp, n, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy)::fp
    COOP_INT n, i
    COOP_SINGLE xstart(n), ystart(n), xend(n), yend(n)
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING lineproperty
    if(n.gt. 9999)then
       write(*, '(A, I10)') "n = ", n
       stop "cannot plot so many arrows"
    endif
    write(fp%unit, "(A)") "ARROWS"
    write(fp%unit, "(I5)") n
    if(present(color))then
       lineproperty=trim(color)
    else
       lineproperty = "black"
    endif
    if(present(linetype))then
       lineproperty = trim(lineproperty)//"_"//trim(linetype)
    else
       lineproperty = trim(lineproperty)//"_solid"
    endif
    if(present(linewidth))then
       lineproperty = trim(lineproperty)//"_"//trim(coop_num2str(linewidth))
    endif
    write(fp%unit, "(A)") trim(lineproperty)
    do i=1, n
       call fp%write_coor(xstart(i), ystart(i), xend(i), yend(i))
    enddo
  end subroutine coop_asy_arrows_s
  
end module coop_asy_mod
