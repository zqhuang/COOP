module coop_asy_mod
  use coop_wrapper_typedef
  use coop_interpolation_mod
  use coop_file_mod
  use coop_list_mod
  use coop_random_mod
  use coop_evalstr_mod
  implicit none
  private

  public::coop_asy, coop_asy_path, coop_asy_error_bar, coop_asy_errorbars, coop_asy_interpolate_curve, coop_asy_gray_color, coop_asy_rgb_color, coop_asy_label, coop_asy_legend, coop_asy_legend_advance, coop_asy_dot, coop_asy_line, coop_asy_labels, coop_asy_dots, coop_asy_lines, coop_asy_contour, coop_asy_curve, coop_asy_density,  coop_asy_topaxis, coop_asy_rightaxis, coop_asy_clip, coop_asy_plot_function, coop_asy_plot_likelihood,  coop_asy_histogram, coop_asy_band,coop_asy_default_width,coop_asy_default_height, coop_asy_color2rgba, coop_asy_rgba2color


#include "constants.h"

  COOP_INT, parameter::sp = coop_single_real_length
  COOP_INT, parameter::dl = coop_real_length
  COOP_SINGLE::coop_asy_default_width = 4.8
  COOP_SINGLE::coop_asy_default_height = 3.9
  integer, parameter::coop_asy_num_line_types = 12

  type, extends(coop_file) :: coop_asy
     COOP_SINGLE::xmin, xmax, ymin, ymax
     COOP_SINGLE::width = 8.
     COOP_SINGLE::height = 6.
     logical::xlog = .false.
     logical::ylog = .false.
     logical::clip = .false.
     logical::adjust_xmin = .true.
     logical::adjust_xmax = .true.
     logical::adjust_ymin = .true.
     logical::adjust_ymax = .true.
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
     procedure::plot_file => coop_asy_plot_file
     procedure::legend => coop_asy_legend_relative
     procedure::compact_legend=>coop_asy_compact_legend_relative
     procedure::legend_advance => coop_asy_legend_advance_relative
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
     procedure::errorbars => coop_asy_errorbars_d
     procedure::expand => coop_asy_expand
     procedure::arrow => coop_asy_arrow_d
     procedure::arrows => coop_asy_arrows_d
     procedure::annotate => coop_asy_annotate_d
     procedure::adjust_virtual_boundary => coop_asy_adjust_virtual_boundary
  end type coop_asy



  COOP_INT , parameter::coop_asy_path_max_nclosed = 4096

  interface coop_asy_arrow
     module procedure coop_asy_arrow_d, coop_asy_arrow_s
  end interface coop_asy_arrow

  interface coop_asy_annotate
     module procedure coop_asy_annotate_d, coop_asy_annotate_s
  end interface coop_asy_annotate


  interface coop_asy_arrows
     module procedure coop_asy_arrows_d, coop_asy_arrows_s
  end interface coop_asy_arrows
  
  
  interface coop_asy_band
     module procedure coop_asy_band_s, coop_asy_band_d
  end interface coop_asy_band

  interface coop_asy_error_bar
     module procedure coop_asy_error_bar_s, coop_asy_error_bar_d
  end interface coop_asy_error_bar


  interface coop_asy_errorbars
     module procedure coop_asy_errorbars_s, coop_asy_errorbars_d
  end interface coop_asy_errorbars
  

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

  interface coop_asy_legend_advance
     module procedure coop_asy_legend_advance_s, coop_asy_legend_advance_d, coop_asy_legend_advance_location     
  end interface coop_asy_legend_advance
  

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
     module procedure coop_asy_dots_s, coop_asy_dots_d, coop_asy_dots_list
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

  interface coop_asy_color2rgba
     module procedure coop_asy_color2rgba_d, coop_asy_color2rgba_s
  end interface coop_asy_color2rgba


  interface coop_asy_rgba2color
     module procedure coop_asy_rgba2color_d, coop_asy_rgba2color_s
  end interface coop_asy_rgba2color


  type coop_asy_path
     COOP_INT::nclosed = 0
     COOP_INT,dimension(:),allocatable:: length
     type(coop_list_realarr) l
   contains
     procedure::append => coop_asy_path_append_s
     procedure::close => coop_asy_path_close
     procedure::free => coop_asy_path_free
     procedure::from_array => coop_asy_path_from_array
     procedure::from_array_gaussianfit => coop_asy_path_from_array_gaussianfit     
     procedure::from_array_center => coop_asy_path_from_array_center     
     procedure::from_function => coop_asy_path_from_function
     procedure::get_perimeter_and_area => coop_asy_path_get_perimeter_and_area
  end type coop_asy_path



contains

  subroutine coop_asy_adjust_virtual_boundary(this, x, y)
    class(coop_asy)::this
    COOP_SINGLE::x, y
    if(this%adjust_xmin .and. x .lt. this%xmin) &
         this%xmin = x
    if(this%adjust_xmax .and. x .gt. this%xmax) &
         this%xmax = x
    if(this%adjust_ymin .and. y .lt. this%ymin) &
         this%ymin = y
    if(this%adjust_ymax .and. y .gt. this%ymax) &
         this%ymax = y
  end subroutine coop_asy_adjust_virtual_boundary

  subroutine coop_asy_init(this,  xmin, xmax, ymin, ymax, width, height, caption, xlabel, ylabel, xlog, ylog, zlog, doclip, nblocks, nxticks, nyticks)
    class(coop_asy) this
    COOP_SINGLE ,optional:: width, height
    COOP_UNKNOWN_STRING,optional:: caption, xlabel, ylabel
    COOP_INT ,optional::nblocks, nxticks, nyticks
    logical,optional::xlog, ylog, zlog, doclip
    COOP_SINGLE ,optional:: xmin, xmax, ymin, ymax
    character(len = 5) tmp
    if(present(width))then
       this%width = width
    else
       this%width = coop_asy_default_width
    endif
    if(present(height))then
       this%height = height
    else
       this%height = coop_asy_default_height
    endif
    write(this%unit, "(2F12.2)") this%width, this%height
    if(present(caption))then
       if(trim(adjustl(caption)) .eq."")then
          write(this%unit, "(A)") "NULL"
       else
          write(this%unit, "(A)") trim(adjustl(caption))
       endif
    else
       write(this%unit, "(A)") "NULL"
    endif
    if(present(xlabel))then
       if(trim(adjustl(xlabel)).eq."")then
          write(this%unit, "(A)") "NULL"
       else
          write(this%unit, "(A)") trim(adjustl(xlabel))
       endif
    else
       write(this%unit, "(A)") "NULL"
    endif
    if(present(ylabel))then
       if(trim(adjustl(ylabel)).eq."")then
          write(this%unit, "(A)") "NULL"
       else
          write(this%unit, "(A)") trim(adjustl(ylabel))
       endif
    else
       write(this%unit, "(A)") "NULL"
    endif
    if(present(xlog))then
       this%xlog  =xlog
       if(xlog)then
          tmp(1:2) = "1 "
       else
          tmp(1:2) = "0 "
       endif
    else
       this%xlog = .false.
       tmp(1:2) = "0 "
    endif
    if(present(ylog))then
       this%ylog = ylog
       if(ylog)then
          tmp(3:4) = "1 "
       else
          tmp(3:4) = "0 "
       endif
    else
       this%ylog = .false.
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
    write(this%unit, "(A)") tmp
    if(present(doclip))then
       this%clip = doclip
    else
       this%clip = .false.
    endif
    if(present(xmin))then
       this%xmin = xmin
       this%adjust_xmin = .false.
    else
       this%xmin = 1.1e30
       this%adjust_xmin = .true.
       this%clip = .false.                  
    endif
    if(present(xmax))then
       this%xmax = xmax
       this%adjust_xmax = .false.
    else
       this%xmax = -1.1e30
       this%adjust_xmax = .true.
       this%clip = .false.                  
    endif
    if(present(ymin))then
       this%ymin = ymin
       this%adjust_ymin = .false.
    else
       this%ymin = 1.1e30
       this%adjust_ymin = .true.
       this%clip = .false.                  
    endif
    if(present(ymax))then
       this%adjust_ymax = .false.
       this%ymax = ymax
    else
       this%adjust_ymax = .true.
       this%ymax = -1.1e30
       this%clip = .false.           
    endif
    if(this%clip)then
       write(this%unit, "(A)") "1"
    else
       write(this%unit, "(A)") "0"
    endif
    if(present(nxticks))then
       write(this%unit,"(2G16.7, I5)") this%xmin, this%xmax, nxticks
    else
       write(this%unit,"(2G16.7, I5)") this%xmin, this%xmax, 5
    endif
    if(present(nyticks))then
       write(this%unit,"(2G16.7, I5)") this%ymin, this%ymax, nyticks
    else
       write(this%unit,"(2G16.7, I5)") this%ymin, this%ymax, 5
    endif
    if(present(nblocks))then
       write(this%unit,"(I5)") nblocks    !!plot only n blocks
    else
       write(this%unit,"(I5)") 0   !!0 means any number of blocks
    endif
    this%color(1) = "black"
    this%color(2) = "red"
    this%color(3) = "blue"
    this%color(4) = "cyan"
    this%color(5) = "violet"
    this%color(6) = "gray"
    this%color(7) = "skyblue"
    this%color(8) = "orange"
    this%color(9) = "green"
    this%color(10) = "brown"
    this%color(11) = "pink"
    this%color(12) = "yellow"
    this%linewidth(1:3) = 2.
    this%linewidth(4:6) = 1.2
    this%linewidth(7:12) = 0.6
    this%linetype(1) = "solid"
    this%linetype(2) = "dotted"
    this%linetype(3) = "dashed"
    this%linetype(4) = "dashdotted"
    this%linetype(5) = "longdashed"
    this%linetype(6) = "longdashdotted"
    this%linetype(7:12) = "solid"
  end subroutine coop_asy_init

  subroutine coop_asy_write_limits(this, xmin, xmax, ymin, ymax)
    class(coop_asy) this
    COOP_SINGLE  xmin, xmax, ymin, ymax
    write(this%unit, "(2G16.7)") xmin, xmax
    write(this%unit, "(2G16.7)") ymin, ymax
    call this%adjust_virtual_boundary(xmin, ymin)
    call this%adjust_virtual_boundary(xmax, ymax)    
  end subroutine coop_asy_write_limits

  subroutine coop_asy_write_coor(this, x, y, x2, y2)
    class(coop_asy) this
    COOP_SINGLE  x, y
    COOP_SINGLE ,optional:: x2, y2
    call this%adjust_virtual_boundary(x, y)
    if(present(x2) .and. present(y2))then
       write(this%unit, "(4G16.7)") x, y, x2, y2
       call this%adjust_virtual_boundary(x2, y2)
    else
       write(this%unit, "(2G16.7)") x, y
    endif
  end subroutine coop_asy_write_coor

  subroutine coop_asy_dot_d(this, x, y, color, symbol)
    class(coop_asy) this
    COOP_INT  n
    COOP_REAL  x,y
    COOP_UNKNOWN_STRING,optional:: color
    COOP_UNKNOWN_STRING,optional::symbol
    write(this%unit, "(A)") "DOTS"
    write(this%unit, "(A)") "1"
    if(present(color))then
       write(this%unit, "(A)") trim(color)
    else
       write(this%unit, "(A)") "black"
    endif
    if(present(symbol))then
       write(this%unit, "(A)") trim(symbol)
    else
       write(this%unit, "(A)") "dot"
    endif
    call this%write_coor(real(x, sp), real(y, sp))
  end subroutine coop_asy_dot_d


  subroutine coop_asy_dot_s(this, x, y, color, symbol)
    class(coop_asy) this
    COOP_INT  n
    COOP_SINGLE  x,y
    COOP_UNKNOWN_STRING,optional:: color
    COOP_UNKNOWN_STRING,optional::symbol
    write(this%unit, "(A)") "DOTS"
    write(this%unit, "(A)") "1"
    if(present(color))then
       write(this%unit, "(A)") trim(color)
    else
       write(this%unit, "(A)") "black"
    endif
    if(present(symbol))then
       write(this%unit, "(A)") trim(symbol)
    else
       write(this%unit, "(A)") "dot"
    endif
    call this%write_coor(x, y)
  end subroutine coop_asy_dot_s


  subroutine coop_asy_dots_d(this, x, y, color, symbol)
    class(coop_asy) this
    COOP_INT  n,i
    COOP_REAL ,dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional:: color
    COOP_UNKNOWN_STRING,optional::symbol
    n = coop_getdim( "coop_asy_dot_block", size(x), size(y))
    write(this%unit, "(A)") "DOTS"
    write(this%unit, "(I8)") n
    if(present(color))then
       write(this%unit, "(A)") trim(color)
    else
       write(this%unit, "(A)") "black"
    endif
    if(present(symbol))then
       write(this%unit, "(A)") trim(symbol)
    else
       write(this%unit, "(A)") "dot"
    endif
    do i=1,n
       call this%write_coor(real(x(i),sp), real(y(i),sp))
    enddo
  end subroutine coop_asy_dots_d

  subroutine coop_asy_dots_s(this, x, y, color, symbol)
    class(coop_asy) this
    COOP_INT  n,i
    COOP_SINGLE ,dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional:: color
    COOP_UNKNOWN_STRING,optional::symbol
    n = coop_getdim( "coop_asy_dot_block", size(x), size(y))
    write(this%unit, "(A)") "DOTS"
    write(this%unit, "(I8)") n
    if(present(color))then
       write(this%unit, "(A)") trim(color)
    else
       write(this%unit, "(A)") "black"
    endif
    if(present(symbol))then
       write(this%unit, "(A)") trim(symbol)
    else
       write(this%unit, "(A)") "dot"
    endif
    do i=1,n
       call this%write_coor(x(i),y(i))
    enddo
  end subroutine coop_asy_dots_s


  subroutine coop_asy_dots_list(this, xylist, xunit, yunit, color, symbol)
    class(coop_asy) this
    COOP_INT  i
    type(coop_list_realarr)::xylist
    COOP_SINGLE,optional::xunit, yunit
    COOP_SINGLE:: xy(2)
    COOP_UNKNOWN_STRING,optional:: color
    COOP_UNKNOWN_STRING,optional::symbol
    if(xylist%n .eq. 0) return
    write(this%unit, "(A)") "DOTS"
    write(this%unit, "(I8)") xylist%n
    if(present(color))then
       write(this%unit, "(A)") trim(color)
    else
       write(this%unit, "(A)") "black"
    endif
    if(present(symbol))then
       write(this%unit, "(A)") trim(symbol)
    else
       write(this%unit, "(A)") "dot"
    endif
    do i=1, xylist%n
       call xylist%get_element(i, xy)
       if(present(xunit))xy(1) = xy(1)/xunit
       if(present(yunit))xy(2) = xy(2)/yunit
       call this%write_coor(xy(1), xy(2))
    enddo
  end subroutine coop_asy_dots_list

  
  subroutine coop_asy_line_d(this, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy) this
    COOP_REAL ::xstart, ystart, xend, yend
    COOP_UNKNOWN_STRING, optional::color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING lineproperty
    write(this%unit, "(A)") "LINES"
    write(this%unit, "(A)") "1"
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
    write(this%unit, "(A)") trim(lineproperty)
    call this%write_coor(real(xstart, sp), real(ystart, sp), real(xend, sp), real(yend, sp))
  end subroutine coop_asy_line_d


  subroutine coop_asy_line_s(this, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy) this
    real::xstart, ystart, xend, yend
    COOP_UNKNOWN_STRING, optional::color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING lineproperty
    write(this%unit, "(A)") "LINES"
    write(this%unit, "(A)") "1"
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
    write(this%unit, "(A)") trim(lineproperty)
    call this%write_coor(xstart, ystart, xend, yend)
  end subroutine coop_asy_line_s
  
  subroutine coop_asy_lines_d(this, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy) this
    COOP_REAL ,dimension(:),intent(IN)::xstart, ystart, xend, yend
    COOP_UNKNOWN_STRING, optional::color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING lineproperty
    COOP_INT  n, i
    write(this%unit, "(A)") "LINES"
    n = coop_getdim("coop_asy_lines", size(xstart), size(ystart), size(xend), size(yend))
    write(this%unit, "(I8)") n
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
    write(this%unit, "(A)") trim(lineproperty)
    do i = 1, n
       call this%write_coor(real(xstart(i),sp), real(ystart(i),sp), real(xend(i), sp), real(yend(i), sp))
    enddo
  end subroutine coop_asy_lines_d


  subroutine coop_asy_lines_s(this, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy) this
    COOP_SINGLE ,dimension(:),intent(IN)::xstart, ystart, xend, yend
    COOP_UNKNOWN_STRING, optional::color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING lineproperty
    COOP_INT  n, i
    write(this%unit, "(A)") "LINES"
    n = coop_getdim("coop_asy_lines", size(xstart), size(ystart), size(xend), size(yend))
    write(this%unit, "(I8)") n
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
    write(this%unit, "(A)") trim(lineproperty)
    do i = 1, n
       call this%write_coor( xstart(i), ystart(i), xend(i), yend(i))
    enddo
  end subroutine coop_asy_lines_s

  subroutine coop_asy_add_legend(this,  legend,  color, linetype, linewidth, box)
    class(coop_asy)::this
    COOP_UNKNOWN_STRING::legend
    COOP_UNKNOWN_STRING, optional :: color, linetype
    COOP_SINGLE, optional::linewidth
    logical, optional::box
    COOP_STRING lineproperty
    if(present(box))then
       if(.not. box)then
          write(this%unit, "(A)") "LEGEND_NOBOX"
       else
          write(this%unit, "(A)") "LEGEND"          
       endif
    else
       write(this%unit, "(A)") "LEGEND"
    endif
    write(this%unit, "(A)") "VIRTUAL"
    if(trim(adjustl(legend)) .ne. "")then
       write(this%unit, "(A)") trim(adjustl(legend))
    else
       write(this%unit, "(A)") "NULL"
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
    write(this%unit, "(A)") trim(lineproperty)
    this%ymax = this%ymax + (this%ymax-this%ymin)*0.01
  end subroutine coop_asy_add_legend

  subroutine coop_asy_interpolate_curve_d(this, xraw, yraw, interpolate, color, linetype, linewidth, legend)
    COOP_INT ,parameter::n = 256
    COOP_REAL  x(n), y(n),  minx, maxx, dx
    COOP_INT  w(n)
    COOP_REAL ,dimension(:),allocatable::xx, yy, yy2
    class(coop_asy) this
    COOP_REAL ,dimension(:),intent(IN)::xraw, yraw
    COOP_UNKNOWN_STRING::interpolate
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
       select case(COOP_LOWER_STR(interpolate))
       case("linearlog", "loglog")
          y(1) = exp(sum(log(yraw))/m)
       case default
          y(1) = sum(yraw)/m
       end select
       goto 100
    endif
    npt = n
    do_raw = .false.
    select case(COOP_LOWER_STR(interpolate))
    case("linearlinear")
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
    case("linearlog")
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
    case("loglinear")
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
    case("loglog")
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

100 write(this%unit, "(A)") "CURVE"
    write(this%unit, "(I8)") npt
    if(present(legend))then
       if(trim(legend).ne."")then
          write(this%unit, "(A)") trim(legend)
       else
          write(this%unit, "(A)") "NULL"
       endif
    else
       write(this%unit, "(A)") "NULL"
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
    write(this%unit, "(A)") trim(lineproperty)
    write(this%unit, "(A)") "0"
    if(do_raw)then
       do i=1,npt
          call this%write_coor( real(xraw(i), sp), real(yraw(i), sp) )
       enddo
    else
       do i=1,npt
          call this%write_coor(real(x(i), sp), real(y(i), sp))
       enddo
    endif
  end subroutine coop_asy_interpolate_curve_d

  subroutine coop_asy_plot_likelihood_d(this, xraw, yraw, color, linetype, linewidth, legend, left_tail, right_tail)
    class(coop_asy) this
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
100 write(this%unit, "(A)") "CURVE"
    write(this%unit, "(I8)") npt
    if(present(legend))then
       if(trim(legend).ne."")then
          write(this%unit, "(A)") trim(legend)
       else
          write(this%unit, "(A)") "NULL"
       endif
    else
       write(this%unit, "(A)") "NULL"
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
    write(this%unit, "(A)") trim(lineproperty)
    write(this%unit, "(A)") "0"
    if(do_left_tail)then
       call this%write_coor(real(xraw(1)*2.-xraw(2), sp), 0._sp)
    endif
    do i=1,imax-1
       call this%write_coor( real(xraw(i), sp), real(yraw(i), sp) )
    enddo
    if(imax.gt.1) &
         call this%write_coor( real((xraw(i)+xraw(i-1))/2., sp), real(yraw(i)*0.75+yraw(i-1)*0.25, sp) )
    call this%write_coor( real(xraw(imax), sp), real(yraw(imax), sp) )
    if(imax .lt. m) &
         call this%write_coor( real((xraw(i)+xraw(i+1))/2., sp), real(yraw(i)*0.75+yraw(i+1)*0.25, sp) )
    do i=imax+1, m
       call this%write_coor( real(xraw(i), sp), real(yraw(i), sp) )
    enddo
    if(do_right_tail)then
       call this%write_coor(real(xraw(m)*2.-xraw(m-1), sp), 0._sp)
    endif
  end subroutine coop_asy_plot_likelihood_d

  subroutine coop_asy_interpolate_curve_s(this, xraw, yraw, interpolate, color, linetype, linewidth, legend)
    COOP_INT ,parameter::n = 256
    COOP_REAL  x(n), y(n),  minx, maxx, dx
    COOP_INT  w(n)
    COOP_REAL ,dimension(:),allocatable::xx, yy, yy2
    class(coop_asy) this
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
       select case(COOP_LOWER_STR(interpolate))
       case("linearlog", "loglog")
          y(1) = exp(sum(log(yraw))/m)
       case default
          y(1) = sum(yraw)/m
       end select
       goto 100
    endif
    npt = n
    do_raw = .false.
    select case(COOP_LOWER_STR(interpolate))
    case("linearlinear")
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
    case("linearlog")
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
    case("loglinear")
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
    case("loglog")
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

100 write(this%unit, "(A)") "CURVE"
    write(this%unit, "(I8)") npt
    if(present(legend))then
       if(trim(legend).ne."")then
          write(this%unit, "(A)") trim(legend)
       else
          write(this%unit, "(A)") "NULL"
       endif
    else
       write(this%unit, "(A)") "NULL"
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
    write(this%unit, "(A)") trim(lineproperty)
    write(this%unit, "(A)") "0"
    if(do_raw)then
       do i=1,npt
          call this%write_coor( xraw(i), yraw(i) )
       enddo
    else
       do i=1,npt
          call this%write_coor( real(x(i), sp), real(y(i), sp) )
       enddo
    endif
  end subroutine coop_asy_interpolate_curve_s



  subroutine coop_asy_plot_likelihood_s(this, xraw, yraw, color, linetype, linewidth, legend, left_tail, right_tail)
    class(coop_asy) this
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
100 write(this%unit, "(A)") "CURVE"
    write(this%unit, "(I8)") npt
    if(present(legend))then
       if(trim(legend).ne."")then
          write(this%unit, "(A)") trim(legend)
       else
          write(this%unit, "(A)") "NULL"
       endif
    else
       write(this%unit, "(A)") "NULL"
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
    write(this%unit, "(A)") trim(lineproperty)
    write(this%unit, "(A)") "0"
    if(do_left_tail) call this%write_coor(xraw(1)*2.-xraw(2), 0.)
    do i=1,imax-1
       call this%write_coor( xraw(i), yraw(i) )
    enddo
    if(imax.gt.1) &
         call this%write_coor( (xraw(i)+xraw(i-1))/2.,  yraw(i)*0.75+yraw(i-1)*0.25 )
    call this%write_coor(xraw(imax), yraw(imax))
    if(imax .lt. m) &
         call this%write_coor( (xraw(i)+xraw(i+1))/2., yraw(i)*0.75+yraw(i+1)*0.25 )
    do i=imax+1, m
       call this%write_coor( xraw(i), yraw(i) )
    enddo
    if(do_right_tail) call this%write_coor(xraw(m)*2.-xraw(m-1), 0.)
  end subroutine coop_asy_plot_likelihood_s



  subroutine coop_asy_plot_file(this, filename, interpolate, xcol, ycol, color, linetype, linewidth, legend, filename2)
    class(coop_asy) this
    COOP_INT ,parameter::n = 256
    COOP_UNKNOWN_STRING,optional::xcol, ycol
    COOP_REAL  x(n), y(n),  minx, maxx, dx
    COOP_INT  w(n)
    COOP_REAL ,dimension(:),allocatable::xx, yy, yy2
    type(coop_file) fp, fp2
    COOP_REAL ,dimension(:),allocatable::xraw, yraw, line
    COOP_UNKNOWN_STRING interpolate, filename
    COOP_UNKNOWN_STRING,optional::filename2
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_UNKNOWN_STRING, optional::legend
    COOP_INT  i, j, m, nn, npt, ncols, xl, yl, ncols2, ncols1
    COOP_STRING lineproperty
    logical do_raw, eval_x, eval_y
   
    if(.not. coop_file_exists(filename))then
       call coop_return_error("coop_asy_curve_from_file", "The data file "//trim(filename)//" does not exist", "stop")
    endif
    m = coop_file_numlines(filename)
    ncols1 = coop_file_numColumns(filename)
    if(ncols1 .eq. 0 .or. m.eq.0)then
       call coop_return_error("coop_asy_plot_file", "The data file "//trim(filename)//" is empty", "stop")
    endif
    if(present(filename2))then
       m = min(coop_file_numlines(filename2),m)
       ncols2  = coop_file_numColumns(filename2) 
       if(ncols2 .eq. 0 .or. m.eq.0)then
          call coop_return_error("coop_asy_plot_file", "The data file "//trim(filename2)//" is empty", "stop")
       endif
    else
       ncols2 = 0
    endif
    ncols = ncols1+ncols2
    eval_x = .false.
    eval_y = .false.
    if(present(ycol))then
       if(scan(ycol, "$").ne.0)then
          eval_y = .true.
          yl = min(ncols,2)
       elseif(trim(ycol).ne."")then
          read(ycol,*)yl
          if(yl .gt. ncols .or. yl.lt.0) stop "coop_asy_plot_file: ycol > ncol"
       else
          yl = min(ncols,2)
       endif
    else
       yl = min(ncols, 2)
    endif
    if(present(xcol))then
       if(scan(xcol, "$").ne.0)then
          eval_x = .true.
          xl = min(ncols-1,1)
       elseif(trim(xcol).ne."")then
          read(xcol, *) xl
          if(xl .gt. ncols .or. xl.lt.0) stop "coop_asy_plot_file: xcol > ncol"
       else
          xl = min(ncols-1,1)
       endif
    else
       xl = min(ncols-1,1)
    endif
    
    allocate(line(0:ncols))
    allocate(xraw(m), yraw(m))
    call fp%open(trim(filename), "r")
    if(present(filename2))then
       call fp2%open(filename2,"r")
       do i = 1, m
          if(fp%read_real_array(line(1:ncols1)) .and. fp2%read_real_array(line(ncols1 + 1:ncols)))then
             line(0) = i
             if(eval_x)then
                call coop_eval_math( xcol, xraw(i), vars = line(1:ncols))
             else
                xraw(i) = line(xl)
             endif
             if(eval_y)then
                call coop_eval_math(ycol, yraw(i), vars = line(1:ncols))
             else
                yraw(i) = line(yl)
             endif
          else
             call coop_return_error("coop_asy_plot_file", trim(filename)//" or "//trim(filename2)//" has a wrong format", "stop")
          endif
       enddo
       call fp2%close()
    else
       do i = 1, m
          if(fp%read_real_array(line(1:ncols)))then
             line(0) = i
             if(eval_x)then
                call coop_eval_math( xcol, xraw(i), vars = line(1:ncols))
             else
                xraw(i) = line(xl)
             endif
             if(eval_y)then
                call coop_eval_math(ycol, yraw(i), vars = line(1:ncols))
             else
                yraw(i) = line(yl)
             endif
          else
             call coop_return_error("coop_asy_plot_file", trim(filename)//" has a wrong format", "stop")
          endif
       enddo
    endif
    call fp%close()
    minx = minval(xraw)
    maxx = maxval(xraw)
    y = 0.
    w = 0
    if(minx .ge. maxx)then
       npt = 1
       x(1) = minx
       select case(COOP_LOWER_STR(interpolate))
       case("linearlog", "loglog")
          y(1) = exp(sum(log(yraw))/m)
       case default
          y(1) = sum(yraw)/m
       end select
       goto 100
    endif
    npt = n
    do_raw = .false.
    select case(COOP_LOWER_STR(interpolate))
    case("linearlinear")
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
    case("linearlog")
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
    case("loglinear")
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
    case("loglog")
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

100 write(this%unit, "(A)") "CURVE"
    write(this%unit, "(I8)") npt
    if(present(legend))then
       if(trim(legend).ne."")then
          write(this%unit, "(A)") trim(legend)
       else
          write(this%unit, "(A)") "NULL"
       endif
    else
       write(this%unit, "(A)") "NULL"
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
    write(this%unit, "(A)") trim(lineproperty)
    write(this%unit, "(A)") "0"
    if(do_raw)then
       do i=1,npt
          call this%write_coor( real(xraw(i), sp), real(yraw(i), sp))
       enddo
    else
       do i=1,npt
          call this%write_coor(real(x(i), sp), real(y(i), sp))
       enddo
    endif
    deallocate(xraw, yraw, line)
  end subroutine coop_asy_plot_file


  subroutine coop_asy_curve_d(this, x, y, marker, color, linetype, linewidth, legend)
    class(coop_asy) this
    COOP_UNKNOWN_STRING,optional::marker
    COOP_REAL ,dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_UNKNOWN_STRING, optional::legend
    COOP_INT  i,n
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_curve", size(x), size(y))
    write(this%unit, "(A)") "CURVE"
    write(this%unit, "(I8)") n
    if(present(legend))then
       if(trim(legend).ne."")then
          write(this%unit, "(A)") trim(legend)
       else
          write(this%unit, "(A)") "NULL"
       endif
    else
       write(this%unit, "(A)") "NULL"
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
    write(this%unit, "(A)") trim(lineproperty)
    if(present(marker))then
       if(trim(marker).eq.'')then
          write(this%unit, "(A)") "NULL" !!no marker
       else
          write(this%unit, "(A)") trim(marker)
       endif
    else
       write(this%unit, "(A)") "NULL" !!no marker
    endif
    do i=1,n
       call this%write_coor(real(x(i), sp), real(y(i), sp))
    enddo
  end subroutine coop_asy_curve_d

  subroutine coop_asy_curve_s(this, x, y, marker, color, linetype, linewidth, legend)
    class(coop_asy) this
    COOP_UNKNOWN_STRING,optional::marker
    COOP_SINGLE ,dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_INT  i,n
    COOP_UNKNOWN_STRING,optional::legend
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_curve", size(x), size(y))
    write(this%unit, "(A)") "CURVE"
    write(this%unit, "(I8)") n
    if(present(legend))then
       if(trim(legend).ne."")then
          write(this%unit, "(A)") trim(legend)
       else
          write(this%unit, "(A)") "NULL"
       endif
    else
       write(this%unit, "(A)") "NULL"
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
    write(this%unit, "(A)") trim(lineproperty)
    if(present(marker))then
       if(trim(marker) == "")then
          write(this%unit, "(A)") "NULL" 
       else
          write(this%unit, "(A)")  trim(marker)
       endif
    else
       write(this%unit, "(A)") "NULL" 
    endif
    do i=1,n
       call this%write_coor( x(i), y(i) )
    enddo
  end subroutine coop_asy_curve_s

  subroutine coop_asy_labels_d(this, labels, x, y, color, alignment)
    COOP_STRING, dimension(:),intent(IN)::labels
    COOP_REAL ,dimension(:),intent(IN)::x, y
    class(coop_asy) this
    COOP_UNKNOWN_STRING,optional::color
    COOP_UNKNOWN_STRING,optional::alignment
    COOP_INT  n, i
    if(present(alignment))then
       select case(trim(adjustl(alignment)))
       case("left", "LEFT", "Left", "l", "L")
          write(this%unit, "(A)") "LEFTLABELS"
       case("right", "RIGHT", "Right", "r", "R")
          write(this%unit, "(A)") "RIGHTLABELS"
       case default
          write(this%unit, "(A)") "LABELS"          
       end select
    else
       write(this%unit, "(A)") "LABELS"
    endif    
    n = coop_getdim("coop_asy_labels", size(x), size(y), size(labels))
    write(this%unit,"(I8)") n
    if(present(color))then
       write(this%unit, "(A)") trim(color)
    else
       write(this%unit, "(A)") "black"
    endif
    do i=1,n
       call this%write_coor(real(x(i), sp), real(y(i), sp))
       if(trim(labels(i)).eq."")then
          write(this%unit, "(A)") "NULLL"
       else
          write(this%unit, "(A)") trim(labels(i))
       endif
    enddo
  end subroutine coop_asy_labels_d

  subroutine coop_asy_labels_s(this, labels, x, y, color, alignment)
    COOP_STRING, dimension(:),intent(IN)::labels
    COOP_SINGLE ,dimension(:),intent(IN)::x, y
    class(coop_asy) this
    COOP_UNKNOWN_STRING,optional::color
    COOP_UNKNOWN_STRING,optional::alignment
    COOP_INT  n, i
    if(present(alignment))then
       select case(trim(alignment))
       case("left","LEFT","Left","l","L")
          write(this%unit, "(A)") "LEFTLABELS"
       case("right","RIGHT","Right","r","R")
          write(this%unit, "(A)") "RIGHTLABELS"
       case default
          write(this%unit, "(A)") "LABELS"          
       end select
    else
       write(this%unit, "(A)") "LABELS"
    endif
    n = coop_getdim("coop_asy_labels", size(x), size(y), size(labels))
    write(this%unit,"(I8)") n
    if(present(color))then
       write(this%unit, "(A)") trim(color)
    else
       write(this%unit, "(A)") "black"
    endif
    do i=1,n
       call this%write_coor( x(i), y(i))
       if(trim(labels(i)).eq."")then
          write(this%unit, "(A)") "NULLL"
       else
          write(this%unit, "(A)") trim(labels(i))
       endif
    enddo
  end subroutine coop_asy_labels_s


  subroutine coop_asy_label_d(this, label, x, y, color, alignment)
    class(coop_asy) this
    COOP_UNKNOWN_STRING label
    COOP_UNKNOWN_STRING,optional::color
    COOP_UNKNOWN_STRING,optional::alignment
    COOP_REAL  x, y
    if(present(alignment))then
       select case(trim(alignment))
       case("left","LEFT","Left","l","L")
          write(this%unit, "(A)") "LEFTLABELS"
       case("right","RIGHT","Right","r","R")
          write(this%unit, "(A)") "RIGHTLABELS"
       case default
          write(this%unit, "(A)") "LABELS"          
       end select
    else
       write(this%unit, "(A)") "LABELS"
    endif
    write(this%unit, "(A)") "1"
    if(present(color))then
       write(this%unit, "(A)") trim(color)
    else
       write(this%unit, "(A)") "black"
    endif
    call this%write_coor( real(x, sp), real(y, sp) )
    if(trim(label).eq."")then
       write(this%unit, "(A)") "NULL"
    else
       write(this%unit, "(A)") trim(label)
    endif
  end subroutine coop_asy_label_d

  subroutine coop_asy_label_relative(this, label, xratio, yratio, color, alignment)
    class(coop_asy) this
    COOP_UNKNOWN_STRING label
    COOP_UNKNOWN_STRING,optional::color
    COOP_UNKNOWN_STRING,optional::alignment
    COOP_SINGLE  xratio, yratio
    if(present(alignment))then
       select case(trim(alignment))
       case("left","LEFT","Left","l","L")
          write(this%unit, "(A)") "LEFTLABELS"
       case("right","RIGHT","Right","r","R")
          write(this%unit, "(A)") "RIGHTLABELS"
       case default
          write(this%unit, "(A)") "LABELS"          
       end select
    else
       write(this%unit, "(A)") "RIGHTLABELS"
    endif
    write(this%unit, "(A)") "1"
    if(present(color))then
       write(this%unit, "(A)") trim(color)
    else
       write(this%unit, "(A)") "black"
    endif
    call this%write_coor(this%xrel(xratio), this%yrel(yratio))
    if(trim(label).eq."")then
       write(this%unit, "(A)") "NULL"
    else
       write(this%unit, "(A)") trim(label)
    endif
  end subroutine coop_asy_label_relative

  
  subroutine coop_asy_label_s(this, label, x, y, color, alignment)
    class(coop_asy) this
    COOP_UNKNOWN_STRING label
    COOP_UNKNOWN_STRING,optional::color
    COOP_UNKNOWN_STRING,optional::alignment
    COOP_SINGLE  x, y
    if(present(alignment))then
       select case(trim(alignment))
       case("left","LEFT","Left","l","L")
          write(this%unit, "(A)") "LEFTLABELS"
       case("right","RIGHT","Right","r","R")
          write(this%unit, "(A)") "RIGHTLABELS"
       case default
          write(this%unit, "(A)") "LABELS"          
       end select
    else
       write(this%unit, "(A)") "LABELS"
    endif    
    write(this%unit, "(A)") "1"
    if(present(color))then
       write(this%unit, "(A)") trim(color)
    else
       write(this%unit, "(A)") "black"
    endif
    call this%write_coor( x, y)
    if(trim(label).eq."")then
       write(this%unit, "(A)") "NULL"
    else
       write(this%unit, "(A)") trim(label)
    endif
  end subroutine coop_asy_label_s



  subroutine coop_asy_contour_d(this, x, y, colorfill, smooth, linecolor, linetype, linewidth, legend)
    class(coop_asy) this
    logical,optional::smooth
    COOP_REAL ,dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype, legend
    COOP_SINGLE ,optional::linewidth
    COOP_INT  i,n
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_contour", size(x), size(y))
    write(this%unit, "(A)") "CONTOUR"
    write(this%unit, "(A)") "1" !!type 1 contour
    if(present(colorfill))then
       write(this%unit, "(A)") trim(colorfill)
    else
       write(this%unit, "(A)") "black"
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
       write(this%unit, "(A)") trim(lineproperty)
    else
       write(this%unit, "(A)") "NULL"
    endif
    if(present(legend))then
       if(trim(legend).eq.'')then
          write(this%unit, '(A)') "NULL"
       else
          write(this%unit, '(A)') trim(legend)
       endif
    else
       write(this%unit, '(A)') "NULL"
    endif
    if(present(smooth))then
       if(smooth)then
          write(this%unit, "(A)") "1"
       else
          write(this%unit, "(A)") "0"
       endif
    else
       write(this%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(this%unit, "(A)") "1"  !!1 path 
    write(this%unit, "(I8)") n
    do i=1,n
       call this%write_coor(real(x(i), sp), real(y(i), sp))
    enddo
  end subroutine coop_asy_contour_d

  subroutine coop_asy_contour_s(this, x, y, colorfill, smooth, linecolor, linetype, linewidth, legend)
    class(coop_asy) this
    logical,optional::smooth
    COOP_SINGLE ,dimension(:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype, legend
    COOP_SINGLE ,optional::linewidth
    COOP_INT  i,n
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_contour", size(x), size(y))
    write(this%unit, "(A)") "CONTOUR"
    write(this%unit, "(A)") "1"   !!type I contour
    if(present(colorfill))then
       write(this%unit, "(A)") trim(colorfill)
    else
       write(this%unit, "(A)") "black"
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
       write(this%unit, "(A)") trim(lineproperty)
    else
       write(this%unit, "(A)") "NULL"
    endif
    if(present(legend))then
       if(trim(legend).eq.'')then
          write(this%unit, '(A)') "NULL"
       else
          write(this%unit, '(A)') trim(legend)
       endif
    else
       write(this%unit, '(A)') "NULL"
    endif
    if(present(smooth))then
       if(smooth)then
          write(this%unit, "(A)") "1"
       else
          write(this%unit, "(A)") "0"
       endif
    else
       write(this%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(this%unit, "(A)") "1"  !!1 path 
    write(this%unit, "(I8)") n
    do i=1,n
       call this%write_coor( x(i), y(i) )
    enddo
  end subroutine coop_asy_contour_s


  subroutine coop_asy_contour_path(this, path, colorfill, smooth, linecolor, linetype, linewidth, legend)
    class(coop_asy) this
    logical,optional::smooth
    type(coop_asy_path) path
    COOP_SINGLE  xy(2)
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype, legend
    COOP_SINGLE ,optional::linewidth
    COOP_INT  i, ipath, pl
    COOP_STRING lineproperty

    if(path%nclosed.eq.0 .or. (.not. path%l%isinit()) .or. path%l%n .eq. 0)return
    if(path%l%dim .ne.2) stop "coop_asy_contour_path: wrong dimension of list in path"
    write(this%unit, "(A)") "CONTOUR"
    write(this%unit, "(A)") "1"   !!type I contour
    if(present(colorfill))then
       write(this%unit, "(A)") trim(colorfill)
    else
       write(this%unit, "(A)") "black"
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
       write(this%unit, "(A)") trim(lineproperty)
    else
       write(this%unit, "(A)") "NULL"
    endif
    if(present(legend))then
       if(trim(legend).eq.'')then
          write(this%unit, '(A)') "NULL"
       else
          write(this%unit, '(A)') trim(legend)
       endif
    else
       write(this%unit, '(A)') "NULL"
    endif

    if(present(smooth))then
       if(smooth)then
          write(this%unit, "(A)") "1"            
       else
          write(this%unit, "(A)") "0"            
       endif
    else
       write(this%unit, "(A)") "0"  
    endif
    write(this%unit, "(I8)") path%nclosed
    i = 1
    do ipath = 1, path%nclosed
       write(this%unit, "(I8)") path%length(ipath)
       do pl = 1, path%length(ipath)
          call path%l%get_element(i, xy)
          call this%write_coor(xy(1), xy(2))
          i = i + 1
       enddo
    enddo
  end subroutine coop_asy_contour_path

  subroutine coop_asy_contour_mult_d(this, x, y, colorfill, smooth, linecolor, linetype, linewidth, legend)
    class(coop_asy) this
    logical,optional::smooth
    COOP_REAL ,dimension(:,:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype, legend
    COOP_SINGLE ,optional::linewidth
    COOP_INT  i,n, m,ipath
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_contour", size(x,1), size(y,1))
    m = coop_getdim("coop_asy_contour", size(x,2), size(y,2))
    write(this%unit, "(A)") "CONTOUR"
    write(this%unit, "(A)") "1"   !!type I contour
    if(present(colorfill))then
       write(this%unit, "(A)") trim(colorfill)
    else
       write(this%unit, "(A)") "black"
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
       write(this%unit, "(A)") trim(lineproperty)
    else
       write(this%unit, "(A)") "NULL"
    endif
    if(present(legend))then
       if(trim(legend).eq.'')then
          write(this%unit, '(A)') "NULL"
       else
          write(this%unit, '(A)') trim(legend)
       endif
    else
       write(this%unit, '(A)') "NULL"
    endif
    if(present(smooth))then
       if(smooth)then
          write(this%unit, "(A)") "1"
       else
          write(this%unit, "(A)") "0"
       endif
    else
       write(this%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(this%unit, "(I8)") m
    do ipath = 1, m
       write(this%unit, "(I8)") n
       do i=1,n
          call this%write_coor(real(x(i, ipath), sp), real(y(i, ipath), sp))
       enddo
    enddo
  end subroutine coop_asy_contour_mult_d

  subroutine coop_asy_contour_mult_s(this, x, y, colorfill, smooth, linecolor, linetype, linewidth, legend)
    class(coop_asy) this
    logical,optional::smooth
    COOP_SINGLE ,dimension(:,:),intent(IN)::x,y
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype, legend
    COOP_SINGLE ,optional::linewidth
    COOP_INT  i,n, m,ipath
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_contour", size(x,1), size(y,1))
    m = coop_getdim("coop_asy_contour", size(x,2), size(y,2))
    write(this%unit, "(A)") "CONTOUR"
    write(this%unit, "(A)") "1"   !!type I contour
    if(present(colorfill))then
       write(this%unit, "(A)") trim(colorfill)
    else
       write(this%unit, "(A)") "black"
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
       write(this%unit, "(A)") trim(lineproperty)
    else
       write(this%unit, "(A)") "NULL"
    endif
    if(present(legend))then
       if(trim(legend).eq.'')then
          write(this%unit, '(A)') "NULL"
       else
          write(this%unit, '(A)') trim(legend)
       endif
    else
       write(this%unit, '(A)') "NULL"
    endif

    if(present(smooth))then
       if(smooth)then
          write(this%unit, "(A)") "1"
       else
          write(this%unit, "(A)") "0"
       endif
    else
       write(this%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(this%unit, "(I8)") m
    do ipath = 1, m
       write(this%unit, "(I8)") n
       do i=1,n
          call this%write_coor( x(i, ipath), y(i, ipath) )
       enddo
    enddo
  end subroutine coop_asy_contour_mult_s


  subroutine coop_asy_contour_arr_d(this, z, xmin, xmax, ymin, ymax, cvals, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) this
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
    write(this%unit, "(A)") "CONTOUR"
    write(this%unit, "(A)") "2" !!type 2 contour
    call this%write_limits( real(xmin, sp), real(xmax, sp), real(ymin, sp), real(ymax, sp))
    if(present(linewidth))then
       nc = coop_getdim("coop_asy_contour_arr_d", size(cvals), size(colorfill), size(linecolor), size(linetype), size(linewidth))
    else
       nc = coop_getdim("coop_asy_contour_arr_d", size(cvals), size(colorfill), size(linecolor), size(linetype))
    endif
    write(this%unit, "(I8)") nc
    write(this%unit, "("//trim(coop_num2str(nc))//"G16.7)") cvals
    do i = 1, nc
       write(this%unit, "(A)") trim(colorfill(i))
       if(present(linewidth))then
          write(this%unit, "(A)") trim(coop_asy_linestyle(linecolor(i), linetype(i), linewidth(i)))
       else
          write(this%unit, "(A)") trim(coop_asy_linestyle(linecolor(i), linetype(i)))
       endif
    enddo
    if(smooth)then
       write(this%unit, "(A)") "1"
    else
       write(this%unit, "(A)") "0"
    endif
    write(this%unit, "(2I8)") nx, ny
    do i=1,nx
       write(this%unit, "("//trim(coop_num2str(ny))//"G14.5)") real(z(i,:))
    enddo
  end subroutine coop_asy_contour_arr_d
  


  subroutine coop_asy_contour_arr_s(this, z, xmin, xmax, ymin, ymax, cvals, colorfill, smooth, linecolor, linetype, linewidth)
    class(coop_asy) this
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
    write(this%unit, "(A)") "CONTOUR"
    write(this%unit, "(A)") "2" !!type 2 contour
    call this%write_limits(xmin, xmax, ymin, ymax)
    if(present(linewidth))then
       nc = coop_getdim("coop_asy_contour_arr_s", size(cvals), size(colorfill), size(linecolor), size(linetype), size(linewidth))
    else
       nc = coop_getdim("coop_asy_contour_arr_s", size(cvals), size(colorfill), size(linecolor), size(linetype))
    endif
    write(this%unit, "(I8)") nc
    write(this%unit, "("//trim(coop_num2str(nc))//"G16.7)") cvals
    do i = 1, nc
       write(this%unit, "(A)") trim(colorfill(i))
       if(present(linewidth))then
          write(this%unit, "(A)") trim(coop_asy_linestyle(linecolor(i), linetype(i), linewidth(i)))
       else
          write(this%unit, "(A)") trim(coop_asy_linestyle(linecolor(i), linetype(i)))
       endif
    enddo
    if(smooth)then
       write(this%unit, "(A)") "1"
    else
       write(this%unit, "(A)") "0"
    endif
    write(this%unit, "(2I8)") nx, ny
    do i=1,nx
       write(this%unit, "("//trim(coop_num2str(ny))//"G14.5)") z(i,:)
    enddo
  end subroutine coop_asy_contour_arr_s


  subroutine coop_asy_density_d(this, z, xmin, xmax, ymin, ymax, label, zmin, zmax, color_table)
    class(coop_asy) this
    COOP_REAL ,dimension(:,:)::z
    COOP_REAL  xmin, xmax, ymin, ymax
    COOP_REAL ,optional::zmin, zmax
    COOP_UNKNOWN_STRING, optional::label
    COOP_REAL  minz, maxz
    COOP_INT  nx, ny, i
    COOP_UNKNOWN_STRING, optional::color_table
    nx = size(z, 1)
    ny = size(z, 2)
    write(this%unit, "(A)") "DENSITY"
    if(present(color_table))then
       if(trim(color_table).ne."")then
          write(this%unit, "(A)") trim(color_table)
       else
          write(this%unit, "(A)") "Rainbow"
       endif
    else
       write(this%unit, "(A)") "Rainbow"
    endif

    if(present(label))then
       if(trim(label).ne."")then
          write(this%unit,"(A)") trim(label)
       else
          write(this%unit,"(A)") "NULL"
       endif
    else
       write(this%unit,"(A)") "NULL"
    endif
    call this%write_limits(real(xmin, sp), real(xmax, sp), real(ymin, sp), real(ymax, sp))
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

    write(this%unit, "(2G16.7)") minz, maxz
    write(this%unit, "(A)") "0"
    write(this%unit, "(2I8)") nx, ny
    do i=1,nx
       write(this%unit, "("//trim(coop_num2str(ny))//"G14.5)") real(z(i,:))
    enddo
  end subroutine coop_asy_density_d

  subroutine coop_asy_density_s(this, z, xmin, xmax, ymin, ymax, label, zmin, zmax, color_table)
    class(coop_asy) this
    COOP_SINGLE ,dimension(:,:)::z
    COOP_SINGLE  xmin, xmax, ymin, ymax
    COOP_SINGLE ,optional::zmin, zmax
    COOP_SINGLE  minz, maxz
    COOP_INT  nx, ny, i
    COOP_UNKNOWN_STRING, optional::label, color_table
    nx = size(z, 1)
    ny = size(z, 2)
    write(this%unit, "(A)") "DENSITY"
    if(present(color_table))then
       if(trim(color_table).ne."")then
          write(this%unit,"(A)") trim(color_table)
       else
          write(this%unit,"(A)") "Rainbow"
       endif
    else
       write(this%unit,"(A)") "Rainbow"
    endif
    if(present(label))then
       if(trim(label).ne."")then
          write(this%unit,"(A)") trim(label)
       else
          write(this%unit,"(A)") "NULL"
       endif
    else
       write(this%unit,"(A)") "NULL"
    endif
    call this%write_limits(xmin, xmax, ymin, ymax)
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
    write(this%unit, "(2G16.7)") minz, maxz
    write(this%unit, "(A)") "0"
    write(this%unit, "(2I8)") nx, ny
    do i=1,nx
       write(this%unit, "("//trim(coop_num2str(ny))//"G14.5)") z(i,:)
    enddo
  end subroutine coop_asy_density_s


 subroutine coop_asy_irregular_density_d(this, x, y, z, label, xmin, xmax, ymin, ymax, zmin, zmax, color_table)
    class(coop_asy) this
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
    write(this%unit,"(A)") "DENSITY"
    if(present(color_table))then
       if(trim(color_table).ne."")then
          write(this%unit,"(A)") trim(color_table)
       else
          write(this%unit,"(A)") "Rainbow"
       endif
    else
       write(this%unit,"(A)") "Rainbow"
    endif
    if(present(label))then
       if(trim(label).ne."")then
          write(this%unit,"(A)") trim(label)
       else
          write(this%unit,"(A)") "NULL"
       endif
    else
       write(this%unit,"(A)") "NULL"
    endif
    call this%write_limits(real(minx, sp), real(maxx, sp), real(miny, sp), real(maxy, sp))
    write(this%unit, "(2G16.7)") minz, maxz
    write(this%unit,"(A)") "1"
    write(this%unit,"(I10)") n
    do i=1,n
       write(this%unit, "(3G14.5)") real(x(i)), real(y(i)), real(z(i))
    enddo
  end subroutine coop_asy_irregular_density_d

  subroutine coop_asy_irregular_density_s(this, x, y, z, label, xmin, xmax, ymin, ymax, zmin, zmax, color_table)
    class(coop_asy) this
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
    write(this%unit,"(A)") "DENSITY"

    if(present(color_table))then
       if(trim(color_table).ne."")then
          write(this%unit,"(A)") trim(color_table)
       else
          write(this%unit,"(A)") "Rainbow"
       endif
    else
       write(this%unit,"(A)") "Rainbow"
    endif

    if(present(label))then
       if(trim(label).ne."")then
          write(this%unit,"(A)") trim(label)
       else
          write(this%unit,"(A)") "NULL"
       endif
    else
       write(this%unit,"(A)") "NULL"
    endif
    call this%write_limits(minx, maxx, miny, maxy) !real(minx, sp), real(maxx, sp), real(miny, sp), real(maxy, sp))
    write(this%unit, "(2G16.7)") minz, maxz
    write(this%unit,"(A)") "1"
    write(this%unit,"(I10)") n
    do i=1,n
       write(this%unit, "(3G14.5)") x(i), y(i), z(i)
    enddo
  end subroutine coop_asy_irregular_density_s

  

  subroutine coop_asy_clip_d(this, x, y, smooth)
    class(coop_asy) this
    logical,optional::smooth
    COOP_REAL ,dimension(:),intent(IN)::x,y
    COOP_INT  i,n
    n = coop_getdim("coop_asy_clip", size(x), size(y))
    write(this%unit, "(A)") "CLIP"
    write(this%unit, "(I8)") n
    if(present(smooth))then
       if(smooth)then
          write(this%unit, "(A)") "1"
       else
          write(this%unit, "(A)") "0"
       endif
    else
       write(this%unit, "(A)") "0"  !!no smoothing by default
    endif
    do i=1,n
       call this%write_coor(real(x(i), sp), real(y(i), sp))
    enddo

  end subroutine coop_asy_clip_d

  subroutine coop_asy_clip_s(this, x, y, smooth)
    class(coop_asy) this
    logical,optional::smooth
    COOP_SINGLE ,dimension(:),intent(IN)::x,y
    COOP_INT  i,n
    n = coop_getdim("coop_asy_clip", size(x), size(y))
    write(this%unit, "(A)") "CLIP"
    write(this%unit, "(I8)") n
    if(present(smooth))then
       if(smooth)then
          write(this%unit, "(A)") "1"
       else
          write(this%unit, "(A)") "0"
       endif
    else
       write(this%unit, "(A)") "0"  !!no smoothing by default
    endif
    do i=1,n
       call this%write_coor( x(i), y(i))
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

  subroutine coop_asy_path_append_s(this, x, y)
    COOP_SINGLE  x, y
    class(coop_asy_path) this
    call this%l%push( (/ x, y /) )
  end subroutine coop_asy_path_append_s

  subroutine coop_asy_path_append_d(this, x, y)
    COOP_REAL  x, y
    class(coop_asy_path) this
    call coop_asy_path_append_s(this, real(x), real(y))
  end subroutine coop_asy_path_append_d

  subroutine coop_asy_path_close(path)
    class(coop_asy_path) path
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

  subroutine coop_asy_path_free(path)
    class(coop_asy_path) path
    if(allocated(path%length))deallocate(path%length)
    call path%l%init()
    path%nclosed = 0
  end subroutine coop_asy_path_free

  !!generate a path that encloses the region of f(x, y) > threshold
  subroutine coop_asy_path_from_function(this, f, xmin, xmax, ymin, ymax, threshold, resolution, countcut, peakcut, rmscut)
    COOP_INT ,parameter::default_resolution = 64
    COOP_INT ,optional:: resolution
    COOP_SINGLE, optional::countcut, peakcut, rmscut
    class(coop_asy_path) this
    external f
    COOP_SINGLE  f
    COOP_SINGLE  xmin, xmax, ymin, ymax, threshold, dx, dy, dxby2, dyby2
    COOP_INT  n
    logical,dimension(:,:),allocatable::above
    COOP_SINGLE, dimension(:,:),allocatable::image
    COOP_INT ,dimension(:,:,:),allocatable::lines
    COOP_INT  i, j, imin, jmin, last
    call this%free()
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
    if(present(countcut) .or. present(peakcut) .or. present(rmscut))then
       allocate(image(n, n))
       do j=1, n
          do i=1,n
             image(i, j) = f(xmin+dx*(i-0.5), ymin+dy*(j-0.5)) 
          enddo
       enddo
       image = image - sum(image)/n**2       
       if(present(countcut))then
          call coop_array_get_threshold(image, countcut, threshold)
       elseif(present(rmscut))then
          threshold = sqrt(sum(image**2)/n**2)*rmscut
       else
          threshold = maxval(image)*peakcut
       endif
       do j = 1, n
          do i=1,n
             above(i, j) = (image(i, j).gt. threshold)
          enddo
       enddo
       deallocate(image)
    else
       do j = 1, n
          do i=1,n
             above(i, j) = (f(xmin+dx*(i-0.5), ymin+dy*(j-0.5)) .gt. threshold)
          enddo
       enddo
    endif
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
          call this%close()
          goto 100 
       else
          call insert_point()
       endif
    enddo

  contains

    subroutine insert_point()
      COOP_SINGLE  x, y
      COOP_SINGLE  t(2), s
      x =  xmin+i*dx
      y = ymin+j*dy
      if(i.eq.0 .or. i.eq.n)then
         t(1) = 0 
      else
         t(1) = (f(x+dxby2, y) - f(x-dxby2, y))
      endif
      if(j.eq.0 .or. j.eq.n)then
         t(2) = 0
      else
         t(2) =  (f(x, y+dyby2) - f(x, y-dyby2))
      endif
      s= sum(t**2)
      if(s.gt.0)then
         s = (threshold - f(x, y))/s
         x = min(xmax, max(x + max(min(t(1)*s, 0.5), -0.5)*dx, xmin))
         y = min(ymax, max(y + max(min(t(2)*s, 0.5), -0.5)*dy, ymin))
      endif
      call this%append(x, y)
    end subroutine insert_point

  end subroutine coop_asy_path_from_function


  subroutine coop_asy_path_from_array_center(this, f, xmin, xmax, ymin, ymax, threshold)
    class(coop_asy_path) this
    COOP_SINGLE  f(:,:)
    COOP_SINGLE  xmin, xmax, ymin, ymax, threshold, dx, dy
    COOP_INT  nx, ny, n
    nx = size(f, 1)
    ny = size(f, 2)
    n = min(max(nx, ny,16)*2, 256) !!higher resolution is useless since we do bilinear interpolation
    dx = (xmax-xmin)/nx
    dy = (ymax- ymin)/ny
    call this%from_function(finterp, xmin, xmax, ymin, ymax, threshold, n)

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



  subroutine coop_asy_path_from_array(this, f, xmin, xmax, ymin, ymax, threshold)
    class(coop_asy_path) this
    COOP_SINGLE  f(:,:)
    COOP_SINGLE  xmin, xmax, ymin, ymax, threshold, dx, dy
    COOP_INT  nx, ny, n
    nx = size(f, 1)
    ny = size(f, 2)
    n = min(max(nx, ny,16)*2, 256) !!higher resolution is useless since we do bilinear interpolation
    dx = (xmax-xmin)/(nx-1)
    dy = (ymax- ymin)/(ny-1)
    call this%from_function(finterp, xmin, xmax, ymin, ymax, threshold, n)
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

  subroutine coop_asy_path_from_array_gaussianfit(this, f, xmin, xmax, ymin, ymax, threshold)
    COOP_REAL, parameter::ss = 2.d0*(1.d0)**2
    COOP_REAL, parameter::n00 = 1.d0
    COOP_REAL, parameter::n10 = exp(-1.d0/ss)
    COOP_REAL, parameter::n11 = exp(-2.d0/ss)
    COOP_REAL, parameter::n20 = exp(-4.d0/ss)
    COOP_REAL, parameter::n21 = exp(-5.d0/ss)
    COOP_REAL, parameter::n22 = exp(-8.d0/ss)
    COOP_REAL, parameter::nsum = n00 + (n10+n11+n20+n22)*4.d0 + n21*8.d0
    COOP_REAL, parameter::norm00 = n00/nsum
    COOP_REAL, parameter::norm10 = n10/nsum
    COOP_REAL, parameter::norm11 = n11/nsum
    COOP_REAL, parameter::norm20 = n20/nsum
    COOP_REAL, parameter::norm21 = n21/nsum
    COOP_REAL, parameter::norm22 = n22/nsum        
    class(coop_asy_path) this
    COOP_SINGLE  f(:,:)
    COOP_SINGLE  xmin, xmax, ymin, ymax, threshold, dx, dy, cov(2,2), invcov(2,2), mean(2), wtot, delta, fbyg_raw(-1:size(f,1)+2, -1:size(f,2)+2), vec(2), fbyg(size(f,1), size(f, 2))
    COOP_INT  nx, ny, n, ix, iy
    nx = size(f, 1)
    ny = size(f, 2)
    n = min(max(nx, ny,16)*6, 256)
    dx = (xmax-xmin)/(nx-1)
    dy = (ymax- ymin)/(ny-1)
    wtot = sum(f)
    mean = 0.d0    
    do ix = 1, nx
       mean(1) = mean(1) + sum(f(ix, :))*(ix-1)*dx
    enddo
    do iy=1, ny
       mean(2) = mean(2) + sum(f(:, iy))*(iy-1)*dy
    enddo
    mean(1) = mean(1)/wtot
    mean(2) = mean(2)/wtot
    cov = 0.d0
    do ix= 1, nx
       do iy = 1, ny
          cov(1, 1) = cov(1, 1) + ((ix-1)*dx-mean(1))**2*f(ix, iy)
          cov(2, 2) = cov(2, 2) + ((iy-1)*dy-mean(2))**2*f(ix, iy)
          cov(1, 2) = cov(1, 2) + ((ix-1)*dx-mean(1))*((iy-1)*dy-mean(2))*f(ix, iy)
       enddo
    enddo
    cov = cov/wtot
    delta = cov(1,1)*cov(2,2)-cov(1,2)**2
    invcov(1,1) = cov(2,2)/delta
    invcov(2,2) = cov(1,1)/delta
    invcov(1,2) = -cov(1,2)/delta
    invcov(2,1) = invcov(1,2)
    do ix=1, nx
       do iy = 1, ny
          vec = (/ (ix-1)*dx, (iy-1)*dy /) - mean
          fbyg_raw(ix, iy) = f(ix, iy)*exp(0.5d0*dot_product(vec, matmul(invcov, vec)))
       enddo
    enddo
    fbyg_raw(0, :) = fbyg_raw(1, :)
    fbyg_raw(-1, :) = fbyg_raw(1, :)    
    fbyg_raw(nx+1, :) = fbyg_raw(nx, :)
    fbyg_raw(nx+2, :) = fbyg_raw(nx, :)    
    fbyg_raw(:, 0) = fbyg_raw(:, 1)
    fbyg_raw(:, -1) = fbyg_raw(:, 1)    
    fbyg_raw(:, ny+1) = fbyg_raw(:, ny)
    fbyg_raw(:, ny+2) = fbyg_raw(:, ny)    
    
    fbyg = fbyg_raw(1:nx, 1:ny)*norm00 &
         + (fbyg_raw(0:nx-1, 1:ny) + fbyg_raw(2:nx+1, 1:ny) + fbyg_raw(1:nx, 0:ny-1) + fbyg_raw(1:nx, 2:ny+1))*norm10 &
         + (fbyg_raw(0:nx-1, 0:ny-1) + fbyg_raw(2:nx+1, 2:ny+1) + fbyg_raw(0:nx-1, 2:ny+1) + fbyg_raw(2:nx+1, 0:ny-1))*norm11 &
         + (fbyg_raw(-1:nx-2, 1:ny) + fbyg_raw(3:nx+2, 1:ny) + fbyg_raw(1:nx, 3:ny+2) + fbyg_raw(1:nx, -1:ny-2))*norm20 &
         + (fbyg_raw(-1:nx-2, 2:ny+1) + fbyg_raw(-1:nx-2, 0:ny-1) + fbyg_raw(3:nx+2, 2:ny+1) + fbyg_raw(3:nx+2, 0:ny-1) + fbyg_raw(0:nx-1, 3:ny+2) + fbyg_raw(2:nx+1, 3:ny+2) + fbyg_raw(2:nx+1, -1:ny-2) + fbyg_raw(0:nx-1, -1:ny-2) )*norm21 &
         + (fbyg_raw(-1:nx-2, -1:ny-2) + fbyg_raw(3:nx+2, -1:ny-2) + fbyg_raw(-1:nx-2, 3:ny+2) + fbyg_raw(3:nx+2, 3:ny+2))*norm22                          
    
    mean(1) = mean(1) + xmin
    mean(2) = mean(2) + ymin
    call this%from_function(finterp, xmin, xmax, ymin, ymax, threshold, n)    
  contains
    
    function  finterp(x, y)
      COOP_SINGLE  x, y, finterp
      COOP_INT  i, j
      COOP_SINGLE  ri, rj, vec(2)
      ri = (x-xmin)/dx + 1.
      rj = (y-ymin)/dy + 1.
      i = min(max(floor(ri), 1), nx-1)
      j = min(max(floor(rj), 1), ny-1)
      ri = ri - i
      rj = rj - j
      vec = (/ x, y /) - mean
      finterp = ((fbyg(i, j)*(1.-ri) + fbyg(i+1, j)*ri)*(1.-rj)+(fbyg(i, j+1)*(1.-ri) + fbyg(i+1, j+1)*ri)*rj)*exp(-0.5d0*dot_product(vec, matmul(invcov, vec)))
    end function finterp
    
  end subroutine coop_asy_path_from_array_gaussianfit


  subroutine coop_asy_legend_default(this)
    class(coop_asy) this
    call coop_asy_legend_location(this, "N")
  end subroutine coop_asy_legend_default

  subroutine coop_asy_legend_location(this, loc, cols, box)
    class(coop_asy) this
    COOP_UNKNOWN_STRING :: loc
    COOP_INT , optional::cols
    logical,optional::box
    if(present(box))then
       if(.not. box)then
          write(this%unit, "(A)") "LEGEND_NOBOX"
       else
          write(this%unit, "(A)") "LEGEND"          
       endif
    else
       write(this%unit, "(A)") "LEGEND"
    endif
    if(trim(adjustl(loc)).ne."")then
       write(this%unit, "(A)") trim(adjustl(loc))
    else
       write(this%unit, "(A)") "N"
    end if
    if(present(cols))then
       write(this%unit, "(I8)")  cols
    else
       write(this%unit, "(I8)")  1
    endif
  end subroutine coop_asy_legend_location

  subroutine coop_asy_legend_d(this, x, y, cols, box)
    class(coop_asy) this
    COOP_REAL x, y
    COOP_INT ,optional::cols
    logical,optional::box
    if(present(box))then
       if(.not. box)then
          write(this%unit, "(A)") "LEGEND_NOBOX"
       else
          write(this%unit, "(A)") "LEGEND"          
       endif
    else
       write(this%unit, "(A)") "LEGEND"
    endif
    write(this%unit, "(A)") "NULL"
    write(this%unit, "(2G15.4)") x, y
    if(present(cols))then
       write(this%unit, "(I8)")  cols
    else
       write(this%unit, "(I8)")  1
    endif
  end subroutine coop_asy_legend_d


  subroutine coop_asy_legend_advance_d(this, x, y, box_color, xmargin, ymargin, linelength, hskip, vskip, cols)
    class(coop_asy) this
    COOP_REAL::x, y
    COOP_INT::cols
    COOP_UNKNOWN_STRING::box_color
    COOP_SINGLE::xmargin, ymargin, linelength, hskip, vskip
    write(this%unit, "(A)") "LEGEND_ADVANCE"
    write(this%unit, "(A)") trim(adjustl(box_color))
    write(this%unit, "(G15.4)") xmargin
    write(this%unit, "(G15.4)") ymargin
    write(this%unit, "(G15.4)") linelength
    write(this%unit, "(G15.4)") hskip
    write(this%unit, "(G15.4)") vskip    
    write(this%unit, "(I6)")  cols    
    write(this%unit, "(A)") "NULL"
    write(this%unit, "(2G15.4)") x, y
  end subroutine coop_asy_legend_advance_d


  subroutine coop_asy_legend_advance_s(this, x, y, box_color, xmargin, ymargin, linelength, hskip, vskip, cols)
    class(coop_asy) this
    COOP_SINGLE::x, y
    COOP_INT::cols
    COOP_UNKNOWN_STRING::box_color
    COOP_SINGLE::xmargin, ymargin, linelength, hskip, vskip
    write(this%unit, "(A)") "LEGEND_ADVANCE"
    write(this%unit, "(A)") trim(adjustl(box_color))
    write(this%unit, "(G15.4)") xmargin
    write(this%unit, "(G15.4)") ymargin
    write(this%unit, "(G15.4)") linelength
    write(this%unit, "(G15.4)") hskip
    write(this%unit, "(G15.4)") vskip    
    write(this%unit, "(I6)")  cols    
    write(this%unit, "(A)") "NULL"
    write(this%unit, "(2G15.4)") x, y
  end subroutine coop_asy_legend_advance_s


  subroutine coop_asy_legend_advance_relative(this, xratio, yratio, box_color, xmargin, ymargin, linelength, hskip, vskip, cols)
    class(coop_asy) this
    COOP_SINGLE::xratio, yratio
    COOP_INT::cols
    COOP_UNKNOWN_STRING::box_color
    COOP_SINGLE::xmargin, ymargin, linelength, hskip, vskip
    write(this%unit, "(A)") "LEGEND_ADVANCE"
    write(this%unit, "(A)") trim(adjustl(box_color))
    write(this%unit, "(G15.4)") xmargin
    write(this%unit, "(G15.4)") ymargin
    write(this%unit, "(G15.4)") linelength
    write(this%unit, "(G15.4)") hskip
    write(this%unit, "(G15.4)") vskip    
    write(this%unit, "(I6)")  cols    
    write(this%unit, "(A)") "NULL"
    write(this%unit, "(2G15.4)") this%xrel(xratio), this%yrel(yratio)
  end subroutine coop_asy_legend_advance_relative
  
  
  subroutine coop_asy_legend_advance_location(this, loc, box_color, xmargin, ymargin, linelength, hskip, vskip, cols)
    class(coop_asy) this
    COOP_UNKNOWN_STRING::loc
    COOP_INT::cols
    COOP_UNKNOWN_STRING::box_color
    COOP_SINGLE::xmargin, ymargin, linelength, hskip, vskip
    write(this%unit, "(A)") "LEGEND_ADVANCE"
    write(this%unit, "(A)") trim(adjustl(box_color))
    write(this%unit, "(G15.4)") xmargin
    write(this%unit, "(G15.4)") ymargin
    write(this%unit, "(G15.4)") linelength
    write(this%unit, "(G15.4)") hskip
    write(this%unit, "(G15.4)") vskip    
    write(this%unit, "(I6)")  cols
    if(trim(adjustl(loc)).ne."")then
       write(this%unit, "(A)") trim(adjustl(loc))
    else
       write(this%unit, "(A)") "N"
    endif
  end subroutine coop_asy_legend_advance_location
  

  subroutine coop_asy_legend_s(this, x, y, cols, box)
    class(coop_asy) this
    COOP_SINGLE  x, y
    COOP_INT ,optional::cols
    logical,optional::box
    if(present(box))then
       if(.not. box)then
          write(this%unit, "(A)") "LEGEND_NOBOX"
       else
          write(this%unit, "(A)") "LEGEND"          
       endif
    else
       write(this%unit, "(A)") "LEGEND"
    endif
    write(this%unit, "(A)") "NULL"
    write(this%unit, "(2G15.4)") x, y
    if(present(cols))then
       write(this%unit, "(I8)")  cols
    else
       write(this%unit, "(I8)")  1
    endif
  end subroutine coop_asy_legend_s

  subroutine coop_asy_legend_relative(this, xratio, yratio, cols, box)
    class(coop_asy) this
    COOP_SINGLE  xratio, yratio
    COOP_INT ,optional::cols
    logical,optional::box
    write(this%unit, "(A)") "LEGEND_ADVANCE"
    if(present(box))then
       if(box)then
          write(this%unit, "(A)") "black"          
       else
          write(this%unit, "(A)") "invisible"
       endif
    else
       write(this%unit, "(A)") "black"       
    endif
    write(this%unit, "(A)") "0.5"
    write(this%unit, "(A)") "0.3"
    write(this%unit, "(A)") "1"
    write(this%unit, "(A)") "0.95"
    write(this%unit, "(A)") "0.95"
    if(present(cols))then
       write(this%unit, "(I6)")  cols
    else
       write(this%unit, "(A)")  "1"       
    endif
    write(this%unit, "(A)") "NULL"
    write(this%unit, "(2G15.4)") this%xrel(xratio), this%yrel(yratio)
  end subroutine coop_asy_legend_relative
  

  subroutine coop_asy_compact_legend_relative(this, xratio, yratio, cols)
    class(coop_asy) this
    COOP_SINGLE  xratio, yratio
    COOP_INT ,optional::cols
    write(this%unit, "(A)") "LEGEND_ADVANCE"
    write(this%unit, "(A)") "invisible"
    write(this%unit, "(A)") "0.05"
    write(this%unit, "(A)") "0.05"
    write(this%unit, "(A)") "0.75"
    write(this%unit, "(A)") "0.85"
    write(this%unit, "(A)") "0.85"
    if(present(cols))then
       write(this%unit, "(I6)")  cols
    else
       write(this%unit, "(A)")  "1"
    endif
    write(this%unit, "(A)") "NULL"
    write(this%unit, "(2G15.4)") this%xrel(xratio), this%yrel(yratio)
  end subroutine coop_asy_compact_legend_relative



  subroutine coop_asy_topaxis_s(this, xmin, xmax, islog,label)
    COOP_SINGLE  xmin, xmax
    logical islog
    class(coop_asy) this
    COOP_UNKNOWN_STRING label
    write(this%unit, "(A)")"EXTRA_AXIS"
    write(this%unit, "(A)") "top"
    if(trim(adjustl(label)).ne."")then
       write(this%unit, "(A)") trim(adjustl(label))
    else
       write(this%unit, "(A)") "NULL"
    endif
    if(islog)then
       write(this%unit, "(A)") "1"
    else
       write(this%unit, "(A)") "0"
    endif
    write(this%unit, "(2G16.7)") xmin, xmax
  end subroutine coop_asy_topaxis_s

  subroutine coop_asy_topaxis_d(this, xmin, xmax, islog,label)
    COOP_REAL  xmin, xmax
    logical islog
    class(coop_asy) this
    COOP_UNKNOWN_STRING label
    call coop_asy_topaxis_s(this, real(xmin), real(xmax), islog,label)
  end subroutine coop_asy_topaxis_d


  subroutine coop_asy_rightaxis_s(this, ymin, ymax, islog, label)
    COOP_SINGLE  ymin, ymax
    logical islog
    class(coop_asy) this
    COOP_UNKNOWN_STRING label
    write(this%unit, "(A)")"EXTRA_AXIS"
    write(this%unit, "(A)") "right"
    if(trim(adjustl(label)).ne."")then
       write(this%unit, "(A)") trim(adjustl(label))
    else
       write(this%unit, "(A)") "NULL"
    endif
    if(islog)then
       write(this%unit, "(A)") "1"
    else
       write(this%unit, "(A)") "0"
    endif
    write(this%unit, "(2G16.7)") ymin, ymax
  end subroutine coop_asy_rightaxis_s

  subroutine coop_asy_rightaxis_d(this, ymin, ymax, islog, label)
    COOP_REAL  ymin, ymax
    logical islog
    COOP_UNKNOWN_STRING label
    class(coop_asy) this
    call coop_asy_rightaxis_s(this, real(ymin), real(ymax), islog,label)
  end subroutine coop_asy_rightaxis_d

  subroutine coop_asy_errorbars_d(this, x, y, dy, color, barsize, dy_minus, dx, dx_minus, center_color, center_symbol)
    class(coop_asy)::this
    COOP_REAL,dimension(:)::x, y, dy
    COOP_SINGLE, optional::barsize
    COOP_UNKNOWN_STRING,optional::color, center_symbol, center_color
    COOP_REAL,dimension(:), optional::dx, dx_minus, dy_minus
    COOP_INT::n, i
    COOP_SINGLE,dimension(:),allocatable::dym, dxm, dxp
    n = coop_getdim("coop_asy_contour", size(x), size(y), size(dy))
    allocate(dym(n), dxm(n), dxp(n))
    if(present(dy_minus))then
       if(size(dy_minus).ne.n) call coop_return_error("asy_errorbars", "dy_minus size discrepancy", "stop")
       dym = dy_minus
    else
       dym = dy
    endif
    if(present(dx))then
       if(size(dx).ne.n) call coop_return_error("asy_errorbars", "dx size discrepancy", "stop")       
       dxp = dx
    else
       dxp = 0.
    endif
    if(present(dx_minus))then
       if(size(dx_minus).ne.n) call coop_return_error("asy_errorbars", "dx_minus size discrepancy", "stop")       
       dxm = dx_minus
    else
       dxm = dxp
    endif
    write(this%unit, "(A)") "ERRORBARS"
    write(this%unit, "(I5)") n
    if(present(color))then
       write(this%unit, "(A)") trim(adjustl(color))
    else
       write(this%unit, "(A)") "black"
    endif
    if(present(barsize))then
       write(this%unit, "(G16.7)") max(barsize, 0.)
    else
       write(this%unit, "(A)") "0."  !!automatically determined by asymptote
    endif
    if(present(center_symbol) )then
       write(this%unit, "(A)") trim(adjustl(center_symbol))
    else
       write(this%unit, "(A)") "DOT"
    endif
    if(present(center_color))then
       write(this%unit, "(A)") trim(adjustl(center_color))
    else
       write(this%unit, "(A)") "black_solid_5"
    endif
    do i= 1,n
       write(this%unit, "(6G16.7)") real(x(i)), real(y(i)), dxp(i), dxm(i), real(dy(i)), dym(i)
    enddo
    deallocate(dym, dxm, dxp)
  end subroutine coop_asy_errorbars_d


  subroutine coop_asy_errorbars_s(this, x, y, dy, color, barsize, dy_minus, dx, dx_minus, center_color, center_symbol)
    class(coop_asy)::this
    COOP_SINGLE,dimension(:)::x, y, dy
    COOP_SINGLE, optional::barsize
    COOP_UNKNOWN_STRING,optional::color, center_symbol, center_color
    COOP_SINGLE,dimension(:), optional::dx, dx_minus, dy_minus
    COOP_INT::n, i
    COOP_SINGLE,dimension(:),allocatable::dym, dxm, dxp
    n = coop_getdim("coop_asy_contour", size(x), size(y), size(dy))
    allocate(dym(n), dxm(n), dxp(n))
    if(present(dy_minus))then
       if(size(dy_minus).ne.n) call coop_return_error("asy_errorbars", "dy_minus size discrepancy", "stop")
       dym = dy_minus
    else
       dym = dy
    endif
    if(present(dx))then
       if(size(dx).ne.n) call coop_return_error("asy_errorbars", "dx size discrepancy", "stop")       
       dxp = dx
    else
       dxp = 0.
    endif
    if(present(dx_minus))then
       if(size(dx_minus).ne.n) call coop_return_error("asy_errorbars", "dx_minus size discrepancy", "stop")       
       dxm = dx_minus
    else
       dxm = dxp
    endif
    write(this%unit, "(A)") "ERRORBARS"
    write(this%unit, "(I5)") n
    if(present(color))then
       write(this%unit, "(A)") trim(adjustl(color))
    else
       write(this%unit, "(A)") "black"
    endif
    if(present(barsize))then
       write(this%unit, "(G16.7)") max(barsize, 0.)
    else
       write(this%unit, "(A)") "0."  !!automatically determined by asymptote
    endif
    if(present(center_symbol) )then
       write(this%unit, "(A)") trim(adjustl(center_symbol))
    else
       write(this%unit, "(A)") "DOT"
    endif
    if(present(center_color))then
       write(this%unit, "(A)") trim(adjustl(center_color))
    else
       write(this%unit, "(A)") "black_solid_5"
    endif
    do i= 1,n
       write(this%unit, "(6G16.7)") real(x(i)), real(y(i)), dxp(i), dxm(i), real(dy(i)), dym(i)
    enddo
    deallocate(dym, dxm, dxp)
  end subroutine coop_asy_errorbars_s  
  
  subroutine coop_asy_error_bar_d(this, x, y, dy_minus, dy_plus, dx_minus, dx_plus, color)
    class(coop_asy) this
    COOP_REAL  x, y
    COOP_REAL ,optional::dy_minus, dy_plus, dx_minus, dx_plus
    COOP_UNKNOWN_STRING, optional::color
    if(present(color))then
       call coop_asy_label(this, "$\bullet$", x, y, color)
!!$       if(present(dy_minus))then
!!$          call coop_asy_line(this, x, y, x, y-dy_minus, linewidth = 1., color=color)
!!$          call coop_asy_label(this, "-", x, y-dy_minus, color=color)
!!$       endif
!!$       if(present(dy_plus))then
!!$          call coop_asy_line(this, x, y, x, y+dy_plus, linewidth = 1., color= color)
!!$          call coop_asy_label(this, "-", x, y+dy_plus, color=color)
!!$       endif
!!$       if(present(dx_minus))then
!!$          call coop_asy_line(this, x, y, x-dx_minus, y, linewidth = 1., color=color)
!!$          call coop_asy_label(this, "{\tiny $|$}", x-dx_minus, y, color=color)
!!$       endif
!!$       if(present(dx_plus))then
!!$          call coop_asy_line(this, x, y, x+dx_plus, y, linewidth = 1., color=color)
!!$          call coop_asy_label(this, "{\tiny $|$}", x+dx_plus, y, color=color)
!!$       endif
    else
       call coop_asy_label(this, "$\bullet$", x, y)
       if(present(dy_minus))then
          call coop_asy_line(this, x, y, x, y-dy_minus, linewidth = 1.)
          call coop_asy_label(this, "-", x, y-dy_minus)
       endif
       if(present(dy_plus))then
          call coop_asy_line(this, x, y, x, y+dy_plus, linewidth = 1.)
          call coop_asy_label(this, "-", x, y+dy_plus)
       endif
       if(present(dx_minus))then
          call coop_asy_line(this, x, y, x-dx_minus, y, linewidth = 1.)
          call coop_asy_label(this, "{\tiny $|$}", x-dx_minus, y)
       endif
       if(present(dx_plus))then
          call coop_asy_line(this, x, y, x+dx_plus, y, linewidth = 1.)
          call coop_asy_label(this, "{\tiny $|$}", x+dx_plus, y)
       endif
    endif
  end subroutine coop_asy_error_bar_d

  subroutine coop_asy_error_bar_s(this, x, y, dy_minus, dy_plus, dx_minus, dx_plus, color)
    class(coop_asy) this
    COOP_SINGLE  x, y
    COOP_SINGLE ,optional::dy_minus, dy_plus, dx_minus, dx_plus
    COOP_UNKNOWN_STRING, optional::color
    if(present(color))then
       call coop_asy_label(this, "$\bullet$", x, y, color)
       if(present(dy_minus))then
          call coop_asy_line(this, x, y, x, y-dy_minus, linewidth = 1., color=color)
          call coop_asy_label(this, "-", x, y-dy_minus, color=color)
       endif
       if(present(dy_plus))then
          call coop_asy_line(this, x, y, x, y+dy_plus, linewidth = 1., color= color)
          call coop_asy_label(this, "-", x, y+dy_plus, color=color)
       endif
       if(present(dx_minus))then
          call coop_asy_line(this, x, y, x-dx_minus, y, linewidth = 1., color=color)
          call coop_asy_label(this, "{\tiny $|$}", x-dx_minus, y, color=color)
       endif
       if(present(dx_plus))then
          call coop_asy_line(this, x, y, x+dx_plus, y, linewidth = 1., color=color)
          call coop_asy_label(this, "{\tiny $|$}", x+dx_plus, y, color=color)
       endif
    else
       call coop_asy_label(this, "$\bullet$", x, y)
       if(present(dy_minus))then
          call coop_asy_line(this, x, y, x, y-dy_minus, linewidth = 1.)
          call coop_asy_label(this, "-", x, y-dy_minus)
       endif
       if(present(dy_plus))then
          call coop_asy_line(this, x, y, x, y+dy_plus, linewidth = 1.)
          call coop_asy_label(this, "-", x, y+dy_plus)
       endif
       if(present(dx_minus))then
          call coop_asy_line(this, x, y, x-dx_minus, y, linewidth = 1.)
          call coop_asy_label(this, "{\tiny $|$}", x-dx_minus, y)
       endif
       if(present(dx_plus))then
          call coop_asy_line(this, x, y, x+dx_plus, y, linewidth = 1.)
          call coop_asy_label(this, "{\tiny $|$}", x+dx_plus, y)
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

  subroutine coop_asy_histogram(x, nbins, filename, xlabel, ylabel, fit_gaussian)
    COOP_INT nbins    
    COOP_UNKNOWN_STRING::filename
    COOP_REAL::x(:)
    COOP_REAL::xcopy(size(x))
    COOP_REAL c(nbins), xb(nbins), xbar, sigma, A
    COOP_INT::i, n, nstep, ntail, j, next, nstep1, nstep2, nstep3
    COOP_UNKNOWN_STRING,optional::xlabel, ylabel
    logical, optional::fit_gaussian
    type(coop_asy)::fig
    n = size(x)
    if(n.lt. 48 .or. nbins.lt. 6) stop "asy_histogram assumes at least 6 bins and 48 samples"
    xcopy = x
    call coop_quicksort(xcopy)
    nstep1 = n/nbins/4
    if(nstep1 .lt. 2)stop "too many bins for asy_histogram: max number of bins is n/8"        
    nstep2 = nstep1*2
    nstep3 = nstep1*3
    if(nbins.gt.6)then
       nstep = (n- (nstep1+nstep2+nstep3)*2)/(nbins-6)
    endif
    ntail = n - nstep*(nbins-6) - (nstep1+nstep2+nstep3)*2
    nstep1 = nstep1 + ntail/6
    nstep2 = nstep2 + ntail/6
    nstep3 = (n - (nstep1+nstep2)*2 - nstep*(nbins-6) ) / 2
    j = nstep1
    xb(1) = sum(xcopy(1:j))/nstep1
    c(1) = nstep1/((xcopy(j+1)+xcopy(j))/2.d0 - xcopy(1) - (xcopy(2)-xcopy(1))/2.d0)

    do i=2, nbins-1
       next = j+1
       if(i.eq.2 .or. i.eq. nbins-1)then
          j = j + nstep2
       elseif(i.eq.3.or.i.eq.nbins-2)then
          j = j + nstep3
       else
          j = j + nstep
       endif
       xb(i) = sum(xcopy(next:j))/(j-next+1)
       c(i) = (j-next+1)/((xcopy(j+1)+xcopy(j))/2.d0-(xcopy(next)+xcopy(next-1))/2.d0)
    enddo
    next = j+1
    xb(nbins) = sum(xcopy(next:n))/(n-j)
    c(nbins) = (n-j)/((xcopy(n)+(xcopy(n)-xcopy(n-1))/2.d0)-(xcopy(j)+xcopy(next))/2.d0)
    call fig%open(trim(adjustl(filename)))
    if(present(xlabel))then
       if(present(ylabel))then
          call fig%init(xlabel = xlabel, ylabel = ylabel, ymin = 0., ymax = 1.12)
       else
          call fig%init(xlabel = xlabel, ylabel = "P", ymin = 0., ymax = 1.12)          
       endif
    else
       if(present(ylabel))then
          call fig%init(xlabel = "x", ylabel = ylabel, ymin = 0., ymax = 1.12)          
       else
          call fig%init(xlabel = "x", ylabel = "P",ymin = 0., ymax = 1.12)
       endif
    endif
    if(present(fit_gaussian))then
       if(fit_gaussian)then
          call coop_fit_gaussian(x, nbins, xbar, sigma, A)
          c = c/n*sigma*sqrt(coop_2pi)
          call coop_asy_curve(fig, xb, c, color="red", linewidth = 1.5)         
          c = A/n*exp(-(xb-xbar)**2/sigma**2/2.d0)
          call coop_asy_curve(fig, xb, c, color="blue", linetype="dotted", linewidth = 1.)                              
       else
          xbar = sum(x)/n
          sigma = sqrt(sum((x-xbar)**2)/n)
          c = c/n*sigma*sqrt(coop_2pi)
          call coop_asy_curve(fig, xb, c)          
       endif
    else
       c = c/maxval(c)
       call coop_asy_curve(fig, xb, c)
    endif
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



  subroutine coop_asy_band_d(this, x, ylower, yupper, colorfill, smooth, linecolor, linetype, linewidth, legend)
    class(coop_asy) this
    logical,optional::smooth
    COOP_REAL ,dimension(:),intent(IN)::x,ylower, yupper
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype, legend
    COOP_SINGLE ,optional::linewidth
    COOP_INT  i,n
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_contour", size(x), size(ylower), size(yupper))
    write(this%unit, "(A)") "CONTOUR"
    write(this%unit, "(A)") "1" !!type 1 contour
    if(present(colorfill))then
       write(this%unit, "(A)") trim(colorfill)
    else
       write(this%unit, "(A)") "lightgray"
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
       write(this%unit, "(A)") trim(lineproperty)
    else
       write(this%unit, "(A)") "NULL"
    endif
    if(present(legend))then
       if(trim(legend).eq.'')then
          write(this%unit, '(A)') "NULL"
       else
          write(this%unit, '(A)') trim(legend)
       endif
    else
       write(this%unit, '(A)') "NULL"
    endif

    if(present(smooth))then
       if(smooth)then
          write(this%unit, "(A)") "1"
       else
          write(this%unit, "(A)") "0"
       endif
    else
       write(this%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(this%unit, "(A)") "1"  !!1 path 
    write(this%unit, "(I8)") n*2
    if(x(2).gt. x(1))then
       do i = 1, n
          call this%write_coor(real(x(i), sp), real(ylower(i), sp))
       enddo
       do i = n, 1, -1
          call this%write_coor(real(x(i), sp), real(yupper(i), sp))
       enddo
    else
       do i = n, 1, -1
          call this%write_coor(real(x(i), sp), real(ylower(i), sp))
       enddo
       do i = 1, n
          call this%write_coor(real(x(i), sp), real(yupper(i), sp))
       enddo
    endif
  end subroutine coop_asy_band_d

  subroutine coop_asy_band_s(this, x, ylower, yupper, colorfill, smooth, linecolor, linetype, linewidth, legend)
    class(coop_asy) this
    logical,optional::smooth
    COOP_SINGLE ,dimension(:),intent(IN)::x,ylower, yupper
    COOP_UNKNOWN_STRING,optional::colorfill
    COOP_UNKNOWN_STRING,optional:: linecolor, linetype, legend
    COOP_SINGLE ,optional::linewidth
    COOP_INT  i,n
    COOP_STRING lineproperty
    n = coop_getdim("coop_asy_contour", size(x), size(ylower), size(yupper))
    write(this%unit, "(A)") "CONTOUR"
    write(this%unit, "(A)") "1" !!type 1 contour
    if(present(colorfill))then
       write(this%unit, "(A)") trim(colorfill)
    else
       write(this%unit, "(A)") "lightgray"
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
       write(this%unit, "(A)") trim(lineproperty)
    else
       write(this%unit, "(A)") "NULL"
    endif
    if(present(legend))then
       if(trim(legend).eq.'')then
          write(this%unit, '(A)') "NULL"
       else
          write(this%unit, '(A)') trim(legend)
       endif
    else
       write(this%unit, '(A)') "NULL"
    endif

    if(present(smooth))then
       if(smooth)then
          write(this%unit, "(A)") "1"
       else
          write(this%unit, "(A)") "0"
       endif
    else
       write(this%unit, "(A)") "0"  !!no smoothing by default
    endif
    write(this%unit, "(A)") "1"  !!1 path 
    write(this%unit, "(I8)") n*2
    if(x(2).gt. x(1))then
       do i = 1, n
          call this%write_coor(real(x(i), sp), real(ylower(i), sp))
       enddo
       do i = n, 1, -1
          call this%write_coor(real(x(i), sp), real(yupper(i), sp))
       enddo
    else
       do i = n, 1, -1
          call this%write_coor(real(x(i), sp), real(ylower(i), sp))
       enddo
       do i = 1, n
          call this%write_coor(real(x(i), sp), real(yupper(i), sp))
       enddo
    endif
  end subroutine coop_asy_band_s

  function coop_asy_xrel(this, xratio) result(xrel)
    class(coop_asy)::this
    COOP_SINGLE  xratio
    COOP_SINGLE  xrel
    if(this%xlog)then
       xrel = this%xmax**xratio*this%xmin**(1.d0-xratio)
    else
       xrel = xratio*this%xmax + (1.d0-xratio)*this%xmin
    endif
  end function coop_asy_xrel

  function coop_asy_yrel(this, yratio) result(yrel)
    class(coop_asy)::this
    COOP_SINGLE  yratio
    COOP_SINGLE  yrel
    if(this%ylog)then
       yrel = this%ymax**yratio*this%ymin**(1.d0-yratio)
    else
       yrel = yratio*this%ymax + (1.d0-yratio)*this%ymin
    endif
  end function coop_asy_yrel


  !!force adjusting boundary
  subroutine coop_asy_expand(this, xl, xr, yl, yr)
    class(coop_asy) this
    COOP_SINGLE xl, xr, yl, yr, dx, dy
    write(this%unit, "(A)") "EXPAND"
    write(this%unit, "(4G16.7)")  xl, xr, yl, yr
    dx = this%xmax - this%xmin
    dy = this%ymax - this%ymin
    this%xmin = this%xmin - dx*xl
    this%xmax = this%xmax + dx*xr
    this%ymin = this%ymin - dy*yl
    this%ymax = this%ymax + dy*yr
  end subroutine coop_asy_expand

  subroutine coop_asy_annotate_d(this, text, x, y, xtext, ytext, color, linetype, linewidth)
    class(coop_asy)::this
    COOP_UNKNOWN_STRING::text
    COOP_UNKNOWN_STRING,optional::color, linetype
    COOP_REAL::x, y, xtext, ytext
    COOP_SINGLE,optional::linewidth
    COOP_STRING::lineproperty
    write(this%unit, "(A)") "ANNOTATE"
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
    write(this%unit, "(A)") trim(lineproperty)
    write(this%unit, "(4G16.7)") x, y, xtext, ytext
    if(trim(text).eq."")then
       write(this%unit, "(A)") "NULL"
    else
       write(this%unit, "(A)") trim(text)
    endif
  end subroutine coop_asy_annotate_d

  subroutine coop_asy_annotate_s(this, text, x, y, xtext, ytext, color, linetype, linewidth)
    class(coop_asy)::this
    COOP_UNKNOWN_STRING::text
    COOP_UNKNOWN_STRING,optional::color, linetype
    COOP_SINGLE::x, y, xtext, ytext
    COOP_SINGLE,optional::linewidth
    COOP_STRING::lineproperty
    write(this%unit, "(A)") "ANNOTATE"
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
    write(this%unit, "(A)") trim(lineproperty)
    write(this%unit, "(4G16.7)") x, y, xtext, ytext
    if(trim(text).eq."")then
       write(this%unit, "(A)") "NULL"
    else
       write(this%unit, "(A)") trim(text)
    endif
  end subroutine coop_asy_annotate_s



  subroutine coop_asy_arrow_d(this, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy)::this
    COOP_REAL xstart, ystart, xend, yend
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING lineproperty
    write(this%unit, "(A)") "ARROWS"
    write(this%unit, "(A)") "1"
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
    write(this%unit, "(A)") trim(lineproperty)
    call this%write_coor(real(xstart, sp), real(ystart, sp), real(xend, sp), real(yend, sp))    
  end subroutine coop_asy_arrow_d

  subroutine coop_asy_arrow_s(this, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy)::this
    COOP_SINGLE xstart, ystart, xend, yend
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING lineproperty
    write(this%unit, "(A)") "ARROWS"
    write(this%unit, "(A)") "1"
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
    write(this%unit, "(A)") trim(lineproperty)
    call this%write_coor(xstart, ystart, xend, yend)
  end subroutine coop_asy_arrow_s

  subroutine coop_asy_arrows_d(this, n, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy)::this
    COOP_INT n, i
    COOP_REAL xstart(n), ystart(n), xend(n), yend(n)
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING lineproperty
    if(n.gt. 9999)then
       write(*, '(A, I10)') "n = ", n
       stop "cannot plot so many arrows"
    endif
    write(this%unit, "(A)") "ARROWS"
    write(this%unit, "(I5)") n
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
    write(this%unit, "(A)") trim(lineproperty)
    do i=1, n
       call this%write_coor(real(xstart(i), sp), real(ystart(i), sp), real(xend(i), sp), real(yend(i), sp))
    enddo
  end subroutine coop_asy_arrows_d
  
  subroutine coop_asy_arrows_s(this, n, xstart, ystart, xend, yend, color, linetype, linewidth)
    class(coop_asy)::this
    COOP_INT n, i
    COOP_SINGLE xstart(n), ystart(n), xend(n), yend(n)
    COOP_UNKNOWN_STRING,optional:: color, linetype
    COOP_SINGLE ,optional::linewidth
    COOP_STRING lineproperty
    if(n.gt. 9999)then
       write(*, '(A, I10)') "n = ", n
       stop "cannot plot so many arrows"
    endif
    write(this%unit, "(A)") "ARROWS"
    write(this%unit, "(I5)") n
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
    write(this%unit, "(A)") trim(lineproperty)
    do i=1, n
       call this%write_coor(xstart(i), ystart(i), xend(i), yend(i))
    enddo
  end subroutine coop_asy_arrows_s

  subroutine coop_asy_path_get_perimeter_and_area(this, ipath, perimeter, area)
    class(coop_asy_path)::this
    COOP_SINGLE:: perimeter, area
    COOP_INT::ipath
    COOP_SINGLE:: xy(2), lastxy(2), beginxy(2)
    COOP_INT::i, j
    if(this%nclosed .lt. ipath) stop "asy_path_perimeter: ipath overflow"    
    if(ipath.gt.1)then
       i = sum(this%length(1:ipath-1))+1
    else
       i = 1
    endif
    call this%l%get_element(i, beginxy)
    lastxy = 0.
    perimeter = 0.
    area = 0.
    do j = 2, this%length(ipath)
       i = i  + 1
       call this%l%get_element(i, xy)
       xy = xy - beginxy
       perimeter = perimeter + sqrt(sum((xy-lastxy)**2))
       area = area + lastxy(1)*xy(2) - lastxy(2)*xy(1)
       lastxy = xy
    enddo
    xy = 0.d0
    perimeter = perimeter + sqrt(sum((xy-lastxy)**2))
    area = area + lastxy(1)*xy(2) - lastxy(2)*xy(1)
    area = area/2. !abs(area)/2.
  end subroutine coop_asy_path_get_perimeter_and_area

  subroutine coop_asy_color2rgba_d(color, rgba)
    COOP_UNKNOWN_STRING::color
    COOP_REAL::rgba(4)
    COOP_STRING::sc
    COOP_INT::i, istart, iend, l
    sc = trim(adjustl(color))
    call coop_str2lower(sc)
    l = len_trim(sc)
    if(sc(1:3) .eq. "rgb")then
       if(sc(4:4) .eq. "a")then
          call coop_str2arr(sc(6:), rgba)
          rgba(1:3) = rgba(1:3)/255.
       elseif(sc(4:4) .eq. ":")then
          call coop_str2arr(sc(5:), rgba(1:3))
          rgba(1:3) = rgba(1:3)/255.
          rgba(4) = 1.
       else
          goto 100
       endif
    elseif(sc(1:5) .eq. "gray:")then
       read(sc(6:), *) rgba(1)
       rgba(1) = rgba(1)/255.
       rgba(2:3) = rgba(1)
       rgba(4) = 1.
    else
       select case(trim(sc))
       case("red","r")
          rgba = (/ 1., 0., 0., 1. /)
       case("green","g")
          rgba = (/ 0., 1., 0., 1. /)
       case("blue","b")
          rgba = (/ 0., 0., 1., 1. /)
       case("yellow", "y")
          rgba = (/ 1., 1., 0., 1. /)
       case("magenta","m")
          rgba = (/ 1., 0., 1., 1. /)
       case("cyan","turquoise", "c")
          rgba = (/ 0., 1., 1., 1. /)
       case("olive")
          rgba = (/ 0.5, 0.5, 0., 1. /)
       case("black","k")
          rgba = (/ 0., 0. ,0., 1. /)
       case("white","w")
          rgba = (/ 1., 1., 1., 1. /)
       case("gray","grey")
          rgba = (/ 0.5, 0.5, 0.5, 1. /)
       case("orange", "o")
          rgba = (/ 1., 0.5, 0.2, 1. /)
       case("violet", "v")
          rgba = (/ 0.55, 0.22, 0.79, 1. /)
       case("brown")
          rgba = (/ 0.5, 0.25, 0., 1. /)
       case("pink")
          rgba = (/ 0.98, 0.69, 0.73, 1. /)
       case("gold")
          rgba = (/ 0.83, 0.63, 0.09, 1. /)
       case("purple","p")
          rgba = (/ 0.56, 0.21, 0.94, 1. /)
       case("maroon")
          rgba = (/ 0.51, 0.02, 0.25, 1. /)
       case("slateblue")
          rgba = (/ 0.21, 0.45, 0.78, 1. /)
       case("skyblue")
          rgba = (/ 0.24, 0.6, 1., 1. /)
       case("lawngreen","grassgreen", "springgreen")
          rgba = (/ 0.53, 0.97, 0.10, 1. /)
       case("darkgray","darkgrey")
          rgba = (/ 0.25, 0.25, 0.25, 1. /)
       case("lightgray","lightgrey")
          rgba = (/ 0.75, 0.75, 0.75, 1. /)          
       case("darkblue")
          rgba = (/ 0., 0., 0.55, 1. /)
       case("darkgreen")
          rgba = (/ 0., 0.55, 0., 1. /)
       case("darkred")
          rgba = (/ 0.55, 0., 0., 1. /)
       case("lightred")
          rgba = (/ 1., 0.1, 0.1, 1. /)
       case("lightblue")
          rgba = (/ 0.1, 1., 0.1, 1. /)
       case("lightgreen")
          rgba = (/ 0.1, 0.1, 1., 1. /)
       case("royalblue")
          rgba = (/ 0.25, 0.41, 0.95, 1./)
       case("darkcyan")
          rgba = (/ 0., 0.55, 0.55, 1. /)
       case("darkmagenta")
          rgba = (/ 0.55, 0., 0.55, 1. /)
       case("darkyellow")
          rgba = (/ 0.55, 0.55, 0., 1. /)
       case("darkbrown")
          rgba = (/ 0.35, 0.15, 0., 1. /)
       case("invisible")
          rgba = 0.d0
       case default
          goto 100
       end select
    endif
    return
100 write(*,*) trim(color)
    stop "Unknown color name."
  end subroutine coop_asy_color2rgba_d

  subroutine coop_asy_rgba2color_d(rgba, color)
    COOP_REAL::rgba(4)
    COOP_STRING::color
    color = "RGBA:"//COOP_STR_OF(nint(rgba(1)*255.))//":"//COOP_STR_OF(nint(rgba(2)*255.))//":"//COOP_STR_OF(nint(rgba(3)*255.))//":"//trim(adjustl(coop_num2str(rgba(4), "(F10.3)")))
  end subroutine coop_asy_rgba2color_d


  subroutine coop_asy_color2rgba_s(color, rgba)
    COOP_UNKNOWN_STRING::color
    COOP_SINGLE::rgba(4)
    COOP_REAL::tmp(4)
    call coop_asy_color2rgba_d(color, tmp)
    rgba = tmp
  end subroutine coop_asy_color2rgba_s

  subroutine coop_asy_rgba2color_s(rgba, color)
    COOP_SINGLE::rgba(4)
    COOP_STRING::color
    color = "RGBA:"//COOP_STR_OF(nint(rgba(1)*255.))//":"//COOP_STR_OF(nint(rgba(2)*255.))//":"//COOP_STR_OF(nint(rgba(3)*255.))//":"//trim(adjustl(coop_num2str(rgba(4), "(F10.3)")))
  end subroutine coop_asy_rgba2color_s

  
end module coop_asy_mod

