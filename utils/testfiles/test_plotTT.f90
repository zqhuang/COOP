program Daubechies
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_asy)::fig
  character(LEN=*),parameter::genre = "TT"
  character(LEN=*),parameter::interp = "loglinear"  
  character(LEN=*),parameter::root = "cls/wavelet0123"  
  integer,parameter::nls = 2508
  COOP_REAL::lcdm_power(2:nls,5)
  COOP_REAL::lcdm_ells(2:nls)
  integer i
  call load_lcdm()
  call fig%open("PlanckDl"//genre//".txt")
  if(interp .eq. 'loglinear')then
     call fig%init(width=6., height = 4.5, xlabel="$\ell$", ylabel="$D_l^{"//genre//"}$", xlog=.true., xmin=1.8, xmax=2600.)
  else
     call fig%init(width=6., height = 4.5, xlabel="$\ell$", ylabel="$D_l^{"//genre//"}$", xlog=.false., xmin=1.8, xmax=2520.)     
  endif
  i=1
  call plot_from_file(root//"_"//COOP_STR_OF(i)//".theory_cl", genre, "skyblue", "dotted", 0.8, "wavelet $1\sigma$ samples")  
  do i=2, 50
     call plot_from_file(root//"_"//COOP_STR_OF(i)//".theory_cl", genre, "skyblue", "dotted", 0.8)
  enddo
  call plot_from_file(root//"_best.theory_cl", genre, "blue", "solid", 1.5, "wavelet best-fit")
  call plot_lcdm(genre, "red", "solid")  
  call plot_cl_data(genre, "darkgray")
  if(genre=="TT")then
     call fig%legend(xratio = 0.05, yratio = 0.92)
     call fig%label(xratio = 0.05, yratio = 0.6, label="wavelet reconstruction")
  else
     call fig%legend(xratio = 0.05, yratio = 0.3)
     call fig%label(xratio = 0.05, yratio = 0.1, label="wavelet reconstruction")
     
  endif
  call fig%close()

contains

  subroutine load_lcdm()
    integer l
    character(LEN=1024)::line
    type(coop_file)::fp
    call fp%open("planckcls/COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt")
    read(fp%unit,"(A)")line    
    do l = 2, nls
       read(fp%unit, *) lcdm_ells(l), lcdm_power(l, :)
       if(l .ne. nint(lcdm_ells(l))) stop "plot_lcdm error: l does not match"
    enddo
    call fp%close()
  end subroutine load_lcdm


  function lcdm_cl(genre, rl)
    character(LEN=*) genre
    COOP_REAL::rl, lcdm_cl
    integer l, ind
    select case(trim(genre))
    case("TT")
       ind = 1
    case("TE")
       ind = 2
    case("EE")
       ind = 3
    case("BB")
       ind = 4
    case("PP")
       ind = 5
    case default
       stop "plot_lcdm error: unknown genre"
    end select
    l = min(max(floor(rl), 2), nls-1)
    lcdm_cl = (rl-l)*lcdm_power(l+1, ind) + (l+1.d0-rl)*lcdm_power(l, ind)
  end function lcdm_cl

  subroutine plot_from_file(filename, genre, color, linetype, linewidth, legend)
    character(LEN=*) filename, genre, color, linetype
    character(LEN=*),optional::legend
    real linewidth
    COOP_INT::ind, l
    COOP_INT,parameter::nls = 2500
    character(LEN=1024)::line      
    type(coop_file)::fp
    COOP_REAL::power(2:nls, 5), ells(2:nls)
    select case(trim(genre))
    case("TT")
       ind = 1
    case("TE")
       ind = 2
    case("EE")
       ind = 3
    case("BB")
       ind = 4
    case("PP")
       ind = 5
    case default
       stop "plot_lcdm error: unknown genre"
    end select
    call fp%open(filename)
    read(fp%unit, "(A)")line
    do l=2, nls
       read(fp%unit, *) ells(l), power(l, :)
    enddo
    call fp%close()
    if(present(legend))then
       call fig%interpolate_curve(xraw=ells, yraw=power(:, ind), interpolate=interp, color=color, linetype=linetype, linewidth=linewidth, legend = legend)
    else
       call fig%interpolate_curve(xraw=ells, yraw=power(:, ind), interpolate=interp, color=color, linetype=linetype, linewidth=linewidth)
    endif
  end subroutine plot_from_file
  

  subroutine plot_cl_data(genre, color)
    character(LEN=*) genre, color
    type(coop_file)::fp    
    COOP_REAL::ell, Dl, err_up, err_down
    COOP_INT::i
    character(LEN=1024)::line  
    call fp%open("planckcls/"//trim(genre)//"_binned.txt")
    read(fp%unit, "(A)")line
    i = 0
    do
       read(fp%unit, *, ERR=100, END=100) ell, Dl, err_down, err_up
       call fig%error_bar(x=ell, y=Dl, dy_plus = err_up, dy_minus=err_down, color=color)
       i = i+1
    enddo
100 call fp%close()
    print*, i , " data points  in file "//genre//"_binned.txt"
  end subroutine plot_cl_data

  subroutine plot_lcdm(genre, color, linetype)
    character(LEN=*) genre, color, linetype
    COOP_INT::ind
    select case(trim(genre))
    case("TT")
       ind = 1
    case("TE")
       ind = 2
    case("EE")
       ind = 3
    case("BB")
       ind = 4
    case("PP")
       ind = 5
    case default
       stop "plot_lcdm error: unknown genre"
    end select
    call fig%interpolate_curve(xraw=lcdm_ells, yraw=lcdm_power(:, ind), interpolate=interp, color=color, linetype=linetype, linewidth=1., legend="$\Lambda$CDM best-fit")
  end subroutine plot_lcdm

  
end program Daubechies
