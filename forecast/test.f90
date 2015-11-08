program test
  use coop_wrapper_firstorder
  use coop_forecast_mod
  implicit none
#include "constants.h"
  type(coop_asy)::fig
  COOP_INT, parameter::lmin = 2, lmax = 2500
  COOP_REAL::Cl(5, lmin:lmax), ells(lmin:lmax), lcdm_Cl(5, lmin:lmax), lcdm_lensedCl(5, lmin:lmax)
  call coop_MPI_init()
  ells = coop_uniform_array(lmax-lmin+1, lmin, 1)
  call fig%open("Cl_alphaM.txt")
  call fig%init(xlabel = "$\ell$", ylabel = "$C_\ell / C_\ell^{\rm \Lambda CDM} - 1$", xlog=.true., ymin = -0.08, ymax = 0.08)
  call loadCls("chains/lcdm_scalCls.txt", 5)
  lcdm_Cl = Cl
  call loadCls("chains/lcdm_lensedCls.txt", 4)
  lcdm_lensedCl = Cl
  
  call loadCls("chains/alphaM_p2_scalCls.txt", 5)
  call fig%interpolate_curve(ells, Cl(1,:)/lcdm_Cl(1,:) - 1.d0,  "LogLinear", color="red", linetype="solid", linewidth=2.2, legend="$\alpha_{M0} = 0.2$, unlensed TT")
  call fig%interpolate_curve(ells, Cl(4,:)/lcdm_Cl(4,:) - 1.d0, "LogLinear", color="red", linetype="dashed", linewidth=2.2, legend="$\alpha_{M0} = 0.2$, $\phi\phi$")

!!$  call loadCls("chains/alphaM_p2_lensedCls.txt", 3)
!!$  call fig%interpolate_curve(ells, Cl(1,:)/lcdm_lensedCl(1,:) - 1.d0,  "LogLinear", color="red", linetype="dotted", linewidth=2.5, legend="$\alpha_{M0} = 0.2$, lensed TT")

  

  call loadCls("chains/alphaM_p4_scalCls.txt", 5)
  call fig%interpolate_curve(ells, Cl(1,:)/lcdm_Cl(1,:) - 1.d0,  "LogLinear", color="blue", linetype="solid", linewidth=1.2, legend="$\alpha_{M0} = 0.4$, unlensed TT")
  call fig%interpolate_curve(ells, Cl(4,:)/lcdm_Cl(4,:) - 1.d0, "LogLinear", color="blue", linetype="dashed", linewidth=1.2, legend="$\alpha_{M0} = 0.4$, $\phi\phi$")

!!$  call loadCls("chains/alphaM_p4_lensedCls.txt", 3)
!!$  call fig%interpolate_curve(ells, Cl(1,:)/lcdm_lensedCl(1,:) - 1.d0,  "LogLinear", color="blue", linetype="dotted", linewidth=1.5, legend="$\alpha_{M0} = 0.4$, lensed TT")
  
  
  call fig%legend(0.38, 0.36, 1)
  call fig%close()
  call coop_MPI_finalize()

contains
  subroutine loadCls(fname, ncols)
    COOP_UNKNOWN_STRING::fname
    type(coop_file)::fp
    COOP_INT::l, il, ncols
    call fp%open(fname, "r")
    do l= lmin, lmax
       read(fp%unit, *)  il, Cl(1:ncols, l)
    enddo
    call fp%close()
  end subroutine loadCls
end program test
