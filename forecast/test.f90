program test
  use coop_wrapper_firstorder
  use coop_forecast_mod
  implicit none
#include "constants.h"
  type(coop_asy)::fig
  COOP_INT, parameter::lmin = 2, lmax = 2500
  COOP_REAL::Cl(5, lmin:lmax), ells(lmin:lmax), lcdm_Cl(5, lmin:lmax)
  call coop_MPI_init()
  ells = coop_uniform_array(lmax-lmin+1, lmin, 1)
  call fig%open("Cl_alphaM.txt")
  call fig%init(xlabel = "$\ell$", ylabel = "$C_\ell / C_\ell^{\rm \Lambda CDM} - 1$", xlog=.true.)
  call loadCls("chains/lcdm_scalCls.txt")
  lcdm_Cl = Cl
  call loadCls("chains/alphaM_p1_scalCls.txt")
  call fig%curve(ells, Cl(1,:)/lcdm_Cl(1,:) - 1.d0, color="red", linetype="solid", linewidth=2., legend="$\alpha_{M0} = 0.1$, unlensed TT")
  call loadCls("chains/alphaM_p2_scalCls.txt")
  call fig%curve(ells, Cl(1,:)/lcdm_Cl(1,:) - 1.d0, color="red", linetype="dashed", linewidth=2., legend="$\alpha_{M0} = 0.2$, unlensed TT")
  call loadCls("chains/alphaM_p3_scalCls.txt")
  call fig%curve(ells, Cl(1,:)/lcdm_Cl(1,:) - 1.d0, color="red", linetype="dotted", linewidth=2., legend="$\alpha_{M0} = 0.3$, unlensed TT")
  call fig%legend(0.3, 0.4)
  call fig%close()
  call coop_MPI_finalize()

contains
  subroutine loadCls(fname)
    COOP_UNKNOWN_STRING::fname
    type(coop_file)::fp
    COOP_INT::l, il
    call fp%open(fname, "r")
    do l= lmin, lmax
       read(fp%unit, *)  il, Cl(:, l)
    enddo
    call fp%close()
  end subroutine loadCls
end program test
