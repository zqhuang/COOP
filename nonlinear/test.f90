program test
  use coop_wrapper_firstorder
  use coop_halofit_mod
  implicit none
#include "constants.h"
  COOP_INT,parameter::nz = 20
  type(coop_asy)::fig, figdiff
  type(coop_dynamic_array_real)::rlQS, rlGR, rlfull
  COOP_REAL::z(nz)
  COOP_INT::i, n_want
  z = (/ 20.d0, 10.d0, 5.d0, 2.d0,  1.d0, 0.9d0, 0.8d0, 0.7d0, 0.6d0, 0.5d0, 0.45d0, 0.4d0, 0.35d0, 0.3d0, 0.25d0, 0.2d0, 0.15d0, 0.1d0, 0.05d0, 0.d0/)


  do i=1, 20
     call fig%open("density"//COOP_STR_OF(i)//".txt")
     call fig%init(xlabel = "$rH_0$", ylabel = "${\rho a^3}/{\rho_0}$", ymin = 0., ymax = 100.)
     call figdiff%open("diff"//COOP_STR_OF(i)//".txt")
     call figdiff%init(xlabel = "$rH_0$", ylabel = "rel. err. of QS approx.", ymin = -0.005, ymax  = 0.005)

     n_want = 250

     call rlQS%load_txt("col_QS_"//trim(coop_ndigits(i,2))//".txt")
     call rlGR%load_txt("col_GR_"//trim(coop_ndigits(i,2))//".txt")
     call rlfull%load_txt("col_acc3_"//trim(coop_ndigits(i,2))//".txt")
     call fig%plot(rlGR%f(1:n_want, 1), exp(rlGR%f(1:n_want, 2)), color = "darkgray", linewidth = 2., linetype="solid", legend = "GR; $z="//COOP_STR_OF(z(i))//"$")
     call fig%plot(rlQS%f(1:n_want, 1), exp(rlQS%f(1:n_want, 2)), color = "blue", linewidth = 1.2, linetype="solid", legend = "QS approx; $z="//COOP_STR_OF(z(i))//"$")
     call fig%plot(rlfull%f(1:n_want, 1), exp(rlfull%f(1:n_want, 2)), color = "red", linewidth = 2., linetype="dotted", legend = "full; $z="//COOP_STR_OF(z(i))//"$")
     call fig%legend(0.46, 0.92)
     call fig%close()



     call figdiff%plot(rlQS%f(1:n_want, 1), (exp(rlQS%f(1:n_want, 2)- rlfull%f(1:n_want,2))-1.d0), color = "blue", linewidth = 1.5, linetype="solid", legend = "$z="//COOP_STR_OF(z(i))//"$")
     call figdiff%legend(0.3, 0.4)
     call figdiff%close()

     call system("../utils/fasy.sh density"//COOP_STR_OF(i)//".txt")
     call system("../utils/fasy.sh diff"//COOP_STR_OF(i)//".txt")
  enddo

end program test
