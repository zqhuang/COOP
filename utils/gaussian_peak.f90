module gaussian_peak_utils
  use specfunc_utils
  implicit none

contains

  function gp_n_maxima(dim, sigma0, sigma1, sigma2, nu) result(nex)
    integer dim
    real(dl) sigma0, sigma1, sigma2, nex, nu
    select case(dim)
    case(1)
       if(present(nu))then
          nex = sigma2/sigma1/const_pi*erfc(nu/const_sqrt2)/4.d0
       else
          nex = sigma2/sigma1/const_pi/2.d0
       endif
    case(2)

    end select
  end function gp_n_maxima


  function gp_n_minima(dim, sigma0, sigma1, sigma2, nu) result(nex)
    integer dim
    real(dl) sigma0, sigma1, sigma2, nex
    real(dl),optional::nu
    select case(dim)
    case(1)
       if(present(nu))then
          nex = sigma2/sigma1/const_pi*erfc(nu/const_sqrt2)/4.d0
       else
          nex = sigma2/sigma1/const_pi/2.d0
       endif
    case(2)

    end select
  end function gp_n_minima

  function gp_n_extrema(dim, sigma0, sigma1, sigma2) result(nex)
    integer dim
    real(dl) sigma0, sigma1, sigma2, nex
    nex = gp_n_maxima(dim, sigma0, sigma1, sigma2) +  gp_n_minima(dim, sigma0, sigma1, sigma2)
  end function gp_n_extrema


end module gaussian_peak_utils
