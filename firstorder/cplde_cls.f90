program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  !!----------------------------------------
  !!output Cls file
  COOP_STRING::output = "cls_coupled_DE.txt"
  !!----------------------------------------
  !! declare other variables
  type(coop_cosmology_firstorder)::cosmology
  type(coop_function)::fwp1, fQ
  COOP_INT, parameter::lmin = 2, lmax = 2500
  COOP_REAL::Cls(coop_num_Cls, lmin:lmax)
  COOP_REAL::norm
  COOP_INT::l
  logical success
  type(coop_file)::fp
  !!----------------------------------------
  !!main code
  !!initialize cosmology
  call fwp1%init_polynomial( (/ 0.d0.d0, 0.2d0 /) )
  call fQ%init_polynomial( (/ 0.d0,  0.d0 /) )
  call cosmology%set_coupled_DE_cosmology(Omega_b=0.049d0, Omega_c=0.265d0, h = 0.68d0, tau_re = 0.06d0, As = 2.21d-9, ns = 0.968d0, fwp1 = fwp1, fQ = fQ)
  if(cosmology%h() .eq. 0.d0) stop "Initialization failed; check the input parameters."
  call cosmology%compute_source(0, success = success)
  if(.not. success) stop "Solution blows up exponentially; Model is ruled out."
  call cosmology%source(0)%get_all_cls(lmin, lmax, Cls)
  call fp%open(output,"w")
  write(fp%unit, "(A8, 5A16)") "# ell ", "   TT  ",  "   EE  ",  "   TE   ", "Phi_lens Phi_lens ", " T Phi_lens  "
  do l = lmin, lmax
     norm = l*(l+1.d0)/coop_2pi*cosmology%Tcmb()**2*1.d12
     write(fp%unit, "(I8, 5E16.7)") l, Cls(coop_index_ClTT, l)*norm,  Cls(coop_index_ClEE, l)*norm,  Cls(coop_index_ClTE, l)*norm,  Cls(coop_index_ClLenLen, l)*(l*(l+1.d0))**2, Cls(coop_index_ClTLen, l)*(l*(l+1.d0))**1.5
  enddo
  call fp%close()
  
end program test
