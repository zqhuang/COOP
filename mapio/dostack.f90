program Stacking_Maps
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX
  logical::do_max = .true.
  COOP_STRING::stack_field_name = "T"
  COOP_STRING::map_file = "planck14/dx11_v2_smica_int_cmb_010a_1024.fits"
  COOP_STRING::stack_options_file = ""
  COOP_STRING::imask_file = "planck14/dx11_v2_common_int_mask_010a_1024.fits"
  COOP_STRING::polmask_file = "planck14/dx11_v2_common_pol_mask_010a_1024.fits"
  COOP_INT,parameter::n = 36
  COOP_REAL,parameter::r_degree  = 2.d0
  COOP_REAL,parameter::dr = 2.d0*sin(r_degree*coop_SI_degree/2.d0)/n
  
  type(coop_stacking_options)::sto
  type(coop_healpix_patch)::patch
  type(coop_healpix_maps)::hgm, mask, pmap
  COOP_STRING::output
  
  call sto%import(stack_options_file)
  call hgm%read(map_file)
  call patch%init(stack_field_name, n, dr)
  if(all(patch%tbs%spin).eq.0)then
     call mask%read(imask_file)
  else
     call mask%read(polmask_file)
  endif
  call hgm%stack_on_peaks(sto, patch, mask)
  if(patch%nmaps .gt. 1)then
     do i=1, patch%nmaps
        call patch%plot(i, output//"_"//COOP_STR_OF(i)//".txt")
     enddo
  else
     call patch%plot(1, output)
  endif
#else
  print*, "You need to install healpix"
#endif  
end program Stacking_Maps
