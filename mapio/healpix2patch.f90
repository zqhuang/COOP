program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map
  type(coop_healpix_disc) disc
  COOP_STRING::fin, fout
  COOP_INT::pix, imap, n, i, j
  COOP_REAL::l_deg, b_deg, r_deg, theta, phi, rot_deg, mean
  type(coop_healpix_patch)::patch
  type(coop_file)::fp
  call coop_MPI_init()
  if(iargc() .lt. 4)then
     print*,"Syntax:"
     print*,"./CutPatch -inp INPUT_FILE -out OUTPUT_FILE -l L_IN_DEGREE -b B_IN_DEGREE -radius RADIUS_IN_DEGREE -sig #_OF_MAP -n RESOLUTION[(2n+1)*(2n+1)] -rot ROTATE_ANGLE_IN_DEGREE"     
     stop
  endif
  call coop_get_command_line_argument(key = "inp", arg = fin)
  call coop_get_command_line_argument(key = "out", arg = fout)
  call coop_get_command_line_argument(key = "l", arg = l_deg, default = 0.d0)
  call coop_get_command_line_argument(key = "b", arg = b_deg, default = 90.d0)
  call coop_get_command_line_argument(key = "radius", arg = r_deg, default = 10.d0)
  call coop_get_command_line_argument(key = "sig", arg = imap, default = 1)
  call coop_get_command_line_argument(key = "n", arg = n, default = 50)
  call coop_get_command_line_argument(key = "rot", arg = rot_deg, default=0.d0)          
  call coop_healpix_lb2ang(l_deg, b_deg, theta, phi)
  coop_healpix_fmissval = -300. !!map%bad_data    
  call map%read(fin)
  if(imap .ne. 1)map%map(:,1) = map%map(:, imap)
  map%spin(1) = 0
  call map%ang2pix(theta, phi, pix)
  call map%get_disc(pix, disc)
  call patch%init("I", n, r_deg/n*coop_SI_degree)
  call map%fetch_patch(disc, rot_deg*coop_SI_degree, patch)
  mean = sum(patch%image, patch%image .gt. 0.95*coop_healpix_fmissval)/count(patch%image .gt. 0.95*coop_healpix_fmissval)
  where(patch%image(-n:n, -n:n,1) .gt. -299.d0)
     patch%image(-n:n,-n:n, 1) = patch%image(-n:n,-n:n,1) - mean
  end where
  call fp%open(fout, "w")
  do i = -n, n
     write(fp%unit, "("//COOP_STR_OF(2*n+1)//"E14.5)") patch%image(i, :, 1)
  enddo
  call fp%close()
  call coop_MPI_Finalize()
end program test
