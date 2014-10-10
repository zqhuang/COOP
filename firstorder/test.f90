program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::fod
  type(coop_file)::fp, fpxe
  COOP_INT::i
  COOP_REAL z, xe
  !!set cosmology
  call fod%Set_Planck_bestfit()

    !if you want extended models
!  call fod%set_standard_cosmology(Omega_b=0.047d0, Omega_c=0.26d0, h = 0.68d0, tau_re = 0.08d0, nu_mass_eV = 0.06d0, As = 2.15d-9, ns = 0.962d0, nrun = -0.01d0, r = 0.2d0, nt = -0.01d0, YHe = 0.25d0, Nnu = 3.d0)
  call fpxe%open('xe.txt')
  call fp%open("xe2.txt")
  do 
     read(fpxe%unit, *, END=100, ERR=100) z, xe
     write(fp%unit, "(3E16.7)") z, xe, fod%xeofa(1.d0/(1.d0+z))
  enddo
100 continue
  call fp%close()
  call fpxe%close()
end program test
