program test
  use coop_wrapper_firstorder
  use coop_coupledDE_collapse_mod
  implicit none
#include "constants.h"
#if DO_COUPLED_DE
  type(coop_dictionary)::dict
  COOP_INT,parameter::n_threads = 8
  COOP_STRING::params_file, output_root, output_format
  type(coop_coupledDE_collapse_params)::params(n_threads)
  COOP_INT::n_fpk, n_e, n_pbye
  COOP_REAL::fpk_max, fpk_min, e_max, e_min, pbye_min, pbye_max, output_zmax
  COOP_REAL,dimension(:,:,:),allocatable::zcol
  COOP_REAL,dimension(:),allocatable::fpk, e, pbye, z, a
  COOP_INT::ifpk, ie, ipbye, output_nz, ithread, i

  type(coop_file)::fp
  if(iargc().lt. 1)then
     write(*,*) "========================================================"
     write(*,*) "Syntax:"
     write(*,*) "./CDExport params.ini"
     write(*,*) "========================================================"
     stop
  endif

  !!=================import the cosmology from an ini file ===============
  call coop_get_Input(1, params_file)
  call coop_load_dictionary(params_file, dict)
  call coop_dictionary_lookup(dict, "output_root", output_root)
  call coop_dictionary_lookup(dict, "output_format", output_format, "BINARY")
  select case(trim(output_format))
  case("BINARY")
     call fp%open(trim(output_root)//"_zcol.dat", "s")
  case("FITS")
     stop "FITS format is to be implemented in future versions."
  case default
     stop "This format has not been implemented."
  end select
  call params(1)%init(dict, update_cosmology = .true.)
  do ithread = 2, n_threads
     params(ithread) = params(1)
  enddo
  !!==================import the halo settings ==========================
  call coop_dictionary_lookup(dict, "fpk_max", fpk_max, 30.d0)
  call coop_dictionary_lookup(dict, "fpk_min", fpk_min, 1.5d0)
  call coop_dictionary_lookup(dict, "e_min", e_min, 0.d0)
  call coop_dictionary_lookup(dict, "e_max", e_max, 0.5d0)
  call coop_dictionary_lookup(dict, "pbye_min", pbye_min, -1.d0)
  call coop_dictionary_lookup(dict, "pbye_max", pbye_max, 1.d0)

  call coop_dictionary_lookup(dict, "n_fpk", n_fpk, 50)
  call coop_dictionary_lookup(dict, "n_e", n_e, 30)
  call coop_dictionary_lookup(dict, "n_pbye", n_pbye, 30)
  allocate(zcol(n_fpk, n_e, n_pbye))
  allocate(fpk(n_fpk), e(n_e), pbye(n_pbye))
  call coop_set_uniform(n_fpk, fpk, fpk_min, fpk_max, logscale = .true.)
  call coop_set_uniform(n_e, e, e_min, e_max)
  call coop_set_uniform(n_pbye, pbye, pbye_min, pbye_max)

  !$omp parallel do
  do ithread = 1, n_threads
     do ipbye = 1, n_pbye
        do ie = 1, n_e
           do ifpk = 1, n_fpk
              call params(ithread)%update_fep(fpk(ifpk), e(ie), pbye(ipbye)*e(ie))
              zcol(ifpk, ie, ipbye) = params(ithread)%zvir1()
           enddo
        enddo
     enddo
  enddo
  !$omp end parallel do
  !!for other background functions
  call coop_dictionary_lookup(dict, "output_zmax", output_zmax, 30.d0)
  call coop_dictionary_lookup(dict, "output_nz", output_nz, 5000)
  allocate(z(output_nz), a(output_nz))
  call coop_set_uniform(output_nz, z, 0.d0, output_zmax)
  a = 1.d0/(1.d0+z)
  select case(trim(output_format))
  case("BINARY")
     write(fp%unit) n_fpk, fpk_min, fpk_max
     write(fp%unit) n_e, e_min, e_max
     write(fp%unit) n_pbye, pbye_min, pbye_max
     write(fp%unit) zcol
     call fp%close()
     write(*,*) "====================================================================="
     write(*,*) "z_collapse(f_pk, e, p/e) is saved to "//trim(output_root)//"_zcol.dat"
     write(*,*) "sample code to read it:"
     write(*,*) "real*8 fpk_min, fpk_max, e_min, e_max, pbye_min, pbye_max"
     write(*,*) "integer n_fpk, n_e, n_pbye"
     write(*,*) "real*8,allocatable::zcol(:,:,:)"
     write(*,*) "open(..., access ='stream')"
     write(*,*) "read(...) n_fpk, fpk_min, fpk_max"
     write(*,*) "read(...) n_e, e_min, e_max"
     write(*,*) "read(...) n_pbye, pbye_min, pbye_max"
     write(*,*) "allocate(zcol(n_fpk, n_e, n_pbye))"
     write(*,*) "read(...) zcol"
     write(*,*) "****************************"
     write(*,*) "fpk is log uniform, e and p/e are both uniform"
     write(*,*) "====================================================================="

     call fp%open(trim(output_root)//"_background.dat", "s")
     write(fp%unit) output_nz, output_zmax
     do i=1, output_nz
        write(fp%unit) z(i), params(1)%cosmology%time(a(i)),  params(1)%dadt(a(i))/a(i), params(1)%cosmology%comoving_distance(a(i)),  params(1)%growth_D(a(i)), params(1)%Growth_H_D(a(i))
     enddo
     call fp%close()
     write(*,*) "====================================================================="     
     write(*,*) "background info is saved to "//trim(output_root)//"_background.dat"
     write(*,*) "to read it:"
     write(*,*) "real*8 zmax"
     write(*,*) "integer nz, iz"
     write(*,*) "real*8, dimension(:),allocatable::z, t, H, chi, D, H_D"
     write(*,*) "open(..., access ='stream')"
     write(*,*) "read(...) nz, zmax"
     write(*,*) "allocate(z(nz), t(nz), H(nz), chi(nz), D(nz), H_D(nz)"
     write(*,*) "do iz = 1, nz"
     write(*,*) "  read(...) z(iz), t(iz), H(iz), chi(iz), D(iz), H_D(iz)"
     write(*,*) "enddo"
     write(*,*) "***************************************"
     write(*,*) "units: "
     write(*,*) "t:  1/H_0"
     write(*,*) "H:  H_0"
     write(*,*) "chi: c/H_0"
     write(*,*) "H_D:  H_0  (H_D is defined as d ln D/ dt)"
     write(*,*) "====================================================================="

  case("FITS")
     stop "FITS format is to be implemented in future versions."
  case default
     stop "This format has not been implemented."
  end select

#else
  stop "to use CDExport you need to compile the code with DARK_ENERGY_MODEL = COUPLED_DE in configure.in"
#endif
end program test
