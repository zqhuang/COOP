program test
  use coop_wrapper_firstorder
  use coop_forecast_mod
  implicit none
#include "constants.h"
  type(coop_dataset_SN_Simple), target::SN(1)
  type(coop_data_pool)::pool  
  type(coop_MCMC_params)::mcmc
  COOP_REAL::z_g = -4.d-5
  COOP_STRING::action = "BEST"
  COOP_INT i
  call coop_random_init()
  call mcmc%init(prefix = "chains/kcdm", paramnames = "paramnames/wcdm.paramnames", ini = "myinis/kcdm.ini")
  if(iargc() .lt. 2) stop "Test action z_g"  
  call coop_get_Input(1, action)
  call coop_get_Input(2, z_g)
  write(*,*) "action = "//trim(action)
  write(*,*) "z_g = "//COOP_STR_OF(z_g)
  select case(trim(action))
  case("SIMU","simu")
     call SN(1)%import("data/sn_z_mu_dmu_plow_union2.1.txt")
     mcmc%params = mcmc%center
     mcmc%fullparams(mcmc%used) = mcmc%center
     write(*,*) "fiducial parameters"
     do i=1, mcmc%fulln
        write(*,*) trim(mcmc%name(i))//" = ", mcmc%fullparams(i)
     enddo
     call SN(1)%simulate(mcmc)
     SN(1)%z = SN(1)%z + z_g
     SN(1)%mu = SN(1)%mu + 5.d0*log10(1.d0+z_g/(1.d0+SN(1)%z))
     call SN(1)%export("data/sn_"//COOP_STR_OF(nint(z_g*1.d5))//".txt")
     write(*,*) "output: data/sn_"//COOP_STR_OF(nint(z_g*1.d5))//".txt"
  case("BEST","best")
     call SN(1)%import("data/sn_"//COOP_STR_OF(nint(z_g*1.d5))//".txt")
     pool%SN_Simple =>  SN
     call mcmc%findbest(pool, temperature  = 1.d-4)
     do i=1, mcmc%n
        print*, trim(mcmc%name(mcmc%used(i)))//" = ", mcmc%bestparams(i), " +/- ", sqrt(mcmc%covmat(i,i))
     enddo
     print*, "derived parameters = ", mcmc%derived()
     print*,"best like = ", mcmc%bestlike
  case("MCMC","mcmc")
     call SN(1)%import("data/sn_"//COOP_STR_OF(nint(z_g*1.d5))//".txt")
     pool%SN_Simple =>  SN
     do i = 1, 400000
        if(mod(i, 5000).eq.0)then
           print*, i
           if(i.lt.20000)call mcmc%update_propose
        endif
        call mcmc%mcmc_step(pool)
     enddo
  case default
     print*, trim(action)
     stop "unknown action"
  end select
  
end program test
