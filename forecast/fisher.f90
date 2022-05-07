module coop_fisher_mod
  use coop_halofit_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"

  COOP_REAL,parameter::coop_fisher_tolerance = 1.d-6
  COOP_INT,parameter::coop_parameter_type_slow = 0
  COOP_INT,parameter::coop_parameter_type_fast = 1
  COOP_INT,parameter::coop_parameter_type_nuis = 2
  COOP_REAL,parameter::coop_H0_unit = 1.d5/coop_SI_c

  COOP_REAL,parameter::coop_fisher_cmb_l_pivot = 2000.d0

  type coop_observation
     COOP_STRING::filename = ""
     type(coop_dictionary)::settings
     COOP_STRING::name = ""
     COOP_STRING::genre = ""
     COOP_INT::n_obs = 0
     COOP_INT::n_window = 0
     COOP_INT::n_bins = 0
     COOP_INT::dim_obs = 0
     COOP_INT::dim_nuis = 0
     COOP_INT::init_level = 0
     type(coop_int_table)::paramnames
     COOP_REAL::h_calibrate = 0.7d0
     COOP_REAL,dimension(:,:),allocatable::obs !!(dim_obs, n_obs)
     COOP_REAL,dimension(:,:,:),allocatable::dobs !!(dim_obs, n_obs, paramnames%n)
     COOP_REAL,dimension(:,:),allocatable::nuis  !!(dim_nuis, n_obs)
     COOP_REAL,dimension(:,:,:),allocatable::invcov !!(dim, dim, n_obs)
     !!for MPK
     COOP_REAL,dimension(:,:),allocatable::window_Wsq !!(n_window, n_bins)
     COOP_REAL,dimension(:,:),allocatable::window_modes !!(n_window, n_bins)
     COOP_INT,dimension(:),allocatable::window_used !!(n_bins)
     !!for CMB
     COOP_INT::lmin, lmax, lmin_pol, lmax_pol
   contains
     procedure::free => coop_observation_free
     procedure::init => coop_observation_init
     procedure::get_dobs => coop_observation_get_dobs
     procedure::get_invcov => coop_observation_get_invcov
  end type coop_observation
  
  type coop_fisher
     type(coop_cosmology_firstorder)::cosmology
     type(coop_dictionary)::settings
     type(coop_real_table)::paramtable
     COOP_INT::n_params = 0
     COOP_INT::n_params_used = 0
     COOP_INT::n_slow
     COOP_INT::n_fast
     COOP_INT::n_nuis
     COOP_INT::n_observations = 0
     COOP_INT::max_init_level
     COOP_REAL,dimension(:),allocatable::params
     COOP_REAL,dimension(:),allocatable::step1
     COOP_REAL,dimension(:),allocatable::step2
     COOP_REAL,dimension(:),allocatable::priors
     COOP_INT,dimension(:),allocatable::ind_slow, ind_fast, ind_nuis, param_types, init_level
     COOP_REAL,dimension(:,:),allocatable::fisher
     COOP_REAL,dimension(:,:),allocatable::Cov
     COOP_INT,dimension(:),allocatable::ind_used
     logical,dimension(:),allocatable::is_used
     type(coop_observation),dimension(:),allocatable::observations
   contains
     procedure::init => coop_fisher_init
     procedure::free => coop_fisher_free
     procedure::get_fisher => coop_fisher_get_fisher
  end type coop_fisher

contains

  subroutine coop_observation_get_dobs(this, dobs, paramtable, cosmology)
    class(coop_observation)::this
    COOP_REAL::dobs(this%dim_obs, this%n_obs), MStar
    COOP_REAL,dimension(:),allocatable::b0, b2
    COOP_REAL::sigma_g, sigma_z, sr2, a, Hz, cmb_A_noise, cmb_n_noise, cmb_run_noise, cmb_A_noise_pol, cmb_n_noise_pol, cmb_run_noise_pol, cmb_A_tSZ, cmb_TE_leakage_eps0, cmb_TE_leakage_eps2, cmb_TE_leakage_eps4, leakage, logl
    COOP_INT::nk, nz, nmu, iz
    type(coop_real_table)::paramtable
    type(coop_cosmology_firstorder)::cosmology
    COOP_INT::i, idata, l
    select case(trim(this%genre))
    case("SN")
       call paramtable%lookup("sn_absolute_m", Mstar)
       !$omp parallel do
       do idata = 1, this%n_obs
          dobs(1, idata) = 5.d0*log10(cosmology%luminosity_distance(1.d0/(1.d0+this%nuis(1, idata)))/cosmology%H0Mpc()) + Mstar  - this%obs(1, idata) 
       enddo
       !$omp end parallel do
    case("MPK")
       call coop_dictionary_lookup(this%settings,"n_z", nz)
       call coop_dictionary_lookup(this%settings,"n_k", nk)
       call coop_dictionary_lookup(this%settings,"n_mu", nmu)
       call paramtable%lookup("mpk_sigma_g", sigma_g, 400.d0)
       sigma_g =sigma_g*1.d3/coop_SI_c
       call coop_dictionary_lookup(this%settings,"sigma_z", sigma_z, 0.001d0) 
       sr2 = sigma_g**2+sigma_z**2
       allocate(b0(nz), b2(nz))
       do iz=1, nz
          call paramtable%lookup("mpk_b0_"//COOP_STR_OF(iz), b0(iz))
          call paramtable%lookup("mpk_b2_"//COOP_STR_OF(iz), b2(iz), 0.d0)
       enddo
       !$omp parallel do private(idata, iz, a, Hz)
       do idata = 1, this%n_obs
          iz = (idata-1)/(nk*nmu)+1
          a = 1.d0/(1.d0+this%nuis(1,idata))
          Hz = cosmology%Hratio(a)
          dobs(1, idata) = (b0(iz) + b2(iz)*this%nuis(2, idata)**2 &
               +  cosmology%fgrowth_of_z(z=this%nuis(1, idata), k=this%nuis(2, idata))*this%nuis(3,idata)**2)**2 &
               * cosmology%smeared_matter_power(z=this%nuis(1, idata), k=this%nuis(2, idata), nw = this%window_used(iz), kw = this%window_modes(1:this%window_used(iz), iz), wsq = this%window_wsq(1:this%window_used(iz), iz) ) * ((coop_pi**2*2.d0)/this%nuis(2, idata)**3) &
               * exp(-sr2*((1.d0+this%nuis(1,idata))/Hz *this%nuis(2,idata)*this%nuis(3,idata))**2)*(this%h_calibrate/cosmology%h())**3 - this%obs(1, idata)
       enddo
       !$omp end parallel do
       deallocate(b0, b2)
    case("CMB_TE")
       call paramtable%lookup("cmb_A_noise", cmb_A_noise, 1.d0)
       call paramtable%lookup("cmb_n_noise", cmb_n_noise, 0.d0)
       call paramtable%lookup("cmb_run_noise", cmb_run_noise, 0.d0)
       call paramtable%lookup("cmb_A_noise_pol", cmb_A_noise_pol, 1.d0)
       call paramtable%lookup("cmb_n_noise_pol", cmb_n_noise_pol, 0.d0)
       call paramtable%lookup("cmb_run_noise_pol", cmb_run_noise_pol, 0.d0)
       call paramtable%lookup("cmb_A_tSZ", cmb_A_tSZ, 1.d0)
       call paramtable%lookup("cmb_TE_leakage_eps0", cmb_TE_leakage_eps0, 0.d0)
       call paramtable%lookup("cmb_TE_leakage_eps2", cmb_TE_leakage_eps2, 0.d0)
       call paramtable%lookup("cmb_TE_leakage_eps4", cmb_TE_leakage_eps4, 0.d0)
       do idata = 1, this%n_obs
          leakage = cmb_TE_leakage_eps0+ (cmb_TE_leakage_eps2 + cmb_TE_leakage_eps4*(l/coop_fisher_cmb_l_pivot)**2)*(l/coop_fisher_cmb_l_pivot)**2
          l = idata - 1 + min(this%lmin, this%lmin_pol)
          logl = log(l/coop_fisher_cmb_l_pivot)
          dobs(1, idata) = cosmology%source(0)%Cls_lensed(coop_index_ClTT, l) + this%nuis(3, idata)*cmb_A_noise*exp(logl*(cmb_n_noise+cmb_run_noise/2.d0*logl)) + this%nuis(5, idata)*cmb_A_tSZ - this%obs(1, idata)
          dobs(2, idata) = cosmology%source(0)%Cls_lensed(coop_index_ClEE, l)+ cosmology%source(0)%Cls_lensed(coop_index_ClTT, l)*leakage**2 + this%nuis(4, idata)*cmb_A_noise_pol*exp(logl*(cmb_n_noise_pol+cmb_run_noise_pol/2.d0*logl))- this%obs(2, idata)
          dobs(3, idata) = cosmology%source(0)%Cls_lensed(coop_index_ClTE, l) + cosmology%source(0)%Cls_lensed(coop_index_ClTT, l)*leakage - this%obs(3, idata)
       enddo
    case("CMB_T")
       call paramtable%lookup("cmb_A_tSZ", cmb_A_tSZ, 1.d0)
       call paramtable%lookup("cmb_A_noise", cmb_A_noise, 1.d0)
       call paramtable%lookup("cmb_n_noise", cmb_n_noise, 0.d0)
       call paramtable%lookup("cmb_run_noise", cmb_run_noise, 0.d0)
       do idata = 1, this%n_obs
          l = idata - 1 + min(this%lmin, this%lmin_pol)
          logl = log(l/coop_fisher_cmb_l_pivot)
          dobs(1, idata) = cosmology%source(0)%Cls_lensed(coop_index_ClTT, l) + this%nuis(3, idata)*cmb_A_noise*exp(logl*(cmb_n_noise+cmb_run_noise/2.d0*logl)) + this%nuis(5, idata)*cmb_A_tSZ - this%obs(1, idata)
       enddo
    case("CMB_E")
       call paramtable%lookup("cmb_A_noise_pol", cmb_A_noise_pol, 1.d0)
       call paramtable%lookup("cmb_n_noise_pol", cmb_n_noise_pol, 0.d0)
       call paramtable%lookup("cmb_run_noise_pol", cmb_run_noise_pol, 0.d0)
       do idata = 1, this%n_obs
          l = idata - 1 + min(this%lmin, this%lmin_pol)
          logl = log(l/coop_fisher_cmb_l_pivot)
          dobs(1, idata) = cosmology%source(0)%Cls_lensed(coop_index_ClEE, l) + this%nuis(4, idata)*cmb_A_noise_pol*exp(logl*(cmb_n_noise_pol+cmb_run_noise_pol/2.d0*logl)) - this%obs(1, idata)
       enddo
    case("CMB_B")
       call paramtable%lookup("cmb_A_noise_pol", cmb_A_noise_pol, 1.d0)
       call paramtable%lookup("cmb_n_noise_pol", cmb_n_noise_pol, 0.d0)
       call paramtable%lookup("cmb_run_noise_pol", cmb_run_noise_pol, 0.d0)
       do idata = 1, this%n_obs
          l = idata - 1 + min(this%lmin, this%lmin_pol)
          logl = log(l/coop_fisher_cmb_l_pivot)
          dobs(1, idata) = cosmology%source(0)%Cls_lensed(coop_index_ClBB, l) + this%nuis(4, idata)*cmb_A_noise_pol*exp(logl*(cmb_n_noise_pol+cmb_run_noise_pol/2.d0*logl)) - this%obs(1, idata)
       enddo
    case default
       write(*,*) trim(this%genre)
       stop "unknown observation genre"
    end select
  end subroutine coop_observation_get_dobs

!!set the fiducial model and compute the inverse covariance matrix
  subroutine coop_observation_get_invcov(this, paramtable, cosmology)
    class(coop_observation)::this
    type(coop_real_table)::paramtable
    type(coop_cosmology_firstorder)::cosmology
    COOP_STRING::mpk_output_root
    logical::do_mpk_output 
    type(coop_file),dimension(:),allocatable::fp
    COOP_INT idata
    COOP_REAL::sn_peculiar_velocity, sn_sigma_int, sn_sigma_meas, sn_sigma_lens, sn_sigma_sys, Mstar, mpk, bias, khMpc
    COOP_INT:: nz, nk, nmu, iz, ik, imu, l, n_per_zbin
    COOP_REAL::sigma_z,sigma_g, sr2, Hz, rz, a
    COOP_REAL,dimension(:),allocatable::b0, b2
    this%h_calibrate = cosmology%h()
    select case(trim(this%genre))
    case("SN")
       call paramtable%lookup("sn_absolute_m", Mstar)
       call coop_dictionary_lookup(this%settings, "sn_sigma_int", sn_sigma_int, 0.1d0)
       call coop_dictionary_lookup(this%settings, "sn_sigma_meas", sn_sigma_meas, 0.01d0)
       call coop_dictionary_lookup(this%settings, "sn_sigma_lens", sn_sigma_lens, 0.07d0)
       call coop_dictionary_lookup(this%settings, "sn_sigma_sys_per_bin", sn_sigma_sys, 0.d0)  !!
       call coop_dictionary_lookup(this%settings, "sn_peculiar_velocity", sn_peculiar_velocity, 400.d0)
       sn_peculiar_velocity = sn_peculiar_velocity*1.d3/coop_SI_c
       !$omp parallel do
       do idata = 1, this%n_obs
          this%obs(1, idata) =  5.d0*log10(cosmology%luminosity_distance(1.d0/(1.d0+this%nuis(1, idata)))/cosmology%H0Mpc()) + Mstar
          this%invcov(1, 1, idata) = 1.d0/ ( (sn_sigma_int**2 + sn_sigma_meas**2 + (sn_sigma_lens*this%nuis(1, idata))**2 + (sn_peculiar_velocity/this%nuis(1, idata))**2)/ this%nuis(2, idata) + (sn_sigma_sys*(1.d0+this%nuis(1, idata)))**2 )
       enddo
       !$omp end parallel do
    case("MPK")
       call coop_dictionary_lookup(this%settings,"n_z", nz)
       call coop_dictionary_lookup(this%settings, "mpk_output_root", mpk_output_root, "")
       if(trim(mpk_output_root).ne. "")then
          do_mpk_output = .true.
          allocate(fp(nz))
          do iz=1, nz
             call fp(iz)%open(trim(mpk_output_root)//"_zbin"//COOP_STR_OF(iz)//".dat")
             write(fp(iz)%unit, "(3A16)") "# k [h Mpc^{-1}]",  "P(k) [h^{-3}Mpc^3]",  "  smeared P(k) "
          enddo
       else 
          do_mpk_output = .false.
       endif
       call coop_dictionary_lookup(this%settings,"n_k", nk)
       call coop_dictionary_lookup(this%settings,"n_mu", nmu)
       n_per_zbin = nk*nmu
       call paramtable%lookup("mpk_sigma_g", sigma_g, 400.d0)
       sigma_g =sigma_g*1.d3/coop_SI_c
       call coop_dictionary_lookup(this%settings,"sigma_z", sigma_z, 0.001d0) 
       sr2 = sigma_g**2+sigma_z**2
       if(this%n_obs .ne. nz*n_per_zbin)then
          write(*,*) "Error in data file "//trim(this%name)//": n_z * n_k*n_mu does not equal to n_obs"
          stop
       endif
       allocate(b0(nz), b2(nz))
       do iz=1, nz
          call paramtable%lookup("mpk_b0_"//COOP_STR_OF(iz), b0(iz))
          call paramtable%lookup("mpk_b2_"//COOP_STR_OF(iz), b2(iz), 0.d0)
       enddo
       do idata = 1, this%n_obs
          iz = (idata-1)/n_per_zbin+1
          a = 1.d0/(1.d0+this%nuis(1,idata))
          Hz = cosmology%Hratio(a)
          bias = (b0(iz) + b2(iz)*this%nuis(2, idata)**2 &
               +  cosmology%fgrowth_of_z(z=this%nuis(1, idata), k=this%nuis(2, idata))*this%nuis(3,idata)**2) **2
          mpk = cosmology%smeared_matter_power(z=this%nuis(1, idata), k=this%nuis(2, idata), nw = this%window_used(iz), kw = this%window_modes(1:this%window_used(iz), iz), wsq = this%window_wsq(1:this%window_used(iz), iz) ) * ((coop_pi**2*2.d0)/this%nuis(2, idata)**3)
          this%obs(1, idata) = bias * mpk * exp(-sr2*((1.d0+this%nuis(1,idata))/Hz *this%nuis(2,idata)*this%nuis(3,idata))**2)
          this%invcov(1,1,idata) = (this%nuis(2,idata)**2*this%nuis(5,idata) &
               * this%nuis(6,idata) * cosmology%comoving_distance(a) **2/Hz*this%nuis(4,idata) &
               * this%nuis(7,idata)/coop_2pi) &
               / (this%obs(1, idata)+1.d0/this%nuis(8,idata))**2
          if(do_mpk_output)then
             imu = mod(idata-1, nmu)+1
             khMpc =  this%nuis(2, idata)*coop_H0_unit*exp(this%nuis(5,idata)/this%nuis(2, idata)*((imu-0.5d0)/nmu-0.5d0))
             write(fp(iz)%unit, "(3E16.7)") khMpc, cosmology%matter_power(z=this%nuis(1, idata), k=khMpc/coop_H0_unit) * ((coop_pi**2*2.d0)/khMpc**3), cosmology%smeared_matter_power(z=this%nuis(1, idata), k=khMpc/coop_H0_unit, nw = this%window_used(iz), kw = this%window_modes(1:this%window_used(iz), iz), wsq = this%window_wsq(1:this%window_used(iz), iz) ) * ((coop_pi**2*2.d0)/khMpc**3)
          endif
       enddo
       if(do_mpk_output)then
          do iz = 1, nz
             call fp(iz)%close()
          enddo
          deallocate(fp)
       endif
       deallocate(b0, b2)
    case("CMB_TE")
       do idata = 1, this%n_obs
          l = idata - 1 + min(this%lmin, this%lmin_pol) 
          this%obs(1, idata) = cosmology%source(0)%Cls_lensed(coop_index_ClTT, l) + this%nuis(3, idata) + this%nuis(5, idata)  !!Cl + Cl^SZ + Nl_TT
          this%obs(2, idata) = cosmology%source(0)%Cls_lensed(coop_index_ClEE, l) + this%nuis(4, idata)   !!Cl + Nl_EE
          this%obs(3, idata) = cosmology%source(0)%Cls_lensed(coop_index_ClTE, l)
          this%invcov(1,1,idata) = 2.d0*this%obs(1, idata)**2
          this%invcov(2,2,idata) = 2.d0*this%obs(2, idata)**2
          this%invcov(3,3,idata) = (this%obs(1, idata)*this%obs(2, idata)+this%obs(3, idata)**2)
          this%invcov(1,2,idata) = 2.d0*this%obs(3, idata)**2
          this%invcov(2,1,idata) = this%invcov(1,2,idata)
          this%invcov(1,3,idata) = 2.d0*this%obs(1, idata)*this%obs(3,idata)
          this%invcov(3,1,idata) = this%invcov(1,3,idata)
          this%invcov(2,3,idata) = 2.d0*this%obs(2, idata)*this%obs(3,idata)
          this%invcov(3,2,idata) = this%invcov(2,3,idata)
          this%invcov(:,:,idata) = this%invcov(:,:,idata)*(dble(l)**4*1.d24)
          call coop_sympos_inverse(3,3,this%invcov(:,:,idata))
          this%invcov(:,:,idata) = this%invcov(:,:,idata)*((2.d0*l+1.d0)*(dble(l)**4*1.d24))
          this%invcov(1, 1, idata) = this%invcov(1, 1, idata)*this%nuis(1, idata)
          this%invcov(2, 2, idata) = this%invcov(2, 2, idata)*this%nuis(2, idata)
          this%invcov(3, 3, idata) = this%invcov(3, 3, idata)*this%nuis(2, idata)
          this%invcov(1, 2, idata) = this%invcov(1, 2, idata)*min(this%nuis(2, idata), this%nuis(1, idata))
          this%invcov(1, 3, idata) = this%invcov(1, 3, idata)*min(this%nuis(2, idata), this%nuis(1, idata))          
          this%invcov(2, 3, idata) = this%invcov(2, 3, idata)*min(this%nuis(2, idata), this%nuis(1, idata))
          this%invcov(2, 1, idata) = this%invcov(1, 2, idata)          
          this%invcov(3, 1, idata) = this%invcov(1, 3, idata)          
          this%invcov(3, 2, idata) = this%invcov(2, 3, idata)          
       enddo
    case("CMB_T")
       do idata = 1, this%n_obs
          l = idata - 1 + min(this%lmin, this%lmin_pol) 
          this%obs(1, idata) = cosmology%source(0)%Cls_lensed(coop_index_ClTT, l) + this%nuis(3, idata) + this%nuis(5, idata)
          this%invcov(1,1,idata) = this%nuis(1, idata)*(l+0.5d0)/this%obs(1, idata)**2
!          if(mod(l, 10).eq.0)print*, l, this%nuis(1, idata), cosmology%source(0)%Cls_lensed(coop_index_ClTT, l), this%nuis(3, idata), this%nuis(5, idata)
       enddo
    case("CMB_E")
       do idata = 1, this%n_obs
          l = idata - 1 + min(this%lmin, this%lmin_pol) 
          this%obs(1, idata) = cosmology%source(0)%Cls_lensed(coop_index_ClEE, l) + this%nuis(4, idata)
          this%invcov(1,1,idata) = this%nuis(2, idata)*(l+0.5d0)/this%obs(1, idata)**2
       enddo
    case("CMB_B")
       do idata = 1, this%n_obs
          l = idata - 1 + min(this%lmin, this%lmin_pol) 
          this%obs(1, idata) = cosmology%source(0)%Cls_lensed(coop_index_ClBB, l) + this%nuis(4, idata)
          this%invcov(1,1,idata) = this%nuis(2, idata)*(l+0.5d0)/this%obs(1, idata)**2
       enddo
    case default
       write(*,*) trim(this%genre)
       stop "unknown observation genre"
    end select
  end subroutine coop_observation_get_invcov

  subroutine coop_observation_init(this, paramtable, filename)
    COOP_REAL,parameter::Nl_infinity = 1.d30
    class(coop_observation)::this
    type(coop_file)::wfp, fpsz
    COOP_UNKNOWN_STRING,optional::filename
    type(coop_real_table)::paramtable
    type(coop_list_string)::ls
    type(coop_list_real)::lr
    COOP_STRING::windowfile, windowall, tSZ_template
    COOP_REAL::kwarr(2)
    COOP_REAL,dimension(:),allocatable::z, kmin, kmax, nobs, dz, mu, k, dk
    COOP_REAL::dmu, kl, sigma_W, fsky, Nl_T, Nl_pol, fg_dust_r, fg_dust_A, fg_dust_alpha, fg_dust_T, fg_dust_beta, fsky_pol, T353, obs_yr, Cl_dust_residual, tsz_template_freq
    COOP_REAL,dimension(:),allocatable::beam_fwhm, sigmaT, sigmapol, freq
    COOP_INT::i, nz,nk,nmu,iz,ik,imu, n_channels, l, iw, junk, lmin_tsz, lmax_tsz
    COOP_REAL,dimension(:),allocatable::cl_tsz
    COOP_REAL::tsz_factor  !!frequency dependence
    if(present(filename))this%filename = trim(adjustl(filename))
    if(trim(this%filename) .eq. "") stop "observation_init: empty file name"
    if(.not. coop_file_exists(this%filename)) then
       write(*,*) "cannot find file "//trim(this%filename)
       stop
    endif
    call coop_load_dictionary(this%filename, this%settings)
    call coop_dictionary_lookup(this%settings, "genre", this%genre)
    this%genre = trim(adjustl(this%genre))
    call coop_str2upper(this%genre)
    call coop_dictionary_lookup(this%settings, "name", this%name, "COOP_OBSERVATION_"//trim(this%genre))
    call coop_dictionary_lookup(this%settings, "n_obs", this%n_obs)
    select case(trim(this%genre))
    case("SN")
       this%init_level = coop_init_level_set_background
       this%dim_obs = 1  !!distance moduli
       this%dim_nuis = 2  !!z, n_samples
    case("MPK")
       this%init_level = coop_init_level_set_pert
       this%dim_obs = 1 !!matter power spectrum
       this%dim_nuis = 8   !!z, k, mu, dz, dk, dmu, fsky, nobs, z_bin_index
    case("CMB_TE")
       this%init_level = coop_init_level_set_Cls
       this%dim_obs = 3 !!TT, TE, EE
       this%dim_nuis = 5 !! fsky, fsky_pol, N_l(TT), N_l(EE), Cl_tSZ(TT)
    case("CMB_T", "CMB_E", "CMB_B")
       this%init_level = coop_init_level_set_Cls
       this%dim_obs = 1
       this%dim_nuis = 5 !!l, fsky, fsky_pol, N_l(TT), N_l(EE), Cl_tSZ(TT)
    case default
       write(*,*) trim(this%genre)
       stop "Error: unknown observation genre"
    end select
    call coop_dictionary_lookup(this%settings, "params", ls)
    do i=1, ls%n
       if(paramtable%index(ls%element(i)).ne. 0)then
          call this%paramnames%insert(ls%element(i), i)  
       else
!!          write(*,*) "skipping parameter "//trim(ls%element(i))//" for observation "//trim(this%name)//"."
       endif
    enddo
    allocate(this%obs(this%dim_obs, this%n_obs), this%nuis(this%dim_nuis, this%n_obs), this%invcov(this%dim_obs, this%dim_obs, this%n_obs), this%dobs(this%dim_obs, this%n_obs, this%paramnames%n))
    call ls%free()

    select case(trim(this%genre))
    case("SN")
       call coop_dictionary_lookup(this%settings, "z", lr)
       if(lr%n .ne. this%n_obs)then
          write(*,*) "Error in "//trim(this%filename)
          write(*,*) "Number of redshift bins does not equal to n_obs"
          stop
       endif
       !$omp parallel do
       do i = 1, this%n_obs
          this%nuis(1, i) = lr%element(i)
       enddo
       !$omp end parallel do
       call lr%free()

       call coop_dictionary_lookup(this%settings, "n_samples", lr)
       if(lr%n .ne. this%n_obs)then
          write(*,*) "Error in "//trim(this%filename)
          write(*,*) "Length of n_samples list does not equal to n_obs"
          stop
       endif
       !$omp parallel do
       do i = 1, this%n_obs
          this%nuis(2, i) = lr%element(i)
       enddo
       !$omp end parallel do
       call lr%free()
    case("MPK")
       call coop_dictionary_lookup(this%settings,"n_z", nz)
       call coop_dictionary_lookup(this%settings,"n_k", nk)
       call coop_dictionary_lookup(this%settings,"n_mu", nmu)
       call coop_dictionary_lookup(this%settings,"fsky", fsky)
       if(fsky .le. 0.d0 .or. fsky .gt. 1.d0)then
          write(*,*) "fsky = ", fsky , ": out of range 0< fsky <=1)"
          stop
       endif
       if(this%n_obs .ne. nz*nk*nmu)then
          write(*,*) "Error in data file "//trim(this%name)//": n_z * n_k * n_mu does not equal to n_obs"
          stop
       endif
       allocate(kmin(nz), kmax(nz), nobs(nz), dz(nz), mu(nmu), k(nk), dk(nk), z(nz))
       call coop_dictionary_lookup(this%settings, "z", lr)
       if(lr%n .ne. nz)then
          write(*,*) "Error in data file "//trim(this%name)//": length of z list does not equal to n_z"
          stop
       endif
       this%n_bins = nz
       do iz = 1, nz
          z(iz) = lr%element(iz)
       enddo
       call lr%free()
       call coop_dictionary_lookup(this%settings, "kmin", lr)
       if(lr%n .ne. nz)then
          write(*,*) "Error in data file "//trim(this%name)//": length of kmin list does not equal to n_z"
          stop
       endif
       do iz = 1, nz
          kmin(iz) = lr%element(iz)/coop_H0_unit
       enddo
       call lr%free()
       COOP_DEALLOC(this%window_used)
       allocate(this%window_used(this%n_bins))
       this%n_window = 0 
       call coop_dictionary_lookup(this%settings, "window_all",  windowall, "")
       if(trim(windowall).ne. "")then
          write(*,*) "Using window: "//trim(windowall)
       endif
       do iz = 1, nz
          if(trim(windowall).ne. "")then
             windowfile = windowall
          else
             call coop_dictionary_lookup(this%settings, "window"//COOP_STR_OF(iz),  windowfile, "")
          endif
          if(trim(windowfile).ne."")then
             if(coop_file_exists(windowfile))then
                this%window_used(iz) = coop_file_NumLines(windowfile)
             elseif(coop_is_number(windowfile))then
                this%window_used(iz) = 1
             else
                write(*,*) "the window file "//trim(windowfile)//" does not exist"
                stop
             endif
          else
             this%window_used(iz) = 1
          endif
          if(this%window_used(iz).gt. this%n_window) this%n_window = this%window_used(iz)
       enddo
       if(this%n_window .gt. 500)then
          write(*,*) "Warning: the window file contain too many k bins and may slow down the fisher calculation significantly."
       endif
       COOP_DEALLOC(this%window_Wsq)
       COOP_DEALLOC(this%window_modes)
       allocate(this%window_Wsq(this%n_window, this%n_bins))
       allocate(this%window_modes(this%n_window, this%n_bins))
       do iz = 1, nz
          if(trim(windowall).ne. "")then
             windowfile = windowall
          else
             call coop_dictionary_lookup(this%settings, "window"//COOP_STR_OF(iz),  windowfile, "")
          endif
          if(trim(windowfile).ne."")then
             if(coop_file_exists(windowfile))then
                call wfp%open(windowfile, 'r')
                write(*,*) "loading window function:"//trim(windowfile)
                iw = 0
                do while(wfp%read_real_array(kwarr))
                   iw = iw + 1
                   this%window_modes(iw, iz) = kwarr(1)/coop_H0_unit
                   this%window_Wsq(iw, iz) = kwarr(2)**2
                enddo
                call wfp%close()
                if(iw .ne. this%window_used(iz))then
                   write(*,*) "Error in window file "//trim(windowfile)
                   stop
                endif
             else
                read(windowfile, *) sigma_W
                sigma_W = sigma_W / coop_H0_unit
                this%window_modes(1, iz) = sigma_W
                this%window_Wsq(1, iz) = exp(-1.d0)
             endif
          else
             sigma_W = kmin(iz)
             this%window_modes(1,iz) = sigma_W
             this%window_Wsq(1,iz) = exp(-1.d0)
          endif
       enddo
       call coop_dictionary_lookup(this%settings, "kmax", lr)
       if(lr%n .ne. nz)then
          write(*,*) "Error in data file "//trim(this%name)//": length of kmax list does not equal to n_z"
          stop
       endif
       do iz = 1, nz
          kmax(iz) = lr%element(iz)/coop_H0_unit
       enddo
       call lr%free()
       call coop_dictionary_lookup(this%settings, "nobs", lr)
       if(lr%n .ne. nz)then
          write(*,*) "Error in data file "//trim(this%name)//": length of nobs list does not equal to n_z"
          stop
       endif
       do iz = 1, nz
          nobs(iz) = lr%element(iz) / coop_H0_unit**3
       enddo
       call lr%free()
       call coop_dictionary_lookup(this%settings, "delta_z", lr)
       if(lr%n .ne. nz)then
          write(*,*) "Error in data file "//trim(this%name)//": length of delta_z list does not equal to n_z"
          stop
       endif
       do iz = 1, nz
          dz(iz) = lr%element(iz)
       enddo
       call lr%free()
       dmu = 2.d0/nmu
       do imu = 1, nmu
          mu(imu) = -1.d0+dmu*(imu-0.5d0)
       enddo
       i = 0
       call coop_dictionary_lookup(this%settings, "k_linear_sampling", kl, 1.d30)  !!default uniform in ln k
       kl = kl/coop_H0_unit
       do iz = 1, nz
          call coop_lss_k_sampling(kmin(iz), kmax(iz), kl, nk, k, dk)
          do ik = 1, nk
             do imu = 1, nmu
                i = i+1
                this%nuis(1, i) = z(iz)
                this%nuis(2, i) = k(ik)
                this%nuis(3, i) = mu(imu)
                this%nuis(4, i) = dz(iz)
                this%nuis(5, i) = dk(ik)
                this%nuis(6, i) = dmu
                this%nuis(7, i) = fsky
                this%nuis(8, i) = nobs(iz)
             enddo
          enddo
       enddo
       deallocate(kmin, kmax, nobs, dz, mu, k,dk,z)
    case("CMB_TE","CMB_T", "CMB_E", "CMB_B")
       call coop_dictionary_lookup(this%settings,"tSZ_template", tSZ_template,"")  
       if(trim(tSZ_template) .ne. "")then
          write(*,*) "reading SZ template "//trim(coop_file_path_of(this%filename))//trim(tSZ_template)
          call fpsz%open(trim(coop_file_path_of(this%filename))//trim(tSZ_template))
          read(fpsz%unit, *) tSZ_template_freq, lmin_tsz, lmax_tsz
          allocate(cl_tsz(lmin_tsz:lmax_tsz))
          do l = lmin_tsz, lmax_tsz
             read(fpsz%unit, *) junk, cl_tsz(l)
             cl_tsz(l) = cl_tsz(l)*(coop_2pi/l/(l+1.d0))/(COOP_DEFAULT_TCMB*1.e6)**2  !!convert to delta T/T power
             if(junk .ne. l)then
                write(*,*) trim(tsz_template)
                stop "tsz template file broken"
             endif
          enddo
       else
          lmax_tsz = -1
          lmin_tsz = 100000
       endif
       call coop_dictionary_lookup(this%settings, "lmin", this%lmin, 2)
       call coop_dictionary_lookup(this%settings, "lmin_pol", this%lmin_pol, this%lmin)
       if(this%lmin .lt. 2 .or. this%lmin_pol .lt. 2) then
          write(*,*) "Error in data file "//trim(this%filename)//": lmin and lmin_pol must be >= 2"
          stop
       endif
       call coop_dictionary_lookup(this%settings, "lmax", this%lmax)
       call coop_dictionary_lookup(this%settings, "lmax_pol", this%lmax_pol, this%lmax)
       coop_Cls_lmax(0) = max(this%lmax, this%lmax_pol, coop_Cls_lmax(0))
       if(max(this%lmax, this%lmax_pol)-min(this%lmin, this%lmin_pol)+1 .ne. this%n_obs)then
          write(*,*) "Error in data file "//trim(this%filename)//": max(lmax, lmax_pol) - min(lmin, lmin_pol) + 1 does not equal to n_obs"
          stop
       endif
       call coop_dictionary_lookup(this%settings,"fsky", fsky)
       call coop_dictionary_lookup(this%settings,"fsky_pol", fsky_pol, fsky)
       call coop_dictionary_lookup(this%settings,"foreground_dust_residual", Fg_dust_r)  
       !!from 0 to 1
       if(fg_dust_r .lt. 0.d0 .or. fg_dust_r .gt. 1.d0)then
          write(*,*) "Error in data file "//trim(this%filename)
          write(*,*) "foreground_dust_residual must be between 0 and 1"
          stop
       endif
       !!dust foreground,  default values from arxiv: 1409.5738
       call coop_dictionary_lookup(this%settings,"foreground_amp80", Fg_dust_A, 100.d0)  !!at 353GHz and l = 80, this depends on the region and fsky, for BICEP2 it is about 13.4; here I use a rough number with fsky = 0.7
      
       call coop_dictionary_lookup(this%settings,"foreground_l_slope", fg_dust_alpha, -2.42d0)  !!l dependence
       call coop_dictionary_lookup(this%settings,"foreground_freq_slope", fg_dust_beta, 1.59d0)  !!frequency dependence
       fg_dust_beta  = fg_dust_beta + 3.d0  !!nu^3 from Blackbody
       call coop_dictionary_lookup(this%settings,"foreground_T353", fg_dust_T,  19.6d0)
       call coop_dictionary_lookup(this%settings,"n_channels", n_channels)
       allocate(beam_fwhm(n_channels), sigmaT(n_channels), sigmapol(n_channels), freq(n_channels))

       tsz_factor = 0.d0
       do i=1, n_channels
          call coop_dictionary_lookup(this%settings,"beam_fwhm_channel"//COOP_STR_OF(i), beam_fwhm(i))
          call coop_dictionary_lookup(this%settings,"frequency_channel"//COOP_STR_OF(i), freq(i))
          tsz_factor = tsz_factor + coop_SZ_frequency_factor(tSZ_template_freq)**2/max(coop_SZ_frequency_factor(freq(i))**2, 1.d-99)
          call coop_dictionary_lookup(this%settings,"i_sensitivity_channel"//COOP_STR_OF(i), sigmaT(i))
          call coop_dictionary_lookup(this%settings,"pol_sensitivity_channel"//COOP_STR_OF(i), sigmapol(i))
       enddo
       tsz_factor = 1.d0/tsz_factor

       call coop_dictionary_lookup(this%settings,"obs_yr", obs_yr)

       beam_fwhm = beam_fwhm*coop_SI_arcmin
       freq = freq*(1.d9*coop_SI_h/coop_SI_kB) !!convert to temperature
       T353 = 353.d0*(1.d9*coop_SI_h/coop_SI_kB)
       sigmaT = sigmaT*1.d-6/COOP_DEFAULT_TCMB/sqrt(obs_yr*3600.*24.*365.2425)
       sigmapol = sigmapol*1.d-6/COOP_DEFAULT_TCMB/sqrt(obs_yr*3600.*24.*365.2425)
       fg_dust_A = fg_dust_A*coop_2pi/80.d0/81.d0 * fg_dust_r/(1.d6*COOP_DEFAULT_TCMB)**2/(T353**fg_dust_beta/(exp(T353/fg_dust_T)-1.d0))**2/80.d0**fg_dust_alpha

       i = 0
       do l = min(this%lmin, this%lmin_pol), max(this%lmax, this%lmax_pol)
          i = i +1
          Cl_dust_residual =  fg_dust_A*dble(l)**fg_dust_alpha/sum(((exp(freq/fg_dust_T)-1.d0)/freq**fg_dust_beta)**2)
          Nl_T = 1.d0/max(sum(exp(-l*(l+1.d0)*(beam_fwhm*coop_sigma_by_fwhm)**2)/sigmaT**2), 1.d-99) + Cl_dust_residual
          Nl_pol = 1.d0/max(sum(exp(-l*(l+1.d0)*(beam_fwhm*coop_sigma_by_fwhm)**2)/sigmapol**2) , 1.d-99) + Cl_dust_residual
          this%nuis(1, i) = fsky
          this%nuis(2, i) =  fsky_pol
          if(l .ge. this%lmin .and. l .le. this%lmax)then
             this%nuis(3, i) = Nl_T
          else
             this%nuis(3, i) = Nl_infinity
          endif
          if(l .ge. this%lmin_pol .and. l .le. this%lmax_pol)then
             this%nuis(4, i) = Nl_pol
          else
             this%nuis(4, i) = Nl_infinity
          endif
          if(l .ge. lmin_tsz .and. l.le.lmax_tsz)then
             this%nuis(5, i) = cl_tsz(l)*tsz_factor
          else
             this%nuis(5, i) = 0.d0
          endif
       enddo
       deallocate(beam_fwhm, sigmaT, sigmapol, freq)
       if(trim(tsz_template).ne."")then
          deallocate(cl_tsz)
          call fpsz%close()
       endif
    case default
       write(*,*) trim(this%genre)
       stop "Error: unknown observation genre"
    end select
    write(*,*) trim(this%filename)//" is sucessfully loaded"

  end subroutine coop_observation_init

  subroutine coop_observation_free(this)
    class(coop_observation)::this
    call this%settings%free()
    call this%paramnames%free()
    COOP_DEALLOC(this%obs)
    COOP_DEALLOC(this%dobs)
    COOP_DEALLOC(this%invcov)
    COOP_DEALLOC(this%nuis)
    COOP_DEALLOC(this%window_Wsq)
    COOP_DEALLOC(this%window_used)
    COOP_DEALLOC(this%window_modes)
    this%n_window = 0
    this%n_bins = 0
    this%n_obs = 0
    this%dim_obs = 0
    this%dim_nuis = 0
    this%genre = ""
    this%name = ""
    this%filename = ""
    this%h_calibrate = 0.7d0
    this%lmin = 0
    this%lmax = -1
    this%lmin_pol = 0
    this%lmax_pol = -1
  end subroutine coop_observation_free


  subroutine coop_fisher_free(this)
    class(coop_fisher)::this
    COOP_INT::i, j
    call this%settings%free()
    call this%cosmology%free()
    call this%paramtable%free()
    COOP_DEALLOC(this%params)
    COOP_DEALLOC(this%step1)
    COOP_DEALLOC(this%step2)
    COOP_DEALLOC(this%priors)
    COOP_DEALLOC(this%param_types)
    COOP_DEALLOC(this%ind_slow)
    COOP_DEALLOC(this%init_level)
    COOP_DEALLOC(this%ind_fast)
    COOP_DEALLOC(this%ind_nuis)
    COOP_DEALLOC(this%fisher)
    COOP_DEALLOC(this%cov)
    COOP_DEALLOC(this%ind_used)
    COOP_DEALLOC(this%is_used)
    if(allocated(this%observations))then
       do i=1, this%n_observations
          call this%observations(i)%free()
       enddo
       deallocate(this%observations)
    endif
    this%n_params = 0
    this%n_params_used = 0
    this%n_observations = 0
  end subroutine coop_fisher_free

  subroutine coop_fisher_init(this, filename)
    class(coop_fisher)::this
    COOP_UNKNOWN_STRING::filename  
    type(coop_list_string)::ls
    type(coop_list_real)::lr
    logical::success
    COOP_INT::i, j, ip
    call this%free()
    call coop_load_dictionary(filename, this%settings)
    call coop_dictionary_lookup(dict = this%settings, key="n_params", val = this%n_params)
    allocate(this%params(this%n_params), this%step1(this%n_params), this%step2(this%n_params), this%priors(this%n_params), this%fisher(this%n_params, this%n_params),  this%cov(this%n_params, this%n_params), this%is_used(this%n_params), this%param_types(this%n_params), this%init_level(this%n_params))
    this%init_level = 0
    this%max_init_level = 0
    this%fisher = 0.d0
    this%cov = 0.d0
    this%is_used = .false.
    call coop_dictionary_lookup(this%settings, "param_names", ls)
    if(ls%n .ne. this%n_params) stop "Error in fisher_init: size of param_names does not equal to n_params"

    do i = 1, this%n_params
       call coop_dictionary_lookup(this%settings,  "param["//trim(ls%element(i))//"]", lr)
       select case(lr%n)
       case(1)
          this%params(i) = lr%element(1)
          this%step1(i) = 0.d0
          this%step2(i) = 0.d0
          this%priors(i) = 0.d0
       case(2)
          this%params(i) = lr%element(1)
          this%step1(i) = lr%element(2)
          this%step2(i) = -this%step1(i)
          this%priors(i) = abs(this%step1(i))*1.d5 !! making this finite rather than infinity helps to beat down the numeric instability from round-off errors.
       case(3)
          this%params(i) = lr%element(1)
          this%step1(i) = lr%element(2)
          this%step2(i) = lr%element(3)
          this%priors(i) = max(abs(this%step1(i)), abs(this%step2(i)))*1.d5 !! making this finite rather than infinity helps to beat down the numeric instability from round-off errors.
       case(4)
          this%params(i) = lr%element(1)
          this%step1(i) = lr%element(2)
          this%step2(i) = lr%element(3)
          this%priors(i) = lr%element(4)
       case default
          write(*,*)  "param["//trim(ls%element(i))//"] seems to be not right, the format is:"
          write(*,*)  "param["//trim(ls%element(i))//"] = fiducial step1 step2 prior"
          stop
       end select
       if(this%priors(i).ne.0.d0 .and. ((this%step1(i).eq.0.d0 .or. this%step2(i).eq.0.d0) .or. abs(this%step1(i) - this%step2(i)) .lt. abs(this%step1(i))*1.d-2))then
          write(*,*)  "param["//trim(ls%element(i))//"] seems to be not right, the format is:"
          write(*,*)  "param["//trim(ls%element(i))//"] = fiducial step1 step2 prior"
          stop
       endif
       call this%paramtable%insert(trim(ls%element(i)), this%params(i), overwrite = .false.)
    enddo
    call ls%free()
    call lr%free()

    this%param_types = coop_parameter_type_nuis
    !!the parameters that require recomputing cosmological perturbations
    call coop_dictionary_lookup(this%settings,"params_slow", ls) !!
    this%n_slow = ls%n
    allocate(this%ind_slow(ls%n))
    do i = 1, ls%n
       j = this%paramtable%index(trim(ls%element(i)))
       if(j.eq.0)then
          write(*,*) trim(ls%element(i))//" appears in param_slow but not in params"
          stop
       endif
       if(this%param_types(j).ne. coop_parameter_type_nuis)then
          write(*,*) "Error: "//trim(ls%element(i))//" duplicated in params_slow"
          stop
       endif
       this%ind_slow(i) = j
       this%param_types(j) = coop_parameter_type_slow
    enddo
    call ls%free()

    !!the parameters that require updating primordial power spectrum
    call coop_dictionary_lookup(this%settings,"params_fast", ls)
    this%n_fast = ls%n
    allocate(this%ind_fast(ls%n))
    do i = 1, ls%n
       j = this%paramtable%index(trim(ls%element(i)))
       if(j.eq.0)then
          write(*,*) trim(ls%element(i))//" appears in param_fast but not in params"
          stop
       endif
       if(this%param_types(j).ne. coop_parameter_type_nuis)then
          write(*,*) "Error: "//trim(ls%element(i))//" duplicated in params_fast"
          stop
       endif
       this%ind_fast(i) = j
       this%param_types(j) = coop_parameter_type_fast
    enddo
    call ls%free()

    j = 0
    this%n_nuis = this%n_params - this%n_slow - this%n_fast
    allocate(this%ind_nuis(this%n_nuis))
    do i=1, this%n_params
       if(this%param_types(i).eq.coop_parameter_type_nuis)then
          j = j + 1
          this%ind_nuis(j) = i
       endif
    enddo

    call coop_dictionary_lookup(dict = this%settings, key="n_observations", val = this%n_observations)
    if(this%n_observations .gt. 0)then
       allocate(this%observations(this%n_observations))
       do i=1, this%n_observations
          call coop_dictionary_lookup(this%settings, "observation"//COOP_STR_OF(i), this%observations(i)%filename)
          call this%observations(i)%init(this%paramtable)
          do j=1, this%observations(i)%paramnames%n
              this%observations(i)%paramnames%val(j) = this%paramtable%index( this%observations(i)%paramnames%key(j) )
             if( this%observations(i)%paramnames%val(j) .eq. 0)then
                !write(*,*) "Error in fisher_init:"
                !write(*,*) "parameter "//trim(this%observations(i)%paramnames%key(j))//" (required by dataset "//trim(this%observations(i)%name)//") is not found."
                cycle
             endif
             this%init_level(this%observations(i)%paramnames%val(j)) = max(  this%init_level(this%observations(i)%paramnames%val(j)), this%observations(i)%init_level)

          enddo
       enddo
    endif

    do i=1, this%n_params
       this%max_init_level = max(this%max_init_level, this%init_level(i))
    enddo
    !!compute the fiducial cosmology
    call coop_dictionary_lookup(this%settings, "w_is_background", this%cosmology%w_is_background, .true.)
    call coop_dictionary_lookup(this%settings, "inflation_consistency", this%cosmology%inflation_consistency, .true.)
    call coop_dictionary_lookup(this%settings, "pp_genre", this%cosmology%pp_genre, COOP_PP_STANDARD)

    call this%cosmology%set_up(this%paramtable, success, level = this%max_init_level)
    if(.not. success)then
       write(*,*) "cannot set up the cosmology, check the parameter range:"
       call this%paramtable%print()
       stop
    endif
    do i = 1, this%n_observations
       call this%observations(i)%get_invcov(this%paramtable, this%cosmology)
    enddo
  end subroutine coop_fisher_init

  subroutine coop_fisher_get_dobs_slow(this, i, cosmology_tmp)
    class(coop_fisher)::this
    COOP_INT::i, iobs, j
    logical::success
    type(coop_real_table)::paramtable_tmp
    type(coop_cosmology_firstorder)::cosmology_tmp
    COOP_REAL,dimension(:,:),allocatable::dobs_tmp
    do iobs = 1, this%n_observations
       j = this%observations(iobs)%paramnames%index(this%paramtable%key(i))
       if(j.ne.0)goto 100
    enddo
    return  !!if no dependence just return
100 continue
    if(this%priors(i).eq. 0.d0)then
       do iobs = 1, this%n_observations
          j = this%observations(iobs)%paramnames%index(this%paramtable%key(i))
          if(j.ne.0)this%observations(iobs)%dobs(:,:,j) = 0.d0
       enddo
       return
    endif
    paramtable_tmp = this%paramtable
    paramtable_tmp%val(i) = this%paramtable%val(i) + this%step1(i)
    !!compute the fiducial cosmology
!    cosmology_tmp = this%cosmology
    call  cosmology_tmp%set_up(paramtable_tmp, success, level = this%init_level(i))
    if(.not. success)then
       write(*,*) "cannot set up the cosmology, check the parameter range:"
       call paramtable_tmp%print()
       stop
    endif
    do iobs = 1, this%n_observations
       j = this%observations(iobs)%paramnames%index(this%paramtable%key(i)) 
       if(j.ne.0)then
          call this%observations(iobs)%get_dobs( this%observations(iobs)%dobs(:,:,j), paramtable_tmp, cosmology_tmp)
       endif
    enddo

    paramtable_tmp%val(i) = this%paramtable%val(i) + this%step2(i)
    call  cosmology_tmp%set_up(paramtable_tmp, success,level = this%init_level(i))
    if(.not. success)then
       write(*,*) "cannot set up the cosmology, check the parameter range:"
       call paramtable_tmp%print()
       stop
    endif
    do iobs = 1, this%n_observations
       j = this%observations(iobs)%paramnames%index(this%paramtable%key(i)) 
       if(j.ne.0)then
          allocate(dobs_tmp(this%observations(iobs)%dim_obs, this%observations(iobs)%n_obs))
          call this%observations(iobs)%get_dobs(dobs_tmp, paramtable_tmp, cosmology_tmp)
          this%observations(iobs)%dobs(:,:,j) = (this%observations(iobs)%dobs(:,:,j) - dobs_tmp * (this%step1(i)/this%step2(i))**2)/(1.d0- this%step1(i)/this%step2(i))
          deallocate(dobs_tmp)
       endif
    enddo
 !   call cosmology_tmp%free()
    call paramtable_tmp%free()
  end subroutine coop_fisher_get_dobs_slow

  subroutine coop_fisher_get_dobs_fast(this, i, cosmology_tmp)
    class(coop_fisher)::this
    COOP_INT::i, iobs, j
    type(coop_real_table)::paramtable_tmp
    type(coop_cosmology_firstorder)::cosmology_tmp
    COOP_REAL,dimension(:,:),allocatable::dobs_tmp
    do iobs = 1, this%n_observations
       j = this%observations(iobs)%paramnames%index(this%paramtable%key(i))
       if(j.ne.0)goto 100
    enddo
    return  !!if no dependence just return
100 continue
    if(this%priors(i).eq. 0.d0)then
       do iobs = 1, this%n_observations
          j = this%observations(iobs)%paramnames%index(this%paramtable%key(i))
          if(j.ne.0)this%observations(iobs)%dobs(:,:,j) = 0.d0
       enddo
       return
    endif
    paramtable_tmp = this%paramtable
!    cosmology_tmp = this%cosmology

    paramtable_tmp%val(i) = this%paramtable%val(i) + this%step1(i)
    if(this%init_level(i).ge. coop_init_level_set_pp)call cosmology_tmp%set_primordial_power(paramtable_tmp) 
    if(this%init_level(i).ge.coop_init_level_set_Cls)call cosmology_tmp%update_Cls(0)

    do iobs = 1, this%n_observations
       j = this%observations(iobs)%paramnames%index(this%paramtable%key(i)) 
       if(j.ne.0)then
          call this%observations(iobs)%get_dobs( this%observations(iobs)%dobs(:,:,j), paramtable_tmp, cosmology_tmp)
       endif
    enddo

    paramtable_tmp%val(i) = this%paramtable%val(i) + this%step2(i)
    if(this%init_level(i).ge. coop_init_level_set_pp)call cosmology_tmp%set_primordial_power(paramtable_tmp)
    if(this%init_level(i).ge.coop_init_level_set_Cls)call cosmology_tmp%update_Cls(0)
    do iobs = 1, this%n_observations
       j = this%observations(iobs)%paramnames%index(this%paramtable%key(i)) 
       if(j.ne.0)then
          allocate(dobs_tmp(this%observations(iobs)%dim_obs, this%observations(iobs)%n_obs))
          call this%observations(iobs)%get_dobs(dobs_tmp, paramtable_tmp, cosmology_tmp)
          this%observations(iobs)%dobs(:,:,j) = (this%observations(iobs)%dobs(:,:,j) - dobs_tmp * (this%step1(i)/this%step2(i))**2)/(1.d0- this%step1(i)/this%step2(i))
          deallocate(dobs_tmp)
       endif
    enddo
 !   call cosmology_tmp%free()
    call paramtable_tmp%free()

  end subroutine coop_fisher_get_dobs_fast

  subroutine coop_fisher_get_dobs_nuis(this, i)
    class(coop_fisher)::this
    COOP_INT::i, iobs, j
    type(coop_real_table)::paramtable_tmp
    COOP_REAL,dimension(:,:),allocatable::dobs_tmp
    do iobs = 1, this%n_observations
       j = this%observations(iobs)%paramnames%index(this%paramtable%key(i))
       if(j.ne.0)goto 100
    enddo
    return  !!if no dependence just return
100 continue
    if(this%priors(i).eq. 0.d0)then
       do iobs = 1, this%n_observations
          j = this%observations(iobs)%paramnames%index(this%paramtable%key(i))
          if(j.ne.0)this%observations(iobs)%dobs(:,:,j) = 0.d0
       enddo
       return
    endif
    paramtable_tmp = this%paramtable
    paramtable_tmp%val(i) = this%paramtable%val(i) + this%step1(i)
    do iobs = 1, this%n_observations
       j = this%observations(iobs)%paramnames%index(this%paramtable%key(i)) 
       if(j.ne.0)then
          call this%observations(iobs)%get_dobs( this%observations(iobs)%dobs(:,:,j), paramtable_tmp, this%cosmology)
       endif
    enddo

    paramtable_tmp%val(i) = this%paramtable%val(i) + this%step2(i)
    do iobs = 1, this%n_observations
       j = this%observations(iobs)%paramnames%index(this%paramtable%key(i)) 
       if(j.ne.0)then
          allocate(dobs_tmp(this%observations(iobs)%dim_obs, this%observations(iobs)%n_obs))
          call this%observations(iobs)%get_dobs(dobs_tmp, paramtable_tmp, this%cosmology)
          this%observations(iobs)%dobs(:,:,j) = (this%observations(iobs)%dobs(:,:,j) - dobs_tmp * (this%step1(i)/this%step2(i))**2)/(1.d0- this%step1(i)/this%step2(i))
          deallocate(dobs_tmp)
       endif
    enddo
    call paramtable_tmp%free()
  end subroutine coop_fisher_get_dobs_nuis

  subroutine coop_fisher_get_fisher(this)
    COOP_INT, parameter::n_threads = 8
    class(coop_fisher)::this
    COOP_INT::i, idata, j, ithread
    type(coop_cosmology_firstorder)::cosmology_tmp(n_threads)
    COOP_REAL, dimension(:,:),allocatable::cov
    this%fisher = 0.d0
    do i = 1, this%n_observations
       this%observations(i)%dobs = 0.d0
    enddo
    write(*,*) "Doing slow parameters."
    !$omp parallel do private(i, ithread)
    do ithread = 1, n_threads
       do i = ithread, this%n_slow, n_threads
          cosmology_tmp(ithread) = this%cosmology
          call coop_fisher_get_dobs_slow(this, this%ind_slow(i), cosmology_tmp(ithread))
       enddo
    enddo
    !$omp end parallel do
    write(*,*) "Slow paramters done."
    write(*,*) "Doing fast parameters."
    !$omp parallel do private(i, ithread)
    do ithread = 1, n_threads
       do i = ithread, this%n_fast, n_threads
          cosmology_tmp(ithread) = this%cosmology
          call coop_fisher_get_dobs_fast(this, this%ind_fast(i), cosmology_tmp(ithread))
       enddo
    enddo
    !$omp end parallel do
    write(*,*) "Fast paramters done."
    write(*,*) "Doing nuisance parameters."
    !$omp parallel do
    do i = 1, this%n_nuis
       call coop_fisher_get_dobs_nuis(this, this%ind_nuis(i))
    enddo
    !$omp end parallel do
    write(*,*) "Nuisance paramters done."
    do i = 1, this%n_observations
       do idata = 1, this%observations(i)%n_obs
          this%fisher(this%observations(i)%paramnames%val(1:this%observations(i)%paramnames%n), this%observations(i)%paramnames%val(1:this%observations(i)%paramnames%n)) &
               = this%fisher(this%observations(i)%paramnames%val(1:this%observations(i)%paramnames%n), this%observations(i)%paramnames%val(1:this%observations(i)%paramnames%n)) &
               + matmul(transpose(this%observations(i)%dobs(:,idata,:)), matmul(this%observations(i)%invcov(:, :, idata), this%observations(i)%dobs(:, idata, :)))
       enddo
    enddo
    !!compute the covariance matrix
    this%n_params_used = 0
    do i=1, this%n_params
       if(this%fisher(i, i) .gt. coop_fisher_tolerance)then
          this%is_used(i) = .true.
          this%n_params_used = this%n_params_used + 1
       else
          if(this%fisher(i,i) .ne. 0.d0)then
             write(*,*) "parameter "//trim(this%paramtable%key(i))//" is unconstrained"
          endif
          this%is_used(i) = .false.
       endif
    enddo


    do i=1, this%n_params
       if(this%is_used(i))this%fisher(i,i) = this%fisher(i,i) + (this%step1(i)/this%priors(i))**2
    enddo

    COOP_DEALLOC(this%ind_used)
    allocate(this%ind_used(this%n_params_used), cov(this%n_params_used, this%n_params_used))
    j = 0
    do i = 1, this%n_params
       if(this%is_used(i))then
          j = j + 1
          this%ind_used(j) = i
       endif
    enddo
    cov = this%fisher(this%ind_used, this%ind_used)
    if(this%n_params_used .gt. 0)then
       call coop_sympos_inverse(this%n_params_used, this%n_params_used, cov)
    endif
    this%cov(this%ind_used, this%ind_used) = cov
    do i = 1, this%n_params
       if(this%is_used(i))then
          this%fisher(i,:) = this%fisher(i,:) /this%step1(i)
          this%fisher(:,i) = this%fisher(:,i) /this%step1(i)
          this%cov(i,:) = this%cov(i,:) *this%step1(i)
          this%cov(:,i) = this%cov(:,i) *this%step1(i)
       endif
    enddo
    deallocate(cov)
  end subroutine coop_fisher_get_fisher



!!for LSS k sampling
!!~ log uniform sampling for k<kl
!!~ uniform sampling for k > kl
!! with smooth transition at k~ kl
  subroutine coop_lss_k_sampling(kmin, kmax, kl, nk, k, dk)
    COOP_INT::nk
    COOP_REAL::kmin, kmax, kl
    COOP_REAL::k(nk), dk(nk)
    COOP_REAL::kop(nk), kopmin, kopmax, weight, dkop
    COOP_INT::i
    if(kl .le. 0.d0)then
       dk = (kmax-kmin)/nk
       call coop_set_uniform(nk, k, kmin + dk(1)/2.d0, kmax-dk(1)/2.d0)
       return
    endif
    if(kl .gt. 0.99d10)then
       dkop = log(kmax/kmin)/nk
       call coop_set_uniform(nk, k, log(kmin)+dkop/2.d0, log(kmax)-dkop/2.d0)
       k = exp(k)
       dk = dkop*k
       return
    endif
    weight = 1.d0/kl
    kopmin = log(kmin) + kmin/kl
    kopmax = log(kmax) + kmax/kl
    dkop = (kopmax-kopmin)/nk
    call coop_set_uniform(nk, kop, kopmin+dkop/2.d0, kopmax-dkop/2.d0)
    do i=1, nk
       call coop_source_kop2k_noindex(kop(i), weight, k(i))
       dk(i) = dkop*k(i)*kl/(k(i)+kl)
    enddo
  end subroutine coop_lss_k_sampling

  !!Delta T^{SZ} (nu) \propto f(h nu /k_B T_CMB) 
  function coop_SZ_frequency_factor(freq_GHz) result(fnu)
    COOP_REAL::Freq_GHz, fnu, xnu
    COOP_REAL::Freq_piv = 1.d-9*coop_SI_kB/coop_SI_h !!=20.83664   
    xnu = Freq_GHz/Freq_piv/COOP_DEFAULT_TCMB
    if(xnu .gt. 0.01d0)then
       fnu =2.d0- (xnu/2.d0/tanh(xnu/2.d0))
    else
       fnu = 1.d0-xnu**2/12.d0*(1.d0-xnu**2/60.d0*(1.d0-xnu**2/42.d0))
    endif
  end function coop_SZ_frequency_factor

end module coop_fisher_mod
