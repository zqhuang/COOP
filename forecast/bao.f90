module coop_bao_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"


  type coop_bao_point
     COOP_REAL::obs
     COOP_REAL::obs_theory
     COOP_SHORT_STRING::genre
     COOP_REAL::z = 0.d0
     COOP_REAL::dA_unit = 1.d0
     COOP_REAL::H_unit = 1.d0
     COOP_REAL::rs_unit = 1.d0
     COOP_REAL::Dv_unit = 1.d0
   contains
     procedure::get_obs_theory => coop_bao_point_get_obs_theory
  end type coop_bao_point
     
  type coop_bao_object
     COOP_STRING::name
     logical::gaussian_like = .true.
     COOP_INT::nobs = 0
     type(coop_bao_point), dimension(:), allocatable::points
     COOP_REAL, dimension(:,:), allocatable::invcov    
     type(coop_nd_prob)::prob
   contains
     procedure::LogLike => coop_bao_object_Loglike
     procedure::init => coop_bao_object_init
     procedure::free => coop_bao_object_free
  end type coop_bao_object

contains


  subroutine coop_bao_object_free(this)
    class(coop_bao_object)::this
    if(allocated(this%points))deallocate(this%points)
    if(allocated(this%invcov))deallocate(this%invcov)
    call this%prob%free()
    this%nobs = 0
  end subroutine coop_bao_object_free

  subroutine coop_bao_object_init(this, filename)
    class(coop_bao_object)::this
    COOP_UNKNOWN_STRING::filename
    COOP_SHORT_STRING::genre
    COOP_STRING::efile
    type(coop_dictionary)::dict
    COOP_INT::i
    COOP_REAL::zeff, obs, uni
    logical success
    call this%free()
    call coop_load_dictionary(COOP_DATAPATH(filename), dict)
    call coop_dictionary_lookup(dict, "name", this%name)
    call coop_dictionary_lookup(dict, "num_data", this%nobs, 1)
    allocate(this%points(this%nobs), this%invcov(this%nobs, this%nobs))
    call coop_dictionary_lookup(dict, "zeff", zeff, -1.d0)
    if(zeff .gt. 0.d0)then
       this%points%z =  zeff
    else
       do i=1, this%nobs
          call coop_dictionary_lookup(dict, "zeff"//COOP_STR_OF(i), zeff, -1.d0)
          if(zeff .le. 0.d0)goto 100
          this%points(i)%z = zeff
       enddo
    endif
    call coop_dictionary_lookup(dict, "genre", genre)
    if(trim(genre).ne."")then
       this%points%genre = genre
    else
       do i=1, this%nobs
          call coop_dictionary_lookup(dict, "genre"//COOP_STR_OF(i), genre)
          if(trim(genre).eq."")goto 100
          this%points(i)%genre = genre
       enddo
    endif
    
    call coop_dictionary_lookup(dict, "rs_unit", uni, coop_logzero)
    if(uni .lt. coop_logzero)then
       this%points%rs_unit = uni
    else
       do i=1, this%nobs
          call coop_dictionary_lookup(dict, "rs_unit"//COOP_STR_OF(i), uni, coop_logzero)
          if(uni .lt. coop_logzero) &
             this%points(i)%rs_unit = uni
       enddo
       
    endif


    call coop_dictionary_lookup(dict, "H_unit", uni, coop_logzero)
    if(uni .lt. coop_logzero)then
       this%points%H_unit = uni
    else
       do i=1, this%nobs
          call coop_dictionary_lookup(dict, "H_unit"//COOP_STR_OF(i), uni, coop_logzero)
          if(uni .lt. coop_logzero) &
             this%points(i)%H_unit = uni
       enddo
       
    endif


    call coop_dictionary_lookup(dict, "DV_unit", uni, coop_logzero)
    if(uni .lt. coop_logzero)then
       this%points%Dv_unit = uni
    else
       do i=1, this%nobs
          call coop_dictionary_lookup(dict, "DV_unit"//COOP_STR_OF(i), uni, coop_logzero)
          if(uni .lt. coop_logzero) &
             this%points(i)%Dv_unit = uni
       enddo
       
    endif
    

    call coop_dictionary_lookup(dict, "DA_unit", uni, coop_logzero)
    if(uni .lt. coop_logzero)then
       this%points%dA_unit = uni
    else
       do i=1, this%nobs
          call coop_dictionary_lookup(dict, "DA_unit"//COOP_STR_OF(i), uni, coop_logzero)
          if(uni .lt. coop_logzero) &
             this%points(i)%dA_unit = uni
       enddo
       
    endif
    
    
    call coop_dictionary_lookup(dict, "prob_file",efile)
    if(trim(efile).ne."")then
       this%Gaussian_like = .false.
       call this%prob%load(COOP_DATAPATH(efile), form = "TABLE")
    else
       this%Gaussian_like = .true.
       if(this%nobs .eq. 1)then  !!for one observable never need covariance matrix
          call coop_dictionary_lookup(dict, "obs", efile)
          if(trim(efile).eq."")then
             call coop_dictionary_lookup(dict, "obs1", efile)
             if(trim(efile).eq."")goto 100
          endif
          read(efile, *, ERR=100, END=100) this%points(1)%obs, this%invcov(1, 1)
          this%invcov(1, 1) = 1.d0/this%invcov(1,1)**2
       else
          call coop_dictionary_lookup(dict, "invcov_file", efile)
          if(trim(efile).ne."")then  !!covariance matrix
             call coop_import_matrix(COOP_DATAPATH(efile), this%invcov, this%nobs, this%nobs, success)
             if(.not. success) goto 100

             do i=1, this%nobs
                call coop_dictionary_lookup(dict, "obs"//COOP_STR_OF(i), this%points(i)%obs)
             enddo
          else  !!diagonal errors
             this%invcov = 0.d0
             do i=1, this%nobs
                call coop_dictionary_lookup(dict, "obs"//COOP_STR_OF(i), efile)
                read(efile, *, ERR=100, END=100) this%points(i)%obs, this%invcov(i, i)
                this%invcov(i,i) = 1.d0/this%invcov(i,i)**2
             enddo
          endif
       endif
    endif
    return
100    write(*,*) "bao data file "//trim(filename)//" is broken"
       stop
  end subroutine coop_bao_object_init

  function coop_bao_object_loglike(this, cosmology) result(loglike)
    class(coop_bao_object)::this
    COOP_REAL::loglike
    type(coop_cosmology_firstorder),optional::cosmology
    COOP_INT::i
    if(this%nobs.le.0)then
       loglike = 0.d0
       return
    endif
    if(present(cosmology))then
       do i=1, this%nobs
          call this%points(i)%get_obs_theory(cosmology)
       enddo
    endif
    if(this%Gaussian_like)then
       loglike = dot_product( this%points%obs_theory - this%points%obs, matmul(this%invcov, this%points%obs_theory - this%points%obs)) /2.d0
    else
       loglike = this%prob%eval(  this%points%obs_theory )
    endif
  end function coop_bao_object_loglike
  

  subroutine coop_bao_point_get_obs_theory(this, cosmology) 
    class(coop_bao_point)::this
    type(coop_cosmology_firstorder)::cosmology
    select case(trim(this%genre))
    case("DV")
       this%obs_theory = (cosmology%Dv_of_z(this%z)/this%Dv_unit)/(cosmology%r_drag / this%rs_unit )
    case("invDV")
       this%obs_theory = (cosmology%r_drag / this%rs_unit)/(cosmology%Dv_of_z(this%z)/this%Dv_unit)
    case("H")
       this%obs_theory = (cosmology%H_of_z(this%z)/this%H_unit)*(cosmology%r_drag/this%rs_unit)
    case("invH")
       this%obs_theory = (this%H_unit/cosmology%H_of_z(this%z))*(this%rs_unit/cosmology%r_drag)
    case("DA")
       this%obs_theory = (cosmology%dA_of_z(this%z)/this%dA_unit)/(cosmology%r_drag/this%rs_unit)
    case default
       write(*,*) trim(this%genre)
       stop "obs_theory: unknown bao genre"
    end select
  end subroutine coop_bao_point_get_obs_theory

  


  
  
  
end module coop_bao_mod
