Module coop_SNlike_evolve_mod
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  private

  public::coop_data_sn

  character, parameter :: coop_data_sn_uplo = 'U' !For LAPACK  


  !Supernova data type
  TYPE coop_Supernova_sample
     COOP_SHORT_STRING::name
     COOP_REAL :: zhel, zcmb, z_var    !The heliocentric and CMB frame redshifts
     COOP_REAL :: mag, mag_var           !The K-corrected peak magnitude
  End Type coop_Supernova_sample

  type coop_data_sn
     COOP_INT::n = 0
     COOP_SHORT_STRING::name="Pantheon"
     COOP_REAL::pecz
     TYPE(coop_Supernova_sample), dimension(:), allocatable::sn
     COOP_REAL, dimension(:),allocatable::pre_vars, lumdists
     logical::has_covmat  = .false.
     COOP_REAL::intrdisp
     COOP_REAL,dimension(:,:),allocatable::covmat, invcov
     type(coop_dictionary)::settings
   contains
     procedure:: free => coop_data_sn_free
     procedure:: read => coop_data_sn_read
     procedure:: read_sn => coop_data_sn_read_sn
     procedure:: prep => coop_data_sn_prep
     procedure:: invertcov => coop_data_sn_invertcov
     procedure:: likelihood => coop_data_sn_likelihood
  end type coop_data_sn



contains


  subroutine coop_data_sn_free(this)
    class(coop_data_sn)::this
    call this%settings%free()
    if(allocated(this%sn))deallocate(this%sn)
    this%n = 0
    if(allocated(this%covmat))deallocate(this%covmat)
    if(allocated(this%invcov))deallocate(this%invcov)
    if(allocated(this%pre_vars))deallocate(this%pre_vars)
    if(allocated(this%lumdists))deallocate(this%lumdists)
  end subroutine coop_data_sn_free

  subroutine coop_data_sn_read(this, filename)
    class(coop_data_sn)::this
    COOP_UNKNOWN_STRING::filename
    COOP_INT::alpha_i, beta_i, ssq, max_betasq, i
    COOP_STRING::covfile
    call this%free()
    call coop_load_dictionary(COOP_DATAPATH(filename), this%settings)
    call this%read_sn(this%settings%value("data_file"))
    call coop_dictionary_lookup(this%settings, "name", this%name)
    call coop_dictionary_lookup(this%settings, "pecz", this%pecz, 1.d-3)
    write(*,"(A)") "Using peculiar velocity: v = "//COOP_STR_OF(nint(this%pecz*3.e5*coop_sqrt3))//" km/s"
    call coop_dictionary_lookup(this%settings, "intrinsicdisp", this%intrdisp, 0.13d0)
    write(*,"(A)") "Using intrinsic dispersion: "//COOP_STR_OF(nint(this%intrdisp*1000)/1000)
    call coop_dictionary_lookup(this%settings, "has_covmat", this%has_covmat, .false.)
    !Now Read in the covariance matricies
    if (this%has_covmat) then
       allocate(this%invcov(this%n, this%n), this%covmat(this%n, this%n))       
       call coop_dictionary_lookup(this%settings, 'covmat_file', covfile)
       CALL coop_data_sn_readcov( covfile, this%covmat, this%n )
    endif
    call this%prep()
  end subroutine coop_data_sn_read

  subroutine coop_data_sn_read_sn(this, filename)
    class(coop_data_sn)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file)::fp
    COOP_LONG_STRING::line
    INTEGER:: i, ntot
    REAL :: dz, dm, ds, dc, dt
    type(coop_list_string)::lstr
    ntot = 0
    call fp%open(COOP_DATAPATH(filename), "r")
    do while(fp%read_string(line))
       ntot = ntot + 1
    enddo
    call fp%close()
    allocate(this%sn(ntot))
    this%n = 0
    call fp%open(COOP_DATAPATH(filename), "r")
    do while(fp%read_string(line))
       this%n  = this%n + 1
       read(line, *) &
            this%sn(this%n)%name, this%sn(this%n)%zcmb, this%sn(this%n)%zhel, dz, &
            this%sn(this%n)%mag, dm
       this%sn(this%n)%z_var = dz**2
       this%sn(this%n)%mag_var = dm**2
    enddo
    call fp%close()
  End subroutine coop_data_sn_read_sn

  subroutine coop_data_sn_prep(this)
    class(coop_data_sn)::this
    COOP_REAL, PARAMETER :: zfacsq = 25.0/(LOG(10.0))**2    
    COOP_INT i, status
    if(this%n .lt. 1) stop "Found no SNe data"
    allocate(this%pre_vars(this%n), this%lumdists(this%n))
    !$omp parallel do
    do i=1, this%n
       this%pre_vars(i) = this%sn(i)%mag_var + this%intrdisp**2 &
            + zfacsq * this%pecz**2 * ( (1.d0 + this%sn(i)%zcmb)/ &
            (this%sn(i)%zcmb*(1+0.5*this%sn(i)%zcmb)) )**2
    enddo
    !$omp end parallel do
    status = 0
    call this%invertcov(this%invcov, status)
  end subroutine coop_data_sn_prep

  subroutine coop_data_sn_readcov(filename, mat, n)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: n
    COOP_REAL, INTENT(OUT) :: mat(n,n)
    type(coop_file)::fp
    INTEGER :: j,k, file_unit, nfile
    COOP_REAL :: tmp
    call fp%open(COOP_DATAPATH(filename), "r")
    READ (fp%unit, '(I5)', end=200, ERR=100) nfile
    if (nfile .ne. n) THEN
       WRITE (*,'("For file ",A," expected size ",I5," got ",I5)') &
            TRIM(filename), n, nfile
       STOP
    endif

    DO j=1,n
       READ (fp%unit,*, end = 200, err=100) mat(j,1:n)
    endDO

    GOTO 120

100 REWIND(fp%unit)  !Try other possible format
    READ (fp%unit, '(I5)', end=200, ERR=100) nfile

    DO j=1,n
       DO k=1,n
          READ (fp%unit,*, end = 200) mat(j,k)
       end DO
    enddo

120 READ (fp%unit,*, err = 150, end =150) tmp
    GOTO 200

150 call fp%close()
    RETURN

200 WRITE (*,*) 'matrix file '//trim(filename)//' is the wrong size'
    WRITE (*,'("Expected: ",I5," by ",I5)') n,n
    STOP

500 WRITE (*,*) 'Failed to open cov matrix file ' // TRIM(filename)
    STOP
  End subroutine coop_data_sn_readcov

  SUBROUTINE coop_data_sn_invertcov(this, invcovmat, status )
    class(coop_data_sn)::this
    INTEGER, INTENT(INOUT) :: status
    COOP_REAL :: invcovmat(:,:)
    INTEGER :: i
    IF (.not. this%has_covmat) return
    invcovmat = this%covmat
    !Update the diagonal terms
    do i =1, this%n
       invcovmat(i, i) = invcovmat(i, i) + this%pre_vars(I)
    enddo
#ifdef HAS_LAPACK
    CALL DPOTRF(coop_data_sn_uplo,this%n,invcovmat,this%n,status)
    IF ( status .NE. 0 ) THEN
       WRITE(*, *) "Error cannot invert the cov matrix"
       stop
    END IF    
    CALL DPOTRI(coop_data_sn_uplo,this%n,invcovmat,this%n,status)
    IF ( status .NE. 0 ) THEN
       WRITE(*, *) "Error cannot invert the cov matrix"
       stop
    END IF    
#else
    stop "For SN likelihood you need to link COOP to Lapack"
#endif    
  End SUBROUTINE coop_data_sn_invertcov


  FUNCTION  coop_data_sn_likelihood(this, lumdists) result(like)
    class(coop_data_sn)::this
    COOP_REAL :: like
    INTEGER :: i, status
    COOP_REAL :: lumdists(this%n)
    COOP_REAL :: alpha, beta
    !We form an estimate for scriptm to improve numerical
    ! accuracy in our marginaliztion
    COOP_REAL :: estimated_scriptm, wtval
    COOP_REAL :: chisq !Utility variables
    COOP_REAL :: alphasq, betasq, alphabeta !More utility variables
    COOP_REAL :: amarg_A, amarg_B, amarg_C
    COOP_REAL :: amarg_D, amarg_E, amarg_F, tempG !Marginalization params
    COOP_REAL :: diffmag(this%n),invvars(this%n)
    invvars = 1.0 / this%pre_vars 
    wtval = SUM( invvars )
    estimated_scriptm= SUM( (this%sn%mag - lumdists)*invvars ) / wtval
    diffmag = this%sn%mag - lumdists - estimated_scriptm

    if ( this%has_covmat) then
#ifdef HAS_LAPACK       
       call DSYMV(coop_data_sn_uplo,this%n,1.0d0,this%invcov,this%n,diffmag,1,0.0d0,invvars,1)
#else
       stop "For SNe likelihood you need to link COOP to Lapack"
#endif    
       amarg_A = DOT_PRODUCT( diffmag, invvars ) ! diffmag*V^-1*diffmag
       amarg_B = SUM( invvars )
       amarg_E = 0.d0
       if ( coop_data_sn_uplo .EQ. 'U' ) then
          do i=1,this%n
             amarg_E = amarg_E + this%invcov(I,I) + 2.0d0*SUM( this%invcov( 1:I-1, I ) )
          enddo
       else
          do I=1,this%n
             amarg_E = amarg_E + this%invcov(I,I) + 2.0d0*SUM( this%invcov( I+1:this%n, I ) )
          enddo
       endif
    else
       amarg_A = SUM( invvars * diffmag**2 )
       amarg_B = SUM( invvars * diffmag )
       amarg_E = wtval
    endif
    chisq = amarg_A + LOG( amarg_E/coop_2pi ) - amarg_B**2/amarg_E
    like = chisq / 2  !Negative log likelihood
  end FUNCTION  coop_data_sn_likelihood


End Module coop_SNlike_evolve_mod
