Module coop_SNlike_JLA_mod
  !Module for handling JLA supernova data
  !Adapted form the SNLS module
  !The differences between this and the supernova_SNLS.f90 module are:
  ! 1) We use the SALT2 X_1 parameter instead of stretch
  ! 2) Intrinsic dispersion is included in the error reported in data files
  ! 3) Dispersion from lensing is also included
  ! 4) Redshift uncertainty is a global 150km/s
  !
  ! Model:
  !  The model for SN magnitudes is
  !   m_obs = 5*log10( D_L ) - alpha*(stretch-1) + beta*colour + scriptm
  !  scriptm is some combination of the absolute luminosity of SNe Ia and
  !  the Hubble constant -- marginalizing over it is equivalent to
  !  marginalizing over the absolute magnitude.  Because the errors on the
  !  individual points don't depend on scriptm, we can do this.  Since they
  !  do depend on alpha and beta, we can't marginalize over them.    Because
  !  cosmomc returns D_L in Mpc, scriptm is M - 25, where M is the magnitude
  !  of a stretch=1, colour=0 SN Ia in whatever band the mag comes in as.
  !  Here we use a flat prior for one or two scriptm values.
  !  Note that the SN data is independend of H_0 -unless- you use
  !   the absdist file.
  ! Covariance matricies:
  !  This code has support for SN-SN covariance matricies.
  !   We write the final covariance
  !  matrix as:
  !     V = D + V_mm + alpha^2 * V_ss + beta^2 + V_cc + 2 alpha V_ms
  !           - 2 beta V_mc - 2 alpha beta V_sc
  !  where, for example, V_ms is the covariance matrix between the
  !  uncorrected magnitude and the stretch of the SN.  D are the diagonal
  !  terms calculated from
  !     D_ii = sigma_m^2 + alpha^2 sigma_s^2 + beta^2 * sigma_c^2
  !                      + 2 alpha cov_m_s - 2 beta cov_m_c
  !                      - 2 alpha beta cov_s_c + intrinsicdisp^2 +
  !                      (5/log 10)^2 sigma_z^2 / z^2
  !  It may seem a little strange that the diagonal term is split off,
  !  but it is convenient in some circumstances, as it allows for
  !  Sherman-Woodbury inversion.  However, we don't implement that here.
  ! Speed:
  !  One might wonder if it is really necessary to explicitly fit for
  !   alpha and beta.  Can't we just internally marginalize over them?
  !  The short answer is, no, you can't, at least not if you want an
  !   unbiased answer.
  !  The long answer is, yes, sure you can internally marginalize over them.
  !   But doing so correctly is actually slower than fitting for them, so
  !   it isn't a great idea.
  !   There are a few things you might consider trying to do the
  !    internal marginalization:
  !     1) Fixing alpha and beta. This is -wrong- and will both underestimate
  !         and bias your results.  This is the way all previous cosmomc
  !         packages work, and so all those papers are -wrong-.
  !     2) Fixing alpha and beta but including some assumed error on them
  !         to make the errors better.  An improvement, but still wrong
  !         because alpha and beta are correlated with the other parameters.
  !         Of course, if other constraints effectively fix the cosmology,
  !         then this works, but that's equivalent to saying that the SN
  !         data is irrelevant to your fit -- so why are you bothering
  !         anyways.
  !     3) Internally minimizing alpha and beta, then plugging these in
  !         to get the chisq.  This is at least interesting, because
  !         this technique usually works, and would make things much
  !         faster.  Sadly, here it doesn't because
  !         this method only applies if the errors are independent of the
  !         parameters you are marginalizing over.  And the errors do depend
  !         on alpha and beta, so this will give you a biased answer.
  !     4) Explicitly making a grid over alpha and beta, computing the
  !         likelihood for each, and then marginalizing.  This finally
  !         actually works.  But, it turns out to be slower than the
  !         alternative.  To get a good result, you really need to have
  !         your grid be 60x60 or larger.  That means inverting the
  !         systematics covariance matrix (which depends on alpha
  !         and beta) > 60^2 times, and it's about
  !         500x500.  Without SNLS, the slowest step in the likelihood
  !         computation is usually the 3000x3000 inversion of the WMAP
  !         pixel space TT cov matrix.  Matrix inversion is N^3, so
  !         that means that this solution for alpha and beta is
  !         60^2*500^3/3000^3 ~ 17 times slower than the WMAP inversion.
  !         For comparison, fitting for alpha and beta explicitly
  !         slows the code down by about 20% for a typical fit.  So,
  !         you can use this method if you want, but it would be kinda
  !         stupid.
  !     5) Comment on the speed topic by MB. The slow step here is
  !         the recomputation of the cholesky decomposition of the
  !         covariance matrix each time alpha or beta changes
  !         (i.e. each step if you are fitting for alpha and
  !         beta). This is about 150MFLOPs for the JLA sample. The
  !         good way to gain a factor 10 is to make this computation
  !         perturbative with respect to a reference alpha,
  !         beta. First order seems accurate enough even for large
  !         departure from my investigation. I did not implement
  !         this here because all those costs are negligible in
  !         regard to the Planck likelihood evaluations. But that
  !         may become interesting for larger SN samples.
  ! Modification History:
  !  Written by Alex Conley, Dec 2006
  !   aconley, Jan 2007: The OpenMP stuff was causing massive slowdowns on
  !      some processors (ones with hyperthreading), so it was removed
  !   aconley, Jul 2009: Added absolute distance support
  !   aconley, May 2010: Added twoscriptm support
  !   aconley, Apr 2011: Fix some non standard F90 usage.  Thanks to
  !                       Zhiqi Huang for catching this.
  !   aconley, April 2011: zhel, zcmb read in wrong order.  Thanks to
  !                       Xiao Dong-Li and Shuang Wang for catching this
  !   mbetoule, Dec 2013: adaptation to the JLA sample
  !   AL, Mar 2014: updates for latest CosmoMC structure
  !   AL, June 2014: updated JLA_marginalize=T handling so it should work (also default JLA.ini)
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  private

  public::coop_data_JLA

  character, parameter :: coop_data_JLA_uplo = 'U' !For LAPACK  


  !Supernova data type
  TYPE coop_Supernova_sample
     COOP_SHORT_STRING::name
     COOP_INT :: dataset       !Subset identifier if subset dependent intrinsic disp is used     
     COOP_REAL :: zhel, zcmb    !The heliocentric and CMB frame redshifts
     COOP_REAL :: z_var         !The variance of the redshift
     COOP_REAL :: mag           !The K-corrected peak magnitude
     COOP_REAL :: mag_var       !The variance of mag
     COOP_REAL :: stretch       !The light-curve fit stretch parameter
     COOP_REAL :: stretch_var   !The variance in the stretch
     COOP_REAL :: colour        !The colour of the SN
     COOP_REAL :: colour_var    !The variance of colour
     COOP_REAL :: thirdvar      !Third variable for scripm split
     COOP_REAL :: thirdvar_var  !Variance in thirdvar
     COOP_REAL :: cov_mag_stretch !Covariance between mag and stretch
     COOP_REAL :: cov_mag_colour  !Covariance between mag and colour
     COOP_REAL :: cov_stretch_colour !Covariance between stretch and colour
     LOGICAL :: has_absdist = .false.   !This SN has an absolute distance
     COOP_REAL::absdist = 0.d0
  End Type coop_Supernova_sample

  type coop_data_JLA
     logical::marginalize = .false.
     COOP_INT::marge_steps = 7
     COOP_REAL::step_width_alpha =0.003d0
     COOP_REAL::step_width_beta = 0.04d0
     COOP_REAL:: alpha_center =  0.14
     COOP_REAL:: beta_center = 3.123
     COOP_INT :: int_points = 0
     COOP_INT::n = 0
     COOP_SHORT_STRING::name="JLA"
     COOP_REAL::pecz
     COOP_REAL,dimension(:), allocatable :: marge_grid, alpha_grid,beta_grid
     TYPE(coop_Supernova_sample), dimension(:), allocatable::sn
     COOP_REAL, dimension(:),allocatable::pre_vars, a1, a2, lumdists
     logical::diag_errors = .false.
     logical::has_mag_covmat  = .false.
     logical::has_stretch_covmat = .false.
     logical::has_colour_covmat = .false.
     logical::has_mag_stretch_covmat = .false.
     logical::has_mag_colour_covmat = .false.
     logical::has_stretch_colour_covmat = .false.
     logical::alphabeta_covmat = .false.
     COOP_REAL,dimension(:,:),allocatable::mag_covmat, stretch_covmat, colour_covmat, mag_stretch_covmat, mag_colour_covmat, stretch_colour_covmat
     COOP_INT::n_absdist = 0     
     COOP_INT,dimension(:),allocatable::ind_absdist
     COOP_INT::n_disp = 0
     COOP_REAL,dimension(:),allocatable :: intrinsicdisp !In magnitudes
     !Variables having to do with optional two-scripmt fit based
     ! on thirdvar cut
     LOGICAL:: twoscriptmfit = .false. !Carry out two scriptm fit
     LOGICAL :: has_thirdvar = .false. !Data has third variable
     COOP_REAL :: scriptmcut = 10.d0 !Cut in thirdvar between two scriptms
     type(coop_dictionary)::settings
   contains
     procedure:: free => coop_data_JLA_free
     procedure:: read => coop_data_JLA_read
     procedure:: read_sn => coop_data_JLA_read_sn
     procedure:: read_absdist => coop_data_JLA_read_absdist
     procedure:: prep => coop_data_JLA_prep
     procedure:: invertcov => coop_data_JLA_invertcov
     procedure:: alpha_beta_like => coop_data_JLA_alpha_beta_like
  end type coop_data_JLA



contains


  subroutine coop_data_JLA_free(this)
    class(coop_data_JLA)::this
    call this%settings%free()
    if(allocated(this%marge_grid))deallocate(this%marge_grid)
    if(allocated(this%alpha_grid))deallocate(this%alpha_grid)
    if(allocated(this%beta_grid))deallocate(this%beta_grid)
    if(allocated(this%sn))deallocate(this%sn)
    this%n = 0
    if(allocated(this%ind_absdist))deallocate(this%ind_absdist)
    this%n_absdist = 0
    if(allocated(this%intrinsicdisp))deallocate(this%intrinsicdisp)
    this%n_disp = 0
    if(allocated(this%mag_covmat))deallocate(this%mag_covmat)
    if(allocated(this%stretch_covmat))deallocate(this%stretch_covmat)
    if(allocated(this%colour_covmat))deallocate(this%colour_covmat)
    if(allocated(this%mag_stretch_covmat))deallocate(this%mag_stretch_covmat)
    if(allocated(this%mag_colour_covmat))deallocate(this%mag_colour_covmat)
    if(allocated(this%stretch_colour_covmat))deallocate(this%stretch_colour_covmat)
    if(allocated(this%pre_vars))deallocate(this%pre_vars)
    if(allocated(this%lumdists))deallocate(this%lumdists)
    if(allocated(this%a1))deallocate(this%a1)
    if(allocated(this%a2))deallocate(this%a2)
  end subroutine coop_data_JLA_free

  subroutine coop_data_JLA_read(this, filename)
    class(coop_data_JLA)::this
    COOP_UNKNOWN_STRING::filename
    COOP_INT::alpha_i, beta_i, ssq, max_betasq, i
    COOP_STRING::absdist_file, covfile
    COOP_REAL idisp_zero
    call this%free()
    call coop_load_dictionary(filename, this%settings)
    call this%read_sn(this%settings%value("data_file"))
    allocate(this%intrinsicdisp(this%n_disp))
    call coop_dictionary_lookup(this%settings, "name", this%name)
    call coop_dictionary_lookup(this%settings, "absdist_file", absdist_file)
    if(trim(absdist_file).ne."")call this%read_absdist(absdist_file)

    call coop_dictionary_lookup(this%settings, "pecz", this%pecz, 1.d-3)
    call coop_dictionary_lookup(this%settings, "twscriptmfit", this%twoscriptmfit, .false.)
    if(this%twoscriptmfit .and. (.not. this%has_thirdvar))then
       stop "JLA: twoscriptmfit was set but thirdvar information not present"
    endif
    call coop_dictionary_lookup(this%settings, "scriptmcut", this%scriptmcut,10.d0)

    call coop_dictionary_lookup(this%settings, "intrinsicdisp", idisp_zero, 0.13d0)
    do i=1, this%n_disp
       call coop_dictionary_lookup(this%settings, "intrinsicdisp"//COOP_STR_OF(i), this%intrinsicdisp(i), idisp_zero)
    enddo
    call coop_dictionary_lookup(this%settings, "has_mag_covmat", this%has_mag_covmat, .false.)
    call coop_dictionary_lookup(this%settings, "has_stretch_covmat", this%has_stretch_covmat, .false.)
    call coop_dictionary_lookup(this%settings, "has_colour_covmat", this%has_colour_covmat, .false.)
    call coop_dictionary_lookup(this%settings, "has_mag_stretch_covmat", this%has_mag_stretch_covmat, .false.)
    call coop_dictionary_lookup(this%settings, "has_mag_colour_covmat", this%has_mag_colour_covmat, .false.)
    call coop_dictionary_lookup(this%settings, "has_stretch_colour_covmat", this%has_stretch_colour_covmat, .false.)        
    this%alphabeta_covmat = ( this%has_stretch_covmat .OR. this%has_colour_covmat .OR. &
         this%has_mag_stretch_covmat .OR. this%has_mag_colour_covmat .OR. &
         this%has_stretch_colour_covmat )


    !First test for covmat
    if ( this%has_mag_covmat .OR. this%has_stretch_covmat .OR. this%has_colour_covmat .OR. &
         this%has_mag_stretch_covmat .OR. this%has_mag_colour_covmat .OR. &
         this%has_stretch_colour_covmat ) then
       this%diag_errors = .false.

       !Now Read in the covariance matricies
       if (this%has_mag_covmat) THEN
          call coop_dictionary_lookup(this%settings, 'mag_covmat_file', covfile)
          allocate( this%mag_covmat( this%n, this%n ) )
          CALL coop_data_JLA_readcov( covfile, this%mag_covmat, this%n )
       endif
       if (this%has_stretch_covmat) THEN
          call coop_dictionary_lookup(this%settings, 'stretch_covmat_file', covfile)          
          allocate( this%stretch_covmat( this%n, this%n ) )
          CALL coop_data_JLA_readcov( covfile, this%stretch_covmat, this%n )
       endif
       if (this%has_colour_covmat) THEN
          call coop_dictionary_lookup(this%settings, 'colour_covmat_file', covfile)                    
          allocate( this%colour_covmat( this%n, this%n ) )
          CALL coop_data_JLA_readcov( covfile, this%colour_covmat, this%n )
       endif
       if (this%has_mag_stretch_covmat) THEN
          call coop_dictionary_lookup(this%settings, 'mag_stretch_covmat_file', covfile)
          allocate( this%mag_stretch_covmat( this%n, this%n ) )
          CALL coop_data_JLA_readcov( covfile, this%mag_stretch_covmat, this%n )
       endif
       if (this%has_mag_colour_covmat) THEN
          call coop_dictionary_lookup(this%settings, 'mag_colour_covmat_file',covfile)
          allocate( this%mag_colour_covmat( this%n, this%n ) )
          CALL coop_data_JLA_readcov( covfile, this%mag_colour_covmat, this%n )
       endif
       if (this%has_stretch_colour_covmat) THEN
          call coop_dictionary_lookup(this%settings,'stretch_colour_covmat_file',covfile)
          allocate(this%stretch_colour_covmat( this%n, this%n ) )
          CALL coop_data_JLA_readcov( covfile, this%stretch_colour_covmat, this%n )
       endif
    else
       this%diag_errors = .true.
    endif



    call coop_dictionary_lookup(this%settings, "JLA_marginalize", this%marginalize, .false.)
    if(this%marginalize)then
       call coop_dictionary_lookup(this%settings, "JLA_marge_steps", this%marge_steps, 7)
       call coop_dictionary_lookup(this%settings, "JLA_step_width_alpha", this%step_width_alpha, 0.003d0)
       call coop_dictionary_lookup(this%settings, "JLA_step_width_beta", this%step_width_beta, 0.04d0)
       allocate(this%alpha_grid((2*this%marge_steps+1)**2), this%beta_grid((2*this%marge_steps+1)**2))
       ssq  = this%marge_steps**2
       this%int_points = 0       
       do alpha_i = -this%marge_steps, this%marge_steps
          max_betasq = ssq - alpha_i**2
          do beta_i = -this%marge_steps, this%marge_steps
             if( beta_i **2 .le. max_betasq)then
                this%int_points = this%int_points + 1
                this%alpha_grid(this%int_points) = this%alpha_center + alpha_i*this%step_width_alpha
                this%beta_grid(this%int_points) = this%beta_center + beta_i*this%step_width_beta               
             endif
          enddo
       enddo
       allocate(this%marge_grid(this%int_points))
    endif
    call this%prep()

  end subroutine coop_data_JLA_read

  subroutine coop_data_jla_read_sn(this, filename)
    class(coop_data_jla)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file)::fp
    COOP_LONG_STRING::line
    INTEGER:: i, ntot
    REAL :: dz, dm, ds, dc, dt
    type(coop_list_string)::lstr
    ntot = 0
    call fp%open(filename, "r")
    do while(fp%read_string(line))
       ntot = ntot + 1
    enddo
    call fp%close()
    allocate(this%sn(ntot))
    this%n = 0
    call fp%open(filename, "r")
    this%has_thirdvar = .false.
    this%sn%has_absdist = .false.
    this%sn%dataset = 1    
    this%n_disp = 1    
    do while(fp%read_string(line))
       !We have a few formats to try.  First, there is the very
       ! long format with thirdvar and dataset.  If that fails,
       ! try without data set.  If that fails, try without
       ! thirdvar but with dataset, and finally with neither

       !A further complication is that if one line has thirdvar,
       ! they had better all have them or else ugliness will probably
       ! result
       this%n = this%n+1
       call coop_string_to_list(line, lstr)
       select case(lstr%n)
       case(16)
          read(line, *) &
               this%sn(this%n)%name, this%sn(this%n)%zcmb, this%sn(this%n)%zhel, dz, &
               this%sn(this%n)%mag, dm, this%sn(this%n)%stretch, ds, &
               this%sn(this%n)%colour,dc,this%sn(this%n)%thirdvar, dt,&
               this%sn(this%n)%cov_mag_stretch, this%sn(this%n)%cov_mag_colour,this%sn(this%n)%cov_stretch_colour, this%sn(this%n)%dataset
          this%has_thirdvar = .true.
       case(15)
          READ (line, *) &
               this%sn(this%n)%name, this%sn(this%n)%zcmb, this%sn(this%n)%zhel, dz, &
               this%sn(this%n)%mag, dm, this%sn(this%n)%stretch, ds, &
               this%sn(this%n)%colour,dc,this%sn(this%n)%thirdvar,dt,&
               this%sn(this%n)%cov_mag_stretch, this%sn(this%n)%cov_mag_colour,this%sn(this%n)%cov_stretch_colour
          this%has_thirdvar = .true.
       case(14)
          read(line, *) &
               this%sn(this%n)%name, this%sn(this%n)%zcmb, this%sn(this%n)%zhel, dz, &
               this%sn(this%n)%mag, dm, this%sn(this%n)%stretch, ds, &
               this%sn(this%n)%colour,dc,this%sn(this%n)%cov_mag_stretch, this%sn(this%n)%cov_mag_colour, &
               this%sn(this%n)%cov_stretch_colour, this%sn(this%n)%dataset
          if(this%n .gt. 1 .and. this%has_thirdvar)then
             stop "problem with third var"
          else
             this%sn(this%n)%thirdvar = 0.d0
             this%sn(this%n)%dataset = 0
             dt = 0.d0
          endif
       case default
          stop "JLA data file has unexpected format"
       end select
       this%sn(this%n)%z_var = dz**2
       this%sn(this%n)%mag_var = dm**2
       this%sn(this%n)%stretch_var = ds**2
       this%sn(this%n)%colour_var = dc**2
       this%sn(this%n)%thirdvar_var = dt**2
       this%n_disp = max(this%sn(this%n)%dataset, this%n_disp)

    Enddo
    call fp%close()
  End subroutine coop_data_jla_read_sn

  subroutine coop_data_JLA_read_absdist(this, filename)
    class(coop_data_JLA)::this
    COOP_UNKNOWN_STRING::filename

  end subroutine coop_data_JLA_read_absdist


  subroutine coop_data_JLA_prep(this)
    class(coop_data_JLA)::this
    COOP_REAL, PARAMETER :: zfacsq = 25.0/(LOG(10.0))**2    
    COOP_INT i
    LOGICAL :: has_A1, has_A2    
    if(this%n .lt. 1) stop "No JLA data read"
    if(minval(this%sn%dataset) .le. 0) stop "dataset number must be > 0"
    allocate(this%pre_vars(this%n))
    !$omp parallel do
    do i=1, this%n
       this%pre_vars(i) = this%sn(i)%mag_var + this%intrinsicdisp(this%sn(i)%dataset)**2       
       if(.not. this%sn(i)%has_absdist)then
          this%pre_vars(i) = this%pre_vars(i) &
               + zfacsq * this%pecz**2 * ( (1.d0 + this%sn(i)%zcmb)/ &
               (this%sn(i)%zcmb*(1+0.5*this%sn(i)%zcmb)) )**2
       endif
    enddo
    !$omp end parallel do
    allocate(this%lumdists(this%n)) 

    IF (this%twoscriptmfit) THEN
       ALLOCATE( this%A1(this%n), this%A2(this%n) )
       has_A1 = .TRUE.
       has_A2 = .FALSE.
       !Assign A1 and A2 as needed
       DO i=1,this%n
          IF (this%sn(i)%thirdvar .LE. this%scriptmcut ) THEN
             this%A1(i) = 1.d0
             this%A2(i) = 0.d0
             has_A1 = .TRUE.
          ELSE
             this%A1(i) = 0.d0
             this%A2(i) = 1.d0
             has_A2 = .TRUE.
          END IF
       END DO

       IF (.NOT. has_A1) THEN
          !Swap
          this%A1 = this%A2
          this%A2(:) = 0.d0
          this%twoscriptmfit = .FALSE.
          has_A1 = .TRUE.
          has_A2 = .FALSE.
       ENDIF

       IF (.NOT. has_A2)  &
            this%twoscriptmfit = .false.
    ENDIF
  end subroutine coop_data_JLA_prep

  subroutine coop_data_JLA_readcov(filename, mat, n)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: n
    COOP_REAL, INTENT(OUT) :: mat(n,n)
    type(coop_file)::fp
    INTEGER :: j,k, file_unit, nfile
    COOP_REAL :: tmp
    call fp%open(filename, "r")
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
  End subroutine coop_data_JLA_readcov

  SUBROUTINE coop_data_jla_invertcov(this, invcovmat, alpha, beta, status )
    class(coop_data_jla)::this
    CHARACTER(LEN=*), PARAMETER :: cholerrfmt = &
         '("Error computing cholesky decomposition for ",F6.3,2X,F6.3)'
    CHARACTER(LEN=*), PARAMETER :: cholinvfmt = &
         '("Error inverting cov matrix for ",F6.3,2X,F6.3)'
    CHARACTER(LEN=*), PARAMETER :: cholsolfmt = &
         '("Error forming inv matrix product for ",F6.3,2X,F6.3)'


    COOP_REAL, INTENT(IN) :: alpha, beta
    INTEGER, INTENT(INOUT) :: status
    COOP_REAL :: invcovmat(:,:)

    INTEGER :: I
    COOP_REAL:: alphasq, betasq, alphabeta

    alphasq = alpha * alpha
    betasq = beta * beta
    alphabeta = alpha * beta

    IF (this%diag_errors) STOP 'Error -- asking to invert with diagonal errors'

    !Build the covariance matrix, then invert it
    IF (this%has_mag_covmat) THEN
       invcovmat = this%mag_covmat
    ELSE
       invcovmat = 0.d0
    END IF
    IF (this%has_stretch_covmat) invcovmat = invcovmat + &
         alphasq * this%stretch_covmat
    IF (this%has_colour_covmat) invcovmat = invcovmat + &
         betasq * this%colour_covmat
    IF (this%has_mag_stretch_covmat) invcovmat = invcovmat + 2.0 * alpha * this%mag_stretch_covmat
    IF (this%has_mag_colour_covmat) invcovmat = invcovmat - 2.0 * beta * this%mag_colour_covmat
    IF (this%has_stretch_colour_covmat) invcovmat = invcovmat - 2.0 * alphabeta * this%stretch_colour_covmat

    !Update the diagonal terms
    do i =1, this%n
       invcovmat(i, i) = invcovmat(i, i) + this%pre_vars(I) &
            + alphasq * this%sn(i)%stretch_var &
            + betasq  * this%sn(i)%colour_var &
            + 2.0 * alpha * this%sn(i)%cov_mag_stretch &
            - 2.0 * beta * this%sn(i)%cov_mag_colour &
            - 2.0 * alphabeta * this%sn(i)%cov_stretch_colour
    enddo

    !Factor into Cholesky form, overwriting the input matrix
    CALL DPOTRF(coop_data_JLA_uplo,this%n,invcovmat,this%n,status)
    IF ( status .NE. 0 ) THEN
       WRITE(*,cholerrfmt) alpha, beta
       RETURN
    END IF

    !Now invert
    !If we could get away with the relative chisquare
    ! this could be done faster and more accurately
    ! by solving the system V*x = diffmag for x to get
    ! V^-1 * diffmag.  But, with the introduction of alpha, beta
    ! this _doesn't_ work, so we need the actual elements of
    ! the inverse covariance matrix.  The point is that the
    ! amarg_E parameter depends on the sum of the elements of
    ! the inverse covariance matrix, and therefore is different
    ! for different values of alpha and beta.
    !Note that DPOTRI only makes half of the matrix correct,
    ! so we have to be careful in what follows
    CALL DPOTRI(coop_data_JLA_uplo,this%n,invcovmat,this%n,status)
    IF ( status .NE. 0 ) THEN
       WRITE(*,cholinvfmt) alpha, beta

       RETURN
    END IF

  End SUBROUTINE coop_data_jla_invertcov


  FUNCTION  coop_data_JLA_alpha_beta_like(this, alpha, beta,  lumdists) result(JLA_alpha_beta_like)
    COOP_REAL, parameter :: coop_inv2pi = 1.d0/coop_2pi    
    class(coop_data_JLA)::this
    COOP_REAL :: JLA_alpha_beta_like
    CHARACTER(LEN=*), PARAMETER :: invfmt = &
         '("Error inverting cov matrix for ",F6.3,2X,F6.3)'

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
    COOP_REAL, allocatable :: invcovmat(:,:)

    allocate(invcovmat(this%n,this%n))

    alphasq   = alpha*alpha
    betasq    = beta*beta
    alphabeta = alpha*beta

    !We want to get a first guess at scriptm to improve the
    ! numerical precision of the results.  We'll do this ignoring
    ! the covariance matrix and ignoring if there are two scriptms
    ! to deal with
    invvars = 1.0 / ( this%pre_vars + alphasq * this%sn%stretch_var &
         + betasq * this%sn%colour_var &
         + 2.0 * alpha * this%sn%cov_mag_stretch &
         - 2.0 * beta * this%sn%cov_mag_colour &
         - 2.0 * alphabeta * this%sn%cov_stretch_colour )

    wtval = SUM( invvars )
    estimated_scriptm= SUM( (this%sn%mag - lumdists)*invvars ) / wtval
    diffmag = this%sn%mag - lumdists + alpha*( this%sn%stretch ) &
         - beta * this%sn%colour - estimated_scriptm

    IF ( this%diag_errors ) THEN
       amarg_A = SUM( invvars * diffmag**2 )
       IF ( this%twoscriptmfit ) THEN
          amarg_B = SUM( invvars * diffmag * this%A1)
          amarg_C = SUM( invvars * diffmag * this%A2)
          amarg_D = 0.0
          amarg_E = DOT_PRODUCT( invvars, this%A1 )
          amarg_F = DOT_PRODUCT( invvars, this%A2 )
       ELSE
          amarg_B = SUM( invvars * diffmag )
          amarg_E = wtval
       ENDIF
    ELSE
       !Unfortunately, we actually need the covariance matrix,
       ! and can't get away with evaluating terms this
       ! V^-1 * x = y by solving V * y = x.  This costs us in performance
       ! and accuracy, but such is life
      call this%invertcov(invcovmat, alpha,beta,status)
       IF (status .NE. 0) THEN
          WRITE (*,invfmt) alpha,beta
          JLA_alpha_beta_like = coop_logZero
          RETURN
       ENDIF

       !Now find the amarg_ parameters
       !We re-use the invvars variable to hold the intermediate product
       !which is sort of naughty
       ! invvars = V^-1 * diffmag (invvars = 1.0*invcovmat*diffmag+0*invvars)
       CALL DSYMV(coop_data_JLA_uplo,this%n,1.0d0,invcovmat,this%n,diffmag,1,0.0d0,invvars,1)

       amarg_A = DOT_PRODUCT( diffmag, invvars ) ! diffmag*V^-1*diffmag

       IF (this%twoscriptmfit) THEN
          amarg_B = DOT_PRODUCT( invvars, this%A1 ) !diffmag*V^-1*A1
          amarg_C = DOT_PRODUCT( invvars, this%A2 ) !diffmag*V^-1*A2

          !Be naughty again and stick V^-1 * A1 in invvars
          CALL DSYMV(coop_data_JLA_uplo,this%n,1.0d0,invcovmat,this%n,this%A1,1,0.0d0,invvars,1)
          amarg_D = DOT_PRODUCT( invvars, this%A2 ) !A2*V^-1*A1
          amarg_E = DOT_PRODUCT( invvars, this%A1 ) !A1*V^-1*A1
          ! now V^-1 * A2
          CALL DSYMV(coop_data_JLA_uplo,this%n,1.0d0,invcovmat,this%n,this%A2,1,0.0d0,invvars,1)
          amarg_F = DOT_PRODUCT( invvars, this%A2 ) !A2*V^-1*A2
       ELSE
          amarg_B = SUM( invvars ) !GB = 1 * V^-1 * diffmag
          !amarg_E requires a little care since only half of the
          !matrix is correct if we used the full covariance matrix
          ! (which half depends on COOP_DATA_JLA_UPLO)
          !GE = 1 * V^-1 * 1
          amarg_C = 0.d0
          amarg_D = 0.d0
          amarg_E = 0.d0
          amarg_F = 0.d0
          IF ( coop_data_JLA_uplo .EQ. 'U' ) THEN
             DO I=1,this%n
                amarg_E = amarg_E + invcovmat(I,I) + 2.0d0*SUM( invcovmat( 1:I-1, I ) )
             END DO
          ELSE
             DO I=1,this%n
                amarg_E = amarg_E + invcovmat(I,I) + 2.0d0*SUM( invcovmat( I+1:this%n, I ) )
             END DO
          END IF
       ENDIF
    END IF

    IF (this%twoscriptmfit) THEN
       !Messy case
       tempG = amarg_F - amarg_D*amarg_D/amarg_E;
       IF (tempG .LE. 0.0) THEN
          WRITE(*,*) "Twoscriptm assumption violation"
          STOP
       ENDIF
       chisq = amarg_A + LOG( amarg_E*coop_inv2pi ) + &
            LOG( tempG * coop_inv2pi ) - amarg_C*amarg_C/tempG - &
            amarg_B*amarg_B*amarg_F / ( amarg_E*tempG ) + 2.0*amarg_B*amarg_C*amarg_D/(amarg_E*tempG )
    ELSE
       chisq = amarg_A + LOG( amarg_E*coop_inv2pi ) - amarg_B**2/amarg_E
    ENDIF
    JLA_alpha_beta_like = chisq / 2  !Negative log likelihood


  end FUNCTION  coop_data_JLA_alpha_beta_like


End Module coop_SNlike_JLA_mod
