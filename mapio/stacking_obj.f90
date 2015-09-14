module coop_stacking_mod
  use coop_wrapper_firstorder
#ifdef HAS_HEALPIX
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  use udgrade_nr
#endif
  implicit none
#include "constants.h"

  COOP_REAL,parameter::coop_stacking_max_threshold = 10.d0
  COOP_INT,parameter::coop_stacking_genre_Null = 0
  COOP_INT,parameter::coop_stacking_genre_Imax = 1
  COOP_INT,parameter::coop_stacking_genre_Imin = 2
  COOP_INT,parameter::coop_stacking_genre_Saddle = 3  
  COOP_INT,parameter::coop_stacking_genre_Imax_Oriented = 4
  COOP_INT,parameter::coop_stacking_genre_Imin_Oriented = 5
  COOP_INT,parameter::coop_stacking_genre_Saddle_Oriented = 6  
  COOP_INT,parameter::coop_stacking_genre_Lmax = 7
  COOP_INT,parameter::coop_stacking_genre_Lmin = 8
  COOP_INT,parameter::coop_stacking_genre_Lmax_Oriented = 9
  COOP_INT,parameter::coop_stacking_genre_Lmin_Oriented = 10
  COOP_INT,parameter::coop_stacking_genre_Pmax_Oriented = 11
  COOP_INT,parameter::coop_stacking_genre_Pmin_Oriented = 12
  COOP_INT,parameter::coop_stacking_genre_Random_Hot  = 13
  COOP_INT,parameter::coop_stacking_genre_Random_Hot_Oriented  = 14
  COOP_INT,parameter::coop_stacking_genre_Random_Cold  = 15
  COOP_INT,parameter::coop_stacking_genre_Random_Cold_Oriented  = 16
  COOP_INT,parameter::coop_stacking_genre_col_oriented = 17

  type coop_stacking_options
     logical::mask_int = .false.
     logical::mask_pol = .false.
     COOP_INT::genre = coop_stacking_genre_Null
     COOP_INT::nmaps = 0
     COOP_INT::nside = 0
     COOP_INT::index_I = 0
     COOP_INT::index_Q = 0
     COOP_INT::index_U = 0
     COOP_INT::index_L = 0
     COOP_STRING::caption = ""
     COOP_REAL::sigma_I = 0.d0
     COOP_REAL::sigma_L = 0.d0
     COOP_REAL::sigma_P = 0.d0
     COOP_SINGLE::I_lower = -1.e20
     COOP_SINGLE::I_upper = 1.e20
     COOP_SINGLE::L_lower = -1.e20
     COOP_SINGLE::L_upper = 1.e20
     COOP_SINGLE::P2_lower = 0.
     COOP_SINGLE::P2_upper = 1.e15
     COOP_SINGLE::I_lower_nu = -1.e20
     COOP_SINGLE::I_upper_nu = 1.e20
     COOP_SINGLE::L_lower_nu = -1.e20
     COOP_SINGLE::L_upper_nu = 1.e20
     COOP_SINGLE::P_lower_nu = 0.
     COOP_SINGLE::P_upper_nu = 1.e20
     COOP_SINGLE::P2byI2_lower = 0.
     COOP_SINGLE::P2byI2_upper = 1.e15
     COOP_SINGLE::P2byL2_lower = 0.
     COOP_SINGLE::P2byL2_upper = 1.e15
     COOP_INT::threshold_option = 0
     logical::addpi = .true.  !!randomly add pi on polarization directions
     logical::angzero = .true.  !!do not rotate
     logical::nested = .true.
     type(coop_list_integer)::peak_pix
     type(coop_list_realarr)::peak_ang
     type(coop_list_realarr)::peak_map
   contains
     procedure::export => coop_stacking_options_export
     procedure::import => coop_stacking_options_import
     procedure::subset =>  coop_stacking_options_subset
     procedure::init => coop_stacking_options_init
     procedure::free => coop_stacking_options_free
     procedure::convert2ring => coop_stacking_options_convert2ring
     procedure::convert2nested => coop_stacking_options_convert2nested     
     procedure::reject=>coop_stacking_options_reject
     procedure::pix => coop_stacking_options_pix
     procedure::export_pix => coop_stacking_options_export_pix
     procedure::export_ang => coop_stacking_options_export_ang     
     procedure::rotate_angle => coop_stacking_options_rotate_angle
     procedure::peak_r => coop_stacking_options_peak_r
     procedure::peak_e => coop_stacking_options_peak_e
     procedure::peak_get_angle_r_e => coop_stacking_options_peak_get_angle_r_e
  end type coop_stacking_options


  type coop_to_be_stacked
     logical::mask_int = .false.
     logical::mask_pol = .false.
     COOP_INT::nmaps = 0
     COOP_INT,dimension(:),allocatable::ind
     COOP_INT,dimension(:),allocatable::spin
     COOP_STRING,dimension(:),allocatable::label
     logical, dimension(:),allocatable::headless_vector !!headless vector
     logical, dimension(:),allocatable::local_rotation !!local rotation
     COOP_REAL, dimension(:),allocatable::zmin !!lowerbound
     COOP_REAL, dimension(:),allocatable::zmax !!upperbound
   contains     
     procedure::init => coop_to_be_stacked_init
     procedure::free => coop_to_be_stacked_free
  end type coop_to_be_stacked


contains


  subroutine coop_stacking_options_subset(this, subset)
    class(coop_stacking_options)::this
    type(coop_stacking_options)::subset
    COOP_INT::i
    call subset%free()
    if(subset%index_I .ne. 0 .and. (abs(subset%I_lower_nu) .lt. coop_stacking_max_threshold .or. abs(subset%I_upper_nu) .lt. coop_stacking_max_threshold ) )then
       if(this%sigma_I .eq. 0.d0) stop "nu != 0 but sigma_I = 0?"
       subset%sigma_I = this%sigma_I
       if(abs(subset%I_lower_nu) .lt. coop_stacking_max_threshold) subset%I_lower = subset%I_lower_nu * subset%sigma_I
       if(abs(subset%I_upper_nu) .lt. coop_stacking_max_threshold) subset%I_upper = subset%I_upper_nu *  subset%sigma_I
    endif
    if(subset%index_L .ne. 0 .and. (abs(subset%L_lower_nu) .lt. coop_stacking_max_threshold .or. abs(subset%L_upper_nu) .lt. coop_stacking_max_threshold ) )then
       if(this%sigma_L .eq. 0.d0)stop "nu !=0 but sigma_L = 0?"
       if(abs(subset%L_lower_nu) .lt. coop_stacking_max_threshold) subset%L_lower = subset%L_lower_nu * subset%sigma_P
       if(abs(subset%L_upper_nu) .lt. coop_stacking_max_threshold) subset%L_upper = subset%L_upper_nu * (this%L_upper/this%L_upper_nu)
    endif
    if(subset%index_Q .ne. 0 .and. subset%index_U .ne. 0 .and. (abs(subset%P_lower_nu) .lt. coop_stacking_max_threshold .or. abs(subset%P_upper_nu) .lt. coop_stacking_max_threshold ) )then
       if(this%sigma_P .eq. 0.d0)stop "nu !=0 but sigma_P = 0?"
       subset%sigma_P = this%sigma_P
       if(abs(subset%P_lower_nu) .lt. coop_stacking_max_threshold) subset%P2_lower = (subset%P_lower_nu*subset%sigma_P)**2 
       if(abs(subset%P_upper_nu) .lt. coop_stacking_max_threshold) subset%P2_upper = (subset%P_upper_nu*subset%sigma_P)**2
    endif
    do i = 1, this%peak_pix%n
       if( .not. subset%reject(this%peak_map%element(i) ))then
          call subset%peak_pix%push(this%peak_pix%element(i))
          call subset%peak_ang%push(this%peak_ang%element(i))
          call subset%peak_map%push(this%peak_map%element(i))
       endif
    enddo
  end subroutine coop_stacking_options_subset


  subroutine coop_stacking_options_convert2ring(this)
    class(coop_stacking_options)::this
    type(coop_list_integer)::pixcopy
    COOP_INT::pix_ring, pix_nested
    COOP_INT  i
#ifdef HAS_HEALPIX
    if(.not.this%nested) return
    pixcopy = this%peak_pix
    call this%peak_pix%init()
    do i=1, pixcopy%n
       call pixcopy%get_element(i, pix_nested)
       call nest2ring(this%nside, pix_nested, pix_ring)
       call this%peak_pix%push(pix_ring)
    enddo
    call pixcopy%init()
#else
    stop "you need to install Healpix"
#endif    
  end subroutine coop_stacking_options_convert2ring


  subroutine coop_stacking_options_convert2nested(this)
    class(coop_stacking_options)::this
    type(coop_list_integer)::pixcopy
    COOP_INT::pix_ring, pix_nested
    COOP_INT i
#ifdef HAS_HEALPIX
    if(this%nested) return
    pixcopy = this%peak_pix
    call this%peak_pix%init()
    do i=1, pixcopy%n
       call pixcopy%get_element(i, pix_nested)
       call ring2nest(this%nside, pix_nested, pix_ring)
       call this%peak_pix%push(pix_ring)
    enddo
    call pixcopy%init()    
#else
    stop "you need to install Healpix"
#endif    
  end subroutine coop_stacking_options_convert2nested
  
  
  subroutine coop_stacking_options_export(this, filename)
    class(coop_stacking_options)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file)fp
    COOP_INT i
    call fp%open(filename, "u")
    write(fp%unit) this%mask_int, this%mask_pol
    write(fp%unit) this%genre, this%nmaps, this%nside,this%index_I, this%index_Q, this%index_U, this%index_L
    write(fp%unit) this%I_lower, this%I_upper, this%L_lower, this%L_upper, this%P2_lower, this%P2_upper, this%I_lower_nu, this%I_upper_nu, this%L_lower_nu, this%L_upper_nu, this%P_lower_nu, this%P_upper_nu, this%P2byI2_lower, this%P2byI2_upper, this%P2byL2_lower, this%P2byL2_upper
    write(fp%unit) this%caption
    write(fp%unit) this%threshold_option, this%addpi, this%nested, this%angzero
    write(fp%unit) this%peak_pix%n
    do i=1, this%peak_pix%n
       write(fp%unit) this%peak_pix%element(i), this%peak_ang%element(i), this%peak_map%element(i)
    enddo
    write(fp%unit) this%sigma_I, this%sigma_L, this%sigma_P
    call fp%close()
  end subroutine coop_stacking_options_export

  subroutine coop_stacking_options_import(this, filename)
    class(coop_stacking_options)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file) fp
    COOP_INT pix, i, n
    COOP_SINGLE::thetaphi(2)
    COOP_SINGLE,dimension(:),allocatable::map
    if(.not. coop_file_exists(filename))then
       write(*,*) "stack option file "//trim(adjustl(filename))//" does not exist"
       stop
    endif
    call this%free()
    call fp%open(filename, "ur")
    read(fp%unit) this%mask_int, this%mask_pol    
    read(fp%unit) this%genre, this%nmaps, this%nside, this%index_I, this%index_Q, this%index_U, this%index_L
    read(fp%unit) this%I_lower, this%I_upper, this%L_lower, this%L_upper, this%P2_lower, this%P2_upper, this%I_lower_nu, this%I_upper_nu, this%L_lower_nu, this%L_upper_nu, this%P_lower_nu, this%P_upper_nu, this%P2byI2_lower, this%P2byI2_upper, this%P2byL2_lower, this%P2byL2_upper
    read(fp%unit) this%caption
    read(fp%unit) this%threshold_option, this%addpi, this%nested, this%angzero
    allocate(map(this%nmaps))
    read(fp%unit) n
    do i=1, n
       read(fp%unit) pix, thetaphi, map
       call this%peak_pix%push(pix)
       call this%peak_ang%push(thetaphi)
       call this%peak_map%push(map)
    enddo
    read(fp%unit) this%sigma_I, this%sigma_L, this%sigma_P    
    call fp%close()
    deallocate(map)
  end subroutine coop_stacking_options_import

  subroutine coop_stacking_options_free(this)
    class(coop_stacking_options)::this
    call this%peak_pix%init()
    call this%peak_ang%init()
    call this%peak_map%init()
  end subroutine coop_stacking_options_free

  subroutine coop_to_be_stacked_free(this, nmaps)
    class(coop_to_be_stacked)::this
    COOP_INT,optional::nmaps
    if(allocated(this%ind))deallocate(this%ind)
    if(allocated(this%spin))deallocate(this%spin)
    if(allocated(this%label))deallocate(this%label)
    if(allocated(this%zmin))deallocate(this%zmin)        
    if(allocated(this%zmax))deallocate(this%zmax)    
    if(allocated(this%headless_vector)) deallocate(this%headless_vector)
    if(allocated(this%local_rotation)) deallocate(this%local_rotation)
    if(present(nmaps))then
       this%nmaps = nmaps
       allocate(this%ind(this%nmaps), this%spin(this%nmaps), this%label(this%nmaps), this%local_rotation(this%nmaps), this%headless_vector(this%nmaps), this%zmin(this%nmaps), this%zmax(this%nmaps))
       this%zmin = 1.1e31
       this%zmax = -1.1e31
       this%headless_vector = .false.
       this%local_rotation = .false.
       this%label = ""
       this%spin = 0
       this%ind = 1
    else
       this%nmaps = 0
    endif
  end subroutine coop_to_be_stacked_free

  subroutine coop_to_be_stacked_init(this, str)
    class(coop_to_be_stacked)::this
    COOP_UNKNOWN_STRING::str
    type(coop_list_string)::l, subl
    COOP_STRING::line
    COOP_INT i, j
    this%mask_int = .false.
    this%mask_pol  = .false.    
    select case(trim(coop_str_numalpha(str)))
    case("LNI")
       call this%free(1)
       this%mask_int = .true.
       this%label = "$\ln I$"
    case("LNQU")
       call this%free(2)
       this%mask_pol = .true.       
       this%ind = (/ 1, 2 /)       
       this%spin = 2
       this%label(1) =  "$\ln Q$"
       this%label(2) =  "$\ln U$"
       this%headless_vector = .true.
    case("I", "T")
       call this%free(1)
       this%mask_int = .true.
       this%label = "$"//trim(adjustl(str))//"(\mu K)$"
    case("LT")
       call this%free(1)
       this%mask_int = .true.
       this%label = "$\nabla^2 T (\mu K)$"
    case("E", "B")
       call this%free(1)
       this%mask_pol = .true.       
       this%label = "$"//trim(adjustl(str))//"(\mu K)$"
    case("zeta", "Z")
       call this%free(1)
       this%mask_int = .true.       
       this%label = "$10^5\zeta$"
    case("QU")
       call this%free(2)
       this%mask_pol = .true.       
       this%ind = (/ 1, 2 /)       
       this%spin = 2
       this%label(1) =  "$Q(\mu K)$"
       this%label(2) =  "$U(\mu K)$"
       this%headless_vector = .true.
    case("QrUr")
       call this%free(2)
       this%mask_pol = .true.       
       this%ind = (/ 1, 2 /)
       this%spin = 2
       this%label(1) =  "$Q_r(\mu K)$"
       this%label(2) =  "$U_r(\mu K)$"
       this%headless_vector = (/ .true., .false. /)
       this%local_rotation = .true.
    case("QTUT")
       call this%free(2)
       this%mask_int = .true.       
       this%ind = (/ 1, 2 /)
       this%spin = 2
       this%label(1) =  "$Q_T(\mu K)$"
       this%label(2) =  "$U_T(\mu K)$"
       this%headless_vector = .true.
    case("QTUZ")
       call this%free(2)
       this%mask_int = .true.       
       this%ind = (/ 1, 2 /)
       this%spin = 2
       this%label(1) =  "$Q_\zeta(\mu K)$"
       this%label(2) =  "$U_\zeta(\mu K)$"
       this%headless_vector = .true.
    case("QZrUZr")
       call this%free(2)
       this%mask_int = .true.       
       this%ind = (/ 1, 2 /)
       this%spin = 2
       this%label(1) =  "$Q_{\zeta,r}(\mu K)$"
       this%label(2) =  "$U_{\zeta,r}(\mu K)$"
       this%headless_vector = (/ .true., .false. /)
       this%local_rotation = .true.
    case("QTrUTr")
       call this%free(2)
       this%mask_int = .true.       
       this%ind = (/ 1, 2 /)
       this%spin = 2
       this%label(1) =  "$Q_{T,r}(\mu K)$"
       this%label(2) =  "$U_{T,r}(\mu K)$"
       this%headless_vector = (/ .true., .false. /)
       this%local_rotation = .true.
    case("QLTULT")
       call this%free(2)
       this%mask_int = .true.       
       this%ind = (/ 1, 2 /)
       this%spin = 2
       this%label(1) =  "$Q_{\nabla^2 T}(\mu K/\mathrm{rad}^2)$"
       this%label(2) =  "$U_{\nabla^2 T}(\mu K/\mathrm{rad}^2)$"
       this%headless_vector = .true.
    case("QLTrULTr")
       call this%free(2)
       this%mask_int = .true.       
       this%ind = (/ 1, 2 /)
       this%spin = 2
       this%label(1) =  "$Q_{\nabla^2 T,r}(\mu K)$"
       this%label(2) =  "$U_{\nabla^2 T,r}(\mu K)$"
       this%headless_vector = (/ .true., .false. /)
       this%local_rotation = .true.       
    case default  !! label1@index1@spin1@headless_vector1@local_roataion1@zmin1@zmax1::label2:index2:spin2@headless_vector2@local_roataion2@zmin2@zmax2
       call coop_string_to_list(str, l, "::")
       call this%free(l%n)
       do i = 1, l%n
          call coop_string_to_list(l%element(i), subl, "@")
          this%label(i) = trim(subl%element(1))
          if(subl%n .ge. 2)then
             call subl%get_element(2, line)
             read(line, *) this%ind(i)
             if(subl%n .ge. 3)then
                call subl%get_element(3, line)
                read(line, *) this%spin(i)
                if(subl%n .ge. 4)then
                   this%headless_vector(i) = (trim(subl%element(4)) .eq. "T" .and. this%spin(i) .eq. 2)
                   if(subl%n .ge. 5)then
                      this%local_rotation(i) = (trim(subl%element(5)) .eq. "T" .and. this%spin(i) .eq. 2)
                      if(subl%n .ge. 6)then
                         call subl%get_element(6, line)
                         read(line, *) this%zmin(i)
                         if(subl%n .ge. 7)then
                            call subl%get_element(7, line)
                            read(line, *) this%zmax(i)
                         endif
                      endif
                   endif
                endif                
             endif
          endif
       enddo
       call l%init()
       call subl%init()
       write(*,*) "Warning: unclassified maps to be stacked -- mask automatic determination can be wrong"
    end select
    !!check spin pairs
    i = 1
    do while(i.le. this%nmaps)
       if(this%spin(i) .ne. 0)then
          if( i .ge. this%nmaps)then
             print*, this%nmaps
             print*, this%spin
             print*, this%ind             
             stop "to_be_stacked_init: nonzero spin must go in pairs"
          endif
          if( this%spin(i+1) .ne. this%spin(i) .or. this%ind(i+1).ne. this%ind(i)+1)then
             print*, this%nmaps
             print*, this%spin
             print*, this%ind
             stop "to_be_stacked_init: nonzero spin must go in pairs"
          endif
          i = i + 2
       else
          i = i + 1
       endif
    enddo
  end subroutine coop_to_be_stacked_init

  subroutine coop_stacking_options_init(this, domax, peak_name, Orient_name, nmaps)
    class(coop_stacking_options)::this
    COOP_UNKNOWN_STRING::peak_name, orient_name
    COOP_SHORT_STRING::p
    logical domax
    COOP_INT::nmaps
    call this%free()
    this%mask_int = .false.
    this%mask_pol  = .false.
    this%nmaps = nmaps    
    p = trim(adjustl(peak_name))
    if(trim(adjustl(orient_name)).eq. "NULL" .or. trim(adjustl(orient_name)) .eq. "RANDOM" .or. trim(adjustl(orient_name)).eq. "NONE" )then !!random orientation
       if(trim(p) .eq. "SADDLE" .or. trim(p).eq."COL")then
          this%genre = coop_stacking_genre_saddle
          this%caption = "saddle points"
       elseif(domax)then
          if(trim(p).eq."RANDOM")then
             this%genre = coop_stacking_genre_random_hot
             this%caption = "hot spots"
          else          
             if(len_trim(p).gt. 9)then
                if(p(1:9).eq."$\nabla^2")then
                   this%genre = coop_stacking_genre_Lmax
                else
                   this%genre = coop_stacking_genre_Imax
                endif
             else
                this%genre = coop_stacking_genre_Imax
             endif
             this%caption = trim(p)//" maxima"
          endif
       else
          if(trim(p).eq."RANDOM")then
             this%genre = coop_stacking_genre_random_cold
             this%caption = "cold spots"             
          else          
             if(len_trim(p) .gt. 9)then
                if(p(1:9) .eq. "$\nabla^2")then
                   this%genre = coop_stacking_genre_Lmin
                else
                   this%genre = coop_stacking_genre_Imin
                endif
             else
                this%genre = coop_stacking_genre_Imin
             endif
             this%caption = trim(p)//" minima"
          endif
       endif
    else
       if(trim(p).eq."SADDLE")then
          this%genre = coop_stacking_genre_saddle_oriented
          this%caption = " saddle points, "//trim(adjustl(Orient_name))//" oriented"
       elseif(trim(p).eq."COL")then
          this%genre = coop_stacking_genre_col_oriented
          this%caption = "cols, "//trim(adjustl(orient_name))//" oriented"
       elseif(domax)then
          if(trim(p).eq."RANDOM")then
             this%genre = coop_stacking_genre_random_hot_Oriented
             this%caption = "hot spots, "//trim(adjustl(Orient_name))//" oriented"
          else
             if(trim(p).eq."P" .or. trim(p).eq."P_T" .or. p(1:2) .eq. "$P")then
                this%genre = coop_stacking_genre_Pmax_Oriented             
             else
                if(len_trim(p).gt. 9)then
                   if(p(1:9).eq."$\nabla^2")then
                      this%genre = coop_stacking_genre_Lmax_Oriented
                   else
                      this%genre = coop_stacking_genre_Imax_Oriented               
                   endif
                else
                   this%genre = coop_stacking_genre_Imax_Oriented    
                endif
             endif
             this%caption = trim(p)//" maxima, "//trim(adjustl(Orient_name))//" oriented"
          endif
       else
          if(trim(p).eq."RANDOM")then
             this%genre = coop_stacking_genre_random_cold_Oriented
             this%caption = "cold spots, "//trim(adjustl(Orient_name))//" oriented"
          else
             if(p(1:2).eq."$P")then
                this%genre = coop_stacking_genre_Pmin_Oriented             
             else
                if(len_trim(p) .gt. 9)then
                   if(p(1:9) .eq. "$\nabla^2")then
                      this%genre = coop_stacking_genre_Lmin_Oriented
                   else
                      this%genre = coop_stacking_genre_Imin_Oriented
                   endif
                else
                   this%genre = coop_stacking_genre_Imin_Oriented
                endif
             endif
             this%caption = trim(adjustl(peak_name))//" minima, "//trim(adjustl(Orient_name))//" oriented"
          endif
       endif
    endif

    !!now set default indices
    select case(nmaps)
    case(1)
       select case(this%genre)
       case(coop_stacking_genre_Imax, coop_stacking_genre_Imin, coop_stacking_genre_random_hot, coop_stacking_genre_random_cold)
          this%index_I = 1
          this%threshold_option = 4          
       case(coop_stacking_genre_Lmax, coop_stacking_genre_Lmin)
          this%index_L = 1
          this%threshold_option = 1
       case default
          write(*,*) this%genre
          write(*,*) nmaps
          stop "nmaps does not match stacking option"
       end select
    case(2)
       select case(this%genre)
       case(coop_stacking_genre_Pmax_Oriented, coop_stacking_genre_Pmin_Oriented, coop_stacking_genre_random_hot, coop_stacking_genre_random_cold, coop_stacking_genre_random_hot_oriented, coop_stacking_genre_random_cold_oriented)
          this%index_Q = 1
          this%index_U = 2
          this%threshold_option = 2          
       case(coop_stacking_genre_Imax, coop_stacking_genre_Imin, coop_stacking_genre_Lmax, coop_stacking_genre_Lmin)
          this%index_I = 1
          this%index_L = 2
          this%threshold_option = 5          
       case default
          write(*,*) this%genre
          write(*,*) nmaps
          stop "nmaps does not match stacking option"          
       end select
    case(3)
       select case(this%genre)
       case(coop_stacking_genre_Pmax_Oriented, coop_stacking_genre_Pmin_Oriented, coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Imin_Oriented, coop_stacking_genre_Imax, coop_stacking_genre_Imin, coop_stacking_genre_random_hot, coop_stacking_genre_random_cold, coop_stacking_genre_random_hot_oriented, coop_stacking_genre_random_cold_oriented)
          this%index_I = 1
          this%index_Q = 2
          this%index_U = 3
          this%threshold_option = 6          
       case(coop_stacking_genre_Lmax, coop_stacking_genre_Lmin, coop_stacking_genre_Lmax_Oriented, coop_stacking_genre_Lmin_Oriented)  
          this%index_L = 1
          this%index_Q = 2
          this%index_U = 3
          this%threshold_option = 3
       end select
    case(4:)
       this%index_I = 1
       this%index_Q = 2
       this%index_U = 3
       this%index_L = 4
       this%threshold_option = 7       
    end select
    select case(trim(coop_str_numUpperalpha(peak_name)))
    case("T", "I", "ZETA", "PT", "PZ", "PZETA", "SADDLE", "COL")
       this%mask_int = .true.
    case("E", "B", "P")
       this%mask_pol = .true.
    case("RANDOM")
       if(this%threshold_option .eq. 2 .and. this%index_Q .ne. 0 .and. this%index_U .ne. 0 )then
          this%mask_pol = .true.
       else
          this%mask_int = .true.
       endif
    case default
       write(*,*) "Unknown class of peaks: cannot automatically determine mask type"
    end select
    select case(trim(coop_str_numUpperalpha(orient_name)))
    case("QU", "QEUE", "QNABLA2EUNABLA2E")
       this%mask_pol = .true.
    case("QTUT", "QNABLA2TUNABLA2T","QNABLA2ZETAUNABLA2ZETA", "QZETAUZETA")
       this%mask_int = .true.
    end select
  end subroutine coop_stacking_options_init


  function coop_stacking_options_reject(this, map) result(rej)
    class(coop_stacking_options)::this
    logical rej
    COOP_SINGLE::map(:)
    COOP_REAL::p2
    select case(this%threshold_option)
    case(7)
       P2 = map(this%index_Q)**2+map(this%index_U)**2 
       rej = (map(this%index_I) > this%I_upper .or. map(this%index_I) < this%I_lower .or. P2> this%P2_upper .or. P2 < this%P2_lower .or.  map(this%index_L) > this%L_upper .or. map(this%index_L) < this%L_lower .or. P2 > this%P2byI2_upper*map(this%index_I)**2 .or. P2 < this%P2byI2_lower * map(this%index_I)**2 .or. P2 > this%P2byL2_upper*map(this%index_L)**2 .or. P2 < this%P2byL2_lower*map(this%index_L)**2)
    case(6)
       P2 = map(this%index_Q)**2+map(this%index_U)**2        
       rej = (map(this%index_I) > this%I_upper .or. map(this%index_I) < this%I_lower .or. P2 > this%P2_upper .or. P2 < this%P2_lower .or. P2 > this%P2byI2_upper*map(this%index_I)**2 .or. P2 < this%P2byI2_lower * map(this%index_I)**2)
    case(5)
       rej = ( map(this%index_I) > this%I_upper .or. map(this%index_I) < this%I_lower  .or.  map(this%index_L) > this%L_upper .or. map(this%index_L) < this%L_lower)
    case(4)
       rej = (map(this%index_I) > this%I_upper .or. map(this%index_I) < this%I_lower)
    case(3)
       P2 = map(this%index_Q)**2+map(this%index_U)**2               
       rej = ( P2 > this%P2_upper .or. P2 < this%P2_lower .or.  map(this%index_L) > this%L_upper .or. map(this%index_L) < this%L_lower .or. P2 > this%P2byL2_upper*map(this%index_L)**2 .or. P2 < this%P2byL2_lower*map(this%index_L)**2)
    case(2)
       P2 = map(this%index_Q)**2+map(this%index_U)**2                      
       rej = (P2 > this%P2_upper .or. P2 < this%P2_lower)
    case(1)
       rej = ( map(this%index_L) > this%L_upper .or.  map(this%index_L) < this%L_lower)
    case default
       rej  = .false.
    end select
  end function coop_stacking_options_reject

  function coop_stacking_options_pix(this, nside, i) result(pix)
    class(coop_stacking_options)::this
    COOP_INT nside, i, pix
    COOP_SINGLE::thetaphi(2)
#ifdef HAS_HEALPIX    
    if(this%nside .eq. nside)then
       pix = this%peak_pix%element(i)
       return
    endif
    call this%peak_ang%get_element(i, thetaphi)
    if(this%nested)then
       call ang2pix_nest(nside, dble(thetaphi(1)), dble(thetaphi(2)), pix)
    else
       call ang2pix_ring(nside, dble(thetaphi(1)), dble(thetaphi(2)), pix)
    endif
#endif    
  end function coop_stacking_options_pix


  subroutine coop_stacking_options_export_pix(this, filename)
    class(coop_stacking_options)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file)::fp
    COOP_INT :: i
    call fp%open(filename, "w")
    do i = 1, this%peak_pix%n
       write(fp%unit, "(I8)") this%peak_pix%element(i)
    enddo
    call fp%close()
  end subroutine coop_stacking_options_export_pix

  subroutine coop_stacking_options_export_ang(this, filename)
    class(coop_stacking_options)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file)::fp
    COOP_INT :: i
    logical dor
    call fp%open(filename, "w")
    do i = 1, this%peak_ang%n
       write(fp%unit, "(3E14.5)") this%peak_ang%element(i),  this%rotate_angle(i)
    enddo
    call fp%close()
  end subroutine coop_stacking_options_export_ang
  

  function coop_stacking_options_rotate_angle(this, i) result(angle)
    class(coop_stacking_options)::this
    COOP_REAL angle
    COOP_SINGLE::map(this%nmaps)
    COOP_INT i
    if(i.gt. this%peak_pix%n) stop "rotate_angle: pix overflow"
    select case(this%genre)
    case(coop_stacking_genre_Imax, coop_stacking_genre_Imin, coop_stacking_genre_Lmax, coop_stacking_genre_Lmin, coop_stacking_genre_Null, coop_stacking_genre_random_hot, coop_stacking_genre_random_cold, coop_stacking_genre_Saddle)
       if(this%angzero)then
          angle = 0.d0
          return
       else
          call random_number(angle)
          angle = angle*coop_2pi
          return
       endif
    case(coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Imin_Oriented, coop_stacking_genre_Lmax_Oriented, coop_stacking_genre_Lmin_Oriented, coop_stacking_genre_Pmax_Oriented, coop_stacking_genre_Pmin_Oriented, coop_stacking_genre_random_hot_oriented, coop_stacking_genre_random_cold_oriented, coop_stacking_genre_Saddle_oriented)
       call this%peak_map%get_element(i, map)
       angle = COOP_POLAR_ANGLE(dble(map(this%index_Q)),dble(map(this%index_U)))/2.d0
       if(this%addpi) angle = angle + coop_rand01()*coop_pi
    case(coop_stacking_genre_Col_Oriented)
       call this%peak_map%get_element(i, map)
       angle = COOP_POLAR_ANGLE(dble(map(this%index_Q)),dble(map(this%index_U)))/2.d0
       if(map(10)*cos(angle) + map(9)*sin(angle) .le. 0.d0)then
          angle = angle + coop_pi
       endif
    case default
       write(*,*) this%genre
       stop "rotate_angle: unknown genre"
    end select
  end function coop_stacking_options_rotate_angle

  function coop_stacking_options_peak_r(this, i) result(r)
    COOP_REAL, parameter::max_radius = coop_pio2
    class(coop_stacking_options)::this
    COOP_INT i
    COOP_REAL r
    COOP_SINGLE::map(this%nmaps)
    call this%peak_map%get_element(i, map)
    r = min(sqrt(abs(map(this%index_I))/max(abs(map(this%index_L)),1.d-20)), max_radius)
  end function coop_stacking_options_peak_r

  function coop_stacking_options_peak_e(this, i) result(e)
    COOP_REAL, parameter::max_radius = coop_pio2
    class(coop_stacking_options)::this
    COOP_INT i
    COOP_REAL e
    COOP_SINGLE::map(this%nmaps)
    call this%peak_map%get_element(i, map)
    e = min(sqrt(map(this%index_Q)**2 + map(this%index_U)**2)/max(abs(map(this%index_L)),1.d-20), 1.d0)
  end function coop_stacking_options_peak_e

  subroutine coop_stacking_options_peak_get_angle_r_e(this, i, angle, r, e)
    COOP_REAL, parameter::max_radius = coop_pio2
    class(coop_stacking_options)::this
    COOP_REAL angle,r,e
    COOP_SINGLE::map(this%nmaps)
    COOP_INT i
    if(i.gt. this%peak_pix%n) stop "rotate_angle: pix overflow"
    call this%peak_map%get_element(i, map)    
    select case(this%genre)
    case(coop_stacking_genre_Imax, coop_stacking_genre_Imin, coop_stacking_genre_random_hot, coop_stacking_genre_random_cold)
       call random_number(angle)
       angle = angle*coop_2pi
       e = 0.d0
       if(this%index_I.gt.0 .and. this%index_L.gt.0)then
          r = min(sqrt(abs(map(this%index_I))/max(abs(map(this%index_L)),1.d-20)), max_radius)
       else
          r = 0.d0
       endif
    case(coop_stacking_genre_Lmax, coop_stacking_genre_Lmin)
       call random_number(angle)
       angle = angle*coop_2pi
       e = 0.d0
       r = 0.d0
    case(coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Imin_Oriented, coop_stacking_genre_random_hot_oriented, coop_stacking_genre_random_cold_oriented)
       angle = COOP_POLAR_ANGLE(map(this%index_Q), map(this%index_U))
       if(this%index_L .gt. 0)then
          e = min(sqrt(map(this%index_Q)**2 + map(this%index_U)**2)/max(abs(map(this%index_L)),1.d-20), 0.99d0)
       else
          e = min(sqrt(map(this%index_Q)**2 + map(this%index_U)**2)/max(abs(map(this%index_I)),1.d-20), 0.99d0)
       endif
       if(this%index_I.gt.0 .and. this%index_L.gt.0)then
          r = min(sqrt(abs(map(this%index_I))/max(abs(map(this%index_L)),1.d-20)), max_radius)
       else
          r = 0.d0
       endif
    case(coop_stacking_genre_Lmax_Oriented, coop_stacking_genre_Lmin_Oriented, coop_stacking_genre_Pmax_Oriented, coop_stacking_genre_Pmin_Oriented)
       angle = COOP_POLAR_ANGLE(map(this%index_Q), map(this%index_U))
       e = 0.d0
       r = 0.d0
    case default
       write(*,*) this%genre
       stop "rotate_angle: unknown genre"
    end select
    r = min(sqrt(abs(map(this%index_I))/max(abs(map(this%index_L)),1.d-20)), max_radius)
  end subroutine coop_stacking_options_peak_get_angle_r_e

  subroutine coop_stacking_options_split_hemispheres(sto, north, south, l_deg, b_deg)
    type(coop_stacking_options)::sto, north, south
    COOP_REAL l_deg, b_deg, theta, phi, cost, sint
    COOP_INT i
    COOP_SINGLE tp(2)
    call coop_healpix_lb2ang(l_deg, b_deg, theta, phi)
    cost = cos(theta)
    sint = sin(theta)
    north = sto
    south = sto
    call north%free()
    call south%free()
    do i=1, sto%peak_ang%n
       call sto%peak_ang%get_element(i, tp)
       if(cos(tp(1))*cost + sin(tp(1))*sint*cos(phi - tp(2)) .ge. 0.d0)then
          call north%peak_pix%push(sto%peak_pix%element(i))
          call north%peak_ang%push(tp)
          call north%peak_map%push(sto%peak_map%element(i))          
       else
          call south%peak_pix%push(sto%peak_pix%element(i))          
          call south%peak_ang%push(tp)
          call south%peak_map%push(sto%peak_map%element(i))          
       endif
    enddo
  end subroutine coop_stacking_options_split_hemispheres


  subroutine coop_healpix_lb2ang(l_deg, b_deg, theta, phi)
    COOP_REAL l_deg, b_deg
    COOP_REAL theta, phi
    if(b_deg .lt. -90.d0 .or. b_deg .gt. 90.d0) stop "b must be between -90 deg to 90 deg"
    theta = (90.d0 - b_deg)*coop_SI_degree
    phi = l_deg * coop_SI_degree
    do while(phi .ge. coop_2pi)
       phi = phi - coop_2pi
    enddo
    do while(phi .lt. 0.d0)
       phi = phi + coop_2pi
    enddo
  end subroutine coop_healpix_lb2ang


  subroutine coop_healpix_ang2lb(theta, phi, l_deg, b_deg)
    COOP_REAL l_deg, b_deg
    COOP_REAL theta, phi
    if(theta .gt. coop_pi+1.d-12 .or. theta .lt. 0.d0) stop "theta must be between 0 and pi"
    b_deg = (coop_pio2 - theta)/coop_SI_degree
    l_deg = phi/coop_SI_degree
    do while(l_deg .gt. 360.d0)
       l_deg = l_deg - 360.d0
    enddo
    do while(l_deg .lt. 0.d0)
       l_deg = l_deg + 360.d0
    enddo
  end subroutine coop_healpix_ang2lb

end module coop_stacking_mod
