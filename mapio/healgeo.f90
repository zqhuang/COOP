module coop_healpix_mod
!!I always assume ring order
  use coop_wrapper_firstorder
  use coop_sphere_mod
  use coop_stacking_mod
#ifdef HAS_HEALPIX
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  use udgrade_nr
  use coord_v_convert,only:coordsys2euler_zyz
#endif
  implicit none

  
#include "constants.h"

  private

  
  !!some stupid conventions...
  !!the direction of headless vector rotate by pi/2 if you use IAU convention
  logical::coop_healpix_IAU_headless_vector = .true.
  COOP_REAL::coop_healpix_QrUrSign = -1.d0  !!WMAP7 invented this...
  
  COOP_SINGLE::coop_healpix_patch_default_figure_width = 4.5
  COOP_SINGLE::coop_healpix_patch_default_figure_height = 4.
  logical::coop_healpix_patch_stack_fft = .false.
  logical::coop_healpix_patch_default_want_caption = .true.
  logical::coop_healpix_patch_default_want_label = .true.
  logical::coop_healpix_patch_default_want_arrow = .false.
  COOP_INT, parameter:: coop_inpaint_nside_start = 8
  COOP_SINGLE,parameter::coop_inpaint_mask_threshold = 0.05
  
  public::coop_fits_to_header, coop_healpix_maps, coop_healpix_disc, coop_healpix_patch, coop_healpix_split,  coop_healpix_output_map, coop_healpix_smooth_mapfile, coop_healpix_patch_get_fr0, coop_healpix_mask_tol,  coop_healpix_mask_hemisphere, coop_healpix_index_TT,  coop_healpix_index_EE,  coop_healpix_index_BB,  coop_healpix_index_TE,  coop_healpix_index_TB,  coop_healpix_index_EB, coop_healpix_flip_mask, coop_healpix_alm_check_done, coop_healpix_want_cls, coop_healpix_default_lmax, coop_planck_TNoise, coop_planck_ENoise, coop_Planck_BNoise, coop_highpass_filter, coop_lowpass_filter, coop_gaussian_filter,coop_healpix_IAU_headless_vector,  coop_healpix_spot_select_mask, coop_healpix_spot_cut_mask, coop_healpix_merge_masks, coop_healpix_patch_default_figure_width, coop_healpix_patch_default_figure_height, coop_healpix_patch_default_want_caption, coop_healpix_patch_default_want_label, coop_healpix_patch_default_want_arrow,  coop_healpix_QrUrSign, coop_ACT_TNoise, coop_ACT_ENoise, coop_healpix_inpaint, coop_healpix_maps_ave_udgrade, coop_healpix_maps_copy_genre, coop_healpix_correlation_function, coop_healpix_mask_reverse, coop_healpix_maps_diffuse, coop_healpix_mask_diffuse, coop_healpix_nside2lmax, coop_healpix_filament, coop_healpix_lmax_by_nside, coop_healpix_patch_stack_fft
  

  logical::coop_healpix_alm_check_done = .false.
  logical::coop_healpix_want_cls = .true.
  COOP_REAL,parameter::coop_healpix_zeta_normalization = 1.d-5
  COOP_UNKNOWN_STRING,parameter::coop_healpix_maps_default_genre = "GENERAL"

  COOP_INT, parameter::coop_healpix_default_lmax=2500
  COOP_REAL, parameter::coop_healpix_lmax_by_nside = 2.5
  COOP_REAL::coop_healpix_mask_tol = 0.d0  !!default mask tolerance
  COOP_INT,parameter::coop_healpix_index_TT = 1
  COOP_INT,parameter::coop_healpix_index_EE = 2
  COOP_INT,parameter::coop_healpix_index_BB = 3
  COOP_INT,parameter::coop_healpix_index_TE = 4
  COOP_INT,parameter::coop_healpix_index_EB = 5
  COOP_INT,parameter::coop_healpix_index_TB = 6

  type, extends(coop_sphere_disc):: coop_healpix_disc
     COOP_INT::nside = 0
     COOP_INT::center
     COOP_INT::ordering = COOP_RING
   contains
     procedure :: pix2ang => coop_healpix_disc_pix2ang
     procedure :: ang2pix => coop_healpix_disc_ang2pix
     procedure :: pix2xy => coop_healpix_disc_pix2xy
     procedure :: xy2pix => coop_healpix_disc_xy2pix
  end type coop_healpix_disc

  type coop_healpix_maps
     COOP_INT::nside = 0
     COOP_INT::nmaps = 0
     COOP_INT:: npix = 0
     COOP_INT::ordering = COOP_RING
     COOP_INT::lmax = -1
     COOP_SHORT_STRING::creator = "COOP"
     COOP_SHORT_STRING::dtype = "MAP"     
     COOP_SHORT_STRING::version = "0.0"
     COOP_SHORT_STRING::coordsys = "G"
     COOP_STRING::filename = ""
     COOP_INT::firstpix = 0
     COOP_INT::lastpix = 0
     COOP_INT::iq = 0
     COOP_INT::iu = 0
     COOP_SINGLE::bad_data =  -1.637500000000e30
     COOP_REAL::fwhm_degree = 0.
     logical::polar = .false.     
     type(coop_dictionary)::header
     COOP_INT,dimension(:),allocatable::spin
     COOP_SHORT_STRING,dimension(:),allocatable::fields
     COOP_SHORT_STRING,dimension(:),allocatable::units     
     COOP_SINGLE, dimension(:,:),allocatable::map
     COOP_SINGLE_COMPLEX, dimension(:,:,:),allocatable::alm
     COOP_SINGLE, dimension(:,:),allocatable::Cl
     logical,dimension(:),allocatable::alm_done
     COOP_SINGLE,dimension(:,:),allocatable::checksum
   contains
     procedure :: fields_to_spins => coop_healpix_maps_fields_to_spins
     procedure :: set_unit => coop_healpix_maps_set_unit
     procedure :: set_units => coop_healpix_maps_set_units     
     procedure :: set_field => coop_healpix_maps_set_field
     procedure :: set_fields => coop_healpix_maps_set_fields     
     procedure :: regularize_in_mask => coop_healpix_maps_regularize_in_mask
     procedure :: ang2pix => coop_healpix_maps_ang2pix
     procedure :: pix2ang => coop_healpix_maps_pix2ang
     procedure :: vec2pix => coop_healpix_maps_vec2pix
     procedure :: pix2vec => coop_healpix_maps_pix2vec
     procedure :: query_disc => coop_healpix_maps_query_disc
     procedure :: query_strip => coop_healpix_maps_query_strip
     procedure :: get_disc => coop_healpix_maps_get_disc
     procedure :: init => coop_healpix_maps_init
     procedure :: free => coop_healpix_maps_free
     procedure :: write => coop_healpix_maps_write
     procedure :: read => coop_healpix_maps_read
     procedure :: extend => coop_healpix_maps_extend
     procedure :: draw_latitude_line => coop_healpix_maps_draw_latitude_line
     procedure :: import => coop_healpix_maps_import
     procedure :: open => coop_healpix_maps_read
     procedure :: allocate_alms => coop_healpix_maps_allocate_alms
     procedure :: get_cls =>   coop_healpix_maps_get_cls
     procedure :: map2alm => coop_healpix_maps_map2alm
     procedure :: alm2map => coop_healpix_maps_alm2map
     procedure :: udgrade => coop_healpix_maps_udgrade
     procedure :: simulate => coop_healpix_maps_simulate
     procedure :: simulate_Tmaps => coop_healpix_maps_simulate_Tmaps
     procedure :: simulate_TQUmaps => coop_healpix_maps_simulate_TQUmaps
     procedure :: qu2EB => coop_healpix_maps_qu2EB
     procedure :: EB2QU => coop_healpix_maps_EB2QU
     procedure :: get_QU => coop_healpix_maps_get_QU
     procedure :: get_QUL => coop_healpix_maps_get_QUL
     procedure :: get_QULDD => coop_healpix_maps_get_QULDD     
     procedure :: get_dervs => coop_healpix_maps_get_dervs
     procedure :: smooth => coop_healpix_maps_smooth
     procedure :: smooth_with_window => coop_healpix_maps_smooth_with_window
     procedure :: T2zeta => coop_healpix_maps_t2zeta
     procedure :: TE2zeta => coop_healpix_maps_te2zeta
     procedure :: E2zeta => coop_healpix_maps_E2zeta
     procedure :: convert2nested => coop_healpix_convert_to_nested
     procedure :: convert2ring => coop_healpix_convert_to_ring
     procedure :: filter_alm =>  coop_healpix_filter_alm
     procedure :: reflect => coop_healpix_maps_reflect
     !!mask
     procedure :: apply_mask => coop_healpix_maps_apply_mask
     procedure :: generate_latcut_mask => coop_healpix_maps_generate_latcut_mask
     procedure :: mask_disc => coop_healpix_maps_mask_disc
     procedure :: mask_strip => coop_healpix_maps_mask_strip
     procedure :: mask_pixel => coop_healpix_maps_mask_pixel
     !!stacking stuff
     procedure :: fetch_patch => coop_healpix_maps_fetch_patch
     procedure :: stack_on_patch => coop_healpix_maps_stack_on_patch
     procedure :: get_peaks => coop_healpix_maps_get_peaks
     procedure :: mark_peaks => coop_healpix_maps_mark_peaks
     procedure :: stack_on_peaks  =>     coop_healpix_maps_stack_on_peaks
     procedure :: stack_filaments_on_peaks  =>     coop_healpix_maps_stack_filaments_on_peaks     
     procedure:: mask_peaks => coop_healpix_mask_peaks
     procedure:: rotate_coor => coop_healpix_maps_rotate_coor
     procedure:: distr_nu_e => coop_healpix_maps_distr_nu_e
     procedure :: fetch_filament => coop_healpix_maps_fetch_filament
     procedure :: stack_on_filament => coop_healpix_maps_stack_on_filament
     procedure :: perimeter_area_list => coop_healpix_maps_perimeter_area_list
     procedure :: zeros => coop_healpix_maps_zeros
     procedure :: local_disk_minkowski0 => coop_healpix_maps_local_disk_minkowski0
     procedure :: scan_local_minkowski0 => coop_healpix_maps_scan_local_minkowski0
     procedure :: local_disk_minkowski1 => coop_healpix_maps_local_disk_minkowski1
     procedure :: scan_local_minkowski1 => coop_healpix_maps_scan_local_minkowski1
     
  end type coop_healpix_maps

  type coop_healpix_patch
     COOP_STRING::caption=""
     type(coop_to_be_stacked):: tbs
     COOP_SHORT_STRING::color_table="Rainbow"
     COOP_SHORT_STRING::genre 
     COOP_INT::n = 0
     COOP_INT::mmax = 0
     COOP_INT::nmaps = 0
     COOP_INT::npix = 0
     COOP_INT::nstack_raw = 0
     COOP_REAL::dr
     COOP_REAL,dimension(:,:,:),allocatable::image
     COOP_COMPLEX,dimension(:,:,:),allocatable::fftimage     
     COOP_REAL,dimension(:),allocatable::r
     COOP_REAL,dimension(:,:,:),allocatable::fr
     COOP_REAL,dimension(:,:,:),allocatable::wcm
     COOP_REAL,dimension(:,:,:),allocatable::wsm
     COOP_INT, dimension(:,:,:),allocatable::icm
     COOP_REAL,dimension(:,:),allocatable::nstack, indisk
     COOP_REAL::num_indisk_tol
   contains
     procedure::free => coop_healpix_patch_free
     procedure::init => coop_healpix_patch_init
     procedure::get_radial_profile => coop_healpix_patch_get_radial_profile
     procedure::get_all_radial_profiles => coop_healpix_patch_get_all_radial_profiles
     procedure::plot => coop_healpix_patch_plot
     procedure::fft => coop_healpix_patch_fft     
     procedure::export => coop_healpix_patch_export
     procedure::import => coop_healpix_patch_import     
     procedure::plot_fft => coop_healpix_patch_plot_fft
  end type coop_healpix_patch


  type, extends(coop_healpix_patch)::coop_healpix_filament
     COOP_INT::ny = 0
   contains
     procedure::plot => coop_healpix_filament_plot
  end type coop_healpix_filament


  type coop_healpix_inpaint
     COOP_INT::lmax = -1
     COOP_INT::ncorr = 0
     COOP_REAL::fsky = 0.d0
     COOP_REAL::dtheta = 0.d0
     COOP_REAL::sigma0 = 0.d0
     COOP_INT::nMT = 0
     COOP_INT::ncT = 0
     COOP_INT::base_nside = 0
     logical::first_realization = .true.
     COOP_REAL,dimension(:),allocatable::als, bls, Cls, sqrtCls, smooth_Cls
     COOP_REAL,dimension(:),allocatable::corr
     COOP_REAL,dimension(:,:),allocatable::vec
     logical,dimension(:),allocatable::lask
     type(coop_healpix_maps),pointer::map => null()
     type(coop_healpix_maps),pointer::mask => null()
     type(coop_healpix_maps)::smooth_mask, smooth_map
     type(coop_healpix_maps)::lMT, lcT, lM, sim
     COOP_REAL, dimension(:,:), allocatable::Fmean
     COOP_REAL, dimension(:,:), allocatable::Ffluc     
     COOP_SINGLE, dimension(:), allocatable::mean
     COOP_INT, dimension(:), allocatable::indMT, indCT, iiMT, iiCT
   contains
     procedure::free => coop_healpix_inpaint_free
     procedure::init => coop_healpix_inpaint_init
     procedure::upgrade => coop_healpix_inpaint_upgrade     
     procedure::set_corr => coop_healpix_inpaint_set_corr
     procedure::eval_cov => coop_healpix_inpaint_eval_cov
     procedure::correlation => coop_healpix_inpaint_correlation
     procedure::lask2mask => coop_healpix_inpaint_lask2mask
  end type coop_healpix_inpaint

  type coop_healpix_correlation_function
     COOP_INT::lmax = -1
     COOP_REAL,dimension(:),allocatable::als, bls, Cls, Cls_smooth
     COOP_REAL::fwhm = 0.d0
   contains
     procedure::free => coop_healpix_correlation_function_free
     procedure::init => coop_healpix_correlation_function_init
     procedure::simulate => coop_healpix_correlation_function_simulate
     procedure::set_beam => coop_healpix_correlation_function_set_beam
     procedure::corr => coop_healpix_correlation_function_corr
     procedure::corr2int => coop_healpix_correlation_function_corr2int
  end type coop_healpix_correlation_function
  


#define COS2RADIUS(cosx) (sqrt(2.d0*(1.d0 - (cosx))))
#define RADIUS2COS(r)  (1.d0-(r)**2/2.d0)

contains

  function coop_unit2muK(unit) result(fac)
    COOP_UNKNOWN_STRING::unit
    COOP_REAL::fac
    if(unit(1:1).eq."K")then
       fac = 1.d6
       return
    endif
    if(unit(1:2).eq."mK")then
       fac = 1.d3
       return
    endif
    fac = 1.d0
  end function coop_unit2muK

  subroutine coop_fits_to_header(filename, header)
    COOP_UNKNOWN_STRING::filename
    COOP_STRING::cfname
    type(coop_dictionary)::header
    COOP_LONG_STRING::str
    COOP_INT nkeys, i, j, istart, iend, ikey
    call header%free()
    cfname  = trim(adjustl(filename))
    str = ""
    call coop_convert_to_C_String(cfname)
    call coop_fits_read_all_headers_to_string(cfname, str, nkeys)
    call coop_convert_to_Fortran_String(str)
    istart = 1
    do i=1, nkeys
       j = scan(str(istart:),"=")
       iend = scan(str(istart:), coop_newline)
       j = j + istart - 1
       iend = iend + istart - 2
       call header%insert(str(istart:j-1), trim(coop_string_strip_quotes(str(j+1:iend))))
       istart = iend+2
    enddo
  end subroutine coop_fits_to_header


  subroutine coop_healpix_maps_set_unit(this, imap, unit)
    class(coop_healpix_maps)::this
    COOP_INT::imap
    COOP_UNKNOWN_STRING::unit
    this%units(imap) = trim(adjustl(unit))
    call this%header%insert("TUNIT"//COOP_STR_OF(imap), trim(unit), overwrite=.true.)
  end subroutine coop_healpix_maps_set_unit


  subroutine coop_healpix_maps_set_units(this, unit)
    class(coop_healpix_maps)::this
    COOP_UNKNOWN_STRING::unit
    COOP_INT i
    do i=1, this%nmaps
       call this%set_unit(i, unit)
    enddo
  end subroutine coop_healpix_maps_set_units
  

  subroutine coop_healpix_maps_set_field(this, imap, field)
    class(coop_healpix_maps)::this
    COOP_INT::imap
    COOP_UNKNOWN_STRING::field
    call this%header%insert("TTYPE"//COOP_STR_OF(imap), trim(adjustl(field)), overwrite=.true.)
    this%fields(imap) = trim(adjustl(field))
    call this%fields_to_spins(imap)
  end subroutine coop_healpix_maps_set_field

  subroutine coop_healpix_maps_set_fields(this, field)
    class(coop_healpix_maps)::this
    COOP_INT::imap
    COOP_UNKNOWN_STRING::field
    do imap = 1, this%nmaps
       call this%set_field(imap, field)
    enddo
  end subroutine coop_healpix_maps_set_fields
  
  

  subroutine coop_healpix_maps_ang2pix(this, theta, phi, pix)
    class(coop_healpix_maps)::this
    COOP_REAL theta, phi
    COOP_INT pix
#ifdef HAS_HEALPIX    
    if(this%ordering .eq. COOP_RING)then
       call ang2pix_ring(this%nside, theta, phi, pix)
    else
       call ang2pix_nest(this%nside, theta, phi, pix)       
    endif
#else
    stop "You need to compile with Healpix to use ang2pix" 
#endif       
  end subroutine coop_healpix_maps_ang2pix

  subroutine coop_healpix_maps_pix2ang(this, pix, theta, phi)
    class(coop_healpix_maps)::this
    COOP_REAL theta, phi
    COOP_INT pix
#ifdef HAS_HEALPIX    
    if(this%ordering .eq. COOP_RING)then
       call pix2ang_ring(this%nside, pix, theta, phi)
    else
       call pix2ang_nest(this%nside, pix, theta, phi)
    endif
#else
    stop "You need to compile with Healpix to use pix2ang" 
#endif       
  end subroutine coop_healpix_maps_pix2ang
  

    subroutine coop_healpix_maps_vec2pix(this, vec, pix)
    class(coop_healpix_maps)::this
    COOP_REAL vec(3)
    COOP_INT pix
#ifdef HAS_HEALPIX    
    if(this%ordering .eq. COOP_RING)then
       call vec2pix_ring(this%nside, vec, pix)
    else
       call vec2pix_nest(this%nside, vec, pix)       
    endif
#else
    stop "You need to compile with Healpix to use ang2pix" 
#endif       
  end subroutine coop_healpix_maps_vec2pix

  subroutine coop_healpix_maps_pix2vec(this, pix, vec)
    class(coop_healpix_maps)::this
    COOP_REAL vec(3)
    COOP_INT pix
#ifdef HAS_HEALPIX    
    if(this%ordering .eq. COOP_RING)then
       call pix2vec_ring(this%nside, pix, vec)
    else
       call pix2vec_nest(this%nside, pix, vec)
    endif
#else
    stop "You need to compile with Healpix to use pix2ang" 
#endif       
  end subroutine coop_healpix_maps_pix2vec
  
  subroutine coop_healpix_spot_cut_mask(nside, l_deg, b_deg, r_deg, filename)
    COOP_INT::nside
    COOP_REAL::l_deg, b_deg, r_deg
    COOP_UNKNOWN_STRING::filename
    type(coop_healpix_maps)::mask
    COOP_INT::listpix(0:nside**2*12-1), nlist, pix
    COOP_REAL::theta, phi
#ifdef HAS_HEALPIX
    call mask%init(nside = nside, nmaps = 1, genre="MASK")
    mask%map(:,1) = 1.
    call coop_healpix_lb2ang(l_deg, b_deg, theta, phi)    
    call mask%ang2pix(theta, phi, pix)
    call mask%query_disc(pix, r_deg*coop_SI_degree, listpix, nlist)
    mask%map(listpix(0:nlist-1), 1) = 0.
    call mask%write(trim(filename))
    call mask%free()
#else
    stop "cannot find HEALPIX"
#endif    
  end subroutine coop_healpix_spot_cut_mask

  subroutine coop_healpix_maps_mask_disc(this, l_deg, b_deg, r_deg)
    class(coop_healpix_maps)::this
    COOP_REAL::l_deg, b_deg, r_deg
    COOP_INT::listpix(0:this%npix-1), nlist, pix
    COOP_REAL::theta, phi
    call coop_healpix_lb2ang(l_deg, b_deg, theta, phi)    
    call this%ang2pix(theta, phi, pix)
    call this%query_disc(pix, r_deg*coop_SI_degree, listpix, nlist)
    this%map(listpix(0:nlist-1), :) = 0.
  end subroutine coop_healpix_maps_mask_disc


  subroutine coop_healpix_maps_draw_latitude_line(this, latitude_deg, width_deg)
    class(coop_healpix_maps)::this
    COOP_REAL::latitude_deg
    COOP_REAL,optional::width_deg
    COOP_INT::listpix(0:this%npix-1), nlist
    call this%query_strip(coop_pio2-(latitude_deg+width_deg/2.d0)*coop_SI_degree,  coop_pio2-(latitude_deg-width_deg/2.d0)*coop_SI_degree, listpix, nlist)
    this%map(listpix(0:nlist-1), :) = this%bad_data
  end subroutine coop_healpix_maps_draw_latitude_line


  subroutine coop_healpix_maps_mask_strip(this, l_deg, b_deg, r1_deg, r2_deg)
    class(coop_healpix_maps)::this
    COOP_REAL::l_deg, b_deg, r1_deg, r2_deg
    COOP_INT::listpix1(0:this%npix-1), nlist1, pix, listpix2(0:this%npix-1), nlist2
    COOP_REAL::theta, phi
    type(coop_healpix_maps)::mask
    call mask%init(nside = this%nside, nmaps = 1, genre = "MASK", nested = (this%ordering .eq. COOP_NESTED) )
    call coop_healpix_lb2ang(l_deg, b_deg, theta, phi)    
    call mask%ang2pix(theta, phi, pix)
    call mask%query_disc(pix, r1_deg*coop_SI_degree, listpix1, nlist1)
    call mask%query_disc(pix, r2_deg*coop_SI_degree, listpix2, nlist2)
    mask%map = 1.
    if(nlist1 .lt. nlist2)then
       mask%map( listpix2(0:nlist2-1), 1) = 0.
       mask%map( listpix1(0:nlist1-1), 1) = 1.
    else
       mask%map( listpix1(0:nlist1-1), 1) = 0.
       mask%map( listpix2(0:nlist2-1), 1) = 1.       
    endif
    call this%apply_mask(mask)
    call mask%free
  end subroutine coop_healpix_maps_mask_strip
  


  subroutine coop_healpix_maps_mask_pixel(this, nside, ipix)
    class(coop_healpix_maps)::this
    COOP_INT::nside, ipix
    if(nside.gt.this%nside)return    
    call this%convert2nested()
    this%map(ipix*(this%nside/nside)**2:(ipix+1)*(this%nside/nside)**2-1, :) = 0.
  end subroutine coop_healpix_maps_mask_pixel
  


  subroutine coop_healpix_spot_select_mask(nside, l_deg, b_deg, r_deg, filename)
    COOP_INT::nside
    COOP_REAL::l_deg, b_deg, r_deg
    COOP_UNKNOWN_STRING::filename
    type(coop_healpix_maps)::mask
    COOP_INT::listpix(0:nside**2*12-1), nlist, pix
    COOP_REAL::theta, phi
#ifdef HAS_HEALPIX
    call coop_healpix_lb2ang(l_deg, b_deg, theta, phi)
    call mask%init(nside = nside, nmaps = 1, genre="MASK")
    mask%map(:,1) = 0.
    call mask%ang2pix(theta, phi, pix)
    call mask%query_disc(pix, r_deg*coop_SI_degree, listpix, nlist)
    mask%map(listpix(0:nlist-1), 1) = 1.
    call mask%write(trim(filename))
    call mask%free()
#else
    stop "cannot find HEALPIX"
#endif    
  end subroutine coop_healpix_spot_select_mask
  

  subroutine coop_healpix_merge_masks(mask1, mask2, maskout)
    COOP_UNKNOWN_STRING::mask1,mask2,maskout
    type(coop_healpix_maps)::m1,m2
    call m1%read(mask1, nmaps_wanted=1)
    call m2%read(mask2, nmaps_wanted=1)
    m1%map = m1%map * m2%map
    call m1%set_unit(1, "1")
    call m1%set_field(1, "MASK")
    call m1%write(maskout)
  end subroutine coop_healpix_merge_masks

  subroutine coop_healpix_maps_generate_latcut_mask(this, nside, latitude_deg, depth_deg)
    class(coop_healpix_maps)::this
    COOP_INT::nside    
    COOP_REAL::latitude_deg
    COOP_REAL::theta, phi, lat, dep
    COOP_INT::i, nlist
    COOP_INT::listpix(0:12*nside**2-1)
    COOP_REAL, optional::depth_deg
    call this%init(nside = nside, nmaps = 1, genre="MASK")    
    this%map(:,1) = 1.
    lat = latitude_deg*coop_SI_degree    
    call query_strip(this%nside, coop_pio2 - lat, coop_pio2 + lat, listpix, nlist, nest = 0, inclusive = 1)
    if(present(depth_deg))then
       dep = depth_deg*coop_SI_degree       
       do i = 0, nlist-1
          call pix2ang_ring(this%nside, listpix(i), theta, phi)
          this%map(listpix(i), 1) = exp(- ((lat - abs(theta - coop_pio2))/dep)**2)
       enddo       
    else
       this%map(listpix(0:nlist-1), 1) = 0.
    endif
  end subroutine coop_healpix_maps_generate_latcut_mask


  subroutine coop_healpix_patch_plot(this, imap, output, use_degree)
    COOP_INT::bgrids
    class(coop_healpix_patch)::this
    COOP_INT imap
    COOP_UNKNOWN_STRING::output
    type(coop_asy)::fig
    COOP_INT nb, i, j, k, ns, iq, iu
    COOP_REAL  xc, yc,  norm, r, theta, minz, maxz
    COOP_REAL,dimension(:),allocatable::xstart, xend, ystart, yend
    logical,optional::use_degree
    logical use_rad
    COOP_SHORT_STRING::xlabel, ylabel
    if(imap .le. 0 .or. imap .gt. this%nmaps) stop "coop_healpix_patch_plot: imap overflow"    
    call fig%open(output)
    if(present(use_degree))then
       if(use_degree)then
          this%r = this%r/coop_SI_degree
       endif
       use_rad = .not. use_degree
    else
       use_rad = .true.
    endif
    if(use_rad)then
       xlabel = "$\varpi\cos\phi$"
       ylabel = "$\varpi\sin\phi$"       
    else
       xlabel = "$\varpi\cos\phi (\mathrm{deg})$"
       ylabel = "$\varpi\sin\phi (\mathrm{deg})$"
    endif
    this%caption = adjustl(this%caption)
    if(coop_healpix_patch_default_want_caption)then
       if(len_trim(this%caption) .gt. 60)then
          call fig%init(caption = "{\tiny "//trim(this%caption)//"}", xlabel =trim(xlabel), ylabel =trim(ylabel), width = coop_healpix_patch_default_figure_width, height = coop_healpix_patch_default_figure_height, xmin = -real(this%r(this%n)), xmax = real(this%r(this%n)), ymin = -real(this%r(this%n)), ymax = real(this%r(this%n)))
       elseif(len_trim(this%caption) .gt. 50)then
          call fig%init(caption = "{\scriptsize "//trim(this%caption)//"}", xlabel =trim(xlabel), ylabel =trim(ylabel), width = coop_healpix_patch_default_figure_width, height = coop_healpix_patch_default_figure_height, xmin = -real(this%r(this%n)), xmax = real(this%r(this%n)), ymin = -real(this%r(this%n)), ymax = real(this%r(this%n)))
       elseif(len_trim(this%caption) .gt. 40)then
          call fig%init(caption = "{\small "//trim(this%caption)//"}", xlabel =trim(xlabel), ylabel =trim(ylabel), width = coop_healpix_patch_default_figure_width, height = coop_healpix_patch_default_figure_height, xmin = -real(this%r(this%n)), xmax = real(this%r(this%n)), ymin = -real(this%r(this%n)), ymax = real(this%r(this%n)))       
       else
          call fig%init(caption = trim(this%caption), xlabel =trim(xlabel), ylabel =trim(ylabel), width = coop_healpix_patch_default_figure_width, height = coop_healpix_patch_default_figure_height, xmin = -real(this%r(this%n)), xmax = real(this%r(this%n)), ymin = -real(this%r(this%n)), ymax = real(this%r(this%n)))
       endif
    else
       call fig%init(xlabel =trim(xlabel), ylabel =trim(ylabel), width = coop_healpix_patch_default_figure_width, height = coop_healpix_patch_default_figure_height, xmin = -real(this%r(this%n)), xmax = real(this%r(this%n)), ymin = -real(this%r(this%n)), ymax = real(this%r(this%n)))          
    endif
    if(this%tbs%zmin(imap) .lt.0.99e30)then
       minz = this%tbs%zmin(imap)
    else
       call coop_array_get_threshold(this%image(:,:,imap), COOP_REAL_OF(0.99), minz)
    endif
    if(this%tbs%zmax(imap) .gt. -0.99e30)then
       maxz = this%tbs%zmax(imap)
    else
       call coop_array_get_threshold(this%image(:,:,imap), COOP_REAL_OF(0.01), maxz)
    endif
    if(coop_healpix_patch_default_want_label)then
       call coop_asy_density(fig, this%image(:,:,imap), -this%r(this%n), this%r(this%n), -this%r(this%n), this%r(this%n), label = trim(this%tbs%label(imap)), zmax = maxz, zmin = minz, color_table = trim(this%color_table))
    else
       call coop_asy_density(fig, this%image(:,:,imap), -this%r(this%n), this%r(this%n), -this%r(this%n), this%r(this%n), zmax = maxz, zmin = minz, color_table = trim(this%color_table))       
    endif
    if(use_rad .and. coop_healpix_patch_default_want_arrow)then
       theta = nint(2.d0*asin(this%r(this%n)/2.d0)/coop_SI_degree*10.d0)/10.d0
       call coop_asy_label(fig, "$\mathbf{-"//COOP_STR_OF(theta)//"}^\circ$", -this%r(this%n)*1.01, -this%r(this%n)*1.23, color="blue")
       call coop_asy_label(fig, "$\mathbf{"//COOP_STR_OF(theta)//"}^\circ$", this%r(this%n)*1.01, -this%r(this%n)*1.23, color="blue")
       call fig%arrow(this%r(this%n),  -this%r(this%n)*1.17, this%r(this%n),  -this%r(this%n)*1.13)
       call fig%arrow(-this%r(this%n),  -this%r(this%n)*1.17, -this%r(this%n),  -this%r(this%n)*1.13)
    endif
    if(this%tbs%spin(imap).eq.2 .and. this%tbs%headless_vector(imap) .and. this%nmaps .ge. 2)then
       if(imap .eq. this%nmaps)then
          iq = imap-1
          iu = imap
       elseif(this%tbs%spin(imap+1) .eq. 2)then
          iq = imap
          iu = imap+1
       elseif(imap .gt. 1)then
          if(this%tbs%spin(imap-1).eq.2)then
             iq = imap-1
             iu = imap
          else
             stop "patch_plot: spin 2 must go in pairs"
          endif
       else
          stop "patch_plot: spin 2 must go in pairs"
       endif
       bgrids = max(this%n/7, 2)
       norm = maxval(this%image(:,:,iq)**2+this%image(:,:,iu)**2)
       if(norm .gt. 0.d0)then
          norm = bgrids*this%dr/2./sqrt(norm)*0.96
       else
          goto 100
       endif
       ns = floor((this%n-0.5d0*bgrids)/bgrids)
       nb = (2*ns+1)**2
       allocate(xstart(nb),  ystart(nb), xend(nb), yend(nb))
       k = 0
       ns = ns*bgrids
       do j = -ns, ns, bgrids
          do i = -ns, ns, bgrids
             xc = i*this%dr
             yc = j*this%dr
             r = sqrt(this%image(i,j,iq)**2+this%image(i,j,iu)**2)*norm
             if(coop_healpix_IAU_headless_vector)then
                theta = 0.5d0*COOP_POLAR_ANGLE(this%image(i,j,iq), this%image(i,j,iu)) + coop_pio2                   
             else
                theta = 0.5d0*COOP_POLAR_ANGLE(this%image(i,j,iq), this%image(i,j,iu))
             endif
             if(this%tbs%local_rotation(imap))then
                theta  = theta + COOP_POLAR_ANGLE(xc, yc) + (1.d0-coop_healpix_QrUrSign)*coop_pio4
             endif
             k = k + 1
             xstart(k) = xc - r*cos(theta)
             ystart(k) = yc - r*sin(theta)
             xend(k) = 2*xc - xstart(k)
             yend(k) = 2*yc - ystart(k)
          enddo
       enddo
       call coop_asy_lines(fig, xstart, ystart, xend, yend, "black", "solid", 2.)
       deallocate(xstart, xend, ystart, yend)
    endif
100 call fig%close()
    if(present(use_degree))then
       if(use_degree)then
          this%r = this%r*coop_SI_degree
       endif
    endif    
  end subroutine coop_healpix_patch_plot

  subroutine coop_healpix_patch_export(this, filename)
    class(coop_healpix_patch)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file)::fp
    call fp%open(trim(filename), "u")
    write(fp%unit) this%color_table
    write(fp%unit) this%genre
    write(fp%unit) this%n
    write(fp%unit) this%dr
    write(fp%unit) this%mmax    
    write(fp%unit) this%nstack
    write(fp%unit) this%nstack_raw
    write(fp%unit) this%image
    call fp%close()
  end subroutine coop_healpix_patch_export

  subroutine coop_healpix_patch_import(this, filename)
    class(coop_healpix_patch)::this
    COOP_UNKNOWN_STRING::filename
    type(coop_file)::fp
    call fp%open(trim(filename), "ur")
    read(fp%unit) this%color_table
    read(fp%unit) this%genre
    read(fp%unit) this%n
    read(fp%unit) this%dr        
    read(fp%unit) this%mmax
    call this%init(genre = this%genre, n = this%n, dr = this%dr, mmax = this%mmax)
    read(fp%unit) this%nstack   
    read(fp%unit) this%nstack_raw
    read(fp%unit) this%image
    call this%get_all_radial_profiles()
    call fp%close()
  end subroutine coop_healpix_patch_import

  subroutine coop_healpix_patch_fft(this, imap)
    class(coop_healpix_patch)::this
    COOP_INT::imap
    COOP_INT::  j, ns
    if(imap .gt. this%nmaps .or. this%n .le. 0) stop "patch_fft: imap overflow"
    ns = this%n*2+1
    !$omp critical
    call coop_fft_forward(ns, ns, this%image(-this%n:this%n,-this%n:this%n,imap), this%fftimage(0:this%n, 0:this%n*2, imap))
    !$omp end critical
    this%image(1:this%n, 0:this%n, imap) = real(this%fftimage(1:this%n, 0:this%n, imap))
    this%image(1:this%n, -this%n:-1, imap) = real(this%fftimage(1:this%n, this%n+1:this%n*2, imap))
    this%image(-this%n:-1, 0:this%n, imap) = aimag(this%fftimage(this%n:1:-1, 0:this%n, imap))
    this%image(-this%n:-1, -this%n:-1, imap) = aimag(this%fftimage(this%n:1:-1, this%n+1:this%n*2, imap))
    this%image(0, 0, imap) = real(this%fftimage(0, 0, imap))
    do j=1, this%n
       this%image(0, j, imap) = real(this%fftimage(0, j, imap))
       this%image(0, -j, imap) = aimag(this%fftimage(0, j, imap))
    enddo
  end subroutine coop_healpix_patch_fft

  subroutine coop_healpix_patch_plot_fft(this, imap, output,label)
    COOP_INT,parameter::lmax=100
    COOP_INT::bgrids
    class(coop_healpix_patch)::this
    COOP_INT imap
    COOP_UNKNOWN_STRING::output
    type(coop_asy)::fig
    COOP_INT i,j,ns
    COOP_REAL   minz, maxz, dk
    COOP_COMPLEX, dimension(:,:),allocatable:: fftmap
    COOP_REAL,dimension(:,:),allocatable::remap,immap
    COOP_SINGLE,parameter::width=6.5, height=5.
    COOP_SHORT_STRING::xlabel, ylabel
    COOP_UNKNOWN_STRING,optional::label    
    
    if(imap .le. 0 .or. imap .gt. this%nmaps) stop "coop_healpix_patch_plot: imap overflow"

    allocate(fftmap(0:this%n, 0:this%n*2))
    !$omp critical
    call coop_fft_forward(this%n*2+1, this%n*2+1, this%image(:,:,imap), fftmap)
    !$omp end critical
    fftmap = fftmap/(2*this%n+1)**2

    dk = coop_2pi/(this%dr*(2*this%n+1))
    ns = min(nint(lmax/dk), nint(this%n/coop_sqrt2))

    allocate(remap(-ns:ns,-ns:ns), immap(-ns:ns,-ns:ns))    
    remap(0:ns, 0:ns) = real(fftmap(0:ns, 0:ns) )
    remap(0:ns, -ns:-1) = real(fftmap(0:ns, this%n*2+1-ns:this%n*2))
    
    immap(0:ns, 0:ns) = aimag(fftmap(0:ns, 0:ns) )
    immap(0:ns, -ns:-1) = aimag(fftmap(0:ns, this%n*2+1-ns:this%n*2))

    remap(-ns:-1,-ns:ns)=remap(ns:1:-1,ns:-ns:-1)
    immap(-ns:-1,-ns:ns)=-immap(ns:1:-1,ns:-ns:-1)

    
    call fig%open(coop_file_add_postfix(output, "_FFTRe"))
    
    xlabel = "$\ell_x$"
    ylabel = "$\ell_y$"
    call fig%init(caption = trim(this%caption)//", Re($f_\ell$)", xlabel =trim(xlabel), ylabel =trim(ylabel), width = width, height = height, xmin = -real(ns*dk), xmax = real(ns*dk), ymin = -real(ns*dk), ymax = real(ns*dk))
    call coop_array_get_threshold(remap, COOP_REAL_OF(0.995), minz)
    call coop_array_get_threshold(remap, COOP_REAL_OF(0.005), maxz)
    if(present(label))then
       call coop_asy_density(fig, remap, -ns*dk, ns*dk,-ns*dk, ns*dk, label = "$\mathrm{Re} "//trim(adjustl(label))//"_\ell$", zmax = maxz, zmin = minz, color_table = trim(this%color_table))
    else
       call coop_asy_density(fig, remap, -ns*dk, ns*dk,-ns*dk, ns*dk, label = trim(adjustl(this%tbs%label(imap))), zmax = maxz, zmin = minz, color_table = trim(this%color_table))       
    endif
    call fig%close()

    call fig%open(coop_file_add_postfix(output, "_FFTIm"))    
    xlabel = "$\ell_x$"
    ylabel = "$\ell_y$"
    call fig%init(caption = trim(this%caption)//", Im($f_\ell$)", xlabel =trim(xlabel), ylabel =trim(ylabel), width = width, height = height, xmin = -real(ns*dk), xmax = real(ns*dk), ymin = -real(ns*dk), ymax = real(ns*dk))
    call coop_array_get_threshold(immap, COOP_REAL_OF(0.995), minz)
    call coop_array_get_threshold(immap, COOP_REAL_OF(0.005), maxz)
    if(present(label))then
       call coop_asy_density(fig, immap, -ns*dk, ns*dk,-ns*dk, ns*dk, zmax = maxz, zmin = minz, label = "$\mathrm{Im}"//trim(adjustl(label))//"_\ell$", color_table = trim(this%color_table))
    else
       call coop_asy_density(fig, immap, -ns*dk, ns*dk,-ns*dk, ns*dk, zmax = maxz, zmin = minz, label =trim(adjustl(this%tbs%label(imap))), color_table = trim(this%color_table))       
    endif
    call fig%close()

    call fig%open(coop_file_add_postfix(output, "_FFTl2Cl"))
    immap = (remap**2+immap**2)*dk**2
    do i=-ns,ns
       do j=-ns,ns
          immap(i,j)=immap(i,j)*dble(i**2+j**2)
       enddo
    enddo
    xlabel = "$\ell_x$"
    ylabel = "$\ell_y$"
    call fig%init(caption = trim(this%caption)//", $\ell^2|f_\ell|^2$", xlabel =trim(xlabel), ylabel =trim(ylabel), width = width, height = height, xmin = -real(ns*dk), xmax = real(ns*dk), ymin = -real(ns*dk), ymax = real(ns*dk))
    call coop_array_get_threshold(immap, COOP_REAL_OF(0.995), minz)
    call coop_array_get_threshold(immap, COOP_REAL_OF(0.005), maxz)
    if(present(label))then
       call coop_asy_density(fig, immap, -ns*dk, ns*dk, -ns*dk, ns*dk, zmax = maxz, zmin = minz, label = "$\ell^2 |"//trim(adjustl(label))//"_\ell|^2 $", color_table = trim(this%color_table))       
    else
       call coop_asy_density(fig, immap, -ns*dk, ns*dk, -ns*dk, ns*dk, zmax = maxz, zmin = minz, label = trim(this%tbs%label(imap)), color_table = trim(this%color_table))
    endif
    call fig%close()    
    
    deallocate(fftmap,remap,immap)
  end subroutine coop_healpix_patch_plot_fft
 

  subroutine coop_healpix_patch_free(this)
    class(coop_healpix_patch) this
    if(allocated(this%image))deallocate(this%image)
    if(allocated(this%fftimage))deallocate(this%fftimage)    
    if(allocated(this%r))deallocate(this%r)
    if(allocated(this%fr))deallocate(this%fr)
    if(allocated(this%wcm))deallocate(this%wcm)
    if(allocated(this%wsm))deallocate(this%wsm)
    if(allocated(this%icm))deallocate(this%icm)
    if(allocated(this%nstack))deallocate(this%nstack)
    if(allocated(this%indisk))deallocate(this%indisk)
    call this%tbs%free()
    this%nmaps = 0
  end subroutine coop_healpix_patch_free

  subroutine coop_healpix_patch_init(this, genre, n, dr, mmax, ny)
    class(coop_healpix_patch) this
    logical,parameter::do_norm = .true.
    logical,parameter::remove_m0 = .true.
    COOP_UNKNOWN_STRING::genre
    COOP_INT n, nsq
    COOP_REAL dr, cosmt, sinmt, theta, rij
    COOP_INT i,j,m
    COOP_INT, optional::mmax, ny
    COOP_REAL sumrc(0:n+1), sumrs(0:n+1), weight(0:n+1)
    call this%free()
    this%genre = trim(genre)
    this%n = n
    select type(this)
    type is(coop_healpix_filament)
       if(present(ny))then       
          this%ny = ny
       else
          stop "filament_init needs ny input"
       endif
    end select
    this%dr = dr
    if(present(mmax))then
       this%mmax = mmax
    else
       this%mmax = 4
    endif
    call this%tbs%init(genre)
    this%nmaps = this%tbs%nmaps    
    this%npix = (2*this%n+1)**2
    if(this%n .lt. 0) return
    select type(this)
    type is(coop_healpix_patch)
       allocate(this%image(-this%n:this%n, -this%n:this%n, this%nmaps))
       allocate(this%fftimage(0:this%n, 0:this%n*2, this%nmaps))       
       allocate(this%nstack(-this%n:this%n, -this%n:this%n))
       allocate(this%indisk(-this%n:this%n, -this%n:this%n))
       allocate(this%r(0:this%n))
       this%image = 0.d0
       this%nstack = 0.d0
       this%indisk = 1.d0
       this%nstack_raw = 0
       !$omp parallel do
       do i=0, this%n
          this%r(i) = this%dr * i
       enddo
       !$omp end parallel do
    type is(coop_healpix_filament)
       allocate(this%image(-this%n:this%n, -this%ny:this%ny, this%nmaps))
       allocate(this%fftimage(0:this%n, 0:this%ny*2, this%nmaps))              
       allocate(this%nstack(-this%n:this%n,  -this%ny:this%ny))
       allocate(this%indisk(-this%n:this%n, -this%ny:this%ny))
       allocate(this%r(0:max(this%n, this%ny)))
       this%image = 0.d0
       this%nstack = 0
       this%indisk = 1.d0
       !$omp parallel do
       do i=0, this%n
          this%r(i) = this%dr * i
       enddo
       !$omp end parallel do
       return
    end select
    allocate(this%fr(0:this%n, 0:this%mmax/2, this%nmaps))
    allocate(this%wcm(-this%n:this%n, -this%n:this%n, 0:this%mmax+1))
    allocate(this%wsm(-this%n:this%n, -this%n:this%n, 0:this%mmax+1))
    allocate(this%icm(-this%n:this%n, -this%n:this%n, 0:1))
    this%wcm = 0.d0
    this%wsm = 0.d0
    this%wcm(0,0,0) = 1.d0
    this%fr = 0.d0   
    nsq = this%n**2
    do j=1, this%n
       i = ceiling(sqrt(dble((this%n-j)*(this%n+j))))
       if(i**2+j**2 .le. nsq) i = i + 1
       this%indisk(i:this%n, j) = 0.d0
       this%indisk(-this%n:-i, j) = 0.d0
       this%indisk(i:this%n, -j) = 0.d0
       this%indisk(-this%n:-i, -j) = 0.d0
    enddo

    this%num_indisk_tol = floor(count(this%indisk .gt. 0.d0)*coop_healpix_mask_tol)+0.5d0
    !!------------------------------------------------------------------
    !!start doing m=0 
    !$omp parallel do private(i, j, rij)
    do j=-this%n, this%n
       do i=-this%n, this%n
          rij =  sqrt(dble(i**2+j**2))
          this%icm(i, j, 0) = floor(rij)
          rij = rij - this%icm(i,j,0)
          this%icm(i, j, 1) = this%icm(i, j, 0) + 1
          if(this%icm(i,j,0) .gt. n) &
               this%icm(i,j,0) = 0
          if(this%icm(i,j,1) .gt. n) &
               this%icm(i,j,1) = 0
          if(this%icm(i,j,0) .ne. 0) &
               this%wcm(i, j, 0) = 1.d0  - rij
          if(this%icm(i, j, 1) .ne. 0) &
               this%wcm(i, j, 1) = rij
       enddo
    enddo
    !$omp end parallel do

    sumrc = 1.d-30
    
    do j=-this%n, this%n
       do i=-this%n, this%n
          sumrc(this%icm(i,j,0)) = sumrc(this%icm(i,j,0)) + this%wcm(i,j,0)
          sumrc(this%icm(i,j,1)) = sumrc(this%icm(i,j,1)) + this%wcm(i,j,1)
       enddo
    enddo
    
    !$omp parallel do private(i, j)
    do j=-this%n, this%n
       do i=-this%n, this%n
          this%wcm(i, j, 0) = this%wcm(i, j, 0)/sumrc(this%icm(i, j, 0))
          this%wcm(i, j, 1) = this%wcm(i, j, 1)/sumrc(this%icm(i, j, 1))
       enddo
    enddo
    !$omp end parallel do

    !!m = 0 is done
    !!------------------------------------------------------------------
    
    do m = 2, this%mmax, 2
       do j=-this%n, this%n
          do i = -this%n, this%n
             if(this%icm(i,j,0).ne.0)then
                theta = atan2(dble(j), dble(i))
                cosmt = cos(m*theta)*2.d0
                sinmt = sin(m*theta)*2.d0
                this%wcm(i, j, m) = this%wcm(i, j, 0)*cosmt
                this%wcm(i, j, m+1) = this%wcm(i, j, 1)*cosmt
                this%wsm(i, j, m) = this%wcm(i, j, 0)*sinmt
                this%wsm(i, j, m+1) = this%wcm(i, j, 1)*sinmt
             endif
          enddo
       enddo

       if(do_norm)then
          !!now normalize
          sumrc = 1.d-6
          sumrs = 1.d-6
          do j=-this%n, this%n
             do i=-this%n, this%n
                if(this%icm(i, j, 0) .ne. 0)then
                   theta = atan2(dble(j), dble(i))
                   cosmt = cos(m*theta)
                   sumrc(this%icm(i,j,0)) = sumrc(this%icm(i,j,0)) + this%wcm(i,j,m)*cosmt
                   sumrc(this%icm(i,j,1)) = sumrc(this%icm(i,j,1)) + this%wcm(i,j,m+1)*cosmt
                   sinmt = sin(m*theta)                
                   sumrs(this%icm(i,j,0)) = sumrs(this%icm(i,j,0)) + this%wsm(i,j,m)*sinmt                

                   sumrs(this%icm(i,j,1)) = sumrs(this%icm(i,j,1)) + this%wsm(i,j,m+1)*sinmt
                endif
             enddo
          enddo
          !$omp parallel do private(i, j)
          do j=-this%n, this%n
             do i=-this%n, this%n
                if(this%icm(i, j, 0).ne. 0)then
                   if(sumrc(this%icm(i, j, 0)) .gt. 0.01d0) this%wcm(i, j, m) = this%wcm(i, j, m)/sumrc(this%icm(i, j, 0))
                   if(sumrs(this%icm(i, j, 0)) .gt. 0.01d0) this%wsm(i, j, m) = this%wsm(i, j, m)/sumrs(this%icm(i, j, 0))
                   if(sumrc(this%icm(i, j, 1)) .gt. 0.01d0) this%wcm(i, j, m+1) = this%wcm(i, j, m+1)/sumrc(this%icm(i, j, 1))
                   if(sumrs(this%icm(i, j, 1)) .gt. 0.01d0) this%wsm(i, j, m+1) = this%wsm(i, j, m+1)/sumrs(this%icm(i, j, 1))
                endif
             enddo
          enddo
          !$omp end parallel do
       end if
       if(remove_m0)then             !!remove m = 0 degeneracy
          sumrc = 0.d0
          sumrs = 0.d0
          do j=-this%n, this%n
             do i = -this%n, this%n
                if(this%icm(i,j,0).ne.0)then
                   theta = atan2(dble(j), dble(i))
                   sumrc(this%icm(i, j, 0)) = sumrc(this%icm(i, j, 0)) + this%wcm(i, j, m)
                   sumrc(this%icm(i, j, 1)) = sumrc(this%icm(i, j, 1)) + this%wcm(i, j, m+1)                
                   sumrs(this%icm(i, j, 0)) = sumrs(this%icm(i, j, 0)) + this%wsm(i, j, m)
                   sumrs(this%icm(i, j, 1)) = sumrs(this%icm(i, j, 1)) + this%wsm(i, j, m+1)
                endif
             enddo
          enddo
          !$omp parallel do private(i, j)
          do j=-this%n, this%n
             do i=-this%n, this%n
                if(this%icm(i, j, 0) .ne. 0)then
                   this%wcm(i, j, m) = this%wcm(i, j, m) - sumrc(this%icm(i, j, 0))*this%wcm(i, j, 0)
                   this%wcm(i, j, m+1) = this%wcm(i, j, m+1) - sumrc(this%icm(i, j, 1))*this%wcm(i, j, 1)
                   this%wsm(i, j, m) = this%wsm(i, j, m) - sumrs(this%icm(i, j, 0))*this%wcm(i, j, 0)
                   this%wsm(i, j, m+1) = this%wsm(i, j, m+1) - sumrs(this%icm(i, j, 1))*this%wcm(i, j, 1)                

                endif
             enddo
          enddo
          !$omp end parallel do
       endif
       
    enddo
    
  end subroutine coop_healpix_patch_init


  subroutine coop_healpix_patch_get_all_radial_profiles(this)
    class(coop_healpix_patch)::this
    COOP_INT i, j, imap, m
    if(this%mmax .lt. 0) return
    do imap = 1, this%nmaps
       do m = 0, this%mmax, 2
          call this%get_radial_profile(imap, m)
       enddo
    enddo
  end subroutine coop_healpix_patch_get_all_radial_profiles

  subroutine coop_healpix_patch_get_radial_profile(this, imap, m)
    class(coop_healpix_patch)::this
    COOP_INT m, imap, halfm
    COOP_INT i,j
    if(m.gt. this%mmax .or. mod(m,2).ne.0 .or. imap.gt.this%nmaps .or. imap.le.0) stop "coop_healpix_patch_get_radial_profile: wrong input arguments"
    halfm = m/2
    this%fr(:, halfm, imap) = 0.d0
    select case(imap)
    case(1)
       do i=-this%n, this%n
          do j=-this%n, this%n
             this%fr(this%icm(i, j, 0), halfm, imap) =  this%fr(this%icm(i, j, 0), halfm, imap) + this%image(i, j, imap) * this%wcm(i, j, m)
             this%fr(this%icm(i, j, 1), halfm, imap) =  this%fr(this%icm(i, j, 1), halfm, imap) + this%image(i, j, imap)* this%wcm(i, j, m+1) 
          enddo
       enddo
    case(2)
       do i=-this%n, this%n
          do j=-this%n, this%n
             this%fr(this%icm(i, j, 0), halfm, imap) =  this%fr(this%icm(i, j, 0), halfm, imap) + this%image(i, j, imap) * this%wsm(i, j, m)
             this%fr(this%icm(i, j, 1), halfm, imap) =  this%fr(this%icm(i, j, 1), halfm, imap) + this%image(i, j, imap)* this%wsm(i, j, m+1) 
          enddo
       enddo
    case default
       stop "Cannot get radial profile for more than 2 maps."
    end select
  end subroutine coop_healpix_patch_get_radial_profile


  subroutine coop_healpix_maps_simulate(this)
    class(coop_healpix_maps) this
    COOP_REAL,dimension(:),allocatable::sqrtCls
    COOP_REAL,dimension(:, :),allocatable::Cls_sqrteig
    COOP_REAL,dimension(:,:,:),allocatable::Cls_rot
    COOP_INT l
    if(this%nmaps.eq.1 .and. this%spin(1).eq.0)then
       allocate(sqrtCls(0:this%lmax))
       !$omp parallel do
       do l = 0, this%lmax
          sqrtCls(l) = sqrt(this%Cl(l,1))
       enddo
       !$omp end parallel do
       call this%simulate_Tmaps(this%nside, this%lmax, sqrtCls)
       deallocate(sqrtCls)
    elseif(this%nmaps.eq.3 .and. this%iq .eq.2)then
       allocate(Cls_sqrteig(3, 0:this%lmax), Cls_rot(3,3,0:this%lmax))
       call coop_healpix_Cls2Rot(this%lmax, this%Cl, Cls_sqrteig, Cls_rot)
       call this%simulate_TQUmaps(this%nside, this%lmax, Cls_sqrteig, Cls_rot)
       deallocate(Cls_sqrteig, Cls_rot)
    else
       stop "unknown coop_healpix_maps_simulate mode"
    endif
  end subroutine coop_healpix_maps_simulate

  subroutine coop_healpix_maps_simulate_Tmaps(this, nside, lmax, sqrtCls, lmin, onlyalm, onlyphase)
    class(coop_healpix_maps) this
    COOP_INT nside
    COOP_INT lmax
    COOP_REAL sqrtCls(0:lmax)
    COOP_INT l,m,  ell_min, ell_max
    COOP_INT,optional::lmin
    logical op
    logical,optional::onlyalm, onlyphase
    COOP_REAL::theta
    if(present(lmin))then
       ell_min = lmin
    else
       ell_min = 0
    endif
    call coop_healpix_nside2lmax(nside, ell_max)
    ell_max = min(lmax, ell_max)
    if(this%nside .ne. nside)then
       call this%init(nside = nside, nmaps = 1, genre = "TEMPERATURE", lmax=ell_max)
    elseif(this%lmax .ne. ell_max)then
       call this%allocate_alms(lmax = ell_max)
    endif

    if(present(onlyphase))then
       op = onlyphase
    else
       op = .false.
    endif
    if(op)then
       !$omp parallel do private(l, m, theta)
       do l=ell_min, ell_max
          this%alm(l, 0, 1) = SqrtCls(l)     
          do m = 1, l
             theta = coop_random_unit()*coop_2pi
             this%alm(l, m, 1) = cmplx(SqrtCls(l)*cos(theta), SqrtCls(l)*sin(theta))
          enddo
       enddo
       !$omp end parallel do       
    else
       !$omp parallel do private(l, m)
       do l=ell_min, ell_max
          this%alm(l, 0, 1) = coop_random_complex_Gaussian(.true.)*SqrtCls(l)     
          do m = 1, l
             this%alm(l, m, 1) = coop_random_complex_Gaussian()*SqrtCls(l)
          enddo
       enddo
       !$omp end parallel do
    endif
    if(present(onlyalm))then
       if(onlyalm)return
    endif
    call this%alm2map( index_list = (/ 1 /) )
  end subroutine coop_healpix_maps_simulate_Tmaps


  subroutine coop_healpix_maps_get_cls(this) !!I assume you have already called    this_map2alm(this)
    class(coop_healpix_maps)this
    COOP_INT l, m, i, j, k
    if(.not.allocated(this%alm)) stop "coop_healpix_maps_get_cls: you have to call coop_healpix_maps_map2alm before calling this subroutine"
    this%cl = 0.
    !$omp parallel do private(i,j,k,l)
    do i=1, this%nmaps
       do j=1, i
          k = coop_matsym_index(this%nmaps, i, j)
          do l = 0, this%lmax
             this%Cl(l, k) = (sum(COOP_MULT_REAL(this%alm(l, 1:l, i), this%alm(l, 1:l, j))) + 0.5d0 * COOP_MULT_REAL(this%alm(l,0,i), this%alm(l,0,j)) )/(l+0.5d0)
          enddo
       enddo
    enddo
    !$omp end parallel do
  end subroutine coop_healpix_maps_get_cls

  subroutine coop_healpix_Cls2Rot(lmax, Cls, Cls_sqrteig, Cls_rot)
    COOP_INT lmax
    COOP_SINGLE,dimension(0:lmax, 6),intent(IN)::Cls !!ordering is TT, EE, BB, TE, EB, TB
    COOP_REAL, dimension(3, 0:lmax),intent(OUT)::Cls_sqrteig
    COOP_REAL, dimension(3, 3, 0:lmax),intent(OUT)::Cls_rot
    COOP_INT l
    COOP_REAL a2(2,2), a3(3,3)
    COOP_REAL psi2(2,2), psi3(3,3)
    if(all(Cls(:,coop_healpix_index_EB).eq.0.d0) .and. all(Cls(:, coop_healpix_index_TB).eq.0.d0))then
       Cls_sqrteig(coop_healpix_index_BB,:) = sqrt(Cls(:, coop_healpix_index_BB))
       Cls_rot(coop_healpix_index_BB, coop_healpix_index_TT, :) = 0
       Cls_rot(coop_healpix_index_BB, coop_healpix_index_EE, :) = 0
       Cls_rot(coop_healpix_index_TT, coop_healpix_index_BB, :) = 0
       Cls_rot(coop_healpix_index_EE, coop_healpix_index_BB, :) = 0
       Cls_rot(coop_healpix_index_BB, coop_healpix_index_BB, :) = 1.
       do l=0, lmax
          a2(1,1) = Cls(l, coop_healpix_index_TT)
          a2(2,2) = Cls(l, coop_healpix_index_EE)
          a2(1,2) = Cls(l, coop_healpix_index_TE)
          a2(2,1) = a2(1,2)
          call coop_matsym_diag(2, 2, a2, psi2)
          Cls_sqrteig(coop_healpix_index_TT, l) = sqrt(a2(1,1))
          Cls_sqrteig(coop_healpix_index_EE, l) = sqrt(a2(2,2))
          Cls_rot( coop_healpix_index_TT:coop_healpix_index_EE, coop_healpix_index_TT:coop_healpix_index_EE, l) = psi2
       end do
    else
       do l=0, lmax
          a3(1,1) = Cls(l,coop_healpix_index_TT)
          a3(2,2) = Cls(l,coop_healpix_index_EE)
          a3(3,3) = Cls(l,coop_healpix_index_BB)
          a3(1,2) = Cls(l,coop_healpix_index_TE)
          a3(2,3) = Cls(l,coop_healpix_index_EB)
          a3(3,1) = Cls(l,coop_healpix_index_TB)
          a3(2,1) = a3(1,2)
          a3(3,2) = a3(2,3)
          a3(1,3) = a3(3,1)
          call coop_matsym_diag(3, 3, a3, psi3)
          Cls_sqrteig(coop_healpix_index_TT, l) = sqrt(a3(1,1))
          Cls_sqrteig(coop_healpix_index_EE, l) = sqrt(a3(2,2))
          Cls_sqrteig(coop_healpix_index_BB, l) = sqrt(a3(3,3))
          Cls_rot(:, :, l) = psi3
       end do
    endif
  end subroutine coop_healpix_Cls2Rot

  subroutine coop_healpix_maps_get_QU(this, idone)
    class(coop_healpix_maps) this
    logical, optional::idone
    COOP_INT::i
    if(this%spin(1).ne.0)then
       stop "get_QU: the first map must be spin 0"
    end if
    if(.not. present(idone))then
       call this%map2alm(index_list = (/ 1 /))
    else
       if(.not. idone)then
          call this%map2alm(index_list = (/ 1 /))
       endif
    endif
    if(this%nmaps.lt.3) call this%extend(3)
    this%alm(:,:,2) = this%alm(:,:,1)
    this%alm(:,:,3) = 0
    do i=2, 3
       call this%set_unit(i, this%units(1))
    enddo    
    select case(trim(coop_str_numUpperalpha(this%fields(1))))
    case("LNI")
       call this%set_field(2, "LNQT")
       call this%set_field(3, "LNUT")       
    case("T","I","TEMPERATURE", "INTENSITY", "ISTOKES")
       call this%set_field(2, "QT")
       call this%set_field(3, "UT")
    case("ZETA", "Z")
       call this%set_field(2, "QZ")
       call this%set_field(3, "UZ")
    case("E", "EPOLARISATION")
       call this%set_field(2, "Q")
       call this%set_field(3, "U")
    case default
       write(*,*) trim(this%fields(1))
       stop "get_QU: only supprt I, ZETA and E maps"
    end select
    call this%alm2map( index_list = (/ 2, 3 /) )
  end subroutine coop_healpix_maps_get_QU

  subroutine coop_healpix_maps_get_QUL(this, idone)
    class(coop_healpix_maps) this
    logical,optional::idone
    COOP_REAL:: resol 
    COOP_INT l, i
    if(this%spin(1).ne.0)then
       stop "get_QUL: the first map must be spin 0"
    end if
    if(this%nmaps.lt.4) call this%extend(4)    
    if(.not.present(idone))then
       call this%map2alm(index_list = (/ 1 /) )
    else
       if(.not. idone)then
          call this%map2alm(index_list = (/ 1 /) )
       endif
    endif
    resol = 2.d0/this%nside/this%nside
    this%alm(:,:, 2:4) = 0.
    do l = 2, min(this%lmax, this%nside*2)
       this%alm(l,0:l,4) = this%alm(l,0:l,1)*(l*(l+1.)*exp(-l*(l+1.)*resol))
       this%alm(l,0:l,2) = this%alm(l,0:l,4)
    enddo
    do i=2, 4
       call this%set_unit(i, this%units(1))
    enddo
    select case(trim(coop_str_numUpperalpha(this%fields(1))))
    case("LNI")
       call this%set_field(2, "LNQLT")
       call this%set_field(3, "LNULT")
       call this%set_field(4, "LNLT")              
    case("T","I","TEMPERATURE", "INTENSITY", "ISTOKES")
       call this%set_field(2, "QLT")
       call this%set_field(3, "ULT")
       call this%set_field(4, "LT")       
    case("ZETA", "Z")
       call this%set_field(2, "QLZ")
       call this%set_field(3, "ULZ")
       call this%set_field(4, "LZ")       
    case("E", "EPOLARISATION")
       call this%set_field(2, "QLE")
       call this%set_field(3, "ULE")
       call this%set_field(4, "LE")       
    case default
       write(*,*) trim(this%fields(1))
       stop "get_QUL: only supprt I, ZETA and E maps"
    end select

    call this%alm2map( index_list = (/ 2, 3, 4 /) )
  end subroutine coop_healpix_maps_get_QUL



  subroutine coop_healpix_maps_get_QULDD(this, idone)
    class(coop_healpix_maps) this
    logical,optional::idone    
    COOP_INT l, lmax, i
    COOP_REAL::resol
#ifdef HAS_HEALPIX
    if(this%spin(1).ne.0)then
       stop "get_QULDD: the first map must be spin 0"
    end if
    if(this%nmaps.lt.6) call this%extend(6)

    !!set alms    
    if(.not.present(idone))then
       call this%map2alm(index_list = (/ 1 /) )
    else
       if(.not. idone)then
          call this%map2alm(index_list = (/ 1 /) )
       endif
    endif
    lmax = min(this%lmax, this%nside*2)
    resol = 2./this%nside/this%nside    
    this%alm(:,:, 2:6) = 0.
    do l = 2, lmax
       this%alm(l,0:l,2) = this%alm(l, 0:l, 1)* (l*(l+1.)*exp(-l*(l+1.)*resol))
       this%alm(l,0:l,4) = this%alm(l, 0:l, 2)
       this%alm(l,0:l,5) = this%alm(l, 0:l, 1)* sqrt(l*(l+1.)*exp(-l*(l+1.)*resol))
    enddo
    
    !!set units and spins
    do i=2, 6
       call this%set_unit(i, this%units(1))
    enddo    
    select case(trim(coop_str_numUpperalpha(this%fields(1))))
    case("LNI")
       call this%set_field(2, "LNQLT")
       call this%set_field(3, "LNULT")
       call this%set_field(4, "LNLT")
       call this%set_field(5, "LNID1")
       call this%set_field(6, "LNID2")       
    case("T","I","TEMPERATURE", "INTENSITY", "ISTOKES")
       call this%set_field(2, "QLT")
       call this%set_field(3, "ULT")
       call this%set_field(4, "LT")
       call this%set_field(5, "ID1")
       call this%set_field(6, "ID2")       
    case("ZETA", "Z")
       call this%set_field(2, "QLZ")
       call this%set_field(3, "ULZ")
       call this%set_field(4, "LZ")
       call this%set_field(5, "ZD1")
       call this%set_field(6, "ZD2")       
    case("E", "EPOLARISATION")
       call this%set_field(2, "QLE")
       call this%set_field(3, "ULE")
       call this%set_field(4, "LE")
       call this%set_field(5, "ED1")
       call this%set_field(6, "ED2")              
    case default
       write(*,*) trim(this%fields(1))
       stop "get_QULDD: only supprt I, ZETA and E maps"
    end select
    call this%alm2map( index_list = (/ 2, 3, 4, 5, 6 /) )
#endif    
  end subroutine coop_healpix_maps_get_QULDD



  subroutine coop_healpix_maps_qu2EB(this)
    class(coop_healpix_maps) this
    if(this%iq.eq.0 .or. this%iu.eq.0 .or. this%iq+1 .ne. this%iu)stop "QU2EB: cannot find Q, U maps"
    call this%map2alm(index_list = (/ this%iq, this%iu /) )
    call this%set_field(this%iq, "E-POLARISATION")
    call this%set_field(this%iu, "B-POLARISATION")
    call this%alm2map( index_list = (/ this%iq, this%iu /) )
    this%iq = 0
    this%iu = 0
  end subroutine coop_healpix_maps_qu2EB

  subroutine coop_healpix_maps_EB2QU(this)
    class(coop_healpix_maps) this
    COOP_INT ie
    ie = 1
    do while(ie.lt.this%nmaps)
       if(trim(this%fields(ie)).eq. "E-POLARISATION" .or. trim(this%fields(ie)).eq."E")then
          if(trim(this%fields(ie+1)).eq."B-POLARISATION" .or. trim(this%fields(ie+1)).eq."B") exit
       endif
       ie  = ie + 1
    enddo
    if(ie.ge.this%nmaps)stop "EB2QU: cannot find E, B maps"
    call this%map2alm(index_list = (/ ie, ie+1 /) )
    call this%set_field(ie, "Q-POLARISATION")
    call this%set_field(ie+1, "U-POLARISATION")
    call this%alm2map(index_list = (/ ie, ie+1 /) )
  end subroutine coop_healpix_maps_EB2QU


  subroutine coop_healpix_maps_simulate_TQUmaps(this, nside, lmax, Cls_sqrteig, Cls_rot)
    class(coop_healpix_maps) this
    COOP_INT lmax, nside
    COOP_REAL,dimension(3, 0:lmax)::Cls_sqrteig
    COOP_REAL,dimension(3, 3, 0:lmax)::Cls_rot
    COOP_INT l, m, lm
    call coop_healpix_nside2lmax(nside, lm)
    lm = min(lmax, lm)
    if(this%nside .ne. nside)then
       call this%init(nside = nside, nmaps = 3, genre="IQU", lmax = lm)
    elseif(this%lmax .ne. lm)then
       call this%allocate_alms(lmax = lm)
    endif
    !$omp parallel do private(l, m)
    do l=0, lm
       this%alm(l, 0, :) = matmul(Cls_rot(:, :, l), Cls_sqrteig(:,l) * (/ coop_random_complex_Gaussian(.true.), coop_random_complex_Gaussian(.true.), coop_random_complex_Gaussian(.true.) /) )
       do m = 1, l
          this%alm(l, m, :) = matmul(Cls_rot(:, :, l), Cls_sqrteig(:,l) * (/ coop_random_complex_Gaussian(), coop_random_complex_Gaussian(), coop_random_complex_Gaussian() /) )
       enddo
    enddo
    !$omp end parallel do
    call this%alm2map()
  end subroutine coop_healpix_maps_simulate_TQUmaps


  subroutine coop_healpix_maps_init(this, nside, nmaps, genre, lmax, nested)
    class(coop_healpix_maps) this
    COOP_INT:: nside, nmaps, i
    COOP_INT, optional::lmax
    COOP_UNKNOWN_STRING,optional::genre
    logical,optional::nested
#ifdef HAS_HEALPIX
    if(allocated(this%map))then
       if(this%nside .eq. nside .and. this%nmaps.eq.nmaps)then
          goto 100
       endif
       deallocate(this%map)
    endif
    if(allocated(this%spin))deallocate(this%spin)
    if(allocated(this%fields))deallocate(this%fields)
    if(allocated(this%units))deallocate(this%units)        
    this%nside = nside
    this%nmaps = nmaps
    this%npix = nside2npix(nside)
    allocate(this%map(0:this%npix - 1, nmaps))
    allocate(this%spin(this%nmaps), this%fields(this%nmaps), this%units(this%nmaps))
100 if(present(lmax)) call this%allocate_alms(lmax)
    this%ordering = COOP_RING !!default ordering
    if(present(genre))then
       select case(trim(adjustl(coop_str_numUpperalpha(genre))))
       case("UNKNOWN")
          call this%set_units("1")
          do i=1, this%nmaps
             call this%set_field(i, coop_healpix_maps_default_genre)
          enddo          
       case("IQU", "TQU")
          if(this%nmaps.lt.3)stop "For IQU map you need at least 3 maps"
          call this%set_units("muK")
          call this%set_field(1, "INTENSITY")
          call this%set_field(2, "Q-POLARISATION")
          call this%set_field(3, "U-POLARISATION")
          do i=4, this%nmaps
             call this%set_field(i, coop_healpix_maps_default_genre)
          enddo
       case("MASK", "M")
          call this%set_units("1")          
          do i=1, this%nmaps
             call this%set_field(i, trim(genre))
          enddo
       case("I", "INTENSITY", "T", "TEMPERATURE", "E", "B", "EPOLARISATION", "BPOLARISATION", "LT", "LZ", "LE", "LB")
          call this%set_units("muK")          
          call this%set_fields(trim(genre))
       case("ZETA", "Z")
          call this%set_units("10^{-5}")
          call this%set_fields(trim(genre))          
       case("QU")
          if(this%nmaps.lt.2) stop "For QU map you need at least 2 maps"
          call this%set_units("muK")
          call this%set_field(1, "Q")
          call this%set_field(2, "U")
          do i=3, this%nmaps
             call this%set_field(i, coop_healpix_maps_default_genre)
          enddo
       case("TEB", "IEB")
          if(this%nmaps.lt.3)stop "For TEB map you need at least 3 maps"        
          call this%set_units("muK")          
          call this%set_field(1, "I")
          call this%set_field(2, "E")
          call this%set_field(3, "B")
          do i=4, this%nmaps
             call this%set_field(i, coop_healpix_maps_default_genre)
          enddo
       case("TE", "IE")
          if(this%nmaps.lt.2) stop "For TE map you need at least 2 maps"
          call this%set_units("muK")
          call this%set_field(1, "T")
          call this%set_field(2, "E")
          do i=3, this%nmaps
             call this%set_field(i, coop_healpix_maps_default_genre)
          enddo
       case("EB")
          if(this%nmaps.lt.2) stop "For EB map you need at least 2 maps"
          call this%set_units("muK")
          call this%set_field(1, "E")
          call this%set_field(2, "B")
          do i=3, this%nmaps
             call this%set_field(i, coop_healpix_maps_default_genre)
          enddo
       case("TQUL")
          if(this%nmaps.lt.4) stop "For TQUL map you need at least 4 maps"
          call this%set_units("muK")
          call this%set_field(1, "I")
          call this%set_field(2, "QLT")
          call this%set_field(3, "ULT")
          call this%set_field(4, "LT")
          do i=5, this%nmaps
             call this%set_field(i, coop_healpix_maps_default_genre)
          enddo
       case("TQULDD")
          if(this%nmaps.lt.6) stop "For TQUL map you need at least 4 maps"
          call this%set_units("muK")
          call this%set_field(1, "I")
          call this%set_field(2, "QLT")
          call this%set_field(3, "ULT")
          call this%set_field(4, "LT")
          call this%set_field(5, "ID1")
          call this%set_field(6, "ID2")          
          do i=7, this%nmaps
             call this%set_field(i, coop_healpix_maps_default_genre)
          enddo                    
       case("ZQUL")
          if(this%nmaps.lt.4) stop "For TQUL map you need at least 4 maps"
          call this%set_units("10^{-5}")
          call this%set_field(1, "ZETA")
          call this%set_field(2, "QLZ")
          call this%set_field(3, "ULZ")
          call this%set_field(4, "LZ")
          do i=5, this%nmaps
             call this%set_field(i,  coop_healpix_maps_default_genre)  
          enddo
       case("ZQULDD")
          if(this%nmaps.lt.6) stop "For TQUL map you need at least 4 maps"
          call this%set_units("10^{-5}")
          call this%set_field(1, "ZETA")
          call this%set_field(2, "QLZ")
          call this%set_field(3, "ULZ")
          call this%set_field(4, "LZ")
          call this%set_field(5, "ZD1")
          call this%set_field(6, "ZD2")          
          do i=7, this%nmaps
             call this%set_field(i,  coop_healpix_maps_default_genre)  
          enddo                              
       case("EQUL")
          if(this%nmaps.lt.4) stop "For TQUL map you need at least 4 maps"
          call this%set_units("muK")
          call this%set_field(1, "E")
          call this%set_field(2, "QLE")
          call this%set_field(3, "ULE")
          call this%set_field(4, "LE")
          do i=5, this%nmaps
             call this%set_field(i,  coop_healpix_maps_default_genre)  
          enddo
       case("EQULDD")
          if(this%nmaps.lt.6) stop "For TQUL map you need at least 4 maps"
          call this%set_units("muK")
          call this%set_field(1, "E")
          call this%set_field(2, "QLE")
          call this%set_field(3, "ULE")
          call this%set_field(4, "LE")
          call this%set_field(5, "ED1")
          call this%set_field(6, "ED2")          
          do i=7, this%nmaps
             call this%set_field(i,  coop_healpix_maps_default_genre)  
          enddo                              
       case("TQTUT")
          if(this%nmaps.lt.3) stop "For TQTUT map you need at least 3 maps"
          call this%set_units("muK")
          call this%set_field(1, "TEMPERATURE")
          call this%set_field(2, "QT")
          call this%set_field(3, "UT")
          do i=4, this%nmaps
             call this%set_field(i, coop_healpix_maps_default_genre)       
          enddo          
       case("DERVS")
          if(this%nmaps .lt. 6) stop "for dervs map you need at least 6 maps"
          call this%set_units("muK")
          call this%set_field(1, "I")
          call this%set_field(2, "QLT")
          call this%set_field(3, "ULT")
          call this%set_field(4, "LT")
          call this%set_field(5, "ID1")
          call this%set_field(6, "ID2")
       case("ZQZUZ")
          if(this%nmaps.lt.3) stop "For ZQZUZ map you need at least 3 maps"
          call this%set_units("10^{-5}")
          call this%set_field(1, "ZETA")
          call this%set_field(2, "QZ")
          call this%set_field(3, "UZ")
          do i=4, this%nmaps
             call this%set_field(i, coop_healpix_maps_default_genre)         
          enddo
       case default
          write(*,"(A)") trim(genre)//": unknown map types"
          stop
       end select
    end if
    if(present(nested))then
       if(nested)this%ordering = COOP_NESTED
    endif
#else
    stop "DID not find healpix"
#endif
  end subroutine coop_healpix_maps_init


  subroutine coop_healpix_maps_extend(this, nmaps)
    class(coop_healpix_maps) this
    COOP_INT nmaps
    COOP_SINGLE,dimension(:,:),allocatable::map
    COOP_INT, dimension(:), allocatable::spin
#ifdef HAS_HEALPIX
    if(nmaps .le. this%nmaps) return
    allocate(map(0:this%npix-1, this%nmaps))
    map = this%map
    allocate(spin(this%nmaps))
    spin  = this%spin
    deallocate(this%map, this%spin)
    allocate(this%map(0:this%npix-1, nmaps), this%spin(nmaps))
    this%map(:, 1:this%nmaps) = map
    this%spin(1:this%nmaps) = spin
    if(allocated(this%alm))then
       deallocate(this%alm)
       call this%allocate_alms(this%lmax)
    endif
    deallocate(map, spin)
#else
    stop "DID not find healpix"
#endif
  end subroutine coop_healpix_maps_extend


  subroutine coop_healpix_maps_allocate_alms(this, lmax)
    class(coop_healpix_maps) this
    COOP_INT lmax
    if(allocated(this%alm))then
       if(this%lmax  .eq. lmax )then
          return
       endif
       deallocate(this%alm, this%alm_done, this%checksum)
    endif
    if(allocated(this%cl))deallocate(this%cl)
    this%lmax = lmax
    allocate(this%alm(0:this%lmax, 0:this%lmax, this%nmaps))
    allocate(this%cl(0:max(this%lmax, coop_healpix_default_lmax), this%nmaps*(this%nmaps+1)/2))
    allocate(this%alm_done(this%nmaps), this%checksum(4, this%nmaps))
    this%alm_done = .false.
    this%checksum = 1.e30
    this%alm = 0.
    this%cl = 0.
  end subroutine coop_healpix_maps_allocate_alms

  subroutine coop_healpix_maps_free(this)
    class(coop_healpix_maps) this
    if(allocated(this%fields))deallocate(this%fields)
    if(allocated(this%units))deallocate(this%units)    
    if(allocated(this%map))deallocate(this%map)
    if(allocated(this%alm))deallocate(this%alm)
    if(allocated(this%cl))deallocate(this%cl)
    if(allocated(this%spin))deallocate(this%spin)
    if(allocated(this%alm_done))deallocate(this%alm_done)
    if(allocated(this%checksum))deallocate(this%checksum)
    call this%header%free()
    this%lmax = -1
    this%nmaps = 0
    this%nside = 0
    this%npix = 0
  end subroutine coop_healpix_maps_free

  subroutine coop_healpix_maps_read(this, filename, nmaps_wanted, nmaps_to_read, known_size, nested)
    class(coop_healpix_maps) this
    COOP_UNKNOWN_STRING filename
    COOP_STRING::lowercase_name
    COOP_INT,optional::nmaps_wanted, nmaps_to_read
    integer(8) npixtot
    COOP_INT nmaps_actual
    logical,optional::known_size, nested
    COOP_REAL::fac, maxabs, rt
    logical::missing_fields
    COOP_SHORT_STRING ordering, the_unit
    COOP_INT i
#ifdef HAS_HEALPIX
    if(.not. coop_file_exists(filename))then
       write(*,*) trim(filename)
       stop "cannot find the file"
    else
       this%filename = trim(adjustl(filename))
    endif
    if(present(known_size))then
       if(known_size)then
          nmaps_actual = this%nmaps
          goto 200
       endif
    else
       call this%free()
    endif
    call coop_fits_to_header(filename, this%header)
    call coop_dictionary_lookup(this%header, "NSIDE", this%nside, 0)
    call coop_dictionary_lookup(this%header, "TFIELDS", nmaps_actual, 0)
    call coop_dictionary_lookup(this%header, "ORDERING", ordering, "")
    call coop_dictionary_lookup(this%header, "BAD_DATA", this%bad_data, this%bad_data)
    call coop_dictionary_lookup(this%header, "COORDSYS", this%coordsys, this%coordsys)
    call coop_dictionary_lookup(this%header, "FIRSTPIX", this%firstpix, 0) 
   
    
    
    if(this%nside .eq. 0 .or. nmaps_actual .eq. 0 .or. trim(ordering).eq."")then !!try robust routine
       npixtot = getsize_fits(trim(filename), nmaps = nmaps_actual, nside = this%nside, ordering = this%ordering)
       if(this%nside .eq. 0 .or. nmaps_actual .eq. 0) stop "coop_healpix_maps_read failed"       
    else
       select case(trim(ordering))
       case("RING", "ring", "Ring")
          this%ordering = COOP_RING
       case("NESTED", "nested", "Nested")
          this%ordering = COOP_NESTED
       case default
          this%ordering = COOP_UNKNOWN_ORDERING
       end select
    endif
    this%npix =nside2npix(this%nside)

    call coop_dictionary_lookup(this%header, "LASTPIX", this%lastpix, this%npix-1)

    if(present(nmaps_wanted))then       
       this%nmaps = nmaps_wanted
       if(nmaps_wanted .lt. nmaps_actual)then          
          nmaps_actual = nmaps_wanted
       endif
    else
       this%nmaps = nmaps_actual
    endif
    if(this%nmaps.eq.0) stop "cannot read 0 maps"
    if(allocated(this%spin))then
       if(size(this%spin) .ne. this%nmaps)then
          deallocate(this%spin)
          allocate(this%spin(this%nmaps))
       endif
    else
       allocate(this%spin(this%nmaps))
    endif
    if(allocated(this%fields))then
       if(size(this%fields) .ne. this%nmaps)then
          deallocate(this%fields)
          allocate(this%fields(this%nmaps))
       endif
    else
       allocate(this%fields(this%nmaps))
    endif
    if(allocated(this%units))then
       if(size(this%units) .ne. this%nmaps)then
          deallocate(this%units)
          allocate(this%units(this%nmaps))
       endif
    else
       allocate(this%units(this%nmaps))
    endif
    
    if(allocated(this%map))then
       if(size(this%map, 1).ne. this%npix .or. size(this%map, 2).ne.this%nmaps)then
          deallocate(this%map)
          allocate(this%map(0:this%npix-1, this%nmaps))
       endif
    else
       allocate(this%map(0:this%npix-1, this%nmaps))
    endif
    
200 if(present(nmaps_to_read))then
       call input_map(trim(filename), this%map, this%npix, min(nmaps_actual, nmaps_to_read), fmissval = 0.)
       if(this%nmaps .gt. nmaps_to_read)then
          do i=nmaps_to_read+1, this%nmaps
             call this%set_unit(i, this%units(1))
          enddo
       endif
    else
       call input_map(trim(filename), this%map, this%npix, nmaps_actual, fmissval = 0.)
    endif
    if(present(nested))then
       if(nested)then
          call this%convert2nested()
       else
          call this%convert2ring()
       endif
    else
       call this%convert2ring()
    endif
    missing_fields = .false.
    lowercase_name = trim(filename)       
    call coop_str2lower(lowercase_name)
    
    do i=1, this%nmaps
       this%units(i) = trim(this%header%value("TUNIT"//COOP_STR_OF(i)))
       if(trim(this%units(i)).eq."")then
          maxabs = maxval(abs(this%map(:,i)))
          rt = sqrt(sum(this%map(:,i)**2)/this%npix)
          if( index(lowercase_name, "mask") .ne. 0 .or. index(lowercase_name, "_ln.fits") .ne. 0 .or. (maxabs .lt. 2. .and. rt .gt. 0.2 .and. this%nmaps .eq. 1) )then
             call this%set_unit(i, "1")
          elseif( maxabs .lt. 100. .and. rt .lt. 0.2)then
             call this%set_unit(i, "K")
          else
             call this%set_unit(i, "muK")             
          endif
          call this%header%update("TUNIT"//COOP_STR_OF(i), trim(this%units(i)))          
       endif
       fac = coop_unit2muK(this%units(i))
       if(abs(fac-1.d0).gt.1.d-6)then
          this%map(:, i) = this%map(:,i)*fac
          this%units(i) = "muK"
          call this%header%update("TUNIT"//COOP_STR_OF(i), "muK")
       endif
       call coop_dictionary_lookup(this%header, "TTYPE"//COOP_STR_OF(i), this%fields(i))
       missing_fields = missing_fields .or. (trim(this%fields(i)).eq."")
    enddo

    if(missing_fields)then
       select case(this%nmaps)
       case(1)
          if(index(lowercase_name, "mask").ne.0 )then
             call this%set_field(1,  "MASK")
          elseif(index(lowercase_name, "zeta").ne.0)then
             call this%set_field(1, "ZETA")
          elseif(index(lowercase_name, "_t").ne. 0 .or. index(lowercase_name, "_i").ne.0)then
             call this%set_field(1, "T")
          elseif(index(lowercase_name, "_e").ne. 0)then
             call this%set_field(1, "E")
          elseif(index(lowercase_name, "_b").ne. 0)then
             call this%set_field(1, "B")
          else
             write(*,*)"Warning: cannot determine the map types for file "//trim(filename)
             call this%set_field(1, "T")
          endif
       case(2)
          if(index(lowercase_name, "qu").ne.0 .or. index(lowercase_name, "pol").ne.0)then
             call this%set_field(1, "Q")
             call this%set_field(2, "U")
          elseif(index(lowercase_name, "eb").ne.0)then
             call this%set_field(1, "E")
             call this%set_field(2, "B")
          elseif(index(lowercase_name, "qtut").ne.0)then
             call this%set_field(1, "QT")
             call this%set_field(2, "UT")
          elseif(index(lowercase_name, "qzuz").ne.0)then
             call this%set_field(1, "QZ")
             call this%set_field(2, "UZ")
          elseif(index(lowercase_name, "qltult").ne.0)then
             call this%set_field(1, "QLT")
             call this%set_field(2, "ULT")
          elseif(index(lowercase_name, "qlzulz").ne.0)then
             call this%set_field(1, "QLZ")
             call this%set_field(2, "ULZ")
          else
             if(.not. present(nmaps_to_read))write(*,*)"WARNING: cannot determine the map types for file "//trim(filename)
             call this%set_field(1, "Q")
             call this%set_field(2, "U")             
          endif
       case(3)
          if(index(lowercase_name, "iqu").ne.0 .or. index(lowercase_name, "tqu").ne.0)then          
             call this%set_field(1, "I")
             call this%set_field(2, "Q")
             call this%set_field(3, "U")
          elseif(index(lowercase_name, "teb").ne.0 .or. index(lowercase_name, "ieb").ne.0)then
             call this%set_field(1, "T")
             call this%set_field(2, "E")
             call this%set_field(3, "B")
          elseif(index(lowercase_name, "tqtut").ne.0)then
             call this%set_field(1, "T")
             call this%set_field(2, "QT")
             call this%set_field(3, "UT")
          elseif(index(lowercase_name, "zqzuz").ne.0)then             
             call this%set_field(1, "ZETA")
             call this%set_field(2, "QZ")
             call this%set_field(3, "UZ")
          else
             if(.not. present(nmaps_to_read))write(*,*)"Warning: cannot determine the map types for file "//trim(filename)
             call this%set_field(1, "I")
             call this%set_field(2, "Q")
             call this%set_field(3, "U")                                       
          endif
       case(4)
          if(index(lowercase_name, "tqul").ne.0)then                       
             call this%set_field(1, "T")
             call this%set_field(2, "QLT")
             call this%set_field(3, "ULT")
             call this%set_field(4, "LT")
          elseif(index(lowercase_name, "zqul").ne.0)then                       
             call this%set_field(1, "ZETA")
             call this%set_field(2, "QLZ")
             call this%set_field(3, "ULZ")
             call this%set_field(4, "LZ")
          elseif(index(lowercase_name, "equl").ne.0)then                       
             call this%set_field(1, "E")
             call this%set_field(2, "QLE")
             call this%set_field(3, "ULE")
             call this%set_field(4, "LE")
          else
             if(.not. present(nmaps_to_read))write(*,*)"Warning: cannot determine the map types for file "//trim(filename)
             call this%set_field(1, "T")
             call this%set_field(2, "QLT")
             call this%set_field(3, "ULT")
             call this%set_field(4, "LT")
          endif
       case(6)
          if(index(lowercase_name, "tquldd") .ne. 0)then
             call this%set_field(1, "T")
             call this%set_field(2, "QLT")
             call this%set_field(3, "ULT")
             call this%set_field(4, "LT")
             call this%set_field(5, "ID1")
             call this%set_field(6, "ID2")
          elseif(index(lowercase_name, "zquldd") .ne. 0)then
             call this%set_field(1, "Z")
             call this%set_field(2, "QLZ")
             call this%set_field(3, "ULZ")
             call this%set_field(4, "LZ")
             call this%set_field(5, "ZD1")
             call this%set_field(6, "ZD2")
          elseif(index(lowercase_name, "equldd") .ne. 0)then
             call this%set_field(1, "E")
             call this%set_field(2, "QLE")
             call this%set_field(3, "ULE")
             call this%set_field(4, "LE")
             call this%set_field(5, "ED1")
             call this%set_field(6, "ED2")
          else
             if(.not. present(nmaps_to_read))write(*,*)"Warning: cannot determine the map types for file "//trim(filename)
             call this%set_field(1, "T")
             call this%set_field(2, "QLT")
             call this%set_field(3, "ULT")
             call this%set_field(4, "LT")
             call this%set_field(5, "ID1")
             call this%set_field(6, "ID2")
          endif
       case default
          if(.not. present(nmaps_to_read))then
             write(*,*) "nmaps = ", this%nmaps
             stop "read: cannot determine the fields for this nmaps"
          endif
       end select
    else
       call this%fields_to_spins()
    endif

#else
    stop "DID NOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_read


  subroutine coop_healpix_maps_fields_to_spins(this, imap)
    class(coop_healpix_maps)::this
    COOP_INT::i
    COOP_INT,optional::imap
    if(present(imap))then
       i = imap
       select case(trim(coop_str_numUpperalpha(this%fields(i))))
       case("INTENSITY", "TEMPERATURE", "MASK", "EPOLARISATION", "BPOLARISATION", "E", "B", "ZETA", "Z", "I", "T", "M", "LT", "LZ", "LE", "ISTOKES", "", "LNI", "GENERAL")
          this%spin(i) = 0
       case("ID1", "ID2", "ZD1", "ZD2", "ED1", "ED2", "BD1", "BD2", "LNID1", "LNID2")
          this%spin(i) = 1
       case("QPOLARISATION", "Q", "QT", "QLT", "QLZ", "QZ", "QLE", "QSTOKES", "LNQ", "LNQT", "LNQLT")
          this%spin(i) = 2
          this%iq = i
       case("UPOLARISATION", "U", "UT", "ULT", "ULZ", "UZ", "ULE", "USTOKES", "LNU", "LNUT", "LNULT")
          this%spin(i) = 2
          this%iu = i
       case default
          write(*,*) "WARNING: map type "//trim(this%fields(i))//" is unknown"
          this%spin(i) = 0
       end select
    else       
       do i = 1, this%nmaps
          select case(trim(coop_str_numUpperalpha(this%fields(i))))
          case("INTENSITY", "TEMPERATURE", "MASK", "EPOLARISATION", "BPOLARISATION", "E", "B", "ZETA", "Z", "I", "T", "M", "LT", "LZ", "LE", "ISTOKES", "", "LNI")
             this%spin(i) = 0
          case("ID1", "ID2", "ZD1", "ZD2", "ED1", "ED2", "BD1", "BD2", "LNID1", "LNID2")
             this%spin(i) = 1
          case("QPOLARISATION", "Q", "QT", "QLT", "QLZ", "QZ", "QLE", "QSTOKES", "LNQ", "LNQT", "LNQLT")
             this%spin(i) = 2
             this%iq = i
          case("UPOLARISATION", "U", "UT", "ULT", "ULZ", "UZ", "ULE", "USTOKES", "LNU", "LNUT", "LNULT")
             this%spin(i) = 2
             this%iu = i
          case default
             write(*,*) "WARNING: map type "//trim(this%fields(i))//" is unknown"
             this%spin(i) = 0
          end select
       enddo
    endif
  end subroutine coop_healpix_maps_fields_to_spins

  subroutine coop_healpix_maps_import(this, filename, index_start, index_end)
    class(coop_healpix_maps) this
    COOP_UNKNOWN_STRING filename
    type(coop_healpix_maps)::tmp
    COOP_INT index_start, index_end, spin(index_end-index_start + 1)
#ifdef HAS_HEALPIX
    if(index_end .gt. this%nmaps) call this%extend(index_end)
    call tmp%read(filename = trim(filename), nmaps_wanted = index_end-index_start + 1, nested = (this%ordering .eq. COOP_NESTED) )
    if(tmp%nside .ne. this%nside) stop "import: nside must be the same"
    this%map(:, index_start:index_end) = tmp%map
    this%spin(index_start:index_end) = tmp%spin
    this%fields(index_start:index_end) = tmp%fields
    this%units(index_start:index_end) = tmp%units
    call tmp%free()
#else
    stop "DID NOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_import


  subroutine coop_healpix_maps_write(this, filename, index_list)
    class(coop_healpix_maps)this
    COOP_UNKNOWN_STRING filename
    COOP_INT,dimension(:),optional::index_list
    character(LEN=80),dimension(:),allocatable::header    
    logical pol
    COOP_INT i, ind
#ifdef HAS_HEALPIX    
    if(present(index_list))then
       if(any(index_list .lt. 1 .or. index_list .gt. this%nmaps)) stop "coop_healpix_write_map: index out of range"
       pol = any(this%spin(index_list).eq.2)
    else
       pol =any(this%spin.eq.2)
    endif
    call coop_delete_file(trim(filename))
    allocate(header(this%header%n+64))
    do i=1, size(header)
       header(i)  = ''
    enddo
    if(allocated(this%alm))then
       call write_minimal_header(header,dtype = trim(this%dtype), nside=this%nside, order = this%ordering, creator=trim(this%creator), version =trim(this%version), nlmax = this%lmax, nmmax = this%lmax, polar= pol, coordsys = trim(this%coordsys), fwhm_degree = this%fwhm_degree)
    else
       call write_minimal_header(header,dtype = trim(this%dtype), nside=this%nside, order = this%ordering, creator='COOP', version = '0.0', polar=pol, coordsys = trim(this%coordsys), fwhm_degree = this%fwhm_degree)
    endif
    do i=1, this%header%n
       select case(trim(coop_str_numUpperAlpha(this%header%key(i))))
       case("SIMPLE", "BITPIX", "NAXIS", "EXTEND", "XTENSION", "NAXIS1", "NAXIS2", "PCOUNT", "GCOUNT", "TFIELDS", "DATE")  !!these will be added by Healpix automatically
          cycle
       case("TTYPE1", "TTYPE2", "TTYPE3", "TTYPE4", "TTYPE5", "TTYPE6", "TUNIT1", "TUNIT2", "TUNIT3", "TUNIT4", "TUNIT5", "TUNIT6")
             call add_card(header, trim(this%header%key(i)), trim(this%header%val(i)), update = .true.)
       case default          
       end select
    enddo
    if(present(index_list))then
       call output_map(this%map(:, index_list), header, trim(filename))
    else
       call output_map(this%map, header, trim(filename))
    endif
    deallocate(header)
#else
       stop "DID NOT FIND HEALPIX"
#endif
   
  end subroutine coop_healpix_maps_write

  subroutine coop_healpix_convert_to_nested(this)
    class(coop_healpix_maps) this
#ifdef HAS_HEALPIX
    if(.not. allocated(this%map)) stop "coop_healpix_convert_to_nested: map is not allocated yet"
    if(this%ordering .eq. COOP_RING) then
       call convert_ring2nest(this%nside, this%map)
       this%ordering = COOP_NESTED
    elseif(this%ordering .ne. COOP_NESTED)then
       write(*,*) "ordering = ", this%ordering
       stop "coop_healpix_convert_to_nested: unknown ordering"
    endif
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_convert_to_nested

  subroutine coop_healpix_convert_to_ring(this)
    class(coop_healpix_maps) this
#ifdef HAS_HEALPIX
    if(.not. allocated(this%map)) stop "coop_healpix_convert_to_ring: map is not allocated yet"
    if(this%ordering .eq. COOP_NESTED)then
       call convert_nest2ring(this%nside, this%map)
       this%ordering = COOP_RING
    elseif(this%ordering .ne. COOP_RING)then
       write(*,*) "ordering = ", this%ordering
       stop "coop_healpix_convert_to_ring: UNKNOWN ordering"
    endif
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_convert_to_ring

  subroutine coop_healpix_maps_map2alm(this, lmax, index_list)
    class(coop_healpix_maps) this
    COOP_INT,optional::lmax
    COOP_INT i, l, j, lm, n
    COOP_INT,dimension(:),optional::index_list
    COOP_SINGLE_COMPLEX, dimension(:,:,:),allocatable::alm
#ifdef HAS_HEALPIX
    call this%convert2ring()
    if(present(lmax))then
       lm = min(lmax, floor(this%nside*coop_healpix_lmax_by_nside))
    else
       lm =  min(coop_healpix_default_lmax, floor(this%nside*coop_healpix_lmax_by_nside))
    endif
    call this%allocate_alms(lm)
    if(present(index_list))then
       j = 1
       n = size(index_list)
       do while(j.le.n)
          i = index_list(j)
          if(i.gt.this%nmaps)then
             j  = j + 1
             cycle
          endif
          if(coop_healpix_alm_check_done)then
             if(this%alm_done(i))then
                if( this%checksum(1, i).eq. this%map(0, i) .and. this%checksum(2, i).eq. this%map(this%npix/2, i) .and. this%checksum(3, i) .eq. this%map(this%npix-1, i) .and. nint(this%checksum(4, i)) .eq. this%spin(i)  )then
                   j  = j + 1
                   cycle
                endif
             endif
             this%alm_done(i) = .true.
             this%checksum(1, i) = this%map(0, i)
             this%checksum(2, i) = this%map(this%npix/2, i)
             this%checksum(3, i) = this%map(this%npix-1, i)
             this%checksum(4, i) = this%spin(i)
          endif
          if(this%spin(i).eq.0)then
             call map2alm(this%nside, this%lmax, this%lmax, this%map(:,i), this%alm(:,:,i:i))
             j = j + 1
          else
             if(i.lt. this%nmaps .and. j.lt. n)then
                if(index_list(j+1).eq.i+1 .and. this%spin(i+1) .eq. this%spin(i))then
                   if(.not. allocated(alm))allocate(alm(2, 0:this%lmax, 0:this%lmax))
                   call map2alm_spin(this%nside, this%lmax, this%lmax, this%spin(i), this%map(:,i:i+1), alm)
                   this%alm(:,:,i) = alm(1, :, :)
                   this%alm(:,:,i+1) = alm(2, :, :)
                   j = j + 2
                   cycle
                endif
             endif
             write(*,*) this%nmaps
             write(*,*) j
             write(*,*) index_list
             write(*,*)  this%spin
             stop "coop_healpix_maps_map2alm: nonzero spin maps must appear in pairs"
          endif
       enddo
    else
       i = 1
       do while(i.le. this%nmaps)
          if(this%spin(i).eq.0)then
             call map2alm(this%nside, this%lmax, this%lmax, this%map(:,i), this%alm(:,:,i:i))
             i = i + 1
          else
             if(.not. allocated(alm))allocate(alm(2, 0:this%lmax, 0:this%lmax))
             if(i.lt. this%nmaps)then
                if(this%spin(i+1) .eq. this%spin(i))then
                   call map2alm_spin(this%nside, this%lmax, this%lmax, this%spin(i), this%map(:,i:i+1), alm)
                   this%alm(:,:,i) = alm(1, :, :)
                   this%alm(:,:,i+1) = alm(2, :, :)
                   i = i + 2
                   cycle
                endif
             endif
             write(*,*) this%spin
             stop "coop_healpix_maps_map2alm: nonzero spin maps must appear in pairs"
          endif
       enddo
    endif
    if(allocated(alm))deallocate(alm)
    if(coop_healpix_want_cls)call this%get_cls
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_map2alm


  subroutine coop_healpix_maps_alm2map(this, index_list)
    class(coop_healpix_maps) this
    COOP_INT i, j, n
    COOP_INT,dimension(:),optional::index_list
    COOP_SINGLE_COMPLEX,dimension(:,:,:),allocatable::alm
#ifdef HAS_HEALPIX
    if(present(index_list))then
       j = 1
       n = size(index_list)
       do while(j.le.n)
          i = index_list(j)
          if(i.gt. this%nmaps)then
             j = j + 1
             cycle
          endif
          if(this%spin(i).eq.0)then
             call alm2map(this%nside, this%lmax, this%lmax, this%alm(:,:,i:i), this%map(:,i))
             j = j + 1
             cycle
          endif
          if(j .lt. n .and. i .lt. this%nmaps)then
             if(index_list(j+1).eq.i+1 .and. this%spin(i+1) .eq. this%spin(i))then
                if(.not.allocated(alm))allocate(alm(2,0:this%lmax, 0:this%lmax))
                alm(1,:,:) = this%alm(:,:,i)
                alm(2,:,:) = this%alm(:,:,i+1)
                call alm2map_spin(this%nside, this%lmax, this%lmax, this%spin(i), alm, this%map(:,i:i+1))
                j = j + 2
                cycle
             endif
          endif
          stop "coop_healpix_maps_alm2map: nonzero spin maps must appear in pairs"
       enddo
    else
       i = 1
       do while(i.le. this%nmaps)
          if(this%spin(i).eq.0)then
             call alm2map(this%nside, this%lmax, this%lmax, this%alm(:,:,i:i), this%map(:,i))
             i = i + 1
          else
             if(i.lt. this%nmaps)then
                if(this%spin(i+1) .eq. this%spin(i))then
                   if(.not.allocated(alm))allocate(alm(2,0:this%lmax, 0:this%lmax))
                   alm(1,:,:) = this%alm(:,:,i)
                   alm(2,:,:) = this%alm(:,:,i+1)
                
                   call alm2map_spin(this%nside, this%lmax, this%lmax, this%spin(i), alm, this%map(:,i:i+1))
                   i = i + 2
                   cycle
                endif
             endif
             stop "coop_healpix_maps_alm2map: nonzero spin maps must appear in pairs"
          endif
       enddo
    endif
    if(allocated(alm))deallocate(alm)
    this%ordering = COOP_RING
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_alm2map

  subroutine coop_healpix_filter_alm(this, fwhm, lpower, window, index_list)
    class(coop_healpix_maps) this
    COOP_REAL,optional::window(0:this%lmax)
    COOP_REAL,optional::fwhm
    COOP_REAL,optional::lpower
    COOP_INT,dimension(:), optional::index_list
    COOP_INT l
    COOP_SINGLE c, w(0:this%lmax)
    w = 1.
    if(present(fwhm))then
       c = sign((coop_sigma_by_fwhm * fwhm)**2/2., fwhm)
       !$omp parallel do
       do l = 0,  this%lmax
          w(l) = w(l)*exp(-l*(l+1.)*c)
       enddo
       !$omp end parallel do
    endif
    if(present(lpower))then
       !$omp parallel do
       do l = 0,  this%lmax
          w(l) = w(l)*(l*(l+1.))**(lpower/2.)
       enddo
       !$omp end parallel do
    endif
    if(present(window))then
       !$omp parallel do
       do l = 0,  this%lmax
          w(l) = w(l)*window(l)
       enddo
       !$omp end parallel do       
    endif
    if(present(index_list))then
       !$omp parallel do
       do l = 0, this%lmax
          this%alm(l,:,index_list) = this%alm(l,:,index_list)*w(l)
       enddo
       !$omp end parallel do
    else
       !$omp parallel do
       do l = 0, this%lmax
          this%alm(l,:,:) = this%alm(l,:,:)*w(l)
       enddo
       !$omp end parallel do
    endif
  end subroutine coop_healpix_filter_alm



  subroutine coop_healpix_maps_get_disc(this, pix, disc)
    class(coop_healpix_maps)::this
    COOP_INT pix
    type(coop_healpix_disc) disc
    COOP_REAL r
#ifdef HAS_HEALPIX
    disc%nside  = this%nside
    disc%center = pix
    disc%ordering = this%ordering
    call this%pix2ang(pix, disc%theta, disc%phi)
    call ang2vec(disc%theta, disc%phi, disc%nz)
    disc%nx = (/  sin(disc%phi) , - cos(disc%phi) , 0.d0 /)
    call coop_vector_cross_product(disc%nz, disc%nx, disc%ny)
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_get_disc

  subroutine coop_healpix_disc_pix2ang(disc, pix, r, phi)
    class(coop_healpix_disc) disc
    COOP_INT pix
    COOP_REAL r, phi, vec(3), x, y
    if(pix .eq. disc%center)then
       r = 0
       phi = 0
       return
    endif
#ifdef HAS_HEALPIX
    if(disc%ordering .eq. COOP_RING)then
       call pix2vec_ring(disc%nside, pix, vec)
    else
       call pix2vec_nest(disc%nside, pix, vec)       
    endif
#else
    stop "CANNOT FIND HEALPIX"
#endif
    r = COS2RADIUS( dot_product(vec, disc%nz) )
    x = dot_product(vec, disc%nx)
    y = dot_product(vec, disc%ny)
    phi = COOP_POLAR_ANGLE(x, y)
  end subroutine coop_healpix_disc_pix2ang

  subroutine coop_healpix_disc_ang2pix(disc, r, phi, pix)
    class(coop_healpix_disc) disc
    COOP_REAL r !!in unit of radian
    COOP_REAL phi, vec(3), cost, sint
    COOP_INT pix
    if(r.le.1.d-8)then
       pix = disc%center
       return
    endif
    cost = RADIUS2COS(abs(r))
    sint = sign(sqrt(1.d0 - cost**2), r)
    vec = sint*cos(phi)* disc%nx + sint*sin(phi)*disc%ny + cost*disc%nz
#ifdef HAS_HEALPIX
    if(disc%ordering .eq. COOP_RING)then
       call vec2pix_ring(disc%nside, vec, pix)
    else
       call vec2pix_nest(disc%nside, vec, pix)
    endif
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_disc_ang2pix


  subroutine coop_healpix_disc_pix2xy(disc, pix, x, y)
    class(coop_healpix_disc) disc
    COOP_INT pix
    COOP_REAL r, phi, vec(3), x, y
    if(pix .eq. disc%center)then
       x = 0
       y = 0
       return
    endif
#ifdef HAS_HEALPIX
    if(disc%ordering .eq. COOP_RING)then
       call pix2vec_ring(disc%nside, pix, vec)
    else
       call pix2vec_nest(disc%nside, pix, vec)
    endif
#else
    stop "CANNOT FIND HEALPIX"
#endif
    r = COS2RADIUS( dot_product(vec, disc%nz) )
    x = dot_product(vec, disc%nx)
    y = dot_product(vec, disc%ny)
    r = r/sqrt(x**2+y**2)
    x = r * x
    y = r * y
  end subroutine coop_healpix_disc_pix2xy


  subroutine coop_healpix_disc_xy2pix(disc, x, y, pix)
    class(coop_healpix_disc) disc
    COOP_REAL x, y !!in unit of radian
    COOP_REAL vec(3), cost, sint, r
    COOP_INT pix
    r = sqrt(x**2+y**2)
    if(r.lt.1.d-8)then
       pix = disc%center
       return
    endif
    cost = RADIUS2COS(r)
    sint = sqrt(1.d0 - cost**2)
    vec = sint*(x/r)* disc%nx + sint*(y/r)*disc%ny + cost*disc%nz
#ifdef HAS_HEALPIX
    if(disc%ordering .eq. COOP_RING)then
       call vec2pix_ring(disc%nside, vec, pix)
    else
       call vec2pix_nest(disc%nside, vec, pix)       
    endif
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_disc_xy2pix

  subroutine coop_healpix_rotate_qu(qu, phi, spin)
    COOP_SINGLE qu(2)
    COOP_INT spin
    COOP_REAL phi, cosp, sinp
    cosp = cos(spin*phi)
    sinp = sin(spin*phi)
    qu = (/ qu(1)*cosp + qu(2)*sinp,  -qu(1)*sinp + qu(2)*cosp /)
  end subroutine coop_healpix_rotate_qu


  subroutine coop_healpix_patch_get_fr0(patch, nvar, var)
    COOP_INT nvar
    type(coop_healpix_patch)::patch
    COOP_REAL var(nvar)
    call patch%get_radial_profile(1, 0)
    var = patch%fr(0:patch%n, 0, 1)
  end subroutine coop_healpix_patch_get_fr0


  subroutine coop_healpix_maps_stack_on_patch(this, disc, angle, patch, tmp_patch, mask)
    class(coop_healpix_maps) this
    type(coop_healpix_disc) disc
    type(coop_healpix_maps),optional::mask
    COOP_REAL angle
    type(coop_healpix_patch) patch, tmp_patch
    if(present(mask))then
       call this%fetch_patch(disc, angle, tmp_patch, mask)
       if(sum(tmp_patch%nstack*tmp_patch%indisk) .lt. patch%num_indisk_tol)return
    else
       call this%fetch_patch(disc, angle, tmp_patch)
    endif
    patch%image = patch%image + tmp_patch%image
    patch%nstack = patch%nstack + tmp_patch%nstack
    patch%nstack_raw = patch%nstack_raw + tmp_patch%nstack_raw
  end subroutine coop_healpix_maps_stack_on_patch

  subroutine coop_healpix_maps_fetch_patch(this, disc, angle, patch, mask)
    class(coop_healpix_maps)::this
    type(coop_healpix_disc) disc
    type(coop_healpix_maps),optional::mask
    COOP_REAL angle
    type(coop_healpix_patch) patch
    COOP_INT i, j, pix,  k
    COOP_REAL x, y, r, phi
    COOP_SINGLE qu(2)
    if(.not. present(mask))patch%nstack = 1.d0
    patch%nstack_raw  = 1
    if(all(patch%tbs%spin .eq. 0))then
       do j = -patch%n, patch%n
          do i = -patch%n, patch%n
             x = patch%dr * i
             y = patch%dr * j
             r = sqrt(x**2+y**2)
             if(r.ge. 2.d0)cycle
             phi = COOP_POLAR_ANGLE(x, y) + angle
             call disc%ang2pix( r, phi, pix)
             if(present(mask))then
                patch%nstack(i, j) = mask%map(pix, 1)
                patch%image(i, j, :) = this%map(pix, patch%tbs%ind)*mask%map(pix,1)
             else 
                patch%image(i, j, :) = this%map(pix,patch%tbs%ind)
             endif
          enddo
       enddo
    elseif(patch%nmaps .eq. 2 .and. all(patch%tbs%spin .ne. 0))then
       if(all(patch%tbs%local_rotation))then
          do j = -patch%n, patch%n
             do i = -patch%n, patch%n
                x = patch%dr * i
                y = patch%dr * j
                r = sqrt(x**2+y**2)
                if(r.ge. 2.d0)cycle                                
                phi = COOP_POLAR_ANGLE(x, y) + angle
                call disc%ang2pix( r, phi, pix)
                qu = this%map(pix, patch%tbs%ind)
                call coop_healpix_rotate_qu(qu, phi, patch%tbs%spin(1))
                qu = qu * coop_healpix_QrUrSign
                if(present(mask))then
                   patch%nstack(i, j) = mask%map(pix, 1)
                   patch%image(i, j, :) = qu*mask%map(pix,1)
                else 
                   patch%image(i, j, :) = qu
                endif
             enddo
          enddo
       else
          do j = -patch%n, patch%n
             do i = -patch%n, patch%n
                x = patch%dr * i
                y = patch%dr * j
                r = sqrt(x**2+y**2)
                if(r.ge. 2.d0)cycle                                
                phi = COOP_POLAR_ANGLE(x, y) + angle
                call disc%ang2pix( r, phi, pix)
                qu = this%map(pix, patch%tbs%ind)
                call coop_healpix_rotate_qu(qu, angle, patch%tbs%spin(1))
                if(present(mask))then
                   patch%nstack(i, j) = mask%map(pix, 1)
                   patch%image(i, j, :) = qu*mask%map(pix,1)
                else 
                   patch%image(i, j, :) = qu
                endif
             enddo
          enddo
       endif
    else
       do j = -patch%n, patch%n
          do i = -patch%n, patch%n
             x = patch%dr * i
             y = patch%dr * j
             r = sqrt(x**2+y**2)
             if(r.ge. 2.d0)cycle                             
             phi = COOP_POLAR_ANGLE(x, y) + angle
             call disc%ang2pix( r, phi, pix)
             k = 1
             do while(k.le. patch%nmaps)
                select case(patch%tbs%spin(k))
                case(0)
                   if(present(mask))then
                      patch%nstack(i, j) = mask%map(pix, 1)
                      patch%image(i, j, k) = this%map(pix, patch%tbs%ind(k))*mask%map(pix,1)
                   else 
                      patch%image(i, j, k) = this%map(pix,patch%tbs%ind(k))
                   endif
                case default
                   qu = this%map(pix, patch%tbs%ind(k:k+1))
                   if(patch%tbs%local_rotation(k))then
                      call coop_healpix_rotate_qu(qu, phi, patch%tbs%spin(k))
                      qu = qu * coop_healpix_QrUrSign
                   else
                      call coop_healpix_rotate_qu(qu, angle,patch%tbs%spin(k))
                   endif
                   if(present(mask))then
                      patch%nstack(i, j) = mask%map(pix, 1)
                      patch%image(i, j, k:k+1) = qu*mask%map(pix,1)
                   else 
                      patch%image(i, j, k:k+1) = qu
                   endif
                end select
                if(patch%tbs%spin(k).ne.0)then
                   k = k + 2
                else
                   k = k + 1
                endif
             enddo
          enddo
       enddo
    endif
    if(coop_healpix_patch_stack_fft)then
       do i = 1, patch%nmaps
          call patch%fft(i)
       enddo
    endif
  end subroutine coop_healpix_maps_fetch_patch

  subroutine coop_healpix_smooth_mapfile(mapfile, fwhm)
    COOP_UNKNOWN_STRING mapfile
    type(coop_healpix_maps) map
    COOP_REAL fwhm
    call map%read(mapfile)
    call map%smooth(fwhm)
    call map%write(trim(coop_file_add_postfix(trim(mapfile),"_"//trim(coop_num2str(nint(fwhm/coop_SI_arcmin)))//"a")))
    write(*,*) "output: "//trim(coop_file_add_postfix(trim(mapfile),"_"//trim(coop_num2str(nint(fwhm/coop_SI_arcmin)))//"a"))
    call map%free()
  end subroutine coop_healpix_smooth_mapfile



  subroutine coop_healpix_maps_smooth(map, fwhm, index_list, l_lower, l_upper)
    class(coop_healpix_maps) map
    COOP_REAL fwhm
    COOP_INT, optional::l_lower, l_upper
    COOP_INT,dimension(:),optional::index_list
    COOP_INT lmax
    if(fwhm .gt. 0.d0)then
       lmax = min(ceiling(3./max(abs(fwhm)*coop_sigma_by_fwhm, 1.d-6)), floor(map%nside*coop_healpix_lmax_by_nside), coop_healpix_default_lmax)
    else
       lmax = min(floor(map%nside*coop_healpix_lmax_by_nside), coop_healpix_default_lmax)
    endif
    if(present(l_upper))then
       lmax = min(lmax, l_upper)
    endif
    if(lmax .lt. 2) stop "Huge smoothing scale cannot be done!"    
    if(lmax*abs(fwhm).lt.0.01)return
    if(present(index_list))then
       if(any(index_list .gt. map%nmaps)) stop "smooth: index_list overflow"
       call map%map2alm(lmax, index_list)
       if(present(l_lower)) map%alm(0:l_lower-1, :, :) = 0.
       call map%filter_alm(fwhm = fwhm, index_list = index_list)
       call map%alm2map(index_list)
    else
       call map%map2alm(lmax)
       if(present(l_lower)) map%alm(0:l_lower-1, :, :) = 0.       
       call map%filter_alm(fwhm =fwhm)
       call map%alm2map()
    endif
  end subroutine coop_healpix_maps_smooth

  subroutine coop_healpix_mask_diffuse(from, to, fwhm)
    class(coop_healpix_maps)::from, to
    type(coop_healpix_maps)::hm
    COOP_REAL::fwhm, thetasq
    COOP_INT::l, lmax
    call hm%init(nside = min(from%nside*2, 2048), nmaps = 1, genre = "MASK", nested = .true.)
    call coop_healpix_maps_ave_udgrade(from, hm)

    thetasq = (fwhm*coop_sigma_by_fwhm)**2/2.d0
    call coop_healpix_nside2lmax(hm%nside, lmax)
    lmax = min(lmax, ceiling(2.5d0/sqrt(thetasq)))
    call hm%map2alm(lmax = lmax )
    do l = 0, hm%lmax
       hm%alm(l, 0:l, 1) = hm%alm(l, 0:l, 1)*exp(-l*(l+1.d0)*thetasq)
    enddo
    call hm%alm2map()
    if(to%nside .ne. from%nside) call to%init(nside = from%nside, nmaps=1, genre = "MASK")
    call coop_healpix_maps_ave_udgrade(hm, to)
    where (from%map(:, 1) .gt. 0.)
       to%map(:, 1) = from%map(:, 1)
    elsewhere
       to%map(:, 1) = min(max(to%map(:, 1), 0.), 1.)
    end where
  end subroutine coop_healpix_mask_diffuse


  subroutine coop_healpix_maps_diffuse(from, to, mask, fwhm)
    class(coop_healpix_maps)::from, to, mask
    type(coop_healpix_maps)::hm
    COOP_REAL::fwhm, thetasq
    COOP_INT::l, lmax
    call mask%convert2nested()
    call hm%init(nside = min(from%nside*2, 2048), nmaps = from%nmaps, genre = "I", nested = .true.)
    call coop_healpix_maps_ave_udgrade(from = from, to = hm, mask = mask)
    thetasq = (fwhm*coop_sigma_by_fwhm)**2/2.d0
    call coop_healpix_nside2lmax(hm%nside, lmax)
    lmax = min(lmax, ceiling(2.5d0/sqrt(thetasq)))
    call hm%map2alm(lmax = lmax )
    do l = 0, hm%lmax
       hm%alm(l, 0:l, 1) = hm%alm(l, 0:l, 1)*exp(-l*(l+1.d0)*thetasq)
    enddo
    call hm%alm2map()
    if(to%nside .ne. from%nside) call to%init(nside = from%nside, nmaps=1, genre = "I")
    call coop_healpix_maps_ave_udgrade(hm, to)
    where (mask%map(:, 1) .gt. 0.)
       to%map(:, 1) = from%map(:, 1)
    end where
  end subroutine coop_healpix_maps_diffuse
  
  

  subroutine coop_healpix_maps_smooth_with_window(map, fwhm, lmax, window, index_list)
    class(coop_healpix_maps) map
    COOP_REAL fwhm
    COOP_INT lmax
    COOP_INT,dimension(:),optional::index_list
    COOP_REAL window(0:lmax)
    if(present(index_list))then
       if(any(index_list .gt. map%nmaps)) stop "smooth: index_list overflow"
       call map%map2alm(lmax, index_list)
       call map%filter_alm(fwhm = fwhm, window = window, index_list = index_list)
       call map%alm2map(index_list)
    else
       call map%map2alm(lmax)
       call map%filter_alm(fwhm = fwhm, window = window)
       call map%alm2map()
    endif
  end subroutine coop_healpix_maps_smooth_with_window


  subroutine coop_healpix_getQU(Emap_file, QUmap_file)
    COOP_UNKNOWN_STRING Emap_file, QUmap_file
    type(coop_healpix_maps) hge, hgqu
    call hge%read(Emap_file,  nmaps_wanted = 1)
    call hge%map2alm()
    call hgqu%init(nside = hge%nside, nmaps = 2, genre="QU", lmax=hge%lmax)
    hgqu%alm(:, :, 1) = hge%alm(:, :, 1)
    hgqu%alm(:, :, 2) = 0
    call hgqu%alm2map()
    call hgqu%write(QUmap_file)
    call hgqu%free()
    call hge%free()
  end subroutine coop_healpix_getQU


  subroutine coop_healpix_output_map(map, header, fname)
    real, dimension(:,:):: map
    character(LEN=80),dimension(:):: header
    COOP_UNKNOWN_STRING fname
#ifdef HAS_HEALPIX
    call coop_delete_file(trim(fname))
    call output_map(map, header, fname)
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_output_map
  
  subroutine coop_healpix_mask_map(mapfile, maskfile, output, index_list)
    type(coop_healpix_maps) map, mask
    COOP_INT i
    COOP_UNKNOWN_STRING mapfile, maskfile, output
    COOP_INT,dimension(:),optional::index_list
    call map%read(mapfile)
    call mask%read(maskfile,  nmaps_wanted = 1)
    if(map%ordering .eq. COOP_RING)then
       call mask%convert2ring()
    elseif(map%ordering .eq. COOP_NESTED)then
       call mask%convert2nested()
    else
       write(*,*) "coop_healpix_mask_map: mask file has unknown ordering."
       stop
    endif
    call coop_delete_file(output)
    if(present(index_list))then
       if(any(index_list .lt. 1 .or. index_list .gt. map%nmaps)) stop "coop_healpix_write_map: index out of range"
       do i=1, size(index_list)
          map%map(:,index_list(i)) = map%map(:, index_list(i))*mask%map(:,1)
       enddo
       call map%write(trim(coop_file_add_postfix(output,"_masked")), index_list)
       call map%write(output)
    else
       do i=1, map%nmaps
          map%map(:,i) = map%map(:, i)*mask%map(:,1)
       enddo
       call map%write(output)
    endif
    call map%free
    call mask%free

  end subroutine coop_healpix_mask_map

  subroutine coop_healpix_maps_query_disc(this, pix, radius, listpix, nlist)
    class(coop_healpix_maps)::this
    COOP_INT pix, listpix(0:), nlist
    COOP_REAL radius, vec(3)
#ifdef HAS_HEALPIX
    call this%pix2vec(pix, vec)
    if(this%ordering .eq. COOP_RING)then
       call query_disc(this%nside, vec, radius, listpix, nlist, nest = 0, inclusive = 0)
    else
       call query_disc(this%nside, vec, radius, listpix, nlist, nest = 1, inclusive = 0)
    endif
#else
      stop "CANNOT FIND HEALPIX"    
#endif    
    end subroutine coop_healpix_maps_query_disc

  subroutine coop_healpix_maps_query_strip(this, theta1, theta2, listpix, nlist)
    class(coop_healpix_maps)::this
    COOP_INT listpix(0:), nlist
    COOP_REAL theta1, theta2
#ifdef HAS_HEALPIX
    if(this%ordering .eq. COOP_RING)then
       call query_strip(this%nside, theta1, theta2, listpix, nlist, nest = 0, inclusive = 0)
    else
       call query_strip(this%nside, theta1, theta2, listpix, nlist, nest = 1, inclusive = 0)
    endif
#else
      stop "CANNOT FIND HEALPIX"    
#endif    
    end subroutine coop_healpix_maps_query_strip
    

  subroutine coop_healpix_mask_hemisphere(mask, theta, phi)
    COOP_REAL vec(3), theta, phi
    COOP_INT nlist
    type(coop_healpix_maps)::mask
    COOP_INT,dimension(:),allocatable::listpix
#ifdef HAS_HEALPIX
    call mask%convert2ring()
    allocate(listpix(0:mask%npix-1))
    call ang2vec(theta, phi, vec)
    call query_disc(mask%nside, vec, coop_pio2, listpix, nlist, nest = 0, inclusive = 0)
    mask%map(listpix(0:nlist-1),:) = 0.
    deallocate(listpix)
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_mask_hemisphere


  subroutine coop_healpix_split(filename)
    COOP_UNKNOWN_STRING:: filename
    type(coop_healpix_maps) this
    COOP_INT i
    call this%read(filename)
    do i=1, this%nmaps
       call this%write(trim(coop_file_add_postfix(filename, "_submap"//trim(coop_ndigits(i, 3)))), (/ i /) )
    enddo
    call this%free()
  end subroutine coop_healpix_split



  subroutine coop_healpix_flip_mask(mask, flip)
    type(coop_healpix_maps)mask, flip
    COOP_INT pix, conjpix
    COOP_REAL theta, phi
#ifdef HAS_HEALPIX
    if(allocated(flip%map))then
       if(flip%nside .ne. mask%nside)then
          call flip%free()
          call flip%init(nside = mask%nside, nmaps = 1, genre="MASK")
       endif
    else
       call flip%free()
       call flip%init(nside = mask%nside, nmaps = 1, genre = "MASK")
    endif
    call mask%convert2ring()
    !$omp parallel do private(pix, conjpix, theta, phi)
    do pix=0, mask%npix-1
       call mask%pix2ang(pix, theta, phi)
       theta = coop_pi-theta
       phi = coop_pi +  phi
       call mask%ang2pix(theta, phi, conjpix)
       flip%map(pix, 1) = mask%map(conjpix, 1)
    enddo
    !$omp end parallel do
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_flip_mask


  subroutine coop_healpix_maps_regularize_in_mask(this, mask, imap)
    class(coop_healpix_maps)::this
    type(coop_healpix_maps)::mask
    COOP_INT::i, imap
    COOP_SINGLE::mmax, mmin, diff
    if(this%nside .ne. mask%nside .or. this%ordering .ne. mask%ordering) stop "regularize_in_mask: maps do not match the mask"
    mmax = maxval(this%map(:,imap)*mask%map(:,1))
    mmin = minval(this%map(:,imap)*mask%map(:,1))
    diff = (mmax-mmin)/10.d0
    if(diff .eq. 0.)then
       this%map(:,imap) = mmax
       return
    endif
    do i=0, this%npix-1
       if(mask%map(i,1).gt.0.5)cycle
       if(this%map(i, imap).gt. mmax)then
          this%map(i, imap) = mmax + log(1.+log(1.+(this%map(i, imap)-mmax)/diff))*diff
       elseif(this%map(i, imap) .lt. mmin)then
          this%map(i, imap) = mmin - log(1.+log(1.+(mmin-this%map(i, imap))/diff))*diff
       endif
    enddo
  end subroutine coop_healpix_maps_regularize_in_mask

 
  subroutine coop_healpix_maps_t2zeta(this, fwhm_arcmin, want_unconstrained)
    class(coop_healpix_maps)::this
    COOP_INT::l, i
    logical,optional::want_unconstrained 
    type(coop_cosmology_firstorder)::fod
    COOP_REAL,dimension(:,:),allocatable::Cls, Cls_lensed
    COOP_REAL::fwhm_arcmin, norm, ucnorm, CTT
    logical douc
#if DO_ZETA_TRANS
    if(present(want_unconstrained))then
       douc = want_unconstrained
       if(douc .and. this%nmaps .lt. 3) &
            stop "maps_t2zeta for unconstrained map you need nmaps>=3"
    else
       douc = .false.
    endif
    call this%map2alm(index_list = (/ 1 /) )
    
    call fod%Set_Planck_bestfit()
    call fod%compute_source(0)
    allocate(Cls(coop_num_cls, 2:this%lmax), Cls_lensed(coop_num_cls, 2:this%lmax))
    call fod%source(0)%get_All_Cls(2, this%lmax, Cls)
    call coop_get_lensing_Cls(2, this%lmax, Cls, Cls_lensed)
    
    norm = 1.d0/coop_healpix_zeta_normalization/ COOP_DEFAULT_TCMB 

    !$omp parallel do private(i, ucnorm, CTT)
    do l = 2, this%lmax
       CTT = ( Cls(coop_index_ClTT, l) + abs(Cls_lensed(coop_index_ClTT,l)) +  coop_Planck_TNoise(l)/coop_healpix_beam(l, fwhm_arcmin) )
       
       this%alm(l, :, 1) = this%alm(l,:,1) *  (Cls(coop_index_ClTzeta,l)/CTT * norm)!!the lensing Cl is added as a "noise"
       if(douc)then

          ucnorm = Cls(coop_index_Clzetazeta, l) -  Cls(coop_index_ClTzeta,l)**2/CTT
          if(ucnorm .lt. 0.d0) stop "clzetazeta - cltzeta**2/cltt negative?"
          ucnorm = sqrt(ucnorm)*norm
          this%alm(l, 0, 2) = coop_random_complex_Gaussian(.true.)*ucnorm
          do i=1, l
             this%alm(l, i, 2) = coop_random_complex_Gaussian()*ucnorm
          enddo
       endif
    enddo
    !$omp end parallel do
    call this%set_units(" 10^{-5}")
    call this%set_field(1, "ZETA")
    
    
    if(douc)then
       this%alm(:,:,3) = this%alm(:,:,2) + this%alm(:,:,1)
       call this%set_field(2, "ZETA")
       call this%set_field(3, "ZETA")       
       this%alm(0:1,:,3) = 0.
       call this%alm2map( index_list = (/ 1, 2, 3 /) )
       this%alm_done(1:3) = .false.
    else
       this%alm(0:1,:,1) = 0.
       call this%alm2map( index_list = (/ 1 /) )
       this%alm_done(1) = .false.
    endif
    deallocate(Cls, Cls_lensed)

#else
    stop "To use t2zeta you have to turn on DO_ZETA_TRANS in include/constants.h"
#endif
  end subroutine coop_healpix_maps_t2zeta


   subroutine coop_healpix_maps_E2zeta(this, fwhm_arcmin, want_unconstrained)
    class(coop_healpix_maps)::this
    COOP_INT::l, i
    logical,optional::want_unconstrained 
    type(coop_cosmology_firstorder)::fod
    COOP_REAL,dimension(:,:),allocatable::Cls, Cls_lensed
    COOP_REAL::fwhm_arcmin, norm, ucnorm, CEE
    logical douc
#if DO_ZETA_TRANS
    if(present(want_unconstrained))then
       douc = want_unconstrained
       if(douc .and. this%nmaps .lt. 3) &
            stop "maps_t2zeta for unconstrained map you need nmaps>=3"
    else
       douc = .false.
    endif
    if(this%spin(1).eq.0)then
       call this%map2alm(index_list = (/ 1 /) )
    else
       call this%map2alm(index_list = (/ 1, 2 /) )
    endif
    call fod%Set_Planck_bestfit()
    call fod%compute_source(0)
    allocate(Cls(coop_num_cls, 2:this%lmax), Cls_lensed(coop_num_cls, 2:this%lmax))
    call fod%source(0)%get_All_Cls(2, this%lmax, Cls)
    call coop_get_lensing_Cls(2, this%lmax, Cls, Cls_lensed)
    
    norm = 1.d0/coop_healpix_zeta_normalization/ COOP_DEFAULT_TCMB 
    !$omp parallel do private(i, ucnorm, CEE, l)
    do l = 2, this%lmax
       CEE = ( Cls(coop_index_ClEE, l) + abs(Cls_lensed(coop_index_ClEE,l)) +  coop_Planck_ENoise(l)/coop_healpix_beam(l, fwhm_arcmin) )
       
       this%alm(l, :, 1) = this%alm(l,:,1) *  (Cls(coop_index_ClEzeta,l)/CEE * norm)!!the lensing Cl is added as a "noise"
       if(douc)then

          ucnorm = Cls(coop_index_Clzetazeta, l) -  Cls(coop_index_ClEzeta,l)**2/CEE
          if(ucnorm .lt. 0.d0) stop "clzetazeta - cltzeta**2/cltt negative?"
          ucnorm = sqrt(ucnorm)*norm
          this%alm(l, 0, 2) = coop_random_complex_Gaussian(.true.)*ucnorm
          do i=1, l
             this%alm(l, i, 2) = coop_random_complex_Gaussian()*ucnorm
          enddo
       endif
    enddo
    !$omp end parallel do

    call this%set_units(" 10^{-5}")
    call this%set_field(1, "ZETA")
    
    if(douc)then
       this%alm(:,:,3) = this%alm(:,:,2) + this%alm(:,:,1)
       this%spin(1:3) = 0
       this%alm(0:1,:,3) = 0.
       call this%alm2map( index_list = (/ 1, 2, 3 /) )
       this%alm_done(1:3) = .false.
    else
       call this%set_field(2, "ZETA")
       call this%set_field(3, "ZETA")       
    
       this%spin(1) = 0
       this%alm(0:1,:,1) = 0.
       call this%alm2map( index_list = (/ 1 /) )
       this%alm_done(1) = .false.
    endif
    deallocate(Cls, Cls_lensed)

#else
    stop "To use t2zeta you have to turn on DO_ZETA_TRANS in include/constants.h"
#endif
  end subroutine coop_healpix_maps_E2zeta



  subroutine coop_healpix_maps_te2zeta(this, fwhm_arcmin, want_unconstrained)
    class(coop_healpix_maps)::this
    COOP_INT::l, i
    logical, optional::want_unconstrained
    type(coop_cosmology_firstorder)::fod
    COOP_REAL,dimension(:,:),allocatable::Cls, Cls_lensed
    COOP_REAL CTT, CTE, CEE, delta, bl
    COOP_REAL::fwhm_arcmin, norm, ucnorm, coef_T, coef_E
    logical douc
#if DO_ZETA_TRANS
    COOP_REAL::knorm
    if(this%nmaps .lt. 3) stop "TE2Zeta not enough maps"
    if(present(want_unconstrained))then
       douc = want_unconstrained
    else
       douc = .false.
    endif
    if(this%spin(1).eq.0 .and. this%spin(2).eq.0)then  !!TEB maps
       call this%map2alm(index_list = (/ 1, 2 /) )
    else  !!IQU maps
       if(this%nmaps .lt. 3) stop "TE2Zeta requires IQU maps or TE maps"
       if(this%spin(1).ne. 0 .or. this%spin(2).ne.2 .or. this%spin(3) .ne. 2) stop "TE2Zeta spin does not match"
       call this%map2alm(index_list = (/ 1, 2, 3 /) )
    endif
    call fod%Set_Planck_bestfit()
    call fod%compute_source(0)
    allocate(Cls(coop_num_cls, 2:this%lmax), Cls_lensed(coop_num_cls, 2:this%lmax))
    call fod%source(0)%get_All_Cls(2, this%lmax, Cls)
    call coop_get_lensing_Cls(2, this%lmax, Cls, Cls_lensed)
    
    norm = 1.d0/coop_healpix_zeta_normalization/COOP_DEFAULT_TCMB

    !$omp parallel do private(CTT, CTE, CEE, delta, bl, ucnorm, coef_T, coef_E)
    do l = 2, this%lmax
       bl = coop_healpix_beam(l, fwhm_arcmin)
       CTT = Cls(coop_index_ClTT, l) + abs(Cls_lensed(coop_index_ClTT,l))+ coop_Planck_TNoise(l)/bl
       CEE = Cls(coop_index_ClEE, l) + abs(Cls_lensed(coop_index_ClEE, l)) + Coop_Planck_Enoise(l)/bl
       CTE = Cls(coop_index_ClTE, l)
       delta =  CTT*CEE-CTE**2
       coef_T = (Cls(coop_index_ClTzeta,l)*CEE - Cls(coop_index_ClEzeta,l) * CTE)/delta
       coef_E = (-Cls(coop_index_ClTzeta,l)*CTE + Cls(coop_index_ClEzeta,l) * CTT)/delta
       this%alm(l, :, 1) = (coef_T * norm) * this%alm(l,:,1) + (coef_E * norm) * this%alm(l,:,2)
       if(douc)then
          ucnorm = Cls(coop_index_Clzetazeta, l) - (CTT*coef_T**2 +CEE*coef_E**2 +  2.d0*coef_E*coef_T*CTE)
          if(ucnorm .lt. 0.d0) stop "clzetazeta - cltzeta**2/cltt negative?"
          ucnorm = sqrt(ucnorm)*norm
          this%alm(l, 0, 2) = coop_random_complex_Gaussian(.true.)*ucnorm
          do i=1, l
             this%alm(l, i, 2) = coop_random_complex_Gaussian()*ucnorm
          enddo
       endif
    enddo

    call this%set_units(" 10^{-5}")
    call this%set_field(1, "ZETA")
    
    if(douc)then
       this%alm(:,:,3) = this%alm(:,:,2) + this%alm(:,:,1)
       this%alm(0:1,:,1:3) = 0
       call this%set_field(2, "ZETA")
       call this%set_field(3, "ZETA")              
       call this%alm2map( index_list = (/ 1, 2, 3 /) )
       this%alm_done(1:3) = .false.
    else
       this%alm(0:1,:,1) = 0
       call this%alm2map( index_list = (/ 1 /) )
       this%alm_done(1) = .false.
    endif
    deallocate(Cls, Cls_lensed)
#else
    stop "To use t2zeta you have to turn on DO_ZETA_TRANS in include/constants.h"
#endif
  end subroutine coop_healpix_maps_te2zeta


  function coop_healpix_beam(l, fwhm_arcmin) result(bl)
    COOP_INT l
    COOP_REAL fwhm_arcmin, bl
    bl = exp(-min(l*(l+1)*(fwhm_arcmin*coop_SI_arcmin*coop_sigma_by_fwhm)**2, 10.d0))
  end function coop_healpix_beam


  function coop_ACT_TNoise(l) result(Nl)
    COOP_INT l
    COOP_REAL Nl
    COOP_REAL,parameter::sigma = 12.d0/2.726e6*coop_SI_arcmin
    Nl = sigma**2
  end function coop_ACT_TNoise


  function coop_ACT_ENoise(l) result(Nl)
    COOP_INT l
    COOP_REAL Nl
    COOP_REAL,parameter::sigma = 12.d0/2.726e6*coop_SI_arcmin*coop_sqrt2
    Nl = sigma**2
  end function coop_ACT_ENoise


  function coop_ACT_BNoise(l) result(Nl)
    COOP_INT l
    COOP_REAL Nl
    COOP_REAL,parameter::sigma = 10.d0/2.726e6*coop_SI_arcmin*coop_sqrt2
    Nl = sigma**2
  end function coop_ACT_BNoise
  

  function coop_Planck_TNoise(l) result(Nl)
    COOP_INT l
    COOP_REAL Nl
    COOP_INT, parameter::fit_n = 10
    COOP_REAL, parameter,dimension(fit_n)::coef =  (/   -34.541512859596317      , &
  -2.3327044684438025      , &
  -5.2434038408357253E-2 , &
  0.10730432003605284      , &
  0.50599237614744652      , &
   4.4555688282625905E-2 , &
 -0.11595894402202822      , &
  -7.9077770071418474E-3 , &
 -0.26911077968221031      , &
  0.16444457464651557  /)
    COOP_REAL, parameter::rmin = coop_ln2, rmax = log(3000.d0)
    call coop_chebeval(fit_n, rmin, rmax, coef, log(dble(l)), Nl)
    Nl = exp(Nl)/COOP_DEFAULT_TCMB **2
  end function Coop_Planck_TNoise

  function coop_Planck_ENoise(l) result(Nl)
    COOP_INT l
    COOP_REAL Nl
    COOP_INT, parameter::fit_n = 10
    COOP_REAL, parameter,dimension(fit_n)::coef = (/  &
         -34.960562790571522      , &
         -0.87917841773973660      , &
         0.80533716121595145      , &
         0.19774868605364659      , &
         6.2456251840707022E-002 , &
         -6.4067689277299111E-002 , &
         6.9826114409022866E-003 , &
         -5.2937498702857466E-002 , &
         -1.3724094895074757E-002 , &
         -3.1217087209044592E-002 /)
    COOP_REAL, parameter::rmin = coop_ln2, rmax = log(3000.d0)
    call coop_chebeval(fit_n, rmin, rmax, coef, log(dble(l)), Nl)
    Nl = exp(Nl)/2.726**2
  end function Coop_Planck_ENoise
  
  function coop_Planck_BNoise(l) result(Nl)
    COOP_INT l
    COOP_REAL Nl
    COOP_INT, parameter::fit_n = 10
    COOP_REAL, parameter,dimension(fit_n)::coef =  (/ &
           -34.981682219713093      , &
           -0.88852205874105650      , &
           0.87526333776073173      , &
           0.13522507399381523      , &
           0.12803135081705186      , &
           -8.4600554185378374E-2 , &
           -1.3976431844054060E-2 , &
           -2.2053627296518163E-2 , &
           -3.6967260171160810E-2 , &
           -1.4665139293283336E-2 /)
    COOP_REAL, parameter::rmin = coop_ln2, rmax = log(3000.d0)
    call coop_chebeval(fit_n, rmin, rmax, coef, log(dble(l)), Nl)
    Nl = exp(Nl)/2.726**2
  end function Coop_Planck_BNoise

  function Coop_highpass_filter(l1, l2, l) result(w)
    COOP_INT l1, l2, l
    COOP_REAL w
    if(l.le. l1)then
       w = 0.d0
       return
    endif
    if(l.ge.l2)then
       w = 1.d0
       return
    endif
    w = sin(dble(l-l1)/dble(l2-l1)*coop_pio2)
  end function Coop_highpass_filter

  function coop_gaussian_filter(fwhm_arcmin, l) result(w)
    COOP_INT l
    COOP_REAL fwhm_arcmin, w
    w = exp(-((coop_sigma_by_fwhm*coop_SI_arcmin/coop_sqrt2)*fwhm_arcmin)**2*l*(l+1.d0))
  end function coop_gaussian_filter

  function coop_lowpass_filter(l1, l2, l) result(w)
    COOP_INT l1, l2, l
    COOP_REAL w
    if(l.ge. l2)then
       w = 0.d0
       return
    endif
    if(l.le.l1)then
       w = 1.d0
       return
    endif
    w = sin(dble(l2-l)/dble(l2-l1)*coop_pio2)
  end function coop_lowpass_filter


  subroutine coop_healpix_maps_get_peaks(this, sto, mask, restore, nside_scan)
    class(coop_healpix_maps)::this  
    type(coop_stacking_options)::sto
    type(coop_healpix_maps),optional::mask
    logical::domask
    COOP_INT i, nneigh, list(8), index_peak, ip, ipr
    COOP_INT::total_weight
    COOP_REAL::thetaphi(2)
    logical, optional::restore
    logical,dimension(:),allocatable::livept
    COOP_INT::nside_rand, npix_rand
    COOP_INT,optional::nside_scan
    COOP_INT, parameter:: num_rand = 500000 !!maximum number of points wanted for random selection of points (sto%genre = coop_stacking_genre_random_hot etc.)
    type(coop_healpix_maps)::zeros1, zeros2
#ifdef HAS_HEALPIX
    if(sto%nmaps .ne. this%nmaps)stop "get_peaks: nmaps mismatch"
    call sto%free()
    sto%nside = this%nside
    call this%convert2nested()
    if(present(mask))then
       if(mask%nside .ne. this%nside)then
          write(*,*) "map nside = ", this%nside
          write(*,*) "mask nside = ", mask%nside
          stop "coop_healpix_get_peaks: mask and map must have the same nside"
       endif
       call mask%convert2nested()
       total_weight = count(mask%map(:,1) .gt. 0.5)
       domask = .true.
    else
       domask = .false.
       total_weight = this%npix
    endif
    if(total_weight .le. 1) stop "get_peaks: no unmasked pixels"
    
    if(sto%index_I .ne. 0 .and. (abs(sto%I_lower_nu).lt.coop_stacking_max_threshold .or. abs(sto%I_upper_nu).lt. coop_stacking_max_threshold))then  !!rescale I
       if(domask)then
          sto%sigma_I  = sqrt(sum(dble(this%map(:,sto%index_I))**2, mask%map(:,1).gt.0.5)/total_weight)
       else
          sto%sigma_I = sqrt(sum(dble(this%map(:,sto%index_I))**2)/total_weight)
       endif
       sto%I_lower = sto%I_lower_nu* sto%sigma_I
       sto%I_upper = sto%I_upper_nu* sto%sigma_I
    endif
    if(sto%index_L .ne. 0 .and. (abs(sto%L_lower_nu).lt.coop_stacking_max_threshold .or. abs(sto%L_upper_nu).lt. coop_stacking_max_threshold) )then  !!rescale L
       if(domask)then
          sto%sigma_L  = sqrt(sum(dble(this%map(:,sto%index_L))**2, mask%map(:,1).gt.0.5)/total_weight)
       else
          sto%sigma_L = sqrt(sum(dble(this%map(:, sto%index_L))**2)/total_weight)
       endif
       sto%L_lower = sto%L_lower_nu* sto%sigma_L
       sto%L_upper = sto%L_upper_nu* sto%sigma_L
    endif
    if(sto%index_Q .ne. 0 .and. sto%index_U .ne. 0 .and. (abs(sto%P_lower_nu).lt.coop_stacking_max_threshold .or. abs(sto%P_upper_nu).lt. coop_stacking_max_threshold) )then  !!rescale P2
       if(domask)then
          sto%sigma_P  = sqrt(sum(dble(this%map(:,sto%index_Q)**2+this%map(:,sto%index_U)**2), mask%map(:,1).gt.0.5)/total_weight)
       else
          sto%sigma_P = sqrt(sum(dble(this%map(:,sto%index_Q)**2+this%map(:,sto%index_U)**2))/total_weight)
       endif
       sto%P2_lower = (sto%P_lower_nu * sto%sigma_P)**2
       sto%P2_upper = (sto%P_upper_nu * sto%sigma_P)**2
    endif
    select case(sto%genre)
    case(coop_stacking_genre_saddle, coop_stacking_genre_saddle_Oriented)
       index_peak = sto%index_I
       if(index_peak .ne. 1 .or. sto%index_L .ne. 4 .or. sto%index_Q .ne. 2 .or. sto%index_U .ne. 3 .or. this%nmaps .lt. 6)then
          stop "get_peaks: wrong configuration for saddle points stacking"
       endif
    case(coop_stacking_genre_Imax, coop_stacking_genre_Imin, coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Imin_Oriented)
       index_peak = sto%index_I
       if(index_peak .le. 0 .or. index_peak .gt. this%nmaps) stop "map index of peak overflow"
    case(coop_stacking_genre_Lmax, coop_stacking_genre_Lmin, coop_stacking_genre_Lmax_Oriented, coop_stacking_genre_Lmin_Oriented)
       index_peak = sto%index_L
       if(index_peak .le. 0 .or. index_peak .gt. this%nmaps) stop "map index of peak overflow"
    case(coop_stacking_genre_random_hot, coop_stacking_genre_random_cold, coop_stacking_genre_random_hot_oriented, coop_stacking_genre_random_cold_oriented)
       if(sto%index_Q .ne. 0 .and. sto%index_U .ne. 0 .and. (abs(sto%P_lower_nu).lt.coop_stacking_max_threshold .or. abs(sto%P_upper_nu).lt. coop_stacking_max_threshold) )then
          index_peak = sto%index_Q
       else
          index_peak = sto%index_I
       endif
    case default
       index_peak = sto%index_Q
    end select
    if(sto%nested)then
       if(domask)then    
          select case(sto%genre)
          case(coop_stacking_genre_saddle, coop_stacking_genre_saddle_Oriented)
             call this%zeros(5, zeros1, mask)
             call this%zeros(6, zeros2, mask)
             zeros1%map = zeros1%map*zeros2%map
             do i=0, this%npix-1
                if(zeros1%map(i, 1) .gt. 0.5 .and. .not. sto%reject(this%map(i, :)) .and. this%map(i, 2)**2+this%map(i, 3)**2 .gt. this%map(i,4)**2 )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
             call zeros1%free()
             call zeros2%free()
          case(coop_stacking_genre_Imax, coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Lmax, coop_stacking_genre_Lmax_Oriented)
             do i=0, this%npix-1
                if( mask%map(i, 1) .le. 0.5  .or. sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), index_peak).lt.this%map(i, index_peak)) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Imin, coop_stacking_genre_Imin_Oriented, coop_stacking_genre_Lmin, coop_stacking_genre_Lmin_Oriented)
             do i=0, this%npix-1
                if( mask%map(i, 1) .le. 0.5  .or. sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), index_peak) .gt. this%map(i, index_peak)) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Pmax_Oriented)
             do i=0, this%npix-1
                if( mask%map(i, 1) .le. 0.5  .or. sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), sto%index_Q)**2 + this%map(list(1:nneigh), sto%index_U)**2 .lt. this%map(i, sto%index_Q)**2 + this%map(i, sto%index_U)**2) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Pmin_Oriented)
             do i=0, this%npix-1
                if( mask%map(i, 1) .le. 0.5  .or. sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), sto%index_Q)**2 + this%map(list(1:nneigh), sto%index_U)**2 .gt. this%map(i, sto%index_Q)**2 + this%map(i, sto%index_U)**2) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_random_hot, coop_stacking_genre_random_cold, coop_stacking_genre_random_hot_oriented, coop_stacking_genre_random_cold_oriented)   !!nested, with mask
             allocate(livept(0:this%npix-1))
             livept = .false.
             do i = 0, this%npix-1
                if(mask%map(i,1).le.0.5 .or. sto%reject(this%map(i,:)))cycle
                livept(i) = .true.
             enddo
             if(present(nside_scan))then
                nside_rand = nside_scan
             else
                npix_rand  = count(livept)
                nside_rand = this%nside/2**max(nint(coop_log2(dble(npix_rand)/num_rand)/2.d0), 0)
             endif
             npix_rand = nside2npix(nside_rand)
                
             do i = 0, npix_rand - 1
                call pix2ang_nest(nside_rand, i, thetaphi(1), thetaphi(2))
                call this%ang2pix(thetaphi(1), thetaphi(2), ip)
                call this%pix2ang(ip, thetaphi(1), thetaphi(2))
                if(livept(ip))then
                   call sto%peak_pix%push(ip)
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(ip,:))
                endif
             enddo
             deallocate(livept)
          case default
             stop "unknown stacking genre"
          end select
       else   !!no mask, nested
          select case(sto%genre)
          case(coop_stacking_genre_saddle, coop_stacking_genre_saddle_oriented)
             call this%zeros(5, zeros1)
             call this%zeros(6, zeros2)
             zeros1%map = zeros1%map*zeros2%map
             do i=0, this%npix-1
                if(zeros1%map(i, 1) .gt. 0.5 .and. .not. sto%reject(this%map(i, :))  .and. this%map(i, 2)**2+this%map(i, 3)**2 .gt. this%map(i,4)**2)then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
             call zeros1%free()
             call zeros2%free()
          case(coop_stacking_genre_Imax, coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Lmax, coop_stacking_genre_Lmax_Oriented)
             do i=0, this%npix-1
                if( sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), index_peak).lt.this%map(i, index_peak)) )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Imin, coop_stacking_genre_Imin_Oriented, coop_stacking_genre_Lmin, coop_stacking_genre_Lmin_Oriented)
             do i=0, this%npix-1
                if( sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), index_peak) .gt. this%map(i, index_peak))  )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Pmax_Oriented)
             do i=0, this%npix-1
                if( sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), sto%index_Q)**2 + this%map(list(1:nneigh), sto%index_U)**2 .lt. this%map(i, sto%index_Q)**2 + this%map(i, sto%index_U)**2)  )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Pmin_Oriented)
             do i=0, this%npix-1
                if(sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), sto%index_Q)**2 + this%map(list(1:nneigh), sto%index_U)**2 .gt. this%map(i, sto%index_Q)**2 + this%map(i, sto%index_U)**2) )then
                   call sto%peak_pix%push(i)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_random_hot, coop_stacking_genre_random_cold, coop_stacking_genre_random_hot_oriented, coop_stacking_genre_random_cold_oriented)             !!nested, no mask

             allocate(livept(0:this%npix-1))
             livept = .false.
             do i = 0, this%npix-1
                if(sto%reject(this%map(i,:)))cycle
                livept(i) = .true.
             enddo
             npix_rand  = count(livept)
             nside_rand = this%nside/2**max(nint(coop_log2(dble(npix_rand)/num_rand)/2.d0), 0)
             npix_rand = nside2npix(nside_rand)
             do i = 0, npix_rand - 1
                call pix2ang_nest(nside_rand, i, thetaphi(1), thetaphi(2))
                call this%ang2pix(thetaphi(1), thetaphi(2), ip)
                call this%pix2ang(ip, thetaphi(1), thetaphi(2))
                if(livept(ip))then
                   call sto%peak_pix%push(ip)
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(ip,:))
                endif
             enddo
             deallocate(livept)
             
          end select
       end if
    else  !!ring ordering
       if(domask)then     !!with mask
          select case(sto%genre)
          case(coop_stacking_genre_saddle, coop_stacking_genre_saddle_Oriented)
             call this%zeros(5, zeros1, mask)
             call this%zeros(6, zeros2, mask)
             zeros1%map = zeros1%map*zeros2%map
             do i=0, this%npix-1
                if(zeros1%map(i, 1) .gt. 0.5 .and. .not. sto%reject(this%map(i, :))  .and. this%map(i, 2)**2+this%map(i, 3)**2 .gt. this%map(i,4)**2 )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
             call zeros1%free()
             call zeros2%free()
          case(coop_stacking_genre_Imax, coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Lmax, coop_stacking_genre_Lmax_Oriented)
             do i=0, this%npix-1
                if( mask%map(i, 1) .le. 0.5  .or. sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), index_peak).lt.this%map(i, index_peak)) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Imin, coop_stacking_genre_Imin_Oriented, coop_stacking_genre_Lmin, coop_stacking_genre_Lmin_Oriented)
             do i=0, this%npix-1
                if( mask%map(i, 1) .le. 0.5  .or. sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), index_peak) .gt. this%map(i, index_peak)) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Pmax_Oriented)
             do i=0, this%npix-1
                if( mask%map(i, 1) .le. 0.5  .or. sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), sto%index_Q)**2 + this%map(list(1:nneigh), sto%index_U)**2 .lt. this%map(i, sto%index_Q)**2 + this%map(i, sto%index_U)**2) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Pmin_Oriented)
             do i=0, this%npix-1
                if( mask%map(i, 1) .le. 0.5  .or. sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), sto%index_Q)**2 + this%map(list(1:nneigh), sto%index_U)**2 .gt. this%map(i, sto%index_Q)**2 + this%map(i, sto%index_U)**2) .and. all(mask%map(list(1:nneigh),1) .gt. 0.5 ) )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_random_hot, coop_stacking_genre_random_cold, coop_stacking_genre_random_hot_oriented, coop_stacking_genre_random_cold_oriented)             !!ring with mask
             
             allocate(livept(0:this%npix-1))
             livept = .false.
             do i = 0, this%npix-1
                if(mask%map(i,1).le.0.5 .or. sto%reject(this%map(i,:)))cycle
                livept(i) = .true.
             enddo
             npix_rand  = count(livept)
             nside_rand = this%nside/2**max(nint(coop_log2(dble(npix_rand)/num_rand)/2.d0), 0)
             npix_rand = nside2npix(nside_rand)
             do i = 0, npix_rand - 1
                call pix2ang_nest(nside_rand, i, thetaphi(1), thetaphi(2))
                call this%ang2pix(thetaphi(1), thetaphi(2), ip)
                call this%pix2ang(ip, thetaphi(1), thetaphi(2))
                if(livept(ip))then
                   call nest2ring(this%nside, ip, ipr)
                   call sto%peak_pix%push(ipr)
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(ip,:))
                endif
             enddo
             deallocate(livept)
          end select
       else  !!ring ordering, no mask
          select case(sto%genre)
          case(coop_stacking_genre_saddle, coop_stacking_genre_saddle_Oriented)
             call this%zeros(5, zeros1)
             call this%zeros(6, zeros2)
             zeros1%map = zeros1%map*zeros2%map
             do i=0, this%npix-1
                if(zeros1%map(i, 1) .gt. 0.5 .and. .not. sto%reject(this%map(i, :)) .and. this%map(i, 2)**2+this%map(i, 3)**2 .gt. this%map(i,4)**2)then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
             call zeros1%free()
             call zeros2%free()             
          case(coop_stacking_genre_Imax, coop_stacking_genre_Imax_Oriented, coop_stacking_genre_Lmax, coop_stacking_genre_Lmax_Oriented)
             do i=0, this%npix-1
                if( sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), index_peak).lt.this%map(i, index_peak)) )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Imin, coop_stacking_genre_Imin_Oriented, coop_stacking_genre_Lmin, coop_stacking_genre_Lmin_Oriented)
             do i=0, this%npix-1
                if( sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), index_peak) .gt. this%map(i, index_peak))  )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Pmax_Oriented)
             do i=0, this%npix-1
                if( sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), sto%index_Q)**2 + this%map(list(1:nneigh), sto%index_U)**2 .lt. this%map(i, sto%index_Q)**2 + this%map(i, sto%index_U)**2) )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_Pmin_Oriented)
             do i=0, this%npix-1
                if(  sto%reject(this%map(i,:)))cycle
                call neighbours_nest(this%nside, i, list, nneigh)  
                if ( all(this%map(list(1:nneigh), sto%index_Q)**2 + this%map(list(1:nneigh), sto%index_U)**2 .gt. this%map(i, sto%index_Q)**2 + this%map(i, sto%index_U)**2) )then
                   call nest2ring(this%nside, i, ip)
                   call sto%peak_pix%push(ip)
                   call this%pix2ang(i, thetaphi(1), thetaphi(2))
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(i, :))
                endif
             enddo
          case(coop_stacking_genre_random_hot, coop_stacking_genre_random_cold, coop_stacking_genre_random_hot_oriented, coop_stacking_genre_random_cold_oriented)             !!ring, no mask
             allocate(livept(0:this%npix-1))
             livept = .false.
             do i = 0, this%npix-1
                if(sto%reject(this%map(i,:)))cycle
                livept(i) = .true.
             enddo
             npix_rand  = count(livept)
             nside_rand = this%nside/2**max(nint(coop_log2(dble(npix_rand)/num_rand)/2.d0), 0)
             npix_rand = nside2npix(nside_rand)
             do i = 0, npix_rand - 1
                call pix2ang_nest(nside_rand, i, thetaphi(1), thetaphi(2))
                call this%ang2pix(thetaphi(1), thetaphi(2), ip)
                call this%pix2ang(ip, thetaphi(1), thetaphi(2))
                if(livept(ip))then
                   call nest2ring(this%nside, ip, ipr)
                   call sto%peak_pix%push(ipr)
                   call sto%peak_ang%push(real(thetaphi))
                   call sto%peak_map%push(this%map(ip,:))
                endif
             enddo
             deallocate(livept)
             
          end select          
       end if
    endif
    if(present(restore))then
       if(restore)then
          call this%convert2ring()
          if(domask) call mask%convert2ring()
       endif
    else  !!by default I do restore
       call this%convert2ring()
       if(domask) call mask%convert2ring()       
    endif
#else
    stop "Cannot find HEALPIX"
#endif    
  end subroutine coop_healpix_maps_get_peaks

  
  subroutine coop_healpix_maps_mark_peaks(this, sto, imap)
    COOP_INT,parameter::max_n = 30
    class(coop_healpix_maps) this
    COOP_INT imap
    type(coop_stacking_options)::sto
    type(coop_healpix_disc) disc    
    COOP_INT i, pix, nx, ny, ix, iy, listpix(0:this%npix-1), nlist
    COOP_SINGLE::thetaphi(2)
    COOP_REAL  angle, e, r, rsq, dr, x, y, cost, sint, eplus, eminus
    dr = sqrt(coop_4pi/this%npix)/3.d0  !!just to be conservative
    this%map(:, imap) = 0.
    do i = 1, sto%peak_ang%n
       call sto%peak_ang%get_element(i, thetaphi)
       call this%ang2pix(dble(thetaphi(1)), dble(thetaphi(2)), pix)
       call this%get_disc(pix, disc)
       call sto%peak_get_angle_r_e(i, angle, r, e)
       cost = cos(angle)
       sint = sin(angle)
       eplus = 1.d0+e
       eminus = 1.d0-e
       nx = min(nint(r/sqrt(eminus)/dr), max_n)
       ny = min(nint(r/sqrt(eplus)/dr), max_n)
       rsq = r**2
       do ix = -nx, nx
          do iy = -ny, ny
             x = ix*dr
             y = iy*dr
             if(x**2*eminus+y**2*eplus .gt. rsq)cycle
             call disc%xy2pix(x*cost - y*sint, x*sint + y*cost, pix)
             this%map(pix, imap) = 1.
          enddo
       enddo
    enddo
  end subroutine coop_healpix_maps_mark_peaks

  subroutine coop_healpix_maps_stack_on_peaks(this, sto, patch, mask, norm)
    COOP_INT,parameter::n_threads = 8
    class(coop_healpix_maps)::this
    type(coop_stacking_options)::sto
    type(coop_healpix_disc),dimension(n_threads)::disc
    type(coop_healpix_patch)::patch
    type(coop_healpix_patch),dimension(n_threads)::p, tmp
    type(coop_healpix_maps),optional::mask
    logical,optional::norm
    COOP_INT ithread, i
#ifdef HAS_HEALPIX
    patch%image = 0.d0
    patch%nstack = 0.d0
    patch%nstack_raw = 0
    if(sto%nested)then
       call this%convert2nested()
       if(present(mask))call mask%convert2nested()
    else
       call this%convert2ring()
       if(present(mask))call mask%convert2ring()
    endif
    !!adjust the index of maps
    i = 1
    do while(i.le.patch%nmaps)
       if(patch%tbs%ind(i) .gt. this%nmaps) stop "stack_on_peaks: index overflow"       
       do while(patch%tbs%spin(i) .ne. this%spin(patch%tbs%ind(i)))
          patch%tbs%ind(i)  = patch%tbs%ind(i) + 1
          if(patch%tbs%ind(i) .gt. this%nmaps) stop "stack_on_peaks: cannot find the map with corresponding spin"
       enddo
       if(patch%tbs%spin(i) .ne. 0)then
          if(i.ge.patch%nmaps .or. patch%tbs%ind(i) .ge. this%nmaps) stop "stack_on_peaks: nonzero spin must go in pairs"
          patch%tbs%ind(i+1) = patch%tbs%ind(i) + 1
          i = i+2
       else
          i = i+1
       endif
    enddo
    do ithread=1, n_threads
       p(ithread) = patch
       tmp(ithread) = patch
    enddo
    !$omp parallel do private(i, ithread)
    do ithread = 1, n_threads
       do i=ithread, sto%peak_pix%n, n_threads
          call this%get_disc(sto%pix(this%nside, i), disc(ithread))
          if(present(mask))then
             call this%stack_on_patch(disc(ithread), sto%rotate_angle(i), p(ithread), tmp(ithread), mask)    
          else
             call this%stack_on_patch(disc(ithread), sto%rotate_angle(i), p(ithread), tmp(ithread))
          endif
       enddo
    enddo
    !$omp end parallel do

    do ithread = 1, n_threads
       patch%image = patch%image + p(ithread)%image
       patch%nstack = patch%nstack + p(ithread)%nstack
       patch%nstack_raw = patch%nstack_raw + p(ithread)%nstack_raw
       call p(ithread)%free()
       call tmp(ithread)%free()
    enddo
    if(present(norm))then
       if(.not.norm)return
    endif
    if(patch%nstack_raw .ne. 0)then
       do i=1, patch%nmaps
          patch%image(:, :, i) = patch%image(:, :, i)/max(patch%nstack, 1.d0)
       enddo
    else
       write(*,*) "warning: no patches has been found"
       patch%image = 0.d0
    endif
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_stack_on_peaks


  
  subroutine coop_healpix_mask_peaks(mask, sto, sto_masked)
    class(coop_healpix_maps)::mask
    type(coop_stacking_options)::sto, sto_masked
    COOP_INT i, pix
    sto_masked%nested = sto%nested
    sto_masked%sigma_I = sto%sigma_I
    sto_masked%sigma_L = sto%sigma_L
    sto_masked%sigma_P = sto%sigma_P    
    if(sto%nested)then
       call mask%convert2nested()
    else
       call mask%convert2ring()
    endif
    call sto_masked%free()
    do i=1, sto%peak_pix%n
       call sto%peak_pix%get_element(i, pix)
       if(mask%map(pix,1).ge. 0.5d0)then
          call sto_masked%peak_pix%push(pix)
          call sto_masked%peak_ang%push(sto%peak_ang%element(i))
          call sto_masked%peak_map%push(sto%peak_map%element(i))
       endif
    enddo
  end subroutine coop_healpix_mask_peaks

  subroutine coop_healpix_maps_rotate_coor(this, output, from, to, zyz_psi, zyz_theta, zyz_phi, l_deg, b_deg)
    !!rotate map between coordinates
    !!from and to can be 'C'/'Q' for Celestial/eQuatorial, 'E' for Ecliptic, and 'G' for Galactic    
    class(coop_healpix_maps)::this
    COOP_UNKNOWN_STRING, optional::from, to
    COOP_REAL,optional::zyz_psi, zyz_theta, zyz_phi, l_deg, b_deg    
    type(coop_healpix_maps),optional::output
#ifdef HAS_HEALPIX
    COOP_REAL psi, theta, phi
    COOP_INT i, istart
    COOP_SINGLE_COMPLEX,dimension(:,:,:),allocatable::alm_TGC
    if(present(from) .and. present(to))then
       call coordsys2euler_zyz(2000.d0, 2000.d0, from, to, psi, theta, phi)
    elseif(present(l_deg) .and. present(b_deg))then
       call coop_healpix_lb2ang(l_deg = l_deg, b_deg = b_deg, theta = theta, phi = phi )
    else
       if(present(zyz_psi))then
          psi = zyz_psi
       else
          psi = 0.d0
       endif
       if(present(zyz_theta))then
          theta = zyz_theta
       else
          theta = 0.d0
       endif
       if(present(zyz_phi))then
          phi = zyz_phi
       else
          phi = 0.d0
       endif
    endif

    if(psi .eq. 0.d0 .and. theta .eq. 0.d0 .and. phi .eq. 0.d0)return
    if(.not. allocated(this%alm))call this%map2alm()
    select case(this%nmaps)
    case(1)
       allocate(alm_TGC(1, 0:this%lmax, 0:this%lmax))
       istart  = 0
    case(2:3)
       allocate(alm_TGC(3, 0:this%lmax, 0:this%lmax))
       istart = 3 - this%nmaps
    case default
       stop "rotate_coor only works for nmaps = 1, 2, 3"
    end select
    if(istart.gt.0) alm_TGC (1:istart, :,:)= 0.
    do i=1, this%nmaps
       alm_TGC(istart + i,:,:)= this%alm(:,:,i)
    enddo
    if(present(l_deg) .and. present(b_deg))then
       call rotate_alm(this%lmax, alm_TGC, 0.d0, 0.d0, -phi)
       call rotate_alm(this%lmax, alm_TGC, 0.d0, -theta, 0.d0)       
    else
       call rotate_alm(this%lmax, alm_TGC, psi, theta, phi)
    endif
    if(present(output))then
       if(output%nmaps .ne. this%nmaps) stop "rotate_coor: different format"
       if(any(output%spin .ne. this%spin)) stop "rotate_coor: different spins"
       if(output%lmax .ne. this%lmax)call output%allocate_alms(lmax = this%lmax)
       
       do i=1, this%nmaps
          output%alm(:,:,i) = alm_TGC(istart + i,:,:)
       enddo
       call output%alm2map()
       if(present(to)) output%coordsys = trim(to)
    else
       do i=1, this%nmaps
          this%alm(:,:,i) = alm_TGC(istart + i,:,:)
       enddo
       call this%alm2map()
       if(present(to)) this%coordsys = trim(to)
    endif
    deallocate(alm_TGC)

#else
    stop "You need to install healpix"
#endif    
  end subroutine coop_healpix_maps_rotate_coor


  subroutine coop_healpix_maps_distr_nu_e(this, mask, numin, numax, emin, emax, nnu, ne, output)
    class(coop_healpix_maps)::this
    type(coop_healpix_maps)::mask
    COOP_REAL:: numin, numax, emin, emax, nu, e, dnu, de, sigmaI
    COOP_INT:: nnu, ne, indexL, indexQ, indexU, indexI,  inu, ie, i
    COOP_REAL::density(nnu, ne)
    COOP_UNKNOWN_STRING::output
    type(coop_asy)::fig
    indexI = 1
    indexQ = 2
    indexU = 3
    if(nnu .le. 1 .or. ne .le. 1)stop "distr_nu_e: nnu and ne must be >1"
    select case(this%nmaps)
    case(3)
       indexL = 1
    case(4)
       indexL = 4
    case default
       stop "distr_nu_e only supports nmaps = 3 or 4"
    end select
    sigmaI = sqrt(sum(dble(this%map(:,indexI)*mask%map(:,1))**2)/sum(dble(mask%map(:,1))))
    dnu = (numax- numin)/(nnu-1)
    de = (emax - emin)/(ne-1)
    density = 0.d0
    do i = 0, this%npix - 1
       nu = this%map(i, indexI)/sigmaI
       e = (this%map(i, indexQ)**2 + this%map(i, indexU)**2)/max(this%map(i, indexL)**2, 1.d-20)
       inu = nint((nu-numin)/dnu + 1.d0)
       if(inu .lt. 1 .or. inu .gt. nnu)then
          cycle
       endif
       ie = nint((e-emin)/de + 1.d0)
       if(ie.lt.1 .or. ie .gt. ne)then
          cycle
       endif
       density(inu, ie) = density(inu, ie) + 1.d0
    enddo
    density = log10(max(density/(dnu*de), 1.d-5))
    call fig%open(output)
    call fig%init(xlabel = "$\nu$", ylabel = "$e$")
    call fig%density(density, numin, numax, emin, emax, "$\log_{10}{\frac{d^2N}{d\nu de}}$")
    
    call fig%close()
  end subroutine coop_healpix_maps_distr_nu_e

  subroutine coop_healpix_maps_ave_udgrade(from, to, mask, imap_from, imap_to)
    class(coop_healpix_maps)::from, to
    class(coop_healpix_maps),optional::mask    
    COOP_INT::div, i, imap, nmaps
    COOP_INT,optional::imap_from, imap_to
    COOP_REAL::summask
    if(present(mask))then
       if(mask%nside .ne. from%nside)stop "ave_udgrade: nside(mask) != nside(from)"
    endif
    nmaps = min(from%nmaps, to%nmaps)
    call coop_healpix_maps_copy_genre(from, to)
    call from%convert2nested()
    call to%convert2nested()
    if(from%nside .eq. to%nside)then
       if(present(imap_from) .and. present(imap_to))then
          if(present(mask))then
             to%map(:, imap_to) = from%map(:, imap_from)*mask%map(:,1)
          else
             to%map(:, imap_to) = from%map(:, imap_from)
          endif
       else
          if(present(mask))then
             do imap = 1, nmaps
                to%map(:, imap) = from%map(:, imap)*mask%map(:,1)
             enddo
          else
             to%map(:, 1:nmaps) = from%map(:, 1:nmaps)                
          endif
       endif
       return
    endif
    if(present(imap_from) .and. present(imap_to))then
       if(present(mask))then
          if(from%nside .gt. to%nside)then !!degrade
             div = (from%nside/to%nside)**2
             do i=0, to%npix-1
                summask =  sum(dble(mask%map(i*div:(i+1)*div-1, 1)))
                if(summask .gt. 0.01d0)then
                   to%map(i, imap_to) = sum(dble(from%map(i*div:(i+1)*div-1, imap_from)*mask%map(i*div:(i+1)*div-1, 1)))/summask
                else
                   to%map(i, imap_to) = 0.
                endif
             enddo
             return
          endif
          !!upgrade
          div = (to%nside/from%nside)**2
          do i=0, from%npix-1
             to%map(i*div:(i+1)*div-1, imap_to) = from%map(i, imap_from)*mask%map(i, 1)
          enddo
       else
          if(from%nside .gt. to%nside)then !!degrade
             div = (from%nside/to%nside)**2
             do i=0, to%npix-1
                to%map(i, imap_to) = sum(dble(from%map(i*div:(i+1)*div-1, imap_from)))/div
             enddo
             return
          endif
          !!upgrade
          div = (to%nside/from%nside)**2
          do i=0, from%npix-1
             to%map(i*div:(i+1)*div-1, imap_to) = from%map(i, imap_from)
          enddo
       endif
    else
       if(present(mask))then
          if(from%nside .gt. to%nside)then !!degrade
             div = (from%nside/to%nside)**2
             do i=0, to%npix-1
                summask =  sum(dble(mask%map(i*div:(i+1)*div-1, 1)))
                if(summask .gt. 0.01d0)then
                   do imap = 1, nmaps
                      to%map(i, imap) = sum(dble(from%map(i*div:(i+1)*div-1, imap)*mask%map(i*div:(i+1)*div-1, 1)))/summask
                   enddo
                else
                   to%map(i, 1:nmaps) = 0.
                endif
             enddo
             return
          endif
          !!upgrade
          div = (to%nside/from%nside)**2
          do imap = 1, nmaps
             do i=0, from%npix-1
                to%map(i*div:(i+1)*div-1, imap) = from%map(i, imap)*mask%map(i, 1)
             enddo
          enddo
       else
          if(from%nside .gt. to%nside)then !!degrade
             div = (from%nside/to%nside)**2
             do imap = 1, nmaps
                do i=0, to%npix-1
                   to%map(i, imap) = sum(dble(from%map(i*div:(i+1)*div-1, imap)))/div
                enddo
             enddo
             return
          endif
          !!upgrade
          div = (to%nside/from%nside)**2
          do imap = 1, nmaps
             do i=0, from%npix-1
                to%map(i*div:(i+1)*div-1, imap) = from%map(i, imap)
             enddo
          enddo
       endif
    endif
  end subroutine coop_healpix_maps_ave_udgrade

  subroutine coop_healpix_maps_udgrade(this, nside)
    class(coop_healpix_maps)::this
    COOP_SINGLE, dimension(:,:),allocatable::newmap
    COOP_INT nside, npix
#ifdef HAS_HEALPIX
    npix = nside2npix(nside)
    allocate(newmap(0:npix-1, this%nmaps))
    if(this%ordering .eq. COOP_NESTED)then
       call udgrade_nest(this%map, this%nside, newmap, nside, 0., .true.)
    else
       call udgrade_ring(this%map, this%nside, newmap, nside, 0., .true.)
    endif
    deallocate(this%map)
    this%nside = nside
    this%npix = npix
    allocate(this%map(0:npix-1, this%nmaps))
    this%map = newmap
    deallocate(newmap)
#else
    stop "HEALPIX library is missing."
#endif    
  end subroutine coop_healpix_maps_udgrade

  
  subroutine coop_healpix_inpaint_free(this)
    class(coop_healpix_inpaint)::this
    nullify(this%map)
    nullify(this%mask)
    if(allocated(this%lask))deallocate(this%lask)
    if(allocated(this%corr))deallocate(this%corr)
    if(allocated(this%als))deallocate(this%als)
    if(allocated(this%bls))deallocate(this%bls)
    if(allocated(this%cls))deallocate(this%cls)
    if(allocated(this%sqrtcls))deallocate(this%sqrtcls)
    if(allocated(this%smooth_cls))deallocate(this%smooth_cls)        
    if(allocated(this%vec))deallocate(this%vec)
    if(allocated(this%Fmean))deallocate(this%Fmean)
    if(allocated(this%Ffluc))deallocate(this%Ffluc)    
    if(allocated(this%mean))deallocate(this%mean)
    if(allocated(this%indMT))deallocate(this%indMT)
    if(allocated(this%indCT))deallocate(this%indCT)
    if(allocated(this%iiMT))deallocate(this%iiMT)
    if(allocated(this%iiCT))deallocate(this%iiCT)            
    call this%lMT%free()
    call this%lcT%free()
    call this%lM%free()
    call this%sim%free()
    call this%smooth_mask%free()
    call this%smooth_map%free()    
    this%lmax = -1
    this%ncorr = 0
    this%dtheta = 0.d0
    this%nMT = 0
    this%ncT = 0
    this%base_nside = 0
    this%first_realization = .true.
  end subroutine coop_healpix_inpaint_free



  subroutine coop_healpix_inpaint_set_corr(this, lmax, Cls)
    class(coop_healpix_inpaint)::this
    COOP_INT, optional::lmax
    COOP_REAL,dimension(:),optional::Cls
    COOP_INT::l, i
    COOP_REAL::thetasq, theta, clr1, clr2, cross, norm
    type(coop_file)::fp
#ifdef HAS_HEALPIX
    if(present(lmax) .and. present(Cls))then
       if(this%lmax .ge. 0) deallocate(this%als, this%bls, this%cls, this%sqrtcls, this%smooth_cls)
       this%lmax = lmax
       allocate(this%als(0:lmax), this%bls(0:lmax), this%cls(0:lmax), this%sqrtcls(0:lmax), this%smooth_cls(0:lmax))
       call coop_sphere_correlation_init(this%lmax, this%als, this%bls)
       this%cls = cls
       this%smooth_cls = cls       
       this%sqrtcls = sqrt(this%cls)
    endif

    
    this%ncorr = max(min(this%lmax*4, this%lM%nside*6), 64)
    this%dtheta = coop_pi/this%ncorr
    if(allocated(this%corr))deallocate(this%corr)
    allocate(this%corr(0:this%ncorr))
    if(this%lM%nside .lt. this%map%nside)then
       theta = 0.3/this%lM%nside  !!0.3 is the universal number that gives the effective pix window function
    else
       theta = 0.d0
    endif
    do l = 0, this%lmax
       this%smooth_Cls(l) = this%cls(l)*exp(-(l*theta)**2)
    enddo
    this%sigma0 = sqrt(coop_sphere_correlation(this%lmax, this%smooth_cls, this%als, this%bls, 1.d0))
    !$omp parallel do
    do i=0, this%ncorr
       this%corr(i) = coop_sphere_correlation(this%lmax, this%smooth_cls, this%als, this%bls, cos(i*this%dtheta))
    enddo
    !$omp end parallel do
    if(allocated(this%vec))then
       if(size(this%vec, 2) .eq. this%lM%npix) return
       deallocate(this%vec)
    endif
    allocate(this%vec(3, 0:this%lM%npix-1))
    !$omp parallel do
    do i=0, this%lM%npix-1
       call pix2vec_nest(this%lM%nside, i, this%vec(:, i))
    enddo
    !$omp end parallel do
#endif    
  end subroutine coop_healpix_inpaint_set_corr
    

  function coop_healpix_inpaint_correlation(this, i, j, costheta) result(c)
    class(coop_healpix_inpaint)::this
    COOP_REAL :: x, c, ri
    COOP_INT,optional :: i, j
    COOP_INT :: loc
    COOP_REAL,optional::costheta
    if(present(i) .and. present(j))then
       x = dot_product(this%vec(:, i), this%vec(:, j))
    elseif(present(costheta))then
       x = costheta
    else
       stop "you need to pass i, j or costheta"
    endif
    if(x .ge. 0.99999999d0)then
       c = this%corr(0)
       return
    endif
    if(x .le. -0.99999999d0)then
       c = this%corr(this%ncorr)
       return
    endif
    ri = acos(x)/this%dtheta
    loc = floor(ri)
    ri = ri - loc
    c = this%corr(loc)*(1.d0-ri) + this%corr(loc+1)*ri
  end function coop_healpix_inpaint_correlation

  subroutine coop_healpix_inpaint_eval_cov(this, i, j, mm, cc, mc, cm, cf,mf, ff) 
    class(coop_healpix_inpaint)::this    
    COOP_REAL,optional::mm, mc, cc,cm, cf, mf, ff
    COOP_REAL::tmp
    logical::do_single
    COOP_INT::i, j, ii, jj
    tmp =  this%correlation(i, j)
    if(present(mm)) mm =   this%lM%map(i,1)*this%lM%map(j, 1)*tmp
    if(present(cc)) cc =  (1.-this%lM%map(i,1))*(1.-this%lM%map(j, 1))*tmp
    if(present(cf)) cf =  (1.-this%lM%map(i, 1))*tmp
    if(present(mf)) mf =  this%lM%map(i, 1)*tmp          
    if(present(cm)) cm =  (1.-this%lM%map(i,1))*this%lM%map(j, 1)*tmp
    if(present(mc)) mc =  this%lM%map(i,1)*(1.-this%lM%map(j, 1))*tmp
    if(present(ff)) ff =  tmp
  end subroutine coop_healpix_inpaint_eval_cov

  subroutine coop_healpix_inpaint_init(this, map, mask, lmax, Cls)
    class(coop_healpix_inpaint)::this    
    type(coop_healpix_maps), target::map, mask
    COOP_INT::lmax, i, iM, ic, j, info
    COOP_REAL::Cls(0:lmax), mindiag
    COOP_REAL, dimension(:, :),allocatable::cov
    COOP_REAL::mm
#ifdef HAS_HEALPIX    
    call this%free()
    if(map%nside .ne. mask%nside) stop "inpaint_init: error"
    this%map => map
    this%mask => mask
    call this%map%convert2nested()
    call this%mask%convert2nested()
    call coop_healpix_mask_diffuse(this%mask, this%smooth_mask, 0.5d0*coop_SI_degree)
    call coop_healpix_maps_diffuse(this%map, this%smooth_map, this%mask, 0.5d0*coop_SI_degree)
    !$omp parallel do
    do i=0, this%smooth_mask%npix-1
       if(this%smooth_mask%map(i, 1) .ge. 0.5)then
          this%smooth_map%map(i,1) = this%smooth_map%map(i,1) / this%smooth_mask%map(i,1)
       else
          this%smooth_map%map(i,1) = this%smooth_map%map(i,1) / max(sqrt(this%smooth_mask%map(i,1)/2.), 0.1)
       endif
    enddo
    !$omp end parallel do
    this%base_nside = min(this%map%nside, coop_inpaint_nside_start)
    call this%lMT%init(nside = this%base_nside, nmaps = 1, genre = "TEMPERATURE", nested = .true.)
    this%lcT = this%lMT
    call this%lM%init(nside = this%base_nside, nmaps = 1, genre = "MASK", nested = .true.)

    
    allocate(this%lask(0:this%mask%npix-1))
    this%lask = (this%mask%map(:,1) .gt. 0.5)
    this%fsky = count(this%lask)/dble(this%mask%npix)
    this%map%map(:,1) = this%map%map(:,1)*this%mask%map(:,1) !!mask the map
    call coop_healpix_maps_ave_udgrade(this%map, this%lMT)  !!M T at low resolution is known;
    call this%lask2mask(this%lM)
    call this%sim%init(nside = this%map%nside, nmaps = 1, genre = "T", nested = .true., lmax = min(lmax, floor(this%map%nside*coop_healpix_lmax_by_nside)))        
    call this%set_corr(lmax, Cls)    !!set up the correlation function calculator

    this%nMT = count(this%lM%map(:,1) .ge. coop_inpaint_mask_threshold)
    this%ncT = count(this%lM%map(:,1) .le. 1.-coop_inpaint_mask_threshold)
    
    allocate(this%mean(this%ncT))
    allocate(this%Ffluc(this%ncT, this%ncT), this%Fmean(this%nCT, this%nMT))
    allocate(this%indMT(this%nMT), this%indcT(this%nCT), this%iiMT(0:this%lM%npix-1), this%iiCT(0:this%lM%npix-1))
    this%iiMT = 0
    this%iiCT = 0
    iM = 0
    ic = 0
    do i=0, this%lM%npix-1
       if(this%lM%map(i,1) .ge. coop_inpaint_mask_threshold)then
          iM = iM + 1
          this%indMT(iM) = i
          this%iiMT(i) = iM
       endif
       if(this%lM%map(i,1) .le. 1.-coop_inpaint_mask_threshold)then
          ic = ic + 1
          this%indCT(ic) = i
          this%iiCT(i) = ic + this%nMT 
       endif
    enddo
    if(iM .ne. this%nMT .or. ic .ne. this%ncT) stop "inpaint_start: unknown error"
    allocate(cov(this%nMT+this%ncT, this%nMT+this%ncT))
    do i = 0, this%lM%npix-1
       do j = 0, i
          if(this%iiMT(i) .eq. 0)then
             if(this%iiMT(j).eq.0)then
                call this%eval_cov(i, j, cc = cov(this%iicT(i), this%iicT(j)))
                cov( this%iicT(j), this%iicT(i)) = cov(this%iicT(i), this%iicT(j))                
             elseif(this%iiCT(j) .eq. 0)then
                call this%eval_cov(i, j, cm = cov(this%iicT(i), this%iiMT(j)))
                cov(this%iiMT(j), this%iicT(i)) = cov(this%iicT(i), this%iiMT(j))                 
             else
                call this%eval_cov(i, j, cc = cov(this%iicT(i), this%iicT(j)),  cm = cov(this%iicT(i), this%iiMT(j)))
                cov(this%iicT(j), this%iicT(i)) =  cov(this%iicT(i), this%iicT(j))
                cov(this%iiMT(j), this%iicT(i)) = cov(this%iicT(i), this%iiMT(j))                                 
             endif
          elseif(this%iiCT(i) .eq. 0)then
             if(this%iiMT(j).eq.0)then
                call this%eval_cov(i, j, mc = cov(this%iiMT(i), this%iicT(j)))
                cov( this%iicT(j), this%iiMT(i) ) = cov( this%iiMT(i), this%iicT(j))                
             elseif(this%iiCT(j) .eq. 0)then
                call this%eval_cov(i, j, mm = cov(this%iiMT(i), this%iiMT(j)))
                cov( this%iiMT(j), this%iiMT(i) ) = cov( this%iiMT(i), this%iiMT(j))
             else
                call this%eval_cov(i, j, mm = cov(this%iiMT(i), this%iiMT(j)), mc = cov(this%iiMT(i), this%iicT(j)))
                cov( this%iiMT(j), this%iiMT(i) ) = cov( this%iiMT(i), this%iiMT(j))                
                cov( this%iicT(j), this%iiMT(i) ) = cov( this%iiMT(i), this%iicT(j))                
             endif
          else
             if(this%iiMT(j).eq.0)then
                call this%eval_cov(i, j, mc = cov(this%iiMT(i), this%iiCT(j)), cc = cov(this%iiCT(i), this%iicT(j)))
                cov( this%iicT(j), this%iiMT(i) ) = cov( this%iiMT(i), this%iicT(j))
                cov( this%iicT(j), this%iicT(i)) = cov(this%iicT(i), this%iicT(j))                                
             elseif(this%iiCT(j) .eq. 0)then
                call this%eval_cov(i, j, mm = cov(this%iiMT(i), this%iiMT(j)), cm = cov(this%iiCT(i), this%iiMT(j)))
                cov( this%iiMT(j), this%iiMT(i) ) = cov( this%iiMT(i), this%iiMT(j))
                cov(this%iiMT(j), this%iicT(i)) = cov(this%iicT(i), this%iiMT(j))                                 
             else
                call this%eval_cov(i, j, mm = cov(this%iiMT(i), this%iiMT(j)), cm = cov(this%iiCT(i), this%iiMT(j)),  mc = cov(this%iiMT(i), this%iiCT(j)), cc = cov(this%iiCT(i), this%iicT(j)))

                cov( this%iicT(j), this%iiMT(i) ) = cov( this%iiMT(i), this%iicT(j))
                cov( this%iicT(j), this%iicT(i)) = cov(this%iicT(i), this%iicT(j))
                cov( this%iiMT(j), this%iiMT(i) ) = cov( this%iiMT(i), this%iiMT(j))
                cov(this%iiMT(j), this%iicT(i)) = cov(this%iicT(i), this%iiMT(j))                                 
             endif
          endif
       enddo
    enddo
    if(this%nMT .gt. 0)then
       call coop_solve_constrained(m = this%nMT + this%nCT, n_known = this%nMT, n_unknown = this%nCT, dim_fmean = this%nCT, dim_ffluc = this%nCT, C = cov, Fmean = this%Fmean, Ffluc = this%Ffluc, epsilon = 1.d-3)
       this%mean = matmul(this%Fmean, this%lMT%map(this%indMT,1))
    else
       this%mean = 0.d0
       this%Ffluc = cov
       call coop_matsym_sqrt(this%FFluc, 1.d-10)
    endif
    this%base_nside = this%lMT%nside
    this%first_realization = .true.
    deallocate(cov)
#endif    
  end subroutine coop_healpix_inpaint_init

  subroutine coop_healpix_inpaint_upgrade(this, reset, nside_want)
    class(coop_healpix_inpaint)::this
    logical, optional::reset
    COOP_INT,optional::nside_want
    COOP_INT::nside_start, nside_target
    COOP_INT::i
#ifdef HAS_HEALPIX
    if(present(nside_want))then
       nside_target = min(max(nside_want, this%base_nside), this%map%nside)
    else
       nside_target = this%map%nside
    endif
    if(present(reset))then
       if(reset)then       
          if(this%lM%nside .ne. this%base_nside)then
             call this%lCT%init(nside = this%base_nside, nmaps = 1, genre = "TEMPERATURE", nested = .true.)
             call this%lMT%init(nside = this%base_nside, nmaps = 1, genre = "TEMPERATURE", nested = .true.)             
             call coop_healpix_maps_ave_udgrade(this%map, this%lMT)
             call this%lM%init(nside = this%base_nside, nmaps = 1, genre = "MASK", nested = .true.)
             call this%lask2mask(this%lM)
          endif
          this%first_realization = .true.
       endif
    endif
    nside_start = this%lCT%nside
    !!step 0: nside = base_nside
    if(nside_start .eq. this%base_nside .and. this%first_realization)then
       do i = 0, this%lM%npix-1
          if(this%lM%map(i,1) .gt.  1.-coop_inpaint_mask_threshold)then
             this%lCT%map(i,1) = this%lMT%map(i,1)* (1.d0-this%lM%map(i,1))/this%lM%map(i,1) 
          endif
       enddo
       this%lCT%map(this%indCT, 1) = this%mean + matmul(this%Ffluc, coop_random_gaussian_vector(this%nCT))
       this%first_realization = .false.
    endif
    if(this%lMT%nside .ge. nside_target)return           
    !! step 1: nside = base_nside * 2
    if(this%lMT%nside .eq. this%base_nside)then
       call do_split()
       call set_sim()
    endif
    if(this%lMT%nside .ge. nside_target)return                  
    call this%lMT%init(nside = nside_target, nmaps = 1, genre = "TEMPERATURE", nested = .true.)
    call coop_healpix_maps_ave_udgrade(this%map, this%lMT)
    call this%lCT%init(nside = nside_target, nmaps = 1, genre = "TEMPERATURE", nested = .true.)    
    call coop_healpix_maps_ave_udgrade(this%sim, this%lCT)
  contains

    subroutine  do_split()
      COOP_INT:: nside_fine, npix_fine, nside_coarse, npix_coarse
      COOP_INT:: i, j, k, l,n_masked, ii, i0, i1, i2, i3, jj
      type(coop_healpix_maps)::llc    
      COOP_INT,dimension(:,:),allocatable::list
      COOP_INT,dimension(:),allocatable::ind_masked, nlist, nbs
      logical,dimension(:),allocatable::is_masked
      COOP_REAL,parameter::radius = 15.d0*coop_SI_degree
      COOP_INT,parameter::nmax = 128
      COOP_REAL::cov(0:nmax-1, 0:nmax-1), Fmean(1, 0:nmax-1), FFluc(1, 1), x(0:nmax-1)
      COOP_INT::listpix(0:nmax-1)
      COOP_REAL::cc(0:3, 0:3)
      
      nside_coarse = this%lMT%nside
      npix_coarse = nside2npix(nside_coarse)
      nside_fine = nside_coarse*2
      npix_fine = nside2npix(nside_fine)
      if(3.d0*(radius*nside_fine)**2 .gt. nmax*1.1)then
         stop "do_split: you need to increase nmax"
      endif
     

      call this%lMT%init(nside = nside_fine, nmaps = 1, genre = "TEMPERATURE", nested = .true.)
      llc = this%lCT
      call this%lCT%init(nside = nside_fine, nmaps = 1, genre = "TEMPERATURE", nested = .true.)
      call coop_healpix_maps_ave_udgrade(llc, this%lCT)
      call coop_healpix_maps_ave_udgrade(this%map, this%lMT)
      call this%lM%init(nside = nside_fine, nmaps = 1, genre = "MASK", nested = .true.)
      call this%lask2mask(this%lM)
      call this%set_corr()
      allocate(is_masked(0:npix_fine-1))
      is_masked = (this%lM%map(:,1) .le.  1.-coop_inpaint_mask_threshold)
      n_masked = count(is_masked)
      allocate(ind_masked(n_masked), list(0:nmax-1, n_masked), nlist(n_masked), nbs(n_masked))
      j = 0
      do i=0, npix_fine-1
         if(is_masked(i))then
            j = j +1
            ind_masked(j) = i
         endif
      enddo
      do j = 1, n_masked
         i0 = ind_masked(j)
         ii = (i0/4)*4         
         select case(i-ii)
         case(0)
            i1 = ii+1
            i2 = ii+2
            i3 = ii+3
         case(1)
            i1 = ii
            i2 = ii+2
            i3 = ii+3
         case(2)
            i1 = ii
            i2 = ii+1
            i3 = ii+3
         case(3)
            i1 = ii
            i2 = ii+1
            i3 = ii+2
         end select
         call query_disc(nside_fine, this%vec(:, i0), radius, list(0:nmax-1, j), nlist(j), nest = 1, inclusive = 0)
         if(nlist(j) .ge. nmax) stop "do split: nlist >= nmax"
         do k=0, nlist(j)-1
            if(list(k,j).eq.i0)then
               list(k,j) = list(0, j)
               list(0, j) = i0               
               exit
            endif
         enddo
         do k=1, nlist(j)-1
            if(list(k,j).eq.i1)then
               list(k,j) = list(1,j)
               list(1,j) = i1
               exit
            endif
         enddo
         do k=2, nlist(j)-1
            if(list(k,j).eq.i2)then
               list(k,j) = list(2,j)
               list(2,j) = i2
               exit
            endif
         enddo
         do k=3, nlist(j)-1
            if(list(k,j).eq.i3)then
               list(k,j) = list(3,j)
               list(3,j) = i3
               exit
            endif
         enddo
         nbs(j) = count(is_masked(list(1:nlist(j)-1, j)))
      enddo
      do
         j = coop_minloc(nbs)
         if(nbs(j) .ge. nmax)exit
         do ii=0,3
            do jj = 0, ii
               call this%eval_cov(list(ii, j), list(jj, j), cc = cc(ii,jj))
            enddo
         enddo
         cov(nlist(j), nlist(j)) = (cc(1,1)+cc(2,2)+cc(3,3)+2.d0*(cc(2,1)+cc(3,1)+cc(3,2)) - 6.d0*sum(cc(1:3,0)) + cc(0,0)*9.d0)/16.d0
         call this%eval_cov(list(0, j), list(0, j), mm = cov(0, 0))
         do i=1, nlist(j)-1
            call this%eval_cov(list(0, j), list(i, j), cf = cc(0,0), mf = cov(0, i))
            cov(i, 0) = cov(0, i)
            do ii=1, 3
               call this%eval_cov(list(ii, j), list(i, j), cf = cc(ii,0))
            enddo
            cov(nlist(j), i) = (cc(0,0)*3.d0-sum(cc(1:3,0)))/4.d0
            cov(i, nlist(j)) = cov(nlist(j), i)
            do k=1, i
               call this%eval_cov(list(i, j), list(k, j), ff = cov(i, k))
               cov(k, i) = cov(i, k)
            enddo
            x(i) = this%lMT%map(list(i,j), 1) + this%lCT%map(list(i,j), 1)
         enddo
         x(0) = this%lMT%map(list(0, j), 1)
         call coop_solve_constrained(m = nmax, n_known = nlist(j), n_unknown = 1, dim_fmean = 1, dim_ffluc = 1, C = cov, Fmean =  Fmean, Ffluc = FFluc,  epsilon = 1.d-2)
         this%lCT%map(list(0, j),1) =  this%lCT%map(list(0, j),1) + dot_product(Fmean(1,0:nlist(j)-1), x(0:nlist(j)-1)) + FFluc(1,1)*coop_random_Gaussian()
         is_masked(ind_masked(j)) = .false.
         i = list(0,j)/4
         if(.not. any(is_masked(4*i:4*i+3)))then
            this%lCT%map(4*i:4*i+3,1) = this%lCT%map(4*i:4*i+3,1) - sum(this%lCT%map(4*i:4*i+3,1))/4.d0 + llc%map(i,1)
         endif
         do i = 1, n_masked
            if(is_masked(ind_masked(i)))then
               nbs(i) = count(is_masked(list(1:nlist(i)-1, i)))
            else
               nbs(i) = nmax
            endif
         enddo
      enddo
      call llc%free()
      deallocate(is_masked, ind_masked, list, nlist, nbs)
    end subroutine do_split


    subroutine set_sim()
      COOP_REAL,parameter::corr_scale = 0.3d0
      COOP_INT::nr, i, l, m
      COOP_REAL::corr, fluc, thiscl, mean
      COOP_REAL::weight(-5:5)
      nr = (this%map%nside/this%lMT%nside)**2      
      do i=0, this%lMT%npix-1
         this%sim%map(i*nr:(i+1)*nr-1,1) = (this%lMT%map(i, 1) + this%lCT%map(i,1))
      enddo
      call this%sim%map2alm(lmax = this%lmax)
      do l= 8, this%sim%lmax
         if(l.le.20)then
            do i = 1, 5
               weight(i) = exp(-(10.d0/l*i)**2)
               weight(-i) = weight(i)
            enddo
            weight(0) = 1.d0
            weight = weight/sum(weight)
         endif
         if(l .lt. this%lMT%nside*coop_healpix_lmax_by_nside)then
            corr = exp(-((l/dble(this%LMT%nside) - l/dble(this%sim%nside))*corr_scale)**2)                     
            thiscl = sum(this%sim%cl(l-5:l+5,1)*weight)
            mean = sqrt(this%cls(l)/thiscl)*corr
            fluc = sqrt(this%cls(l)*(1.-corr**2))
            
         else
            mean = 0.d0
            fluc = sqrt(this%Cls(l))
         endif
         this%sim%alm(l, 0, 1) = this%sim%alm(l, 0, 1)*mean + fluc*coop_random_complex_Gaussian(.true.)
         do m=1, l
            this%sim%alm(l, m, 1) = this%sim%alm(l, m, 1)*mean + fluc*coop_random_complex_Gaussian()            
         enddo
      enddo
      call this%sim%alm2map()
      call this%sim%convert2nested()
      this%sim%map = (this%sim%map*sqrt(1.-this%smooth_mask%map**2)+this%smooth_map%map*this%smooth_mask%map)*(1.-this%mask%map)
    end subroutine set_sim
    
#endif    
  end subroutine coop_healpix_inpaint_upgrade

  subroutine  coop_healpix_inpaint_lask2mask(this, mask)
    class(coop_healpix_inpaint)::this
    type(coop_healpix_maps)::mask
    COOP_INT::i, nfac
    COOP_REAL::rfac
    if(mask%nside .gt. this%mask%nside) stop "last2mask only supports downgrade"
    nfac = (this%mask%nside/mask%nside)**2
    rfac = dble(nfac)
    do i = 0, mask%npix-1
       mask%map(i,1) = count(this%lask(i*nfac:(i+1)*nfac-1))/rfac
    enddo
  end subroutine coop_healpix_inpaint_lask2mask

  subroutine coop_healpix_maps_eval_corr(this, n, corr, xmin, xmax, mask)
    class(coop_healpix_maps)::this
    COOP_INT::n
    COOP_REAL::corr(n), weight(n)
    COOP_REAL::xmin, xmax
    type(coop_healpix_maps),optional::mask
    COOP_REAL,dimension(:,:),allocatable::vec
    COOP_REAL::cost, dx
    COOP_INT::i, j, loc
    dx = (xmax-xmin)/n
    if(this%nside .le. 16)then  !!brute-force
       allocate(vec(3, 0:this%npix-1))
       !$omp parallel do
       do i=0, this%npix-1
          call this%pix2vec(i, vec(:, i))
       enddo
       !$omp end parallel do
       if(present(mask))then
          do i=0, this%npix-1
             if(mask%map(i,1).lt. 0.5)cycle
             do j=1, i-1
                if(mask%map(j,1).lt. 0.5)cycle
                cost = dot_product(vec(:, i) , vec(:, j))
                loc = nint((cost - xmin)/dx+(0.5d0+1.d-10))
                if(loc .ge. 1 .and. loc.le.n)then
                   weight(loc) = weight(loc) + 1.d0
                   corr(loc) = corr(loc) + this%map(i,1)*this%map(j,1)
                endif
             enddo
             !!j = i
             if(xmax .eq. 1.d0)then
                weight(n) = weight(n) + 1
                corr(n) = corr(n) + this%map(i,1)**2
             endif
          enddo
       else
          do i=0, this%npix-1
             do j=1, i-1
                cost = dot_product(vec(:, i) , vec(:, j))
                loc = nint((cost - xmin)/dx+(0.5d0+1.d-10))
                if(loc .ge. 1 .and. loc.le.n)then
                   weight(loc) = weight(loc) + 1.d0
                   corr(loc) = corr(loc) + this%map(i,1)*this%map(j,1)
                endif
             enddo
             !!j = i
             if(xmax .eq. 1.d0)then
                weight(n) = weight(n) + 1
                corr(n) = corr(n) + this%map(i,1)**2
             endif
          enddo
       endif
    else
       
    endif
    where (weight .gt. 0.5d0)
       corr = corr/weight
    elsewhere
       corr = 0.d0
    end where
  end subroutine coop_healpix_maps_eval_corr


  subroutine coop_healpix_maps_copy_genre(from, to)
    class(coop_healpix_maps)::from, to
    COOP_INT::i
    if(to%nmaps .ne. from%nmaps) stop "copy genre only works for same nmaps"
    do i=1, min(from%nmaps, to%nmaps)
       call to%set_field(i, trim(adjustl(from%fields(i))))
       call to%set_unit(i, trim(adjustl(from%units(i))))
    enddo    
    to%spin = from%spin
    to%polar = from%polar
  end subroutine coop_healpix_maps_copy_genre

  subroutine coop_healpix_maps_apply_mask(this, mask, remove_monopole, bad_data)
    class(coop_healpix_maps)::this
    type(coop_healpix_maps)::mask
    logical,optional::remove_monopole
    logical,optional::bad_data
    COOP_INT::i
    COOP_REAL::summask
    if(mask%nside .ne. this%nside) stop "apply_mask:mask and map must have the same nside"
    if(mask%ordering .ne. this%ordering)then
       call mask%convert2nested()
       call this%convert2nested()
    endif
    if(present(remove_monopole))then
       if(remove_monopole)then
          summask =  sum(dble(mask%map(:,1)))
       else
          summask = 0.d0
       endif
    else
       summask = 0.d0
    endif
    do i=1, this%nmaps
       this%map(:,i) = this%map(:,i)*mask%map(:,1)
       if(summask .ge. 1.d0)then
          this%map(:, i) = this%map(:, i) - sum(dble(this%map(:, i)))/summask
       endif
    enddo
    if(present(bad_data))then
       if(bad_data)then
          do i=1, this%nmaps
             where(mask%map(:,1) .lt. 0.5)
                this%map(:, i) = this%bad_data
             end where
          enddo
       endif
    endif
  end subroutine coop_healpix_maps_apply_mask

  subroutine coop_healpix_correlation_function_free(this)
    class( coop_healpix_correlation_function ):: this
    if(allocated(this%als))deallocate(this%als)
    if(allocated(this%bls))deallocate(this%bls)
    if(allocated(this%cls))deallocate(this%cls)
    if(allocated(this%cls_smooth))deallocate(this%cls_smooth)    
    this%lmax = -1
    this%fwhm = 0.d0
  end subroutine coop_healpix_correlation_function_free

  subroutine coop_healpix_correlation_function_init(this, lmax, Cls)
    class( coop_healpix_correlation_function ):: this
    COOP_INT::lmax
    COOP_REAL,optional::Cls(0:lmax)
    if(this%lmax .ne. lmax)then
       call this%free()
       this%lmax = lmax
       allocate(this%als(0:lmax), this%bls(0:lmax), this%cls(0:lmax), this%cls_smooth(0:lmax))
       call coop_sphere_correlation_init(lmax, this%als, this%bls)
    endif
    if(present(Cls))then
       this%cls = cls
       this%cls_smooth = cls
    endif
  end subroutine coop_healpix_correlation_function_init

  subroutine  coop_healpix_correlation_function_simulate(this, fwhm)
    class(coop_healpix_correlation_function):: this        
    COOP_REAL::fwhm, c
    COOP_INT::l
    this%fwhm = fwhm
    c  = -sign((coop_sigma_by_fwhm * fwhm)**2, fwhm)
    do l = 0, this%lmax
       this%Cls_smooth(l) = this%Cls(l) * sum(coop_random_Gaussian_vector(2*l+1))**2/(2*l+1)* exp(c*l*(l+1.d0))
    enddo
  end subroutine coop_healpix_correlation_function_simulate

  subroutine  coop_healpix_correlation_function_set_beam(this, fwhm)
    class(coop_healpix_correlation_function):: this    
    COOP_REAL::fwhm, c
    COOP_INT::l
    this%fwhm = fwhm
    c = -sign((coop_sigma_by_fwhm * fwhm)**2, fwhm)
    !$omp parallel do
    do l = 0, this%lmax
       this%Cls_smooth(l) = this%Cls(l) * exp(c*l*(l+1.d0))
    enddo
    !$omp end parallel do
  end subroutine coop_healpix_correlation_function_set_beam

  function coop_healpix_correlation_function_corr(this, x) result(c)
    class( coop_healpix_correlation_function ):: this        
    COOP_REAL::x, c
    if(this%fwhm .ne. 0.d0)then
       c = coop_sphere_correlation(this%lmax, this%cls_smooth, this%als, this%bls, x)
    else
       c = coop_sphere_correlation(this%lmax, this%cls, this%als, this%bls, x)       
    endif
  end function coop_healpix_correlation_function_corr

  
  function coop_healpix_correlation_function_corr2int(this, xmin, xmax, mask_corr) result(S)
    class( coop_healpix_correlation_function ):: this
    type( coop_healpix_correlation_function ), optional::mask_corr
    
    COOP_REAL::xmin, xmax, S
    S = coop_integrate(c2, max(xmin, -1.d0), min(xmax, 1.d0), 1.d-5)
  contains
    function c2(x)
      COOP_REAL::x, c2
      if(present(mask_corr))then
         c2 = (this%corr(x)/mask_corr%corr(x))**2         
      else
         c2 = this%corr(x)**2
      endif
    end function c2
  end function coop_healpix_correlation_function_corr2int


  subroutine coop_healpix_mask_reverse(fin, fout)
    COOP_UNKNOWN_STRING::fin, fout
    type(coop_healpix_maps)::mask
    call mask%read(fin)
    mask%map =  1. - mask%map
    call mask%write(fout)
    call mask%free()
  end subroutine coop_healpix_mask_reverse

  subroutine coop_healpix_nside2lmax(nside,lmax)
    COOP_INT::nside, lmax
    lmax = floor(nside*coop_healpix_lmax_by_nside)
  end subroutine coop_healpix_nside2lmax



  subroutine coop_healpix_maps_stack_on_filament(this, disc, filament, tmp_filament, mask)
    class(coop_healpix_maps) this
    type(coop_healpix_disc) disc
    type(coop_healpix_maps),optional::mask
    type(coop_healpix_filament) filament, tmp_filament
    if(present(mask))then
       call this%fetch_filament(disc, tmp_filament, mask)
       if(sum(tmp_filament%nstack*tmp_filament%indisk) .lt. filament%num_indisk_tol)return
    else
       call this%fetch_filament(disc, tmp_filament)
    endif
    filament%image = filament%image + tmp_filament%image
    filament%nstack = filament%nstack + tmp_filament%nstack
    filament%nstack_raw = filament%nstack_raw + tmp_filament%nstack_raw
  end subroutine coop_healpix_maps_stack_on_filament  

  subroutine coop_healpix_maps_fetch_filament(this, disc, filament, mask)
    COOP_REAL,parameter::b_min = 0.1
    COOP_REAL,parameter::e2_min = 0.1
    COOP_REAL,parameter::e2_max = 3.
    class(coop_healpix_maps)::this
    type(coop_healpix_disc) disc
    type(coop_healpix_maps),optional::mask
    type(coop_healpix_filament) filament
    COOP_INT i, j, pix,  k, pix1, pix2
    COOP_REAL x, y, x1, y1, x2, y2, cos1, sin1, cos2, sin2, phi1, phi2, p1, p2, phi1_new, phi2_new, const1, const2, fpeak
    logical::left, right
    COOP_SINGLE qu(2)
    filament%image = 0.
    filament%nstack = 1.d0
    filament%nstack_raw  = 1
    pix1 = disc%center
    x1 = 0.d0
    y1 = 0.d0
    x2 = 0.d0
    y2 = 0.d0
    if(this%map(pix1, 2)**2 + this%map(pix1, 3)**2 .le. e2_min * this%map(pix1, 4)**2 .or. this%map(pix1, 2)**2 + this%map(pix1, 3)**2 .gt. e2_max * this%map(pix1,4)**2)then
       return
    endif
    fpeak = abs(this%map(pix1, 1))*b_min
    phi1 = COOP_POLAR_ANGLE(this%map(pix1,2), this%map(pix1,3))/2.d0
    if(coop_random_unit().gt.0.5)then
       const1 = 0.d0
       const2 = coop_pi
    else
       const1 = coop_pi
       const2 = 0.d0
    endif
    if(this%map(pix1, 1) .lt. 0.) then
       const1 = const1 + coop_pio2
       const2 = const2 + coop_pio2
    endif
    phi2 = phi1 + const2
    phi1 = phi1 + const1
    cos1 = cos(phi1)
    sin1 = sin(phi1)
    cos2 = - cos1
    sin2 = - sin1
    do j=-filament%ny, filament%ny
       x = x2 - (filament%dr*j)*sin2
       y = y2 + (filament%dr*j)*cos2
       call disc%xy2pix(x, y, pix)
       filament%image(0, j, 1) = this%map(pix, 1)
       filament%nstack(0, j) = 1.d0
    enddo
    left = .true.
    right = .true.
    do i=1, filament%n
       if(left)then
          x1 = x1 + filament%dr*cos1
          y1 = y1 + filament%dr*sin1
          call disc%xy2pix(x1, y1, pix1)
          phi1_new = COOP_POLAR_ANGLE(this%map(pix1,2), this%map(pix1,3))/2.d0 + const1
          if( abs(this%map(pix1, 1)).lt. fpeak )then
             left = .false.
          else
             phi1 = phi1_new
             cos1 = cos(phi1)
             sin1 = sin(phi1)
             do j=-filament%ny, filament%ny
                x = x1 + (filament%dr*j)*sin1
                y = y1 - (filament%dr*j)*cos1
                call disc%xy2pix(x, y, pix)
                filament%image(-i, j, 1) = this%map(pix, 1)
                filament%nstack(-i, j) = 1.d0
             enddo
          endif
       endif
       if(right)then
          x2 = x2 + filament%dr*cos2
          y2 = y2 + filament%dr*sin2       
          call disc%xy2pix(x2, y2, pix2)
          phi2_new = COOP_POLAR_ANGLE(this%map(pix2,2), this%map(pix2,3))/2.d0 + const2
          if( abs(this%map(pix2, 1)).lt. fpeak)then
             right = .false.
          else
             phi2 = phi2_new
             cos2 = cos(phi2)
             sin2 = sin(phi2)
             do j=-filament%ny, filament%ny
                x = x2 - (filament%dr*j)*sin2
                y = y2 + (filament%dr*j)*cos2
                call disc%xy2pix(x, y, pix)
                filament%image(i, j, 1) = this%map(pix, 1)
                filament%nstack(i, j) = 1.d0
             enddo
          endif
       endif
    enddo
  end subroutine coop_healpix_maps_fetch_filament


  subroutine coop_healpix_filament_plot(this, imap, output, use_degree)
    COOP_INT::bgrids
    class(coop_healpix_filament)::this
    COOP_INT imap
    COOP_UNKNOWN_STRING::output
    type(coop_asy)::fig
    COOP_INT nb, i, j, k, ns, iq, iu
    COOP_REAL  xc, yc,  norm, r, theta, minz, maxz
    COOP_REAL,dimension(:),allocatable::xstart, xend, ystart, yend
    logical,optional::use_degree
    logical use_rad
    COOP_SHORT_STRING::xlabel, ylabel
    if(imap .le. 0 .or. imap .gt. this%nmaps)then
       write(*,*) this%nmaps, imap
       stop "coop_healpix_filament_plot: imap overflow"
    endif
    call fig%open(output)
    if(present(use_degree))then
       if(use_degree)then
          this%r = this%r/coop_SI_degree
       endif
       use_rad = .not. use_degree
    else
       use_rad = .true.
    endif
    if(use_rad)then
       xlabel = "$\varpi\cos\phi$"
       ylabel = "$\varpi\sin\phi$"       
    else
       xlabel = "$\varpi\cos\phi (\mathrm{deg})$"
       ylabel = "$\varpi\sin\phi (\mathrm{deg})$"
    endif
    this%caption = adjustl(this%caption)
    if(coop_healpix_patch_default_want_caption)then
       if(len_trim(this%caption) .gt. 60)then
          call fig%init(caption = "{\tiny "//trim(this%caption)//"}", xlabel =trim(xlabel), ylabel =trim(ylabel), width = coop_healpix_patch_default_figure_width, height = coop_healpix_patch_default_figure_height, xmin = -real(this%r(this%n)), xmax = real(this%r(this%n)), ymin = -real(this%r(this%ny)), ymax = real(this%r(this%ny)))
       elseif(len_trim(this%caption) .gt. 50)then
          call fig%init(caption = "{\scriptsize "//trim(this%caption)//"}", xlabel =trim(xlabel), ylabel =trim(ylabel), width = coop_healpix_patch_default_figure_width, height = coop_healpix_patch_default_figure_height, xmin = -real(this%r(this%n)), xmax = real(this%r(this%n)), ymin = -real(this%r(this%ny)), ymax = real(this%r(this%ny)))
       elseif(len_trim(this%caption) .gt. 40)then
          call fig%init(caption = "{\small "//trim(this%caption)//"}", xlabel =trim(xlabel), ylabel =trim(ylabel), width = coop_healpix_patch_default_figure_width, height = coop_healpix_patch_default_figure_height, xmin = -real(this%r(this%n)), xmax = real(this%r(this%n)), ymin = -real(this%r(this%ny)), ymax = real(this%r(this%ny)))       
       else
          call fig%init(caption = trim(this%caption), xlabel =trim(xlabel), ylabel =trim(ylabel), width = coop_healpix_patch_default_figure_width, height = coop_healpix_patch_default_figure_height, xmin = -real(this%r(this%n)), xmax = real(this%r(this%n)), ymin = -real(this%r(this%ny)), ymax = real(this%r(this%ny)))
       endif
    else
       call fig%init(xlabel =trim(xlabel), ylabel =trim(ylabel), width = coop_healpix_patch_default_figure_width, height = coop_healpix_patch_default_figure_height, xmin = -real(this%r(this%n)), xmax = real(this%r(this%n)), ymin = -real(this%r(this%ny)), ymax = real(this%r(this%ny)))          
    endif
    if(this%tbs%zmin(imap) .lt.0.99e30)then
       minz = this%tbs%zmin(imap)
    else
       call coop_array_get_threshold(this%image(:,:,imap), COOP_REAL_OF(0.99), minz)
    endif
    if(this%tbs%zmax(imap) .gt. -0.99e30)then
       maxz = this%tbs%zmax(imap)
    else
       call coop_array_get_threshold(this%image(:,:,imap), COOP_REAL_OF(0.01), maxz)
    endif
    if(coop_healpix_patch_default_want_label)then
       call coop_asy_density(fig, this%image(:,:,imap), -this%r(this%n), this%r(this%n), -this%r(this%ny), this%r(this%ny), label = trim(this%tbs%label(imap)), zmax = maxz, zmin = minz, color_table = trim(this%color_table))
    else
       call coop_asy_density(fig, this%image(:,:,imap), -this%r(this%n), this%r(this%n), -this%r(this%ny), this%r(this%ny), zmax = maxz, zmin = minz, color_table = trim(this%color_table))       
    endif
    if(use_rad .and. coop_healpix_patch_default_want_arrow)then
       theta = nint(2.d0*asin(this%r(this%n)/2.d0)/coop_SI_degree*10.d0)/10.d0
       call coop_asy_label(fig, "$\mathbf{-"//COOP_STR_OF(theta)//"}^\circ$", -this%r(this%n)*1.01, -this%r(this%n)*1.23, color="blue")
       call coop_asy_label(fig, "$\mathbf{"//COOP_STR_OF(theta)//"}^\circ$", this%r(this%n)*1.01, -this%r(this%n)*1.23, color="blue")
       call fig%arrow(this%r(this%n),  -this%r(this%ny)*1.17, this%r(this%n),  -this%r(this%ny)*1.13)
       call fig%arrow(-this%r(this%n),  -this%r(this%ny)*1.17, -this%r(this%n),  -this%r(this%ny)*1.13)
    endif
100 call fig%close()
    if(present(use_degree))then
       if(use_degree)then
          this%r = this%r*coop_SI_degree
       endif
    endif    
  end subroutine coop_healpix_filament_plot


  subroutine coop_healpix_maps_stack_filaments_on_peaks(this, sto, filament, mask, norm)
    COOP_INT,parameter::n_threads = 8
    class(coop_healpix_maps)::this
    type(coop_stacking_options)::sto
    type(coop_healpix_disc),dimension(n_threads)::disc
    type(coop_healpix_filament)::filament
    type(coop_healpix_filament),dimension(n_threads)::p, tmp
    type(coop_healpix_maps),optional::mask
    logical,optional::norm
    COOP_INT ithread, i
#ifdef HAS_HEALPIX
    filament%image = 0.d0
    filament%nstack = 0.d0
    filament%nstack_raw = 0
    if(sto%nested)then
       call this%convert2nested()
       if(present(mask))call mask%convert2nested()
    else
       call this%convert2ring()
       if(present(mask))call mask%convert2ring()
    endif
    !!adjust the index of maps
    i = 1
    do while(i.le.filament%nmaps)
       if(filament%tbs%ind(i) .gt. this%nmaps) stop "stack_on_peaks: index overflow"       
       do while(filament%tbs%spin(i) .ne. this%spin(filament%tbs%ind(i)))
          filament%tbs%ind(i)  = filament%tbs%ind(i) + 1
          if(filament%tbs%ind(i) .gt. this%nmaps) stop "stack_on_peaks: cannot find the map with corresponding spin"
       enddo
       if(filament%tbs%spin(i) .ne. 0)then
          if(i.ge.filament%nmaps .or. filament%tbs%ind(i) .ge. this%nmaps) stop "stack_on_peaks: nonzero spin must go in pairs"
          filament%tbs%ind(i+1) = filament%tbs%ind(i) + 1
          i = i+2
       else
          i = i+1
       endif
    enddo
    do ithread=1, n_threads
       p(ithread) = filament
       tmp(ithread) = filament
    enddo
    !$omp parallel do private(i, ithread)
    do ithread = 1, n_threads
       do i=ithread, sto%peak_pix%n, n_threads
          call this%get_disc(sto%pix(this%nside, i), disc(ithread))
          if(present(mask))then
             call this%stack_on_filament( disc(ithread),  p(ithread), tmp(ithread), mask)    
          else
             call this%stack_on_filament( disc(ithread),  p(ithread), tmp(ithread))
          endif
       enddo
    enddo
    !$omp end parallel do

    do ithread = 1, n_threads
       filament%image = filament%image + p(ithread)%image
       filament%nstack = filament%nstack + p(ithread)%nstack
       filament%nstack_raw = filament%nstack_raw + p(ithread)%nstack_raw
       call p(ithread)%free()
       call tmp(ithread)%free()
    enddo
    if(present(norm))then
       if(.not.norm)return
    endif
    if(filament%nstack_raw .ne. 0)then
       do i=1, filament%nmaps
          filament%image(:, :, i) = filament%image(:, :, i)/max(filament%nstack, 1.d0)
       enddo
    else
       write(*,*) "warning: no filamentes has been found"
       filament%image = 0.d0
    endif
#else
    stop "CANNOT FIND HEALPIX"
#endif
  end subroutine coop_healpix_maps_stack_filaments_on_peaks


  subroutine coop_healpix_maps_perimeter_area_list(this, pix, r, palist, countcut, peakcut, rmscut, imap, plot, sum_pa, want_abs_area, plot_min_perimeter, area_truncation)
    class(coop_healpix_maps)::this
    COOP_INT::pix
    COOP_SINGLE::r, threshold
    COOP_SINGLE,optional::peakcut, countcut, rmscut
    COOP_INT,optional::imap
    COOP_INT::id, n, i
    type(coop_asy_path) path
    type(coop_healpix_disc)::disc
    type(coop_list_realarr)::palist
    COOP_SINGLE::pa(2), plot_p
    COOP_UNKNOWN_STRING, optional::plot
    COOP_SINGLE, optional::plot_min_perimeter
    type(coop_asy)::fig
    logical::doplot, wabs
    logical,optional::want_abs_area
    COOP_REAL::smallest_area
    COOP_SINGLE, optional::area_truncation
    COOP_SINGLE,optional::sum_pa(2)
    if(present(imap))then
       if(imap .gt. this%nmaps) stop "filament_width: imap overflow"
       id = imap
    else
       id = 1
    endif
    if(present(area_truncation))then
       smallest_area  = area_truncation
    else
       smallest_area = coop_SI_arcmin**2
    endif
    if(present(want_abs_area))then
       wabs = want_abs_area
    else
       wabs = .false.
    endif    
    if(pix .lt. 0 .or. pix .ge. this%npix) stop "filament_width: pix overflow"
    call this%get_disc(pix, disc)
    n = floor(r*this%nside)
    if(present(countcut))then
       call path%from_function(fxy, -r, r, -r, r, threshold, n, countcut = countcut)
    elseif(present(peakcut))then
       call path%from_function(fxy, -r, r, -r, r, threshold, n, peakcut = peakcut)
    elseif(present(rmscut))then
       call path%from_function(fxy, -r, r, -r, r, threshold, n, rmscut = rmscut)       
    else
       call path%from_function(fxy, -r, r, -r, r, threshold, n,  rmscut = 1.)       
    endif
    doplot = .false.
    if(present(plot_min_perimeter))then
       plot_p = plot_min_perimeter
    else
       plot_p = 4.d0*r
    endif
    if(present(sum_pa))then
       sum_pa = 0.
       do i=1, path%nclosed
          call path%get_perimeter_and_area(i, pa(1), pa(2))
          if(abs(pa(2)) .gt. smallest_area)then
             if(wabs) pa(2) = abs(pa(2))
             call palist%push(pa)
          endif
          if(pa(1) .gt. plot_p)doplot = .true.
          sum_pa = sum_pa + pa
       enddo
    else
       do i=1, path%nclosed
          call path%get_perimeter_and_area(i, pa(1), pa(2))
          if(abs(pa(2)) .gt. smallest_area)then
             if(wabs) pa(2) = abs(pa(2))
             call palist%push(pa)
          endif
          if(pa(1) .gt. plot_p)doplot = .true.
       enddo
    endif
    if(present(plot) .and. doplot)then
       call fig%open(plot)
       call fig%init(xlabel = "$x$", ylabel = "$y$", xmin = -r, xmax = r, ymin = -r, ymax = r)
       call coop_asy_contour(fig, path, colorfill = "skyblue", smooth = .false., linecolor = "blue", linetype = "solid", linewidth = 0.5)
       call fig%close()
    endif
    call path%free()    
  contains

    function fxy(x, y)
      COOP_SINGLE::x, y, fxy
      COOP_INT::ip
      call disc%xy2pix(dble(x), dble(y), ip)
      fxy = this%map(ip, id)
    end function fxy
    
  end subroutine coop_healpix_maps_perimeter_area_list



  subroutine coop_healpix_maps_zeros(this, imap, zeros, mask)
    COOP_REAL, parameter::rms_cut = 0.05
    class(coop_healpix_maps)::this
    type(coop_healpix_maps)::zeros
    COOP_INT::imap, l
    type(coop_healpix_maps),optional::mask
    COOP_INT i, k, nn(8)
    COOP_REAL::cut
#ifdef HAS_HEALPIX

    call this%convert2nested()
    if(zeros%nside .ne. this%nside)then
       call zeros%init(nside = this%nside, nmaps = 1, genre = "MASK")
       zeros%ordering = this%ordering
    elseif(zeros%nmaps .ne. 1)then
       call zeros%convert2nested()
    else
       zeros%ordering = COOP_NESTED
    endif
    zeros%map(:,1) = 0.
    if(present(mask))then

       if(mask%nside .ne. this%nside) stop "zeros: mask nside and map nside must be the same"
       call mask%convert2nested()
       cut = sqrt(sum(this%map(:, imap)**2*mask%map(:,1))/max(sum(mask%map(:,1)), 1.)) * rms_cut      
       !$omp parallel do private(nn, k, i)
       do i = 0, this%npix-1
          if(mask%map(i, 1) .gt. 0.5)then
             if(abs(this%map(i, imap)) .lt. cut)then
                call neighbours_nest(this%nside, i, nn, k)
                if(all(mask%map(nn(1:k), 1) .gt. 0.5))then
                   if( any( this%map(nn(1:k), imap) * this%map(i, imap) .le. 0. ) .and. abs(this%map(i, imap)) .lt. cut)then
                      zeros%map(i, 1) = 1.
                   endif
                endif
             endif
          endif
       enddo
       !$omp end parallel do
    else
       cut = sqrt(sum(this%map(:, imap)**2)/real(this%npix)) * rms_cut       
       !$omp parallel do private(nn, k, i)       
       do i = 0, this%npix-1
          if( abs(this%map(i, imap)) .lt. cut)then
             call neighbours_nest(this%nside, i, nn, k)
             if( any( this%map(nn(1:k), imap) * this%map(i, imap) .le. 0. ))then
                zeros%map(i, 1) = 1.
             endif
          endif
       enddo
       !$omp end parallel do
    endif
#endif    
  end subroutine coop_healpix_maps_zeros

  subroutine coop_healpix_maps_get_dervs(this, imap, dervs, alms_done)
    class(coop_healpix_maps)::this
    type(coop_healpix_maps)::dervs
    logical, optional::alms_done
    COOP_INT::imap
#ifdef HAS_HEALPIX
    if(imap.gt. this%nmaps) stop "get_dervs: imap over flow"
    if(dervs%nmaps .ne. 6 .or. dervs%nside .ne. this%nside) &
         call dervs%init(nside = this%nside, nmaps = 6, genre = "DERVS")
    if(.not. present(alms_done))then
       call this%map2alm(index_list = (/ imap /) )
    elseif(.not. alms_done)then
       call this%map2alm(index_list = (/ imap /) )       
    endif
    call alm2map_der(nsmax = this%nside, nlmax = this%lmax, nmmax = this%lmax, alm = reshape(this%alm(0:this%lmax, 0:this%lmax, imap), (/ 1, this%lmax+1, this%lmax+1 /) ), map = dervs%map(:, 1), der1 = dervs%map(:, 5:6), der2 = dervs%map(:, 2:4))
    call dervs%set_units(this%units(imap))
#endif    
  end subroutine coop_healpix_maps_get_dervs


  subroutine coop_healpix_maps_local_disk_minkowski0(this, pix, r_deg, nu, imap, mean, rms, V0, gaussian_fit)
    class(coop_healpix_maps)::this
    COOP_INT::pix, imap, i
    COOP_REAL::r_deg
    COOP_REAL::nu(:)
    COOP_SINGLE:: V0(:), rms, mean
    COOP_REAL::xbar, xrms, Amp
    COOP_INT::listpix(0:this%npix-1)
    COOP_INT::nlist
    logical::gaussian_fit
    if(size(V0) .ne. size(nu)) stop "local_disk_minkowsk: nu and V0 must have the same size"
    call this%query_disc(pix, r_deg, listpix, nlist)
    if(gaussian_fit)then
       call coop_fit_gaussian(dble(this%map(listpix(0:nlist-1), imap)), max(nlist/200, 20), xbar, xrms, Amp)
       mean = xbar
       rms = xrms
    else
       mean = sum(this%map(listpix(0:nlist-1), imap))/nlist
       rms = sqrt(sum((this%map(listpix(0:nlist-1), imap) - mean)**2)/nlist)
    endif
    do i = 1, size(V0)
       V0(i) = count(this%map(listpix(0:nlist-1), imap) .ge. mean + nu(i)*rms)/real(nlist)
    enddo
  end subroutine coop_healpix_maps_local_disk_minkowski0


  subroutine coop_healpix_maps_scan_local_minkowski0(this, imap, nu, meanmap, rmsmap, V0map, r_deg, do_gaussian_fit)
    class(coop_healpix_maps)::this
    COOP_INT::imap
    COOP_REAL::nu(:)
    type(coop_healpix_maps)::meanmap, rmsmap, V0map
    COOP_REAL, optional::r_deg
    COOP_REAL::mean, rms, Amp
    COOP_INT::ipix, space, i
    COOP_REAL::theta, phi
    logical,optional::do_gaussian_fit
    logical::gf
    if(present(do_gaussian_fit))then
       gf = do_gaussian_fit
    else
       gf = .false.
    endif
    if(V0map%nside .ge. this%nside .or. meanmap%nside .ne. rmsmap%nside .or. meanmap%nside .ne. V0map%nside .or. V0map%nmaps.ne.size(nu)) stop "wrong nside for scan_local_minkowski"
    call meanmap%convert2nested()
    call rmsmap%convert2nested()
    call V0map%convert2nested()
    call this%convert2nested()
    if(present(r_deg))then
       if(r_deg .gt. 0.d0)then
          do ipix=0, V0map%npix-1
             call V0map%pix2ang(ipix, theta, phi)
             call this%ang2pix(theta, phi, i)
             call this%local_disk_minkowski0(i, r_deg*coop_SI_degree, nu, imap, meanmap%map(ipix, 1), rmsmap%map(ipix,1), V0map%map(ipix,1:V0map%nmaps), gf)
          enddo
          return
       endif
    endif
    space = (this%nside/V0map%nside)**2       
    do ipix=0, V0map%npix-1
       if(gf)then
          call coop_fit_gaussian(dble(this%map(ipix*space:(ipix+1)*space-1, imap)), max(space/200, 20), mean, rms, Amp)
          meanmap%map(ipix, 1) = mean
          rmsmap%map(ipix,1) = rms
       else
          meanmap%map(ipix,1) = sum(this%map(ipix*space:(ipix+1)*space-1, imap))/real(space)
          rmsmap%map(ipix,1) = sqrt(sum( (this%map(ipix*space:(ipix+1)*space-1, imap) - meanmap%map(ipix,1))**2 )/real(space))
       endif
       do i = 1, size(nu)
          V0map%map(ipix,i) = count(this%map(ipix*space:(ipix+1)*space-1, imap) .ge.  meanmap%map(ipix,1) + nu(i)*rmsmap%map(ipix,1))/real(space)
       enddo
    enddo
  end subroutine coop_healpix_maps_scan_local_minkowski0
  
  subroutine coop_healpix_maps_local_disk_minkowski1(this, imap, sourcemap, grad2map, pix, r_deg, nu, mean, rms, sigma1, V1, gaussian_fit)
    class(coop_healpix_maps)::this
    COOP_INT::pix, imap, i, sourcemap, grad2map
    COOP_REAL::r_deg
    COOP_REAL::nu(:)
    COOP_SINGLE:: V1(:), rms, mean, sigma1
    COOP_REAL::xmean, xrms, Amp
    COOP_INT::listpix(0:this%npix-1)
    COOP_INT::nlist
    logical::gaussian_fit
    if(size(V1) .ne. size(nu)) stop "local_disk_minkowsk: nu and V1 must have the same size"
    call this%query_disc(pix, r_deg, listpix, nlist)
    if(gaussian_fit)then
       call coop_fit_gaussian(dble(this%map(listpix(0:nlist-1), imap)), max(nlist/200, 20), xmean, xrms, Amp)
       mean = xmean
       rms = xrms
    else
       mean = sum(this%map(listpix(0:nlist-1), imap))/nlist
       rms = sqrt(sum((this%map(listpix(0:nlist-1), imap) - mean)**2)/nlist)
    endif
    sigma1 = sqrt(sum(this%map(listpix(0:nlist-1), grad2map))/nlist)
    do i = 1, size(V1)
       V1(i) = sum(this%map(listpix(0:nlist-1), sourcemap), mask  = this%map(listpix(0:nlist-1), imap) .ge. mean + nu(i)*rms)/nlist/sigma1*rms
    enddo
  end subroutine coop_healpix_maps_local_disk_minkowski1


  subroutine coop_healpix_maps_scan_local_minkowski1(this, imap, sourcemap, grad2map, nu, meanmap, rmsmap, sigma1map, V1map, r_deg, do_gaussian_fit)
    class(coop_healpix_maps)::this
    COOP_INT::imap, sourcemap, grad2map
    COOP_REAL::nu(:)
    type(coop_healpix_maps)::meanmap, rmsmap, V1map, sigma1map
    COOP_REAL, optional::r_deg
    COOP_INT::ipix, space, i
    COOP_REAL::theta, phi
    COOP_REAL::xmean, xrms, Amp
    logical,optional::do_gaussian_fit
    logical::gf
    if(present(do_gaussian_fit))then
       gf = do_gaussian_fit
    else
       gf = .false.
    endif    
    if(V1map%nside .ge. this%nside .or. meanmap%nside .ne. rmsmap%nside .or. meanmap%nside .ne. V1map%nside .or. V1map%nmaps.ne.size(nu)) stop "wrong nside for scan_local_minkowski"
    call meanmap%convert2nested()
    call rmsmap%convert2nested()
    call V1map%convert2nested()
    call this%convert2nested()
    if(present(r_deg))then
       if(r_deg .gt. 0.d0)then
          do ipix=0, V1map%npix-1
             call V1map%pix2ang(ipix, theta, phi)
             call this%ang2pix(theta, phi, i)
             call this%local_disk_minkowski1(imap, sourcemap, grad2map, i, r_deg*coop_SI_degree, nu, meanmap%map(ipix, 1), rmsmap%map(ipix,1), sigma1map%map(ipix,1), V1map%map(ipix,1:V1map%nmaps), gf)
          enddo
          return
       endif
    endif
    space = (this%nside/V1map%nside)**2       
    do ipix=0, V1map%npix-1
       if(gf)then
          call coop_fit_Gaussian(dble(this%map(ipix*space:(ipix+1)*space-1, imap)), max(space/200, 20), xmean, xrms, Amp)
          meanmap%map(ipix,1) = xmean
          rmsmap%map(ipix, 1) = xrms
       else
          meanmap%map(ipix,1) = sum(this%map(ipix*space:(ipix+1)*space-1, imap))/real(space)
          rmsmap%map(ipix,1) = sqrt(sum( (this%map(ipix*space:(ipix+1)*space-1, imap) - meanmap%map(ipix,1))**2 )/real(space))
       endif
       sigma1map%map(ipix,1) = sqrt(sum(this%map(ipix*space:(ipix+1)*space-1, grad2map))/real(space))
       do i = 1, size(nu)
          V1map%map(ipix,i) =  sum(this%map(ipix*space:(ipix+1)*space-1, sourcemap), mask = this%map(ipix*space:(ipix+1)*space-1, imap) .ge.  meanmap%map(ipix,1) + nu(i)*rmsmap%map(ipix,1))/space/sigma1map%map(ipix,1)*rmsmap%map(ipix,1)
       enddo
    enddo
  end subroutine coop_healpix_maps_scan_local_minkowski1

  subroutine coop_healpix_maps_reflect(this)
    class(coop_healpix_maps)::this
    COOP_INT::inds(0:this%npix-1)
    COOP_INT::i
    COOP_REAL::theta, phi
    do i=0, this%npix-1
       call this%pix2ang(i, theta, phi)
       call this%ang2pix(coop_pi-theta, coop_pi+phi, inds(i))
    enddo
    this%map(:,:) = this%map(inds, :)
  end subroutine coop_healpix_maps_reflect


end module coop_healpix_mod









