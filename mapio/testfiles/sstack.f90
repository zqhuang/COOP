program massive_stack
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX
  logical::remove_mono = .false.
  logical::read_only = .false.
  COOP_UNKNOWN_STRING,parameter::cal_file_prefix = "rprof/"
  COOP_UNKNOWN_STRING,parameter::data_theory_cl_file = "planck14best_lensedCls.dat"
  COOP_UNKNOWN_STRING,parameter::sim_theory_cl_file = "planck13best_lensedCls.dat"  
  COOP_INT::show = 3
  COOP_STRING::cc_method = "smica"
  COOP_STRING::filter = "tophat"
  COOP_STRING::fr_genre = "cold0"
  COOP_STRING::polcase = "pol_case1"
  COOP_STRING::stack_field_name = "T"
  COOP_INT::resol = 1024
  COOP_INT::fwhm

  COOP_STRING::output,  cal_data_file, cal_sim_file
  COOP_REAL::threshold
  COOP_INT::n_sim 
  COOP_STRING,parameter::peak_name = "$T$"
  COOP_STRING::orient_name 
  COOP_INT,parameter::n = 36
  COOP_REAL,parameter::r_degree  = 2.d0
  COOP_REAL,parameter::dr = 2.d0*sin(r_degree*coop_SI_degree/2.d0)/n
  COOP_STRING::postfix
  COOP_REAL::r(0:n), pfr(0:n), pfr_theory_data(0:n), pfr_theory_sim(0:n), norm_theory_data(0:n), norm_theory_sim(0:n), Wfil(0:n), Smean, S_lower(3), S_upper(3)

  COOP_STRING::outputdir
  COOP_UNKNOWN_STRING,parameter::mapdir = "massffp8/"
  logical, parameter::do_nest = .true.
  logical::do_self = .false.
  COOP_STRING::imap_file, polmap_file, imask_file, polmask_file
  COOP_REAL,dimension(:),allocatable::S_m
  type(coop_stacking_options)::sto_max, sto_min
  type(coop_healpix_patch)::patch_max, patch_min, patch_max_data, patch_max_sim, patch_min_data, patch_min_sim
  type(coop_healpix_maps)::imap, imask, polmask, inoise, polnoise, polmap, imask_copy
  logical,dimension(:),allocatable::imask_logic
  logical::iloaded = .false.
  logical::polloaded  = .false.
  type(coop_file)::fp, fpcheck 
  type(coop_asy)::fig
  COOP_INT i, iredo, j
  COOP_INT ind, ind_done
  COOP_REAL sumimask, sumpolmask
  
  if(iargc() .lt. 7)then
     write(*,*) "Syntax:"
     write(*,*) "./SST cc_method resolution stack_field nu n_sim filter fr_type [Orient] [cut_cold_spot] [remove_mono] [ReadOnly]"
     write(*,*) "Examples:"     
     write(*,*) "./SST smica     1024 T  0.5 1000 tophat    cold0"
     write(*,*) "./SST smica     1024 T  0.5 1000 self    cold0"     
     write(*,*) "./SST nilc      512  T  0.  100  linear    cold1 T"
     write(*,*) "./SST commander 1024 QU 2.  500  gaussian  hot0  T F T"
     write(*,*) "./SST sevem     1024 T  1.  300  self    hc0  T S T T"
     stop
  endif
  cc_method = trim(coop_inputArgs(1))
  select case(trim(cc_method))
  case("smica","sevem","nilc")
     polcase= "pol_case1"
  case("commander")
     polcase = "pol_case5"
  case default
     write(*,*) "unknown cc_method = "//trim(cc_method)
     stop
  end select
  call coop_get_Input(2, resol)
  if(resol .ne. 1024 .and. resol .ne. 512) stop "only support resolution 512 and 1024"
  stack_field_name = trim(coop_inputArgs(3))
  call coop_get_Input(4, threshold)
  call coop_get_Input(5, n_sim)
  allocate(S_m(0:n_sim))
  filter = trim(coop_inputArgs(6))
  fr_genre =  trim(coop_inputArgs(7))
  if(trim(coop_inputArgs(8)).eq."T")then
     orient_name = "$(Q_T, U_T)$"
  else
     orient_name = "NULL"
  endif
  remove_mono = (trim(coop_inputArgs(10)) .eq. "T")
  read_only = (trim(coop_inputArgs(11)) .eq. "T")
  
  fwhm = 10240/resol
  outputdir = "st_"//trim(coop_str_numalpha(stack_field_name))//"_"//trim(coop_ndigits(resol,4))//"/"
  if(remove_mono) outputdir = "m"//trim(outputdir)  
  postfix = "_"//trim(coop_ndigits(fwhm,3))//"a_"//trim(coop_ndigits(resol,4))//".fits"
  
  imap_file = "planck14/dx11_v2_"//trim(cc_method)//"_int_cmb"//trim(postfix)
  polmap_file = "planck14/dx11_v2_"//trim(cc_method)//"_"//trim(polcase)//"_cmb_hp_20_40"//trim(postfix)
  

  select case(trim(coop_inputArgs(9)))
  case("N")
     imask_file = "planck14/HemAsym_north_int_mask"//trim(postfix)
     polmask_file = "planck14/HemAsym_north_pol_mask"//trim(postfix)          
     output = trim(outputdir)//trim(cc_method)//"_nu"//trim(coop_num2goodstr(threshold,"-","pt"))//"_"//trim(coop_str_numalpha(peak_name))//"_Orient"//trim(coop_str_numalpha(orient_name))//"_HAnorth.dat"                 
  case("S")
     imask_file = "planck14/HemAsym_south_int_mask"//trim(postfix)
     polmask_file = "planck14/HemAsym_south_pol_mask"//trim(postfix)          
     output = trim(outputdir)//trim(cc_method)//"_nu"//trim(coop_num2goodstr(threshold,"-","pt"))//"_"//trim(coop_str_numalpha(peak_name))//"_Orient"//trim(coop_str_numalpha(orient_name))//"_HAsouth.dat"                 
  case("T")
     imask_file = "planck14/ColdSpotCut_int_mask"//trim(postfix)
     polmask_file = "planck14/ColdSpotCut_pol_mask"//trim(postfix)          
     output = trim(outputdir)//trim(cc_method)//"_nu"//trim(coop_num2goodstr(threshold,"-","pt"))//"_"//trim(coop_str_numalpha(peak_name))//"_Orient"//trim(coop_str_numalpha(orient_name))//"_CutColdSpot.dat"            
  case default
     imask_file = "planck14/dx11_v2_common_int_mask"//trim(postfix)
     polmask_file = "planck14/dx11_v2_common_pol_mask"//trim(postfix)          
     output = trim(outputdir)//trim(cc_method)//"_nu"//trim(coop_num2goodstr(threshold,"-","pt"))//"_"//trim(coop_str_numalpha(peak_name))//"_Orient"//trim(coop_str_numalpha(orient_name))//".dat"
  end select

  call patch_max%init(trim(stack_field_name), n, dr)
  call patch_min%init(trim(stack_field_name), n, dr)  

  ind_done = -1  
  if(coop_file_exists(output))then
     call fp%open(output, "ru")
     do i = 0, n_sim
        read(fp%unit, ERR = 100, END = 100)  ind, patch_max%fr(0:patch_max%n, 0:patch_max%mmax/2, 1:patch_max%nmaps), patch_min%fr(0:patch_min%n, 0:patch_min%mmax/2, 1:patch_min%nmaps)
        ind_done = ind
     enddo
100  call fp%close()
  endif

  coop_healpix_mask_tol = 0.9d0

  if(trim(orient_name).eq."NULL")then
     call sto_max%init(.true., peak_name, orient_name, nmaps = 1)
     call sto_min%init(.false., peak_name, orient_name, nmaps = 1)     
  else
     call sto_max%init(.true., peak_name, orient_name, nmaps = 3)          
     call sto_min%init(.false., peak_name, orient_name, nmaps = 3)     
  endif
!  sto_max%addpi = .false.
!  sto_min%addpi = .false.
  sto_max%I_lower_nu = threshold
  sto_min%I_upper_nu = -threshold
  sto_max%threshold_option = 4
  sto_min%threshold_option = 4
  sto_max%nested = do_nest
  sto_min%nested = do_nest

  do i=0, n
     r(i) = i*dr
  enddo
  

  call load_theory()
  do_self = .false.
  select case(trim(filter))
  case("tophat","t")
     wfil = 1./(n+1.d0)
  case("linear","l")
     wfil = r
     wfil = wfil/sum(wfil)  
  case("gaussian","g")
     wfil = exp(-(r/coop_SI_degree)**2)
     wfil = wfil/sum(wfil)  
  case("self","s")
     wfil = norm_theory_data/sum(norm_theory_data**2)
     do_self = .true.
  case default
     stop "Unknown filter"
  end select
  


  if(ind_done .lt. n_sim)then 
     call imask%read(imask_file, nmaps_wanted = 1)
     allocate(imask_logic(0:imask%npix-1))
     imask_logic = (imask%map(0:imask%npix-1, 1).gt. 0.5)
     imask_copy = imask
     if(do_nest) call imask%convert2nested()
     sumimask =   count(imask_logic)

     if(trim(stack_field_name).eq."QU")then
        call polmask%read(polmask_file, nmaps_wanted = 1)
        if(do_nest) call polmask%convert2nested()
     endif
  endif
     
  if(read_only)then
     call fp%open(output, "ur")     
  else
    call fp%open(output, "u")
  endif
  if(read_only)then
     call fig%open("tmpfig.txt")
  else
     call fig%open(coop_str_replace(output, ".dat", ".txt"))
  endif
  call fig%init(xlabel = "$\varpi$", ylabel = "$\delta "//trim(stack_field_name)//"_0$")
  if(ind_done .ge. 0 .and. n_sim .gt. 0 .and. .not. read_only)then
     print*, "loaded "//COOP_STR_OF(ind_done+1)//" stacked maps"
  endif
  call fpcheck%open(trim(output)//".chk", "w")
  do ind = 0, n_sim
     if(ind.gt.ind_done)then
        if(read_only)then
           write(*,*) "Not enough maps "//COOP_STR_OF(ind_done)//"/"//COOP_STR_OF(n_sim)//"(program terminated in ReadOnly mode)"
           goto 500
        endif
        if(n_sim .gt. 0 )print*, "stacking map#"//COOP_STR_OF(ind)
        iloaded = .false.
        polloaded = .false.
        call find_peaks()
        call stack_map()
        call compute_fr()
        if(ind .eq. 0 )then
           write(fpcheck%unit, "(2I9,20E16.7)")  patch_max%nstack_raw, sto_max%peak_pix%n, sum(abs(patch_max%image))
           write(fpcheck%unit, "(2I9,20E16.7)")  patch_min%nstack_raw, sto_min%peak_pix%n, sum(abs(patch_min%image))
           call fpcheck%close()  
        endif
        write(fp%unit) ind, patch_max%fr(0:patch_max%n, 0:patch_max%mmax/2, 1:patch_max%nmaps), patch_min%fr(0:patch_min%n, 0:patch_min%mmax/2, 1:patch_min%nmaps)
     else
        read(fp%unit) i, patch_max%fr(0:patch_max%n, 0:patch_max%mmax/2, 1:patch_max%nmaps), patch_min%fr(0:patch_min%n, 0:patch_min%mmax/2, 1:patch_min%nmaps)
        if(i.ne.ind) stop "data file is broken"
     endif
     call get_radial_f(pfr, patch_max, patch_min)
     if(ind.eq.0 )then
        S_m(ind) = sum((pfr - pfr_theory_data)*wfil)
        if(do_self) wfil = norm_theory_sim/sum(norm_theory_sim**2)
        call fig%curve(r, pfr, color = "red", linetype= "solid", linewidth = 1.5)        
     else
        S_m(ind) = sum((pfr - pfr_theory_sim)*wfil)
        call fig%curve(r, pfr, color = "blue", linetype= "dotted", linewidth = 1.)
     endif
  enddo

  call fp%close()
  call fig%close()

#define FSTR(x) trim(coop_num2str(x, "(F13.2)"))  
  if(n_sim .ge. 1)then
     call coop_quicksort(S_m(1:n_sim))
     Smean = (S_m(n_sim/2+1)+S_m(n_sim/2))/2.d0
     S_lower(1) = S_m(floor(n_sim*0.1585d0 + 1.d0))  
     S_lower(2) = S_m(floor(n_sim*0.023d0 + 1.d0))
     S_lower(3) = S_m(floor(n_sim*0.0015 + 1.d0))
     S_upper(1) = S_m(ceiling(n_sim*0.8415d0))  
     S_upper(2) = S_m(ceiling(n_sim*0.977d0))
     S_upper(3) = S_m(ceiling(n_sim*0.9985d0))

     if(.not. read_only)then
        write(*,"(A12, F13.2, A12, F13.2)") "S_m median = ", Smean, " S_m data = ", S_m(0)
        do i = 1, 3
           write(*,"(I5,  G13.4, A9, G13.4)") i, S_lower(i), " < S_m < ", S_upper(i)
        enddo
     endif
     write(*,"(A)") "& $"//FSTR(S_m(0))//"("//FSTR(Smean)//"^{+"//FSTR(S_upper(1)-Smean)//"+"//FSTR(S_upper(2) - Smean)//"}_{-"//FSTR(Smean - S_lower(1))//"-"//FSTR(Smean - S_lower(2))//"})$"  
  else
     write(*,"(A, F10.2)")  "S_m = ", S_m(0)
  endif
  deallocate(S_m)
500  continue
  
contains

  subroutine find_peaks()
    select case(sto_max%genre)
    case(coop_stacking_genre_Imax, coop_stacking_genre_Imax_Oriented)
       call load_imap(ind)
    case default
       stop "so far SST only support temperature peaks (oriented or nonoriented)"
    end select
    call imap%get_peaks(sto_max, mask = imask, restore = .not. do_nest)
    call imap%get_peaks(sto_min, mask = imask, restore = .not. do_nest)
  end subroutine find_peaks

  subroutine compute_fr()
    call patch_max%get_all_radial_profiles()
    call patch_min%get_all_radial_profiles()
    patch_max%fr = patch_max%fr*1.d6
    patch_min%fr = patch_min%fr*1.d6    
  end subroutine compute_fr

  subroutine stack_map()
    select case(trim(stack_field_name))
    case("T")
       call load_imap(ind)
       call imap%stack_on_peaks(sto_max, patch_max, imask)
       call imap%stack_on_peaks(sto_min, patch_min, imask)
    case("QU")
       call load_polmap(ind)
       call polmap%stack_on_peaks(sto_max, patch_max, polmask)
       call polmap%stack_on_peaks(sto_min, patch_min, polmask)
    case default
       stop "so far SST only support T and QU stacking"
    end select
  end subroutine stack_map

  subroutine load_imap(i)
    COOP_INT i
    if(iloaded) return
    if(i.eq.0)then
       call imap%read(imap_file, nmaps_wanted = sto_max%nmaps, nmaps_to_read = 1)
    else
       call imap%read(trim(mapdir)//"cmb/int/"//trim(cc_method)//"/dx11_v2_"//trim(cc_method)//"_int_cmb_mc_"//trim(coop_Ndigits(i-1, 5))//trim(postfix), nmaps_to_read = 1, nmaps_wanted = sto_max%nmaps)
       call inoise%read(trim(mapdir)//"noise/int/"//trim(cc_method)//"/dx11_v2_"//trim(cc_method)//"_int_noise_mc_"//trim(coop_Ndigits(i-1, 5))//trim(postfix), nmaps_to_read = 1, nmaps_wanted = sto_max%nmaps)
       imap%map(:, 1) = imap%map(:, 1) + inoise%map(:, 1)       
    endif
    if(imap%nmaps .eq. 3)then
       call imap%regularize_in_mask(imask_copy, 1)
       call imap%get_QU()
    endif
    if(remove_mono)imap%map(:, 1) = imap%map(:, 1) - sum(dble(imap%map(:, 1)), mask = imask_logic)/sumimask
    if(do_nest) call imap%convert2nested()
    iloaded = .true.
  end subroutine load_imap

  subroutine load_polmap(i)
    COOP_INT i
    if(polloaded) return
    if(i.eq.0)then
       call polmap%read(polmap_file, nmaps_wanted = 2)
    else
       call polmap%read(trim(mapdir)//"cmb/pol/"//trim(cc_method)//"/dx11_v2_"//trim(cc_method)//"_"//trim(polcase)//"_cmb_mc_"//trim(coop_Ndigits(i-1, 5))//"_hp_20_40"//trim(postfix), nmaps_wanted = 2)
       call polnoise%read(trim(mapdir)//"noise/pol/"//trim(cc_method)//"/dx11_v2_"//trim(cc_method)//"_"//trim(polcase)//"_noise_mc_"//trim(coop_Ndigits(i-1, 5))//"_hp_20_40"//trim(postfix), nmaps_wanted = 2)
       polmap%map = polmap%map + polnoise%map
    endif
    if(do_nest) call polmap%convert2nested()
    polloaded = .true.
  end subroutine load_polmap


  subroutine load_theory()
    COOP_STRING::data_head, sim_head
    call patch_max_data%init(trim(stack_field_name), n, dr)
    call patch_min_data%init(trim(stack_field_name), n, dr)
    call patch_max_sim%init(trim(stack_field_name), n, dr)
    call patch_min_sim%init(trim(stack_field_name), n, dr)  
    
    if(trim(stack_field_name).eq."T")then
        if(trim(orient_name).eq. "NULL")then
           data_head = "pl14fwhm"//COOP_STR_OF(fwhm)//"nu"//trim(coop_num2goodstr(threshold,"-","pt"))//"int"
           sim_head ="pl13fwhm"//COOP_STR_OF(fwhm)//"nu"//trim(coop_num2goodstr(threshold,"-","pt"))//"int"
        else
           data_head = "pl14fwhm"//COOP_STR_OF(fwhm)//"nu"//trim(coop_num2goodstr(threshold,"-","pt"))//"Orientedintint"
           sim_head = "pl13fwhm"//COOP_STR_OF(fwhm)//"nu"//trim(coop_num2goodstr(threshold,"-","pt"))//"Orientedint"
        endif
        cal_data_file = trim(cal_file_prefix)//trim(data_head)//"_frI.txt"
        cal_sim_file = trim(cal_file_prefix)//trim(sim_head)//"_frI.txt"           
     else
        if(trim(Orient_name).eq."NULL")then
           data_head  = "pl14fwhm"//COOP_STR_OF(fwhm)//"nu"//trim(coop_num2goodstr(threshold,"-","pt"))//"pol"
           sim_head = "pl13fwhm"//COOP_STR_OF(fwhm)//"nu"//trim(coop_num2goodstr(threshold,"-","pt"))//"pol"
        else
           data_head = "pl14fwhm"//COOP_STR_OF(fwhm)//"nu"//trim(coop_num2goodstr(threshold,"-","pt"))//"Orientedpol"
           sim_head = "pl13fwhm"//COOP_STR_OF(fwhm)//"nu"//trim(coop_num2goodstr(threshold,"-","pt"))//"Orientedpol"
        endif
        cal_data_file =  trim(cal_file_prefix)//trim(data_head)//"_frQU.txt"
        cal_sim_file =  trim(cal_file_prefix)//trim(sim_head)//"_frQU.txt"                   
     endif
     if(.not. coop_file_exists(cal_data_file))then
        if(trim(orient_name) .eq. "NULL")then
           call system("./GetTheo "//trim(data_theory_cl_file)//" Tmax "//trim(stack_field_name)//" "//COOP_STR_OF(threshold)//" "//COOP_STR_OF(fwhm)//" "//trim(data_head) )        
        else
           call system("./GetTheo "//trim(data_theory_cl_file)//" Tmax_QTUTOrient "//trim(stack_field_name)//" "//COOP_STR_OF(threshold)//" "//COOP_STR_OF(fwhm)//" "//trim(data_head) )
        endif
     endif
     if(.not. coop_file_exists(cal_sim_file))then
        if(trim(orient_name).eq."NULL")then
           call system("./GetTheo "//trim(sim_theory_cl_file)//" Tmax "//trim(stack_field_name)//" "//COOP_STR_OF(threshold)//" "//COOP_STR_OF(fwhm)//" "//trim(sim_head))
        else
           call system("./GetTheo "//trim(sim_theory_cl_file)//" Tmax_QTUTOrient "//trim(stack_field_name)//" "//COOP_STR_OF(threshold)//" "//COOP_STR_OF(fwhm)//" "//trim(sim_head))
        endif
     endif
     call fp%open(cal_data_file)
     read(fp%unit, *) patch_max_data%fr
     call fp%close()
     call fp%open(cal_sim_file)
     read(fp%unit, *) patch_max_sim%fr
     call fp%close()
     
     patch_min_sim%fr = patch_max_sim%fr
     patch_min_data%fr = patch_max_data%fr
     
     select case(trim(stack_field_name))
     case("T")
        patch_min_data%fr(:, 0, 1) = - patch_min_data%fr(:, 0, 1)
        patch_min_sim%fr(:, 0, 1) = - patch_min_sim%fr(:, 0, 1)
     case("QU")
        patch_min_data%fr(:, 1, :) = - patch_min_data%fr(:, 1, :)
        patch_min_sim%fr(:, 1, :) = - patch_min_sim%fr(:, 1, :)
     case default
        write(*,*) "So far SST only supports T and QU stacking"
        stop 
     end select
     call get_radial_f(pfr_theory_data, patch_max_data, patch_min_data, norm_theory_data)
     call get_radial_f(pfr_theory_sim, patch_max_sim, patch_min_sim, norm_theory_sim)     
   end subroutine load_theory

   subroutine get_radial_f(pfr, patch_max, patch_min, norm)
     type(coop_healpix_patch)::patch_min, patch_max
     COOP_REAL pfr(0:n)
     COOP_REAL,optional::norm(0:n)
     select case(trim(stack_field_name))
     case("T")
        select case(trim(fr_genre))
        case("cold0")
           pfr =  patch_min%fr(:,0,1)
           if(present(norm))norm = pfr
        case("cold1")
           pfr = patch_min%fr(:, 1, 1)
           if(present(norm))norm = pfr           
        case("cold2")
           pfr = patch_min%fr(:, 2, 1)
           if(present(norm))norm = pfr           
        case("hot0")
           pfr = patch_max%fr(:,0,1)
           if(present(norm))norm = pfr           
        case("hot1")
           pfr = patch_max%fr(:,1,1)
           if(present(norm))norm = pfr           
        case("hot2")
           pfr = patch_max%fr(:,2,1)
           if(present(norm))norm = pfr           
        case("hc0")
           pfr = (patch_max%fr(:,0,1) + patch_min%fr(:,0,1))/2.d0
           if(present(norm))norm = patch_max%fr(:,0,1)          
        case default
           print*,  "Unknown fr type (cold0, cold1, hot0, hot1)"
           stop
        end select
     case("QU")
        select case(trim(fr_genre))
        case("cold0")
           pfr = ( patch_min%fr(:,0,1) + patch_min%fr(:,0,2))/2.d0
           if(present(norm))norm = pfr           
        case("cold1")
           pfr = (patch_min%fr(:, 1, 1) + patch_min%fr(:, 1, 2))/2.d0
           if(present(norm))norm = pfr           
        case("cold2")
           pfr = (patch_min%fr(:, 2, 1) + patch_min%fr(:, 2, 2))/2.d0
           if(present(norm))norm = pfr           
        case("hot0")
           pfr = (patch_max%fr(:,0,1) + patch_max%fr(:, 0, 2))/2.d0
           if(present(norm))norm = pfr           
        case("hot1")
           pfr = (patch_max%fr(:,1,1) + patch_max%fr(:, 1, 2))/2.d0
           if(present(norm))norm = pfr           
        case("hot2")
           pfr = (patch_max%fr(:,2,1) + patch_max%fr(:, 2, 2))/2.d0
           if(present(norm))norm = pfr           
        case("hc0")
           pfr = ((patch_max%fr(:,0,1) + patch_max%fr(:, 0, 2)) - (patch_min%fr(:,0,1) + patch_min%fr(:,0,2)))/4.d0
           if(present(norm))norm = (patch_max%fr(:,0,1) + patch_max%fr(:, 0, 2))/2.d0          
        case default
           print*,  "Unknown fr type (cold0, cold1, cold2, hot0, hot1, hot2)"
           stop
        end select
     case default
        stop "so far SST only supports T and QU stacking"
     end select
   end subroutine get_radial_f
   
#else
  print*, "You need to install healpix"
#endif  
end program massive_stack


