module coop_statchains_mod
  use coop_wrapper
  use coop_latex_mod

  implicit none

#include "constants.h"

  integer,parameter::mcmc_stat_num_cls = 3  !!minimum 3

  real,dimension(mcmc_stat_num_cls)::mcmc_stat_cls = (/ 0.683, 0.954, 0.997 /)

  COOP_STRING :: measured_cltt_file = ""
  COOP_STRING :: bestfit_cl_file = ""
  COOP_STRING :: bestfit_run_file = ""
  COOP_STRING :: bestfit_varytau_file = ""
  logical::coop_postprocess_do_cls = .false.
  integer::coop_postprocess_nbins = 0
  integer::coop_postprocess_num_contours = 2

  type coop_mcmc_chain
     COOP_STRING prefix, output
     logical::do_preprocess = .false.
     COOP_INT::np_used = 0
     COOP_INT::np = 0
     COOP_SINGLE:: bestlike, worstlike, totalmult, likecut(mcmc_stat_num_cls)
     COOP_INT ibest, iworst
     COOP_SHORT_STRING,dimension(:),allocatable::name, simplename
     COOP_STRING,dimension(:),allocatable::label
     COOP_SINGLE,dimension(:),allocatable::lower, upper, mean, std, plotlower, plotupper, dx, base !!plotlower and plotupper are used to make plots, with ~0 long tail truncated
     real,dimension(:,:),allocatable::lowsig, upsig
     logical,dimension(:),allocatable::vary, left_is_tail, right_is_tail
     integer,dimension(:),allocatable::map2used
     real,dimension(:,:),allocatable::covmat, corrmat, cov_used
     integer,dimension(:),allocatable::used


     integer np_pca
     integer,dimension(:),allocatable::pca  

     integer n
     real,dimension(:),allocatable::mult, like  !! n
     real,dimension(:,:),allocatable::params !! n, np

     integer nb !!number of bins
     real,dimension(:,:),allocatable::c1d  !!1d counts: nb, np_used
     real,dimension(:,:,:),allocatable::c2d !!2d counts: nb, nb, np_used*(np_used+1)/2
     real,dimension(:,:),allocatable::cut2d
     logical,dimension(:,:),allocatable::want_2d_output
     logical,dimension(:),allocatable::want_1d_output
     logical::do_extensions = .false.
     COOP_SHORT_STRING, dimension(mcmc_stat_num_cls)::color2d
     type(coop_dictionary) inputparams
     type(coop_dictionary) allparams
     type(coop_dictionary) usedparams

     COOP_INT::index_epss = 0     
     COOP_INT::index_ombh2 = 0
     COOP_INT::index_omch2 = 0          
     COOP_INT::index_theta = 0
     COOP_INT::index_tau = 0
     COOP_INT::index_mnu = 0
     COOP_INT::index_logA = 0
     COOP_INT::index_ns = 0
     COOP_INT::index_nrun = 0     
     COOP_INT::index_r = 0
     COOP_INT::index_nt = 0
     COOP_INT::index_omegam = 0
     COOP_INT::index_omegab = 0     
     COOP_INT::index_omegak = 0     
     COOP_INT::index_de_w = 0
     COOP_INT::index_de_wa = 0
     COOP_INT::index_de_Q = 0
     COOP_INT::index_de_tracking_n = 0
     COOP_INT::index_de_dUdphi = 0
     COOP_INT::index_de_epsv = 0     
     COOP_INT::index_de_dlnQdphi = 0
     COOP_INT::index_de_d2Udphi2 = 0
     COOP_INT::index_h = 0
     logical::index_set = .false.
   contains
     procedure:: export_stats => coop_mcmc_chain_export_stats
     procedure::load => coop_mcmc_chain_load
     procedure::analyze => coop_mcmc_chain_analyze
     procedure::index_of => coop_mcmc_chain_index_of
     procedure::all_index_of => coop_mcmc_chain_all_index_of     
     procedure::set_indices => coop_mcmc_chain_set_indices
     procedure::name2value => coop_mcmc_chain_name2value
     procedure::getcosmomcparams =>  coop_mcmc_chain_getCosmoMCParams
     procedure::preprocess => coop_mcmc_chain_preprocess
  end type coop_mcmc_chain

contains

  subroutine coop_mcmc_chain_set_indices(this)
    class(coop_mcmc_chain)::this
    this%index_set = .false.
    this%index_epss = this%index_of("epss")
    this%index_ombh2 = this%index_of("ombh2")
    if(this%index_ombh2.eq.0)then
       this%index_ombh2 = this%index_of("omegabh2")
    endif
    this%index_omch2 = this%index_of("omch2")
    if(this%index_omch2 .eq.0 )then
       this%index_omch2 = this%index_of("omegach2")
    endif
    this%index_theta = this%index_of("theta")
    this%index_tau = this%index_of("tau")
    this%index_logA = this%index_of("logA")
    this%index_ns = this%index_of("ns")
    this%index_mnu = this%index_of("mnu")                        
    this%index_nrun = this%index_of("nrun")
    this%index_r = this%index_of("r")
    this%index_nt =      this%index_of("nt")
    this%index_de_w = this%index_of("de_w")
    this%index_de_wa = this%index_of("de_wa")
    this%index_de_Q = this%index_of("de_Q")
    this%index_de_tracking_n = this%index_of("de_tracking_n")
    this%index_de_dUdphi = this%index_of("de_dUdphi")
    this%index_de_epsv = this%index_of("de_epsv")    
    this%index_de_dlnQdphi = this%index_of("de_dlnQdphi")
    this%index_de_d2Udphi2 = this%index_of("de_d2Udphi2")
    this%index_h = this%index_of("h")
    this%index_omegam = this%index_of("omegam")
    this%index_omegab = this%index_of("omegab")    
    this%index_omegak = this%index_of("omegak")
    this%index_set = .true.
  end subroutine coop_mcmc_chain_set_indices

  subroutine coop_mcmc_chain_load(mc, prefix, ignore_percent)
    class(coop_mcmc_chain) mc
    COOP_UNKNOWN_STRING prefix
    integer,optional::ignore_percent
    integer nfiles, i, j, lens, k
    integer,dimension(:),allocatable::nskip, nlines
    COOP_LONG_STRING fname
    COOP_LONG_STRING inline, tmp
    integer ispace, ind
    type(coop_file) fp
    mc%prefix = trim(prefix)
    nfiles = 0
    do 
       fname = trim(prefix)//"_"//trim(coop_num2str(nfiles+1))//".txt"
       if(coop_file_exists(fname))then
          nfiles = nfiles+1
          if(nfiles.eq.1)then
             mc%np = coop_file_numcolumns(fname) -2 
             call coop_feedback(trim(coop_num2str(mc%np))//" parameters")
             if(allocated(mc%label))deallocate(mc%name, mc%simplename, mc%label, mc%std, mc%mean, mc%lower, mc%upper, mc%vary, mc%covmat,  mc%plotlower, mc%plotupper, mc%dx, mc%map2used, mc%lowsig, mc%upsig, mc%left_is_tail, mc%right_is_tail, mc%base)
             allocate(mc%label(mc%np), mc%simplename(mc%np), mc%std(mc%np), mc%mean(mc%np), mc%lower(mc%np), mc%upper(mc%np),  mc%plotlower(mc%np), mc%plotupper(mc%np), mc%vary(mc%np), mc%covmat(mc%np, mc%np), mc%dx(mc%np), mc%map2used(mc%np), mc%name(mc%np), mc%lowsig(mcmc_stat_num_cls, mc%np), mc%upsig(mcmc_stat_num_cls, mc%np), mc%left_is_tail(mc%np), mc%right_is_tail(mc%np), mc%base(mc%np))
          endif
       else
          exit
       endif
    enddo
    if(nfiles.eq.0)then
       call coop_return_error("coop_mcmc_chain_load", trim(prefix)//"_1.txt is not found on the disk", "stop")
    else
       call coop_feedback( "found "//trim(coop_num2str(nfiles))//" chain files on the disk")
    endif
    fname = trim(prefix)//".inputparams"
    if(coop_file_exists(fname))then
       call coop_load_dictionary(trim(fname), mc%inputparams)
       call coop_dictionary_lookup(mc%inputparams, "pp_model", cosmomc_pp_model, COOP_PP_STANDARD)
       call coop_dictionary_lookup(mc%inputparams, "de_model", cosmomc_de_model, COOP_DE_COSMOLOGICAL_CONSTANT)
       call coop_dictionary_lookup(mc%inputparams, "de_num_params", cosmomc_de_num_params, 2)
       call coop_dictionary_lookup(mc%inputparams, "pp_num_params", cosmomc_pp_num_params, 8)       
    else
       call coop_feedback( "Warning: inputparams file not found;" )
       cosmomc_pp_model = COOP_PP_STANDARD
       cosmomc_de_model = COOP_DE_COSMOLOGICAL_CONSTANT
    endif
    fname = trim(prefix)//".ranges"
    if(coop_file_exists(fname))then
       call coop_load_dictionary(trim(fname), mc%allparams, col_key = 1, col_value = 2)
    else
       call coop_feedback("ranges file not found;")
       if(cosmomc_pp_model .ne. COOP_PP_STANDARD .or. cosmomc_de_model .ne. COOP_DE_COSMOLOGICAL_CONSTANT) stop
    endif
    
    fname = trim(prefix)//".paramnames"
    if(coop_file_exists(fname))then
       call fp%open(fname, "r")
       do i =1, mc%np
          if(fp%read_string(inline))then
             inline = trim(adjustl(inline))
             lens = len_trim(inline)
             if(lens .eq. 0)then
                mc%name(i) = "param"//trim(coop_num2str(i))
                mc%label(i) = mc%name(i)
             else
                ispace = scan(inline,  " "//coop_tab)
                if(ispace .eq. 0)then
                   mc%name(i) = trim(inline)
                   if(trim(mc%name(i)) .eq. "") mc%name(i) = "param"//trim(coop_num2str(i))
                   mc%label(i) = mc%name(i)
                else
                   mc%name(i) = trim(inline(1:ispace-1))
                   mc%label(i) = "$"//trim(adjustl(inline(ispace+1:)))//"$"
                   if(trim(mc%label(i)).eq. "$$") mc%label(i) = mc%name(i)
                endif
             endif
          else
             mc%name(i) = "param"//trim(coop_num2str(i))
             mc%label(i) = mc%name(i)
          endif

       enddo
       call fp%close()
       call coop_load_dictionary(trim(fname), mc%usedparams, col_key = 1)
       call mc%set_indices()
    else
       write(*,*) "paramnames file not found"
       stop
    endif
    do i = 1, mc%np
       mc%simplename(i) = trim(coop_str_numalpha(mc%name(i)))
    enddo
    if(mc%do_extensions)then
       cosmomc_de_index = max(coop_mcmc_chain_all_index_of(mc, "meffsterile"), coop_mcmc_chain_all_index_of(mc, "mnu"), coop_mcmc_chain_all_index_of(mc, "omegak")) + 1
       ind = min(coop_mcmc_chain_all_index_of(mc, "logA"), coop_mcmc_chain_all_index_of(mc, "ns"))
       cosmomc_de2pp_num_params = ind - cosmomc_de_index - cosmomc_de_num_params
       cosmomc_pp_num_origin = coop_mcmc_chain_all_index_of(mc, "Aphiphi") - ind + 1

       call coop_dictionary_lookup(mc%inputparams, "inflation_consistency", cosmomc_pp_inflation_consistency, .true.)

       call coop_feedback( trim(coop_num2str(ind-1))//" hard parameters ")
       call coop_feedback(" dark energy index from "//trim(coop_num2str(cosmomc_de_index)))
       call coop_feedback(trim(coop_num2str(cosmomc_de_num_params))//" dark energy parameters ")
       call coop_feedback(trim(coop_num2str(cosmomc_pp_num_origin))//" cosmomc default primordial power parameters ")
       call coop_feedback("totally "//trim(coop_num2str(cosmomc_pp_num_params))//" primordial power parameters ")
    endif

    allocate(nskip(nfiles),nlines(nfiles))
    do i = 1, nfiles
       fname = trim(prefix)//"_"//trim(coop_num2str(i))//".txt"
       nlines(i)= coop_file_numlines(fname)          
       call coop_feedback( "found "//trim(coop_num2str(nlines(i)))//" lines in "//trim(fname))
    enddo
    if(present(ignore_percent))then
       nskip = nlines * ignore_percent/100
    else
       nskip = nlines/4  !!discard 25%
    endif
    mc%n = sum(nlines-nskip)
    if(allocated(mc%like)) deallocate(mc%like, mc%mult, mc%params)
    allocate(mc%like(mc%n), mc%mult(mc%n), mc%params(mc%n, mc%np))
    k = 0
    do i=1, nfiles
       if(nlines(i)-nskip(i) .gt. 0)then
          fname = trim(prefix)//"_"//trim(coop_num2str(i))//".txt"
          call fp%open(fname, "r")
          call fp%skip_lines(nskip(i))
          do j=nskip(i)+1, nlines(i)
             k = k + 1
             read(fp%unit, "(A)", ERR = 120, END=300) inline
             read(inline,  *, ERR = 50)  mc%mult(k), mc%like(k), mc%params(k,1:mc%np)
          enddo
50        call fp%close()
       endif
    enddo
    if( k .eq. mc%n)then
       call coop_feedback( "Totally "//trim(coop_num2str(k))//" samples are used." )
    else
       mc%n = k
       call coop_feedback(" Warning: some lines seem to be broken. This could be caused by inconsistent ini files with checkpoint. Check your MCMC chains.")
    endif
    deallocate(nskip, nlines)
    if(mc%do_preprocess)then
       call coop_feedback( "Preprocessing ... ")
       call mc%preprocess()
    endif
    call coop_feedback( "Analysing the chain ... ")
    call mc%analyze()
    call coop_feedback( "Chain analysed.")
    return
120 call coop_return_error("coop_mcmc_chain_load", "bad line #"//trim(coop_num2str(j))//" in file #"//trim(coop_num2str(i)), "stop")
300 call coop_return_error("load chain", "End of file during read. Line #"//trim(coop_num2str(j))//" in file #"//trim(coop_num2str(i))//". This should not happen unless you have manually added comment lines in to the chain file", "stop")
  end subroutine coop_mcmc_chain_load

  subroutine coop_mcmc_chain_analyze(mc)
    integer,parameter::n_fine_bins = 2048
    real c(n_fine_bins), dx,  multcut, acc, maxc
    real,dimension(:),allocatable::c2dlist
    class(coop_mcmc_chain) mc
    integer ip, i, loc, j, j2, ip2, k, loc2, icl
    mc%totalmult = sum(mc%mult)
    if(mc%totalmult .le. 0 .or. mc%n.eq.0) stop "coop_mcmc_chain_analyzes: found no samples"

    mc%ibest = 1
    mc%iworst = 1
    mc%worstlike = mc%like(1)
    mc%bestlike = mc%like(1)

    do i = 2, mc%n
       if(mc%like(i).lt.mc%bestlike)then
          mc%bestlike = mc%like(i)
          mc%ibest = i
       elseif(mc%like(i) .gt. mc%worstlike)then
          mc%worstlike = mc%like(i)
          mc%iworst = i
       endif
    enddo

    c = 0
    dx = (mc%worstlike - mc%bestlike)/(n_fine_bins-1.)
    do i=1, mc%n
       loc = nint((mc%like(i) - mc%bestlike)/dx+1.)
       c(loc) = c(loc) + mc%mult(i)
    enddo
    i = 1
    acc = c(i)
    do j = 1, mcmc_stat_num_cls
       multcut = mc%totalmult*mcmc_stat_cls(j)
       do while(acc .lt. multcut)
          i = i + 1
          acc  = acc + c(i)
       enddo
       mc%likecut(j) = mc%bestlike + min(max((i-(acc-multcut)/c(i)-0.5), 0.), n_fine_bins-1.) *dx
    enddo



    do ip =1 ,mc%np
       mc%lower(ip) = minval(mc%params(:, ip))
       mc%upper(ip) = maxval(mc%params(:, ip))
       dx = (mc%upper(ip) - mc%lower(ip))/n_fine_bins
       if (dx .lt. 1.d-13 )then
          mc%mean(ip) = mc%upper(ip)
          mc%std(ip) = 0.
          mc%vary(ip) = .false.
          mc%plotupper(ip) = mc%upper(ip)
          mc%plotlower(ip) = mc%lower(ip)
          mc%dx(ip) = 0.
          mc%lowsig(:, ip) = 0. !mc%upper(ip)
          mc%upsig(:, ip) = 0. ! mc%upper(ip)
          mc%left_is_tail(ip) = .false.
          mc%right_is_tail(ip) = .false.
       else
          mc%mean(ip) = sum(mc%params(:,ip)*mc%mult)/mc%totalmult
          mc%std(ip) = sqrt(sum((mc%params(:,ip)-mc%mean(ip))**2*mc%mult)/mc%totalmult)
          mc%vary(ip) = .true.

          mc%upper(ip) = mc%upper(ip) + dx
          mc%lower(ip) = mc%lower(ip) - dx
          dx = (mc%upper(ip) - mc%lower(ip))/n_fine_bins          
          c = 0
          do i=1, mc%n
             loc = nint((mc%params(i, ip)-mc%lower(ip))/dx + 0.5)
             c(loc) = c(loc) + mc%mult(i)
          enddo
          i = 1
          acc = c(i)
          multcut = mc%totalmult*max(min(sqrt(0.01/mc%n), 0.002), 0.0005)
          do while(acc + c(i+1).lt.multcut)
             i = i + 1
             acc = acc + c(i)
             if(i.gt. n_fine_bins/4) exit
          enddo
          call coop_array_get_threshold(c, 0.1, maxc)
          if(sum(c(1:i+n_fine_bins/200))/(i+n_fine_bins/200) .ge. maxc )then  
             mc%left_is_tail(ip) = .false.
             mc%plotlower(ip) = mc%lower(ip)
          else
             mc%left_is_tail(ip) = .true.
             mc%plotlower(ip) = mc%lower(ip) + dx*(i-1)
          endif
          i = n_fine_bins
          acc = c(i)
          do while(acc .lt. multcut)
             i = i - 1
             acc = acc + c(i)
          enddo
          if(i.lt.n_fine_bins)then
             acc = acc-c(i)
             i = i + 1
          endif
          if( sum(c(i-n_fine_bins/200:n_fine_bins))/(n_fine_bins+n_fine_bins/200+1-i) .ge. maxc )then
             mc%right_is_tail(ip) = .false.
             mc%plotupper(ip) = mc%upper(ip)
          else
             mc%right_is_tail(ip) = .true.
             mc%plotupper(ip) = mc%upper(ip) - dx*(n_fine_bins - i)
          endif
          if(mc%right_is_tail(ip) .eqv. mc%left_is_tail(ip))then
             i = 1
             acc = c(i)
             do j = mcmc_stat_num_cls, 1, -1
                multcut = mc%totalmult * (1. - mcmc_stat_cls(j))/2.
                do while(acc .lt. multcut)
                   i = i + 1
                   acc = acc + c(i)
                enddo
                mc%lowsig(j, ip) = dx*(i-(acc-multcut)/c(i))
             enddo
             multcut = mc%totalmult/2.
             do while(acc .lt. multcut)
                i = i + 1
                acc = acc + c(i)
             enddo
             mc%base(ip) = dx*(i-(acc-multcut)/c(i))
             mc%lowsig(:, ip) = mc%base(ip) - mc%lowsig(:, ip)
             do j = 1, mcmc_stat_num_cls
                multcut = mc%totalmult * (1.+mcmc_stat_cls(j))/2.
                do while(acc .lt. multcut )
                   i = i + 1
                   acc = acc + c(i)
                enddo
                mc%upsig(j, ip) = dx * (i -(acc-multcut)/c(i))
             enddo
             mc%upsig(:, ip) = mc%upsig(:,ip) - mc%base(ip)
             mc%base(ip) = mc%base(ip)+mc%lower(ip)
          elseif(mc%right_is_tail(ip))then
             mc%base(ip) = mc%lower(ip)
             mc%lowsig(:, ip) = 0
             i = 1
             acc = c(i)
             do j = 1, mcmc_stat_num_cls
                multcut = mc%totalmult * mcmc_stat_cls(j)
                do while(acc .lt. multcut )
                   i = i + 1
                   acc = acc + c(i)
                enddo
                mc%upsig(j, ip) = dx * (i -(acc-multcut)/c(i))
             enddo
          else
             mc%base(ip) = mc%upper(ip)
             mc%upsig(:, ip) = 0
             i = 1
             acc = c(i)
             do j = mcmc_stat_num_cls, 1, -1
                multcut = mc%totalmult * (1. - mcmc_stat_cls(j))
                do while(acc .lt. multcut)
                   i = i + 1
                   acc = acc + c(i)
                enddo
                mc%lowsig(j, ip) = dx*(n_fine_bins - i + (acc - multcut)/c(i))
             enddo
          endif
       endif

    enddo
    mc%np_used = count(mc%vary)
    mc%np_pca = 0 !!default
    if(allocated(mc%used))deallocate(mc%used)
    if(allocated(mc%corrmat))deallocate(mc%corrmat, mc%cov_used)
    allocate(mc%used(mc%np_used), mc%corrmat(mc%np_used, mc%np_used), mc%cov_used(mc%np_used, mc%np_used))
    i = 1
    do ip=1, mc%np
       if(mc%vary(ip))then
          mc%used(i) = ip
          mc%map2used(ip) = i
          i = i +1
       else
          mc%map2used(ip) = 0
       endif
    enddo
    mc%covmat = 0.
    do i=1, mc%np_used
       do j= 1, i-1
          mc%cov_used(i, j)  = sum((mc%params(:, mc%used(i))-mc%mean(mc%used(i)))*(mc%params(:, mc%used(j))-mc%mean(mc%used(j)))*mc%mult)/mc%totalmult
          mc%cov_used(j, i) = mc%cov_used(i, j)
          mc%covmat(mc%used(i), mc%used(j)) = mc%cov_used(i, j) 
          mc%covmat(mc%used(j),mc%used(i)) =  mc%cov_used(i, j) 
          mc%corrmat(i, j) = mc%cov_used(i, j) / (mc%std(mc%used(i))*mc%std(mc%used(j)))
          mc%corrmat(j, i) = mc%corrmat(i, j)
       enddo
       mc%cov_used(i, i) = mc%std(mc%used(i))**2
       mc%corrmat(i, i) =  1.
       mc%covmat(mc%used(i), mc%used(i)) =  mc%cov_used(i, i) 
    enddo
    if(coop_postprocess_nbins .gt. 0)then
       mc%nb = coop_postprocess_nbins 
    else
       mc%nb = min(max(11, ceiling(sqrt(mc%n/100.))), 22)
    endif
    if(allocated(mc%c1d))deallocate(mc%c1d, mc%c2d, mc%cut2d, mc%want_1d_output, mc%want_2d_output)
    allocate(mc%c1d(mc%nb, mc%np_used), mc%c2d(mc%nb, mc%nb, mc%np_used*(mc%np_used+1)/2), mc%cut2d(mcmc_stat_num_cls, mc%np_used*(mc%np_used+1)/2), mc%want_2d_output(mc%np_used, mc%np_used), mc%want_1d_output(mc%np_used) )
    allocate(c2dlist(mc%nb*mc%nb))
    mc%c1d = 0.
    mc%c2d = 0.
    mc%want_2d_output = .false.
    do j = 1, mc%np_used
       ip = mc%used(j)
       mc%dx(ip) = (mc%plotupper(ip) - mc%plotlower(ip))/mc%nb
       do  i = 1, mc%n
          loc = nint((mc%params(i, ip) - mc%plotlower(ip))/mc%dx(ip)+0.5)
          if(loc .ge. 1 .and. loc .le. mc%nb)then
             mc%c1d(loc, j) = mc%c1d(loc, j) + mc%mult(i)
          endif
       enddo
    enddo
    mc%c1d = mc%c1d/mc%totalmult
    do j = 1, mc%np_used
       ip = mc%used(j)
       do j2 = 1, j
          ip2 = mc%used(j2)
          k = j*(j-1)/2 + j2
          do i = 1, mc%n
             loc =  nint((mc%params(i, ip) - mc%plotlower(ip))/mc%dx(ip)+0.5)
             loc2 =  nint((mc%params(i, ip2) - mc%plotlower(ip2))/mc%dx(ip2)+0.5)
             if(loc.ge.1 .and. loc.le.mc%nb .and. loc2.ge.1 .and. loc2.le.mc%nb)then
                mc%c2d(loc, loc2, k) = mc%c2d(loc, loc2, k) + mc%mult(i)
             endif
          enddo
          mc%c2d(:,:,k) = mc%c2d(:,:,k)/mc%totalmult
          do i=1, mc%nb
             c2dlist((i-1)*mc%nb+1:i*mc%nb) = mc%c2d(:,i,k)
          enddo
          call coop_quicksort(c2dlist)
          i = mc%nb*mc%nb
          acc = c2dlist(i)
          do icl = 1, mcmc_stat_num_cls
             do while(acc .lt. mcmc_stat_cls(icl) .and. i.gt.1)
                i = i - 1
                acc = acc + c2dlist(i)
             enddo
             if(i.lt.mc%nb*mc%nb .and. i.gt.1 .and. c2dlist(i).gt.0.)then
                mc%cut2d(icl, k) = (c2dlist(i) + c2dlist(i-1))/2. + (c2dlist(i+1)-c2dlist(i-1))/2.*(acc - mcmc_stat_cls(icl))/c2dlist(i)
             else
                mc%cut2d(icl, k) = c2dlist(i)
             endif
          enddo
       enddo
    enddo
    deallocate(c2dlist)
  end subroutine coop_mcmc_chain_analyze


  subroutine coop_mcmc_chain_export_stats(mc, output)
    class(coop_mcmc_chain) mc    
    integer,parameter::nk = 71
    logical, parameter::doAsMarg = .true.
    logical,parameter::do_traj_cov = .true.
    integer,parameter::distlss = 13893.
    integer, parameter::num_1sigma_trajs = 50
    integer, parameter::num_samples_to_get_mean = 2500
    integer,parameter::lmin = coop_pp_lmin, lmax = coop_pp_lmax, num_cls_samples = 50
    COOP_REAL,parameter::standard_ns = 0.967d0
    COOP_REAL,parameter::low_ell_cut = 50
    COOP_REAL,parameter::low_k_cut = low_ell_cut/distlss
    COOP_STRING::allnames
    COOP_UNKNOWN_STRING output
    COOP_SHORT_STRING::rval=""
    type(coop_file)::fcl
    type(coop_asy) fp, fig_spec, fig_pot, fig_eps, fig_cls, fig_dcls
    type(coop_asy_path) path    
    COOP_INT i , ip, j,  j2, k, ik, ik1,  ik2, ndof, l, junk, num_params, index_pp, pp_location, icontour, numpp, num_trajs, ind_lowk, isam, index_H, ind_highk
    COOP_SINGLE total_mult, cltraj_weight, x(mc%nb), ytop
    logical first_1sigma, inflation_consistency, do_dcl
    COOP_REAL  norm, lnkmin, lnkmax, cltt, errup, errdown, mean_lnAs, hubble, dns_trial
    COOP_REAL  lnk(nk), kMpc(nk), ps(nk), pt(nk), lnpsmean(nk), lnptmean(nk), lnpscov(nk, nk), lnps(nk), lnpt(nk),  mineig, clnps(coop_pp_n), clnpt(coop_pp_n), lnps_bounds(-2:2,nk), standard_lnps(nk), lnps_samples(num_samples_to_get_mean, nk), lnpt_samples(num_samples_to_get_mean, nk), mult_samples(num_samples_to_get_mean), Cls_samples(lmin:lmax, num_cls_samples), cls_mean(lmin:lmax), cls_best(lmin:lmax)
    COOP_REAL, dimension(:,:),allocatable::pcamat
    COOP_REAL, dimension(:),allocatable::eig, ipca, cosmomcParams
    COOP_REAL:: ps_trajs(nk, num_1sigma_trajs), pt_trajs(nk, num_1sigma_trajs) 
!!knots statistics    
    COOP_REAL,dimension(:),allocatable::lnk_knots, lnps_knots, k_knots, lnps_mean_knots, lnps_standard_knots
    COOP_REAL,dimension(:,:),allocatable::cov_knots, cov_lowk, cov_highk, cov_all
    COOP_REAL,dimension(:),allocatable::shift_knots
    COOP_STRING:: inline
    COOP_REAL::best_ns, best_ns_chisq, this_ns_chisq, best_ns_lowk_chisq, best_ns_highk_chisq
    
    mc%output = trim(adjustl(output))//trim(coop_file_name_of(mc%prefix))
    !! =================================================================!!
    index_pp = cosmomc_de_index + cosmomc_de_num_params + cosmomc_de2pp_num_params
    num_params = index_pp + cosmomc_pp_num_params - 1

    !!------------------------------------------------------------------!!
    call fp%open(trim(mc%output)//"_1sig.samples", "w")
    write(fp%unit, "("//trim(coop_num2str(mc%np))//"G14.5)") mc%params(mc%ibest, :)

    if(mc%do_extensions)then
       allocate(cosmomcParams(num_params))       
       mean_lnAs = mc%mean(mc%index_of("logA"))       
       call coop_feedback("Generating primordial power spectra trajectories")
       call fig_spec%open(trim(mc%output)//"_power_trajs.txt", "w")
       if(coop_postprocess_do_cls)then
          call fig_cls%open(trim(mc%output)//"_cl_trajs.txt", "w")
          do_dcl  = (trim(bestfit_cl_file).ne."")
          if(do_dcl)then
             call fig_dcls%open(trim(mc%output)//"_dcl_trajs.txt", "w")          
             call fcl%open(trim(bestfit_cl_file), 'r')
             do l = lmin, lmax
                read(fcl%unit, *) ik, cls_best(l)
                if(ik.ne. l) stop "Error in bestfit cl file"
             enddo
             call fcl%close()
          endif
       endif
       call fig_pot%open(trim(mc%output)//"_potential_trajs.txt", "w")
       call fig_eps%open(trim(mc%output)//"_eps_trajs.txt", "w")
       lnpsmean = 0
       lnptmean = 0
       clnps = 0
       clnpt = 0
       lnpscov = 0
       call coop_mcmc_chain_getCosmoMCParams(mc, 1, CosmomcParams)
       call coop_setup_cosmology_from_cosmomc(Cosmomcparams)
       call coop_setup_pp()
       numpp = cosmomc_pp_num_params - cosmomc_pp_num_origin + 1
       if(numpp .gt. 4)then
          allocate(lnk_knots(numpp), lnps_knots(numpp), cov_knots(numpp, numpp), k_knots(numpp), lnps_mean_knots(numpp), lnps_standard_knots(numpp))
          call coop_set_uniform(numpp, lnk_knots(1:numpp), coop_pp_lnkmin, coop_pp_lnkmax)
          k_knots = exp(lnk_knots)
          cov_knots = 0.d0
          lnps_mean_knots = 0.d0
          lnps_standard_knots = mean_lnAs+(standard_ns-1.)*(lnk_knots-coop_pp_scalar_lnkpivot)-log(1.d10)          
       endif
       lnkmin = coop_pp_lnkmin
       lnkmax = coop_pp_lnkmax
       call coop_set_uniform(nk, lnk, lnkmin, lnkmax)
       kMpc = exp(lnk)
       standard_lnps = mean_lnAs+(standard_ns -1.)*(lnk-coop_pp_scalar_lnkpivot)

       call fig_spec%init(xlabel="$ k [{\rm Mpc}^{-1}]$", ylabel = "$10^{10}\mathcal{P}_{{\cal R},\mathrm{t}}$", xlog=.true., ylog = .true., xmin = real(exp(coop_pp_lnkmin-0.08)), xmax = real(exp(coop_pp_lnkmax + 0.08)), ymin = 1., ymax = 250., doclip = .true.)
       if(coop_postprocess_do_cls)then
          call fig_cls%init(xlabel = "$\ell$", ylabel ="$\mathcal{D}_\ell (\mu K ^2)$",  xlog = .true., ylog = .false., xmin = 1., xmax = 2000., ymin = 0., ymax = 6000., doclip = .true.)
          if(do_dcl) call fig_dcls%init(xlabel = "$\ell$", ylabel ="$\Delta \mathcal{D}_\ell (\mu K^2)$",  xlog = .true., ylog = .false., xmin = 1.8, xmax = 300., ymin = -500., ymax = 500., doclip = .true.)       
       endif
       call coop_asy_topaxis(fig_spec, xmin = real(exp(coop_pp_lnkmin-0.08))*distlss,  xmax = real(exp(coop_pp_lnkmax + 0.08))*distlss, islog = .true. , label = "$\ell_k\equiv  k D_{\rm rec}$")
       call fig_pot%init(xlabel="$(\phi - \phi_{\rm pivot})/M_p$", ylabel = "$\ln (V/V_{\rm pivot})$", xmin = -1.5, xmax = 0.5, ymin = -0.2, ymax = 0.6, doclip = .true.)
       call fig_eps%init(xlabel = "$ k [{\rm Mpc}^{-1}]$", ylabel = "$\epsilon$", xlog = .true. ,  xmin = real(exp(coop_pp_lnkmin-0.08)), xmax = real(exp(coop_pp_lnkmax + 0.08)), ymin = -0.005, ymax = 0.145, doclip = .true.)
       call coop_asy_topaxis(fig_eps, xmin = real(exp(coop_pp_lnkmin-0.08))*distlss,  xmax = real(exp(coop_pp_lnkmax + 0.08))*distlss, islog = .true. , label = "$\ell_k\equiv  k D_{\rm rec}$")             

       num_trajs = 0
       first_1sigma = .true.
       index_H = mc%usedparams%index("H0*")
       isam = 0
       do while(isam .lt. num_samples_to_get_mean)
          isam = isam + 1
          if(isam .eq. 1)then
             j = mc%ibest
          else
             j = coop_random_index(mc%n)
          endif
          call coop_mcmc_chain_getCosmoMCParams(mc, j, CosmomcParams)
          if(isam .le. num_cls_samples .and. coop_postprocess_do_cls)then
             hubble = mc%params(j, index_H)/100.
             write(*,"(A)") "Computing Cls #"//COOP_STR_OF(isam)//" / "//COOP_STR_OF(num_cls_samples)                 
             call coop_setup_cosmology_from_cosmomc(Cosmomcparams, hubble, want_firstorder = .true.)
             if(isam.eq.1)then
                call fcl%open(trim(mc%output)//"_bestcls.txt", "w")
                do l = 2, lmax
                   norm = (l*(l+1.d0)/coop_2pi * 2.72558e6**2)
                   write(fcl%unit, "(I5, 4E16.7)") l, coop_pp_total_cls(coop_index_clTT, l)*norm, coop_pp_total_cls(coop_index_clTE, l)*norm, coop_pp_total_cls(coop_index_clEE, l)*norm, coop_pp_total_cls(coop_index_clBB, l)*norm
                enddo
                call fcl%close()           
             endif
             write(*,"(A)")"        cosmology setup done"
             do l = lmin, lmax
                cls_samples(l, isam) = coop_pp_total_cls(coop_index_clTT, l)*(l*(l+1.d0)/coop_2pi * 2.72558e6**2)
             enddo
             write(*,"(A)")"        h = "//COOP_STR_OF(hubble)//", C_220 = "//COOP_STR_OF(cls_samples(220, isam))          
          else
             call coop_setup_cosmology_from_cosmomc(Cosmomcparams)
             call coop_setup_pp()
          endif
          !$omp parallel do 
          do ik = 1, nk
             lnps_samples(isam, ik) = coop_primordial_lnps(kMpc(ik))
             lnpt_samples(isam, ik) = coop_primordial_lnpt(kMpc(ik))
          enddo
          !$omp end parallel do
          mult_samples(isam) = mc%mult(j)
          clnps = clnps + coop_pp_lnps*mc%mult(j)
          clnpt = clnpt + coop_pp_lnpt*mc%mult(j)
          if(numpp.gt.4)then
             do ik = 1, numpp
                lnps_knots(ik) = coop_primordial_lnps(k_knots(ik))
             enddo
             lnps_mean_knots = lnps_mean_knots + lnps_knots*mc%mult(j)
             do ik=1, numpp
                do ik2 = ik, numpp
                   cov_knots(ik, ik2) = cov_knots(ik, ik2) + lnps_knots(ik) * lnps_knots(ik2) * mc%mult(j)
                enddo
             enddo
          endif
          if(num_trajs .lt. num_1sigma_trajs .and. mc%like(j) .lt. mc%likecut(1))then
             call coop_pp_get_potential()          
             write(fp%unit, "("//trim(coop_num2str(mc%np))//"G14.5)") mc%params(j, :)
             num_trajs = num_trajs+1
             ps = 1.e10*exp(lnps_samples(isam, :))
             pt = 1.e10*exp(lnpt_samples(isam, :))
             ps_trajs(:, num_trajs) = ps
             pt_trajs(:, num_trajs) = pt
             if(first_1sigma)then
                first_1sigma = .false.
                if(isam .le. num_cls_samples .and. coop_postprocess_do_cls)then
                   call fig_cls%interpolate_curve(xraw = coop_pp_ells, yraw = Cls_samples(:, isam), interpolate = "LogLinear", color = "blue", linetype = "dotted", legend = "$1\sigma$ samples")
                   if(do_dcl)then
                      cls_mean = 0.
                      call fig_dcls%interpolate_curve(xraw = coop_pp_ells, yraw = cls_mean, interpolate = "LogLinear", color="HEX:011010", linewidth = 0.5)

                      call fig_dcls%interpolate_curve(xraw = coop_pp_ells, yraw = Cls_samples(:, isam) - cls_best, interpolate = "LogLinear", color = "blue", linetype = "dotted", legend = "$1\sigma$ samples")
                   endif
                endif
                call fig_pot%interpolate_curve(xraw = coop_pp_phi, yraw = coop_pp_lnV-coop_pp_lnV(coop_pp_ipivot), interpolate="LinearLinear", color = "blue", linetype = "dotted", legend="$1\sigma$, samples")
                call fig_eps%interpolate_curve(xraw = exp(coop_pp_lnkMpc), yraw = exp(coop_pp_lneps), interpolate = "LogLinear", color = "blue", linetype = "dotted", legend="1-$\sigma$ trajs.")
             else
                if(isam .le. num_cls_samples .and. coop_postprocess_do_cls)then
                   call fig_cls%interpolate_curve(xraw = coop_pp_ells, yraw = Cls_samples(:, isam), interpolate = "LogLinear", color = "blue", linetype = "dotted")
                   if(do_dcl) call fig_dcls%interpolate_curve(xraw = coop_pp_ells, yraw = Cls_samples(:, isam)-cls_best, interpolate = "LogLinear", color = "blue", linetype = "dotted")                
                endif
                call fig_pot%interpolate_curve(xraw = coop_pp_phi, yraw = coop_pp_lnV-coop_pp_lnV(coop_pp_ipivot), interpolate="LinearLinear", color = "blue", linetype = "dotted")
                call fig_eps%interpolate_curve(xraw = exp(coop_pp_lnkMpc), yraw = exp(coop_pp_lneps), interpolate = "LogLinear", color = "blue", linetype = "dotted")
             endif
          endif
       enddo
       call fp%close()

       !!now plot the mean
       if(coop_postprocess_do_cls)then    
          total_mult = sum(mult_samples(1:num_cls_samples))
          do l = lmin, lmax
             cls_mean(l) = sum(cls_samples(l, 1:num_cls_samples)*mult_samples(1:num_cls_samples))/total_mult
          enddo
          call fig_cls%interpolate_curve(xraw = coop_pp_ells, yraw = cls_mean, interpolate = "LogLinear", color="HEX:E01010", linewidth = 1.8, legend="mean")
          if(do_dcl) call fig_dcls%interpolate_curve(xraw = coop_pp_ells, yraw = cls_mean - cls_best, interpolate = "LogLinear", color="HEX:E01010", linewidth = 1.9, legend="mean")      

          if(do_dcl)then
             call fig_cls%interpolate_curve(xraw = coop_pp_ells, yraw = cls_best, interpolate = "LogLinear", color="HEX:011010", linewidth = 1.4, legend="$\Lambda$CDM bestfit")
          endif
          if(trim(bestfit_run_file).ne."")then
             call fcl%open(trim(bestfit_run_file), 'r')
             do l = lmin, lmax
                read(fcl%unit, *) ik, cls_mean(l)
                if(ik.ne. l) stop "Error in bestfit cl file"
             enddo
             call fcl%close()
             call fig_cls%interpolate_curve(xraw = coop_pp_ells, yraw = cls_mean, interpolate = "LogLinear", color="HEX:20DA30", linetype="dashed", linewidth = 1., legend="$n_{\rm run}$ bestfit")
             if(do_dcl) call fig_dcls%interpolate_curve(xraw = coop_pp_ells, yraw = cls_mean - cls_best, interpolate = "LogLinear", color="HEX:20DA30", linetype="dashed", linewidth = 1.8, legend="$n_{\rm run}$ bestfit")          
          endif

          if(trim(bestfit_varytau_file).ne."")then
             call fcl%open(trim(bestfit_varytau_file), 'r')
             do l = lmin, lmax
                read(fcl%unit, *) ik, cls_mean(l)
                if(ik.ne. l) stop "Error in bestfit cl file"
             enddo
             call fcl%close()
             call fig_cls%interpolate_curve(xraw = coop_pp_ells, yraw = cls_mean, interpolate = "LogLinear", color="orange", linetype="solid", linewidth = 1.8, legend="$\tau = 0.04$")
             if(do_dcl)call fig_dcls%interpolate_curve(xraw = coop_pp_ells, yraw = cls_mean - cls_best, interpolate = "LogLinear", color="orange", linetype="dotdashed", linewidth = 1., legend="$\tau = 0.04$")          
          endif


          if(trim(measured_cltt_file).ne."")then
             call fcl%open(measured_cltt_file, "r")
             do 
                if( fcl%read_string(inline) )then
                   read(inline, *) l, junk, junk, cltt, errup, errdown
                   if(l.gt. lmax)exit
                else
                   exit
                endif
                call coop_asy_error_bar(fig_cls, x = dble(l), y = cltt, dy_minus = errdown, dy_plus = errup, color="HEX:303030")
                call coop_asy_error_bar(fig_dcls, x = dble(l), y = cltt - cls_best(l), dy_minus = errdown, dy_plus = errup, color="HEX:6F806F")             
             enddo
             call fcl%close()
          endif
       endif
       total_mult = sum(mult_samples)
       do ik = 1, nk
          lnpsmean(ik) = sum(lnps_samples(:, ik)*mult_samples)/total_mult
          lnptmean(ik) = sum(lnpt_samples(:, ik)*mult_samples)/total_mult
       enddo
       ps = 1.e10 * exp(lnpsmean)
       pt = 1.e10 * exp(lnptmean)

       clnps = clnps/total_mult
       clnpt = clnpt/total_mult
       do ik1=1, nk
          do ik2 = 1, ik1
             lnpscov(ik1, ik2) = sum((lnps_samples(:, ik1) - lnpsmean(ik1))*(lnps_samples(:, ik2) - lnpsmean(ik2))*mult_samples)/total_mult
             lnpscov(ik2, ik1) = lnpscov(ik1, ik2)           
          enddo
       enddo
       coop_pp_lnps = clnps
       coop_pp_lnpt = clnpt
       call coop_pp_get_potential()
       do ik = 1, nk
          call coop_get_bounds(lnps_samples(:, ik), (/ 0.023d0, 0.1585d0, 0.5d0, 0.8415d0, 0.977d0 /), lnps_bounds(-2:2, ik), mult_samples)
       enddo
       call fig_spec%band(kmpc, 1.d10*exp(lnps_bounds(-2,:)), 1.d10*exp(lnps_bounds(2,:)), colorfill = trim(coop_asy_gray_color(0.65)), linecolor="invisible")
       call fig_spec%band(kmpc, 1.d10*exp(lnps_bounds(-1,:)), 1.d10*exp(lnps_bounds(1,:)), colorfill = trim(coop_asy_gray_color(0.4)), linecolor="invisible")
       call fig_spec%curve(kmpc, ps_trajs(:,1), color="HEX:006FED", linetype="dashed", linewidth=0.5,legend="$\mathcal{P}_{\cal R}$ samples")
       call fig_spec%curve(kmpc, pt_trajs(:, 1), color="HEX:8CD3F5", linetype="dotted", linewidth=0.5, legend="$\mathcal{P}_{\mathrm{t}}$ samples")       
       do j=2, num_trajs
          call fig_spec%curve(kmpc, ps_trajs(:,j), color="HEX:006FED", linetype="dashed", linewidth=0.5)
          call fig_spec%curve(kmpc, pt_trajs(:, j), color="HEX:8CD3F5", linetype="dotted", linewidth=0.5)
       enddo
       call fig_spec%curve(kmpc, ps, color = "red", linetype = "solid", linewidth = 1.5, legend="mean $\mathcal{P}_{\cal R}$")
       call fig_spec%curve(kmpc, pt, color = "violet", linetype = "solid", linewidth = 1.2, legend="mean $\mathcal{P}_{\mathrm{t}}$")

       call fig_pot%interpolate_curve(xraw = coop_pp_phi, yraw = coop_pp_lnV-coop_pp_lnV(coop_pp_ipivot), interpolate="LinearLinear", color = "red", linetype = "solid", linewidth = 1.5, legend="mean traj")
       call fig_eps%interpolate_curve(xraw = exp(coop_pp_lnkMpc), yraw = exp(coop_pp_lneps), interpolate = "LogLinear", color = "red", linetype = "solid", linewidth = 1.5, legend="mean traj")



       call fig_spec%curve(kMpc, exp(standard_lnps), color = "black", linewidth=1.2, legend="$m^2\phi^2$ model $\mathcal{P}_{\cal R}$")
       call fig_spec%curve(kMpc, exp(mean_lnAs - 0.01625*lnk)*0.13, color = "cyan", linewidth=1.2, legend="$m^2\phi^2$ model $\mathcal{P}_{\mathrm{t}}$")
       if(numpp .gt. 4)then
          ps(1:numpp) = 1.3
          call coop_asy_dots(fig_spec, k_knots, ps(1:numpp), "black", "$\Delta$")
          ps(1:numpp) = 0.005
          call coop_asy_dots(fig_eps, k_knots, ps(1:numpp), "black", "$\Delta$")
       endif
       if(coop_postprocess_do_cls) call coop_asy_legend(fig_cls, 4., 5000., 1, box = .false.)
       if(do_dcl)then
          call coop_asy_legend(fig_dcls, 45., 420., 1, box = .false.)
          call fig_dcls%close()
       endif
       if(mc%index_of("r") .ne. 0)then
          call coop_asy_label(fig_spec,  "free $r$", 0.012, 8., "black")
       else
          rval = trim(mc%inputparams%value("param[r]"))
          if(trim(rval) .ne. "")then
             call coop_asy_label(fig_spec, "fixed $r="//COOP_STR_OF(coop_str2real(rval))//"$", 0.012, 8., "black")
          endif
       endif
       call coop_asy_legend_advance(fig_spec, real(exp(lnkmin + 1.)), 170., "invisible", 0., 0., 0.8, 0.9, 0.9, 2)
       call coop_asy_legend_advance(fig_eps, real(exp(coop_pp_lnkmin +4.)), 0.115,  "invisible", 0., 0., 0.8, 0.9, 0.9, 1)
       call coop_asy_legend_advance(fig_pot, -0.2, 0.35, "invisible", 0., 0., 0.8, 0.9, 0.9, 1)
       call fig_spec%close()
       call fig_cls%close()
       call fig_eps%close()
       call fig_pot%close()


       if(numpp .gt. 4)then
          lnps_mean_knots = lnps_mean_knots /total_mult
          cov_knots = cov_knots/total_mult
          do ik=1, numpp
             do ik2 = ik, numpp
                cov_knots(ik, ik2) = cov_knots(ik, ik2) - lnps_mean_knots(ik)*lnps_mean_knots(ik2)
                cov_knots(ik2, ik) = cov_knots(ik, ik2)
             enddo
          enddo
          ind_lowk = 1
          do while(k_knots(ind_lowk+1).lt. low_k_cut)
             ind_lowk = ind_lowk + 1
             if(ind_lowk .ge. numpp) stop "low_k_cut not in proper range"
          enddo
          if(doAsMarg)then
             ind_highk = numpp - 1 - ind_lowk
             allocate(cov_lowk(ind_lowk, ind_lowk), cov_highk(numpp - ind_lowk-1, numpp - ind_lowk-1), shift_knots(numpp-1), cov_all(numpp-1, numpp-1))
             index_pp = mc%index_of("pp1")
             print*, "pp1 index = ", index_pp
             ind_highk = numpp-1-ind_lowk
             if(index_pp .ne. 0)then
                cov_lowk = mc%covmat(index_pp:index_pp+ind_lowk-1,index_pp:index_pp+ind_lowk-1)
                cov_highk = mc%covmat(index_pp+ind_lowk:index_pp+numpp-2, index_pp+ind_lowk:index_pp+numpp-2)
                cov_all = mc%covmat(index_pp:index_pp+numpp-2, index_pp:index_pp+numpp-2)
                call coop_matsym_inverse(cov_lowk)
                call coop_matsym_inverse(cov_highk)
                call coop_matsym_inverse(cov_all)
                write(*,*) "number of lowk knots =", ind_lowk
                write(*,*) "number of highk knots =", numpp - ind_lowk-1
                do i = 1, coop_pp_nleft
                   shift_knots(i) = coop_pp_lnk_per_knot * (i-coop_pp_nleft-1)
                enddo
                do i=coop_pp_nleft + 1, numpp - 1
                   shift_knots(i) = coop_pp_lnk_per_knot * (i - coop_pp_nleft)
                enddo
                best_ns_chisq = 1.e30
                best_ns = -1.d0
                do i = -10, 10
                   dns_trial = i*0.0002
                   this_ns_chisq = dot_product(mc%mean(index_pp:index_pp+numpp-2) - shift_knots*dns_trial, matmul(cov_all, mc%mean(index_pp:index_pp+numpp-2)- shift_knots*dns_trial))/(numpp-1)
                   if(this_ns_chisq .lt. best_ns_chisq)then
                      best_ns_chisq = this_ns_chisq
                      best_ns = standard_ns + dns_trial
                      if(ind_lowk.gt.0)best_ns_lowk_chisq = dot_product(mc%mean(index_pp:index_pp+ind_lowk-1) - shift_knots(1:ind_lowk)*dns_trial, matmul(cov_lowk, mc%mean(index_pp:index_pp+ind_lowk-1)- shift_knots(1:ind_lowk)*dns_trial))/ind_lowk
                      if(numpp-ind_lowk-1.gt.0)best_ns_highk_chisq = dot_product(mc%mean(index_pp+ind_lowk:index_pp+numpp-2)- shift_knots(ind_lowk+1:numpp-1)*dns_trial, matmul(cov_highk, mc%mean(index_pp+ind_lowk:index_pp+numpp-2)- shift_knots(ind_lowk+1:numpp-1)*dns_trial )) /(numpp - ind_lowk-1)
                   endif
                enddo
                write(*,"(A, F10.4)") "bestfit n_s = ", best_ns
                write(*,*) "Full chi^2(LCDM) per dof = ", best_ns_chisq
                if(ind_lowk.gt.0) write(*,*) "chi_LCDM^2(low k) per dof = ", best_ns_lowk_chisq, "p-value = ", coop_IncompleteGamma(ind_lowk/2.d0, best_ns_lowk_chisq*ind_lowk/2.d0)/gamma(ind_lowk/2.d0)
                if(numpp-ind_lowk-1.gt.0) write(*,*) "chi_LCDM^2(high k) per dof = ", best_ns_highk_chisq, "p-value = ",  coop_IncompleteGamma((numpp-ind_lowk-1)/2.d0, best_ns_highk_chisq*(numpp-ind_lowk-1)/2.d0)/gamma((numpp-ind_lowk-1)/2.d0)

             else
                write(*,*) "cannot find pp1, skipping chi^2 calculation"
             endif
             deallocate(cov_lowk, cov_highk, cov_all, shift_knots)
          else
             allocate(cov_lowk(ind_lowk, ind_lowk), cov_highk(numpp - ind_lowk, numpp - ind_lowk))
             cov_lowk = cov_knots(1:ind_lowk, 1:ind_lowk)
             cov_highk = cov_knots(ind_lowk+1:numpp, ind_lowk+1:numpp)
             call coop_matsym_inverse(cov_lowk)
             call coop_matsym_inverse(cov_highk)
             lnps_standard_knots = lnps_standard_knots - lnps_mean_knots
             write(*,*) "number of lowk knots =", ind_lowk
             write(*,*) "number of highk knots =", numpp - ind_lowk 
             if(ind_lowk.gt.0)write(*,*) "chi_LCDM^2(low k) per dof = ", dot_product(lnps_standard_knots(1:ind_lowk), matmul(cov_lowk, lnps_standard_knots(1:ind_lowk)))/ind_lowk
             if(numpp-ind_lowk.gt.0)write(*,*) "chi_LCDM^2(high k) per dof = ", dot_product(lnps_standard_knots(ind_lowk+1:numpp), matmul(cov_highk, lnps_standard_knots(ind_lowk+1:numpp))) /(numpp - ind_lowk)
             deallocate(cov_lowk, cov_highk)
          endif
       endif

       !!now do eigen modes    
       do ik=1, nk
          lnpscov(ik,ik) = lnpscov(ik,ik) + 1.d-6
       enddo
       call fp%open(trim(mc%output)//"_pwtraj_eig.txt","w")
       call fp%init(xlabel="$ k [{\rm Mpc}^{-1}]$", ylabel = "$\delta \ln \Delta_S^2$", xlog = .true. , xmin = real(exp(lnkmin - 0.01)), xmax = real(exp(lnkmax + 0.01)), width = 7.2, height = 6.)
       call coop_matsym_diagonalize(lnpscov, lnps, sort = .true.)
       mineig = max(lnps(1), 1.d-5)
       ytop = 0.
       j = 1
       if(lnps(2) .le. mineig*1.001)then
          do while(lnps(j) .lt. 2.d0*mineig)
             j = j + 1
             if(j.gt. nk) exit
          enddo
       endif
       ndof = nk-j+1
       norm = sqrt(dble(nk)/(lnkmax-lnkmin))
       call coop_feedback( "found "//trim(coop_num2str(ndof))//" degrees of freedom in scalar power spectrum")
       print*, "power traj eigen values: ", sqrt(lnps(j:nk))/norm
       if(ndof .ge. 1)then
          if(maxval(lnpscov(:,j)).gt. - minval(lnpscov(:, j)))then
             call coop_asy_curve(fp, kmpc, lnpscov(:,j)*norm, color = "red", linewidth = 1.8, legend = "$\sigma_1 = "//trim(coop_num2str(sqrt(lnps(j))/norm, "(F10.4)"))//"$")
             ytop = max(ytop, maxval(lnpscov(:,j)*norm))
          else
             call coop_asy_curve(fp, kmpc, -lnpscov(:,j)*norm, color = "red", linewidth = 1.8, legend = "$\sigma_1 = "//trim(coop_num2str(sqrt(lnps(j))/norm, "(F10.4)"))//"$")
             ytop = max(ytop, maxval(-lnpscov(:,j)*norm))
          endif
       endif
       if(ndof .ge.2)then
          j = j + 1
          if(maxval(lnpscov(:,j)).gt. - minval(lnpscov(:, j)))then
             call coop_asy_curve(fp, kmpc, -lnpscov(:, j)*norm, color = "blue",  linetype = "dotted", linewidth = 1.8, legend = "$\sigma_2 = "//trim(coop_num2str(sqrt(lnps(j))/norm,  "(F10.4)"))//"$")
             ytop = max(ytop, maxval(-lnpscov(:,j)*norm))
          else
             call coop_asy_curve(fp, kmpc, lnpscov(:, j)*norm, color = "blue",  linetype = "dotted", linewidth = 1.8, legend = "$\sigma_2 = "//trim(coop_num2str(sqrt(lnps(j))/norm,  "(F10.4)"))//"$")
             ytop = max(ytop, maxval(lnpscov(:,j)*norm))
          endif
       endif
       if(ndof.ge.3)then
          j = j + 1
          if(maxval(lnpscov(:,j)).gt. - minval(lnpscov(:, j)))then
             call coop_asy_curve(fp, kmpc, lnpscov(:, j)*norm, color = "black", linewidth=1.2, linetype="dashed", legend = "$\sigma_3 = "//trim(coop_num2str(sqrt(lnps(j))/norm, "(F10.4)"))//"$")
             ytop = max(ytop, maxval(lnpscov(:,j)*norm))
          else
             call coop_asy_curve(fp, kmpc, -lnpscov(:, j)*norm, color = "black", linewidth=1.2, linetype="dashed", legend = "$\sigma_3 = "//trim(coop_num2str(sqrt(lnps(j))/norm, "(F10.4)"))//"$")
             ytop = max(ytop, maxval(-lnpscov(:,j)*norm))
          endif
       endif
       if(ndof.ge.4)then
          j = j + 1
          if(maxval(lnpscov(:,j)).gt. - minval(lnpscov(:, j)))then
             call coop_asy_curve(fp, kmpc, -lnpscov(:, j)*norm, color = "violet", linewidth = 1., linetype="dashdotted", legend = "$\sigma_4 = "//trim(coop_num2str(sqrt(lnps(j))/norm, "(F10.4)"))//"$")
             ytop = max(ytop, maxval(-lnpscov(:,j)*norm))
          else
             call coop_asy_curve(fp, kmpc, lnpscov(:, j)*norm, color = "violet", linewidth = 1., linetype="dashdotted", legend = "$\sigma_4 = "//trim(coop_num2str(sqrt(lnps(j))/norm, "(F10.4)"))//"$")
             ytop = max(ytop, maxval(lnpscov(:,j)*norm))
          endif
       endif
       if(ndof .ge.5)then
          j = j + 1
          if(maxval(lnpscov(:,j)).gt. - minval(lnpscov(:, j)))then
             call coop_asy_curve(fp, kmpc, lnpscov(:, j)*norm, color = "green", linewidth=0.9, linetype = "longdashdotted", legend = "$\sigma_5 = "//trim(coop_num2str(sqrt(lnps(j))/norm, "(F10.4)"))//"$")
             ytop = max(ytop, maxval(lnpscov(:,j)*norm))
          else
             call coop_asy_curve(fp, kmpc, -lnpscov(:, j)*norm, color = "green", linewidth=0.9, linetype = "longdashdotted", legend = "$\sigma_5 = "//trim(coop_num2str(sqrt(lnps(j))/norm, "(F10.4)"))//"$")
             ytop = max(ytop, maxval(-lnpscov(:,j)*norm))
          endif
       endif
       if(ndof .ge.6)then
          j = j + 1
          if(maxval(lnpscov(:,j)).gt. - minval(lnpscov(:, j)))then
             call coop_asy_curve(fp, kmpc, -lnpscov(:, j)*norm, color = "gray", linewidth=0.8, linetype = "longdashed", legend = "$\sigma_6 = "//trim(coop_num2str(sqrt(lnps(j))/norm, "(F10.4)"))//"$")
             ytop = max(ytop, maxval(-lnpscov(:,j)*norm))
          else
             call coop_asy_curve(fp, kmpc, lnpscov(:, j)*norm, color = "gray", linewidth=0.8, linetype = "longdashed", legend = "$\sigma_6 = "//trim(coop_num2str(sqrt(lnps(j))/norm, "(F10.4)"))//"$")
             ytop = max(ytop, maxval(lnpscov(:,j)*norm))
          endif
       endif
       if(ndof .ge.7)then
          j = j + 1
          if(maxval(lnpscov(:,j)).gt. - minval(lnpscov(:, j)))then
             call coop_asy_curve(fp, kmpc, lnpscov(:, j)*norm, color = "cyan", linewidth=0.5, legend = "$\sigma_7 = "//trim(coop_num2str(sqrt(lnps(j))/norm, "(F10.4)"))//"$")
             ytop = max(ytop, maxval(lnpscov(:,j)*norm))
          else
             call coop_asy_curve(fp, kmpc, -lnpscov(:, j)*norm, color = "cyan", linewidth=0.5, legend = "$\sigma_7 = "//trim(coop_num2str(sqrt(lnps(j))/norm, "(F10.4)"))//"$")
             ytop = max(ytop, maxval(-lnpscov(:,j)*norm))
          endif
       endif
       if(ndof .ge.8)then
          j = j + 1
          if(maxval(lnpscov(:,j)).gt. - minval(lnpscov(:, j)))then
             call coop_asy_curve(fp, kmpc, -lnpscov(:, j)*norm, color = "skyblue", linewidth=0.5, legend = "$\sigma_8 = "//trim(coop_num2str(sqrt(lnps(j))/norm, "(F10.4)"))//"$")
             ytop = max(ytop, maxval(-lnpscov(:,j)*norm))
          else
             call coop_asy_curve(fp, kmpc, lnpscov(:, j)*norm, color = "skyblue", linewidth=0.5, legend = "$\sigma_8 = "//trim(coop_num2str(sqrt(lnps(j))/norm, "(F10.4)"))//"$")
             ytop = max(ytop, maxval(lnpscov(:,j)*norm))
          endif
       endif
       call coop_asy_legend(fp, 0.0001, ytop*0.9, box = .false.)
       call fp%close()


       deallocate(CosmoMcParams)
       if(allocated(lnk_knots))deallocate(lnk_knots, k_knots, cov_knots, lnps_knots, lnps_mean_knots, lnps_standard_knots)
    endif
    !!likes file
    call fp%open(trim(mc%output)//".likes", "w")
    write(fp%unit, "(A, G14.5)") "Best -lnlike = ", mc%bestlike
    write(fp%unit, "(A, G14.5)") "Worst -lnlike = ", mc%worstlike
    do j = 1, mcmc_stat_num_cls
       write(fp%unit, "(A, F14.3, A, G14.5)") "prob = ", mcmc_stat_cls(j)," truncation like = ", mc%likecut(j)
    enddo
    write(fp%unit,'(A)') "Best params:"    
    do i=1, mc%np
       write(fp%unit,'(A)') 'param['//trim(mc%name(i))//'] = '//COOP_STR_OF(mc%params(mc%ibest, i))
    enddo
    call fp%close()        

    !!margestats    
    call fp%open(trim(mc%output)//".margestats", "w")
    select case(mcmc_stat_num_cls)
    case(3:)
       write(fp%unit, "(A16, A24, 10A15)") "name", "label", "best", "mean", "std", "base", "-1sigma", "+1sigma", "-2sigma", "+2sigma", "-3sigma", "+3sigma"
       do i=1, mc%np
          write(fp%unit, "(A16, A24, 10G15.6)") trim(mc%name(i)), trim(mc%label(i)), mc%params(mc%ibest, i), mc%mean(i), mc%std(i), mc%base(i), mc%lowsig(1, i), mc%upsig(1, i), mc%lowsig(2, i), mc%upsig(2, i),mc%lowsig(3, i), mc%upsig(3, i)
       enddo
    case(2)
       write(fp%unit, "(2A32, 8A16)") "name", "label", "best", "mean", "std", "base", "-1sigma", "+1sigma", "-2sigma", "+2sigma"
       do i = 1, mc%np
          write(fp%unit, "(2A16, 8G16.7)") trim(mc%name(i)), trim(mc%label(i)),  mc%params(mc%ibest, i), mc%mean(i), mc%std(i), mc%base(i), mc%lowsig(1, i), mc%upsig(1, i), mc%lowsig(2, i), mc%upsig(2, i)
       enddo
    case(1)
       write(fp%unit, "(2A32, 6A16)") "name", "label", "best", "mean", "std", "base", "-1sigma", "+1sigma"
       do i = 1, mc%np
          write(fp%unit, "(2A16, 6G16.7)") trim(mc%name(i)), trim(mc%label(i)),  mc%params(mc%ibest, i), mc%mean(i), mc%std(i), mc%base(i), mc%lowsig(1, i), mc%upsig(1, i)
       enddo
    case default
       write(fp%unit, "(2A32, 3A16)") "name", "label", "best", "mean", "std"
       do i = 1, mc%np
          write(fp%unit, "(2A16, 3G16.7)") trim(mc%name(i)), trim(mc%label(i)), mc%params(mc%ibest, i), mc%mean(i), mc%std(i)
       enddo
    end select
    call fp%close()

!!tex    
    call fp%open(trim(mc%output)//".tex", "w")
    Write(fp%unit,"(A)")"\documentclass[12pt]{article}"
    Write(fp%unit,"(a)")"\usepackage[left=2cm, top=2cm, right=2cm,bottom=3cm,nohead]{geometry}"
    Write(fp%unit,"(a)")"\renewcommand{\arraystretch}{1.25}"
    Write(fp%unit,"(A)")"\begin{document}"
    Write(fp%unit,"(A)")"\begin{table}"
    Write(fp%unit,"(A)")"\begin{center} % or you can use \centering to remove the blank line between caption and table"
    Write(fp%unit,"(A)")"\caption{Cosmological Parameters}"
    Write(fp%unit,"(A)")"\label{tbl:cosmology_parameters}"
    Write(fp%unit,"(A)")"\begin{tabular}{cc}"
    Write(fp%unit,"(A)")"\hline"
    Write(fp%unit,"(a)")"\hline"
    do ip = 1, mc%np_used
       if(mc%want_1d_output(ip))then
          write(fp%unit, "(A)") trim(mc%label(mc%used(ip)))//" & "//trim(latex_range(mc%base(mc%used(ip)), mc%lowsig(:, mc%used(ip)), mc%upsig(:, mc%used(ip))))//" \\ "
          Write(fp%unit,"(a)")"\hline"
       endif
    end do
    Write(fp%unit,"(a)")"\end{tabular}"
    Write(fp%unit,"(a)")"\end{center}"
    Write(fp%unit,"(a)")"\end{table}"
    Write(fp%unit,"(a)")"\end{document}"
    call fp%close()

    !!covmat
    allnames = ""
    do ip=1, mc%np_used
       allnames = trim(allnames)//" "//trim(mc%name(mc%used(ip)))
    enddo    
    call fp%open(trim(mc%output)//".covmat", "w")    
    write(fp%unit, "(A)") trim(allnames)
    call coop_write_matrix(fp%unit, mc%cov_used, mc%np_used, mc%np_used)
    call fp%close()
    call fp%open(trim(mc%output)//".corr", "w")
    call coop_write_matrix(fp%unit, mc%corrmat, mc%np_used, mc%np_used)
    call fp%close()

    call fp%open(trim(mc%output)//".mcorr", "w")
    do i = 1, mc%np_used
       do j=1, mc%np_used
          if(j.ne. i .and. mc%corrmat(i, j) .gt. 0.8)then
             write(fp%unit, "(A, F10.3)") trim(mc%name(mc%used(i)))//"  "// trim(mc%name(mc%used(j))), mc%corrmat(i, j)
          endif
       enddo
    enddo
    call fp%close()

!!$    call fp%open(trim(mc%output)//".covused", "w")
!!$    call coop_write_matrix(fp%unit, mc%cov_used, mc%np_used, mc%np_used)
!!$    call fp%close()
    if(mc%np_pca .gt. 0)then
       call coop_feedback( "Doing PCA")
       allocate(pcamat(mc%np_pca, mc%np_pca),eig(mc%np_pca), ipca(mc%np_pca))
       pcamat = mc%covmat(mc%pca, mc%pca)
       call coop_set_uniform(mc%np_pca, ipca, 1.d0, 1.d0*mc%np_pca)
       call coop_matsym_diagonalize(pcamat, eig, sort=.true.) !!sort eigen values
       eig = sqrt(eig)
       call fp%open(trim(mc%output)//".pcamat", "w")
       write(fp%unit, "(A)") "# format is  i, sigma_i (newline) eigen vector (i = 2, 3, ...)"
       do i = 1, mc%np_pca
          write(fp%unit, "(I8, G14.5)") i, eig(i)
          write(fp%unit, "("//trim(coop_num2str(mc%np_pca))//"G14.5)") pcamat(:, i)
       enddo
       call fp%close()
       call fp%open(trim(mc%output)//"_pcafig.txt", "w")
       call fp%init( xlabel="PCA index", ylabel="eigen modes")
       call coop_asy_curve(fp, ipca, pcamat(:,1), smooth = .false., color = "red", linetype = "solid", linewidth = 2., legend="$\sigma_1="//trim(coop_num2str(eig(1),"(G11.2)"))//"$")
       if(mc%np_pca .ge. 2)then
          call coop_asy_curve(fp, ipca, pcamat(:,2), smooth = .false., color = "blue", linetype = "dashed", linewidth = 1.5, legend = "$\sigma_2="//trim(coop_num2str(eig(2), "(G11.2)"))//"$")
       endif
       if(mc%np_pca .ge. 3)then
          call coop_asy_curve(fp, ipca, pcamat(:,3), smooth = .false., color = "black", linetype = "dotted", linewidth = 1., legend = "$\sigma_3="//trim(coop_num2str(eig(3), "(G11.2)"))//"$")
       endif
       if(mc%np_pca .ge. 4)then
          call coop_asy_curve(fp, ipca, pcamat(:,4), smooth = .false., color = "green", linetype = "dashdotted", linewidth = 0.8, legend =  "$\sigma_4="//trim(coop_num2str(eig(4), "(G11.2)"))//"$")
       endif
       if(mc%np_pca .ge. 5)then
          call coop_asy_curve(fp, ipca, pcamat(:,5), smooth = .false., color = "gray", linetype = "longdashdotted", linewidth = 0.8, legend =  "$\sigma_5="//trim(coop_num2str(eig(5), "(G11.2)"))//"$")
       endif
       call coop_asy_legend(fp)
       call fp%close()
       deallocate(pcamat , eig , ipca)
    endif
    do ip = 1, mc%np_used
       if(mc%want_1d_output(ip))then
          call fp%open(trim(mc%output)//"_"//trim(mc%simplename(mc%used(ip)))//"_1D.txt", "w")
          do i = 1, mc%nb
             x(i) = mc%plotlower(mc%used(ip)) + mc%dx(mc%used(ip))*(i-1)
          enddo
          call fp%init( xlabel = trim(mc%label(mc%used(ip))), ylabel = "P", width=3., height=2.5, ymin=0., ymax=1.05)
          call coop_asy_plot_likelihood(fp, x, mc%c1d(:, ip)/maxval(mc%c1d(:, ip)), left_tail = .true., right_tail = .true., linewidth=1.5)
          call fp%close()
       endif
    enddo
    if(coop_postprocess_num_contours .gt. 0)then
       do j = 1, mc%np_used
          do j2 = 1, j
             k = j*(j-1)/2 + j2
             if(mc%want_2d_output(j, j2))then
                call fp%open(trim(mc%output)//"_"//trim(mc%simplename(mc%used(j)))//"_"//trim(mc%simplename(mc%used(j2)))//"_2D.txt", "w")
                call fp%init( xlabel = trim(mc%label(mc%used(j))), ylabel = trim(mc%label(mc%used(j2))), xmin=mc%plotlower(mc%used(j)), xmax = mc%plotupper(mc%used(j)), ymin=mc%plotlower(mc%used(j2)), ymax = mc%plotupper(mc%used(j2)), width=3., height=2.5 )
                do icontour = coop_postprocess_num_contours, 1, -1
                   call coop_asy_path_from_array(path, mc%c2d(:, :, k), mc%plotlower(mc%used(j)), mc%plotupper(mc%used(j)), mc%plotlower(mc%used(j2)), mc%plotupper(mc%used(j2)), mc%cut2d(icontour, k))

                   call coop_asy_contour(fp, path, colorfill = trim(mc%color2d(icontour)), smooth = .false., linecolor = "black", linetype = "solid")
                enddo
                call fp%close()
             endif
             if(j2 .ne. j .and. mc%want_2d_output(j2, j))then
                call fp%open(trim(mc%output)//"_"//trim(mc%simplename(mc%used(j2)))//"_"//trim(mc%simplename(mc%used(j)))//"_2D.txt", "w")
                call fp%init( xlabel = trim(mc%label(mc%used(j2))), ylabel = trim(mc%label(mc%used(j))), xmin=mc%plotlower(mc%used(j2)), xmax = mc%plotupper(mc%used(j2)), ymin=mc%plotlower(mc%used(j)), ymax = mc%plotupper(mc%used(j)) )
                do icontour = coop_postprocess_num_contours, 1, -1
                   call coop_asy_path_from_array(path, transpose(mc%c2d(:, :, k)), mc%plotlower(mc%used(j2)), mc%plotupper(mc%used(j2)),  mc%plotlower(mc%used(j)), mc%plotupper(mc%used(j)), mc%cut2d(icontour, k))

                   call coop_asy_contour(fp, path, colorfill = trim(mc%color2d(icontour)), smooth = .false., linecolor = "black", linetype = "solid")
                enddo
                call fp%close()

             endif
          enddo
       enddo
    endif
  end subroutine coop_mcmc_chain_export_stats


  subroutine coop_mcmc_chain_getCosmoMCParams(this, ind, Params)
    class(coop_mcmc_chain) this
    COOP_REAL  Params(:)
    integer nparams    
    integer i, ind
    nparams = size(Params)
    do i = 1, nparams
       params(i) = this%name2value(ind,  trim(this%allparams%key(i)))
    enddo
  end subroutine coop_mcmc_chain_getCosmoMCParams


  function coop_mcmc_chain_all_index_of(mc, name) result(ind)
    class(coop_mcmc_chain) mc
    integer ind
    COOP_UNKNOWN_STRING name
    ind = mc%allparams%index(trim(name))
  end function coop_mcmc_chain_all_index_of


  function coop_mcmc_chain_index_of(this, name) result(ind)
    class(coop_mcmc_chain) this
    COOP_UNKNOWN_STRING name
    COOP_INT ind
    ind = this%usedparams%index(trim(adjustl(name)))   
    return
  end function coop_mcmc_chain_index_of

  function coop_mcmc_chain_name2value(this, isample, name) result(val)
    class(coop_mcmc_chain) this
    integer isample
    COOP_UNKNOWN_STRING name
    real val
    integer ind
    COOP_SHORT_STRING str
    if(isample .le. 0 .or. isample .gt. this%n)  stop "coop_mcmc_chain_name2value: index overflow"
    ind = this%index_of(name)
    if(ind.eq.0)then
       str = this%inputparams%value("param["//trim(name)//"]")
       if(trim(str).eq."") stop "coop_mcmc_chain_name2value: cannot find the name"
       read(str, *) val
    else
       val = this%params( isample, ind )
    end if
  end function coop_mcmc_chain_name2value

  subroutine coop_mcmc_chain_preprocess(this)
    class(coop_mcmc_chain) this
    COOP_INT::i
    if(.not. this%do_preprocess) return
    call coop_feedback("rescaling probability with sqrt(epsilon_s)")
    if(this%index_epss .ne. 0 .and. this%index_de_dUdphi .ne. 0)then
       do i=1, this%n
          this%mult(i) = this%mult(i)*sqrt(this%params(i, this%index_epss))
       enddo
    endif
       
  end subroutine coop_mcmc_chain_preprocess





End Module Coop_statchains_mod
