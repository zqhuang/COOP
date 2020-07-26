module coop_statchains_mod
  use coop_wrapper
  use coop_latex_mod
  use camb_mypp

  implicit none

#include "constants.h"
#define PK_SPLINE 1
#define PK_WAVELET 2
  
  COOP_INT,parameter::mcmc_stat_num_cls = 3  !!minimum 3

  COOP_SINGLE,dimension(mcmc_stat_num_cls)::mcmc_stat_cls = (/ 0.683, 0.954, 0.997 /)

  COOP_STRING :: measured_cltt_file = ""
  COOP_STRING :: bestfit_cl_file = ""
  COOP_STRING :: bestfit_run_file = ""
  COOP_STRING :: bestfit_varytau_file = ""
  logical::coop_postprocess_do_cls = .false.
  COOP_INT::coop_postprocess_nbins = 0
  COOP_INT::coop_postprocess_num_contours = 2

  type coop_mcmc_chain
     COOP_STRING prefix, output
     logical::do_preprocess = .false.
     COOP_INT::np_used = 0
     COOP_INT::np = 0
     COOP_INT::nprim = 0
     COOP_INT::nwant = 0
     COOP_SINGLE:: bestlike, worstlike, totalmult, likecut(mcmc_stat_num_cls)
     COOP_INT ibest, iworst
     COOP_SHORT_STRING,dimension(:),allocatable::name, simplename
     logical,dimension(:),allocatable::derived
     COOP_STRING,dimension(:),allocatable::label
     COOP_SINGLE,dimension(:),allocatable::lower, upper, mean, skewness, kurtosis, std, plotlower, plotupper, dx, base !!plotlower and plotupper are used to make plots, with ~0 long tail truncated
     COOP_SINGLE,dimension(:,:),allocatable::lowsig, upsig
     logical,dimension(:),allocatable::vary, left_is_tail, right_is_tail
     COOP_INT,dimension(:),allocatable::map2used
     COOP_SINGLE,dimension(:,:),allocatable::covmat, corrmat
     COOP_INT,dimension(:),allocatable::used, prim, want  !!used(i): full-index of the i-th used parameter;  prim(i): used-index of the i-th used primary parameter


     COOP_INT np_pca
     COOP_INT,dimension(:),allocatable::pca  

     COOP_INT::n = 0
     COOP_SINGLE,dimension(:),allocatable::mult, like  !! n
     COOP_SINGLE,dimension(:,:),allocatable::params !! n, np

     COOP_INT::nb = 10 !!number of bins
     COOP_SINGLE,dimension(:,:),allocatable::c1d  !!1d counts: nb, np_used
     COOP_SINGLE,dimension(:,:,:),allocatable::c2d !!2d counts: nb, nb, np_used*(np_used+1)/2
     COOP_SINGLE,dimension(:,:),allocatable::cut2d
     logical,dimension(:,:),allocatable::want_2d_output
     logical,dimension(:),allocatable::want_1d_output
     COOP_INT::extmode = 0
     COOP_SHORT_STRING, dimension(mcmc_stat_num_cls)::color2d
     COOP_STRING::datasets 
     type(coop_dictionary) settings
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
     COOP_INT::index_pp = 0
     COOP_INT::index_de_w = 0
     COOP_INT::index_de_wa = 0
     COOP_INT::index_de_Q = 0
     COOP_INT::index_de_tracking_n = 0
     COOP_INT::index_de_dUdphi = 0
     COOP_INT::index_de_epsv = 0     
     COOP_INT::index_de_dlnQdphi = 0
     COOP_INT::index_de_d2Udphi2 = 0
     COOP_INT::index_h = 0
     COOP_INT::index_h0 = 0     
     logical::index_set = .false.
     COOP_REAL::default_ns = 0.967
     COOP_REAL::default_r = 0.
     COOP_REAL::fit_skewness_threshold = 0.1d0
     COOP_REAL::fit_kurtosis_threshold = 0.02d0     
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
    this%index_pp = this%index_of("pp1")
    if(this%index_pp .eq. 0)then
       this%index_pp = this%index_of("A0m2")
    endif
    this%index_nt =      this%index_of("nt")
    this%index_de_w = this%index_of("de_w")
    if(this%index_de_w .eq. 0)this%index_de_w = this%index_of("w")
    this%index_de_wa = this%index_of("de_wa")
    if(this%index_de_wa .eq. 0)this%index_de_wa = this%index_of("wa")    
    this%index_de_Q = this%index_of("de_Q")
    this%index_de_tracking_n = this%index_of("de_tracking_n")
    this%index_de_dUdphi = this%index_of("de_dUdphi")
    this%index_de_epsv = this%index_of("de_epsv")    
    this%index_de_dlnQdphi = this%index_of("de_dlnQdphi")
    this%index_de_d2Udphi2 = this%index_of("de_d2Udphi2")
    this%index_h = this%index_of("h")
    this%index_H0 = this%index_of("H0*")    
    this%index_omegam = this%index_of("omegam")
    this%index_omegab = this%index_of("omegab")    
    this%index_omegak = this%index_of("omegak")
    this%index_set = .true.
  end subroutine coop_mcmc_chain_set_indices

  subroutine coop_mcmc_chain_load(mc, prefix, ignore_percent)
    class(coop_mcmc_chain) mc
    COOP_UNKNOWN_STRING prefix
    COOP_INT,optional::ignore_percent
    COOP_INT nfiles, i, j, lens, k
    COOP_INT,dimension(:),allocatable::nskip, nlines
    COOP_LONG_STRING fname
    COOP_LONG_STRING inline, tmp
    COOP_INT ispace, ind
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
             if(allocated(mc%label))deallocate(mc%name, mc%derived, mc%simplename, mc%label, mc%std, mc%mean,mc%skewness, mc%kurtosis, mc%lower, mc%upper, mc%vary,  mc%plotlower, mc%plotupper, mc%dx, mc%map2used, mc%lowsig, mc%upsig, mc%left_is_tail, mc%right_is_tail, mc%base, mc%covmat)
             allocate(mc%label(mc%np), mc%simplename(mc%np), mc%std(mc%np), mc%mean(mc%np), mc%skewness(mc%np),mc%kurtosis(mc%np),  mc%lower(mc%np), mc%upper(mc%np), mc%covmat(mc%np, mc%np), mc%plotlower(mc%np), mc%plotupper(mc%np), mc%vary(mc%np), mc%dx(mc%np), mc%map2used(mc%np), mc%name(mc%np),mc%derived(mc%np), mc%lowsig(mcmc_stat_num_cls, mc%np), mc%upsig(mcmc_stat_num_cls, mc%np), mc%left_is_tail(mc%np), mc%right_is_tail(mc%np), mc%base(mc%np))
          else
             i = coop_file_numcolumns(fname) -2              
             if(i .ne. mc%np)then
                write(*,*) "Inconsistent number of columns in  "//trim(fname), ": ", i
                stop
             endif
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
       call coop_dictionary_lookup(mc%inputparams, "param[ns]", mc%default_ns, 0.967d0)
       call coop_dictionary_lookup(mc%inputparams, "param[r]", mc%default_r, 0.d0)       
       if(mc%extmode .eq. PK_WAVELET) then
          mypp_nknots = 64  !!wavelet mode
       else
          call coop_dictionary_lookup(mc%inputparams, "mypp_nknots", mypp_nknots, 0)          
       endif
       print*, "number of knots  = ", mypp_nknots
    else
       call coop_feedback( "Warning: inputparams file not found;" )
    endif
    fname = trim(prefix)//".ranges"
    if(coop_file_exists(fname))then
       call coop_load_dictionary(trim(fname), mc%allparams, col_key = 1, col_value = 2)
    else
       call coop_feedback("ranges file not found;")
       stop
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
                mc%derived(i) = .false.                
             else
                ispace = scan(inline,  " "//coop_tab)
                if(ispace .eq. 0)then
                   mc%name(i) = trim(inline)
                   mc%derived(i) = .false.                                   
                   if(trim(mc%name(i)) .eq. "") mc%name(i) = "param"//trim(coop_num2str(i))
                   mc%label(i) = mc%name(i)
                else
                   mc%name(i) = trim(inline(1:ispace-1))
                   mc%derived(i) = (inline(ispace-1:ispace-1) .eq. "*")
                   mc%label(i) = "$"//trim(adjustl(inline(ispace+1:)))//"$"
                   if(trim(mc%label(i)).eq. "$$") mc%label(i) = mc%name(i)
                endif
             endif
          else
             mc%name(i) = "param"//trim(coop_num2str(i))
             mc%label(i) = mc%name(i)
             mc%derived(i) = .false.                             
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
    allocate(nskip(nfiles),nlines(nfiles))
    do i = 1, nfiles
       fname = trim(prefix)//"_"//trim(coop_num2str(i))//".txt"
       nlines(i)= coop_file_numlines(fname)          
       call coop_feedback( "found "//trim(coop_num2str(nlines(i)))//" lines in "//trim(fname))
    enddo
    if(present(ignore_percent))then
       if(ignore_percent.lt.100)then
          nskip = nlines * ignore_percent/100
       else
          nskip = ignore_percent
       endif
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
       print*, k, mc%n
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
    COOP_INT,parameter::n_coarse_bins = 32
    COOP_INT,parameter::n_cc = 64    
    COOP_INT,parameter::n_fine_bins = n_coarse_bins*n_cc
    COOP_SINGLE:: c(n_fine_bins), dx,  multcut, acc, maxc, cc(n_coarse_bins)
    COOP_SINGLE,dimension(:),allocatable::c2dlist
    class(coop_mcmc_chain) mc
    COOP_INT ip, i, loc, j, j2, ip2, k, loc2, icl
    mc%totalmult = sum(mc%mult)
    if(mc%totalmult .le. 0 .or. mc%n.eq.0)then
       stop "coop_mcmc_chain_analyzes: found no samples"
    endif

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

    write(*,*) "Like cuts done"

    do ip =1 ,mc%np
       mc%lower(ip) = minval(mc%params(:, ip))
       mc%upper(ip) = maxval(mc%params(:, ip))
       dx = (mc%upper(ip) - mc%lower(ip))/n_fine_bins
       if (dx .lt. 1.d-13 )then
          mc%mean(ip) = mc%upper(ip)
          mc%std(ip) = 0.
          mc%skewness(ip) = 0.
          mc%kurtosis(ip) = 0.
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
          mc%skewness(ip) = sum((mc%params(:,ip)-mc%mean(ip))**3*mc%mult)/mc%totalmult /mc%std(ip)**3
          mc%kurtosis(ip) = sum((mc%params(:,ip)-mc%mean(ip))**4*mc%mult)/mc%totalmult/mc%std(ip)**4/3.d0 - 1.d0
          mc%vary(ip) = .true.

!!$          mc%upper(ip) = mc%upper(ip) + dx/5.d0
!!$          mc%lower(ip) = mc%lower(ip) - dx/5.d0
          dx = (mc%upper(ip) - mc%lower(ip))/n_fine_bins          
          c = 0
          do i=1, mc%n
             loc = nint((mc%params(i, ip)-mc%lower(ip))/dx + 0.5)
             c(loc) = c(loc) + mc%mult(i)
          enddo
          do i=1, n_coarse_bins
             cc(i) = sum(c((i-1)*n_cc+1:i*n_cc))/n_cc
          enddo
          i = 1
          acc = c(i)
          multcut = mc%totalmult*max(min(sqrt(0.01/mc%n), 5.e-3), 1.e-5)
          do while(acc + c(i+1).lt.multcut)
             i = i + 1
             acc = acc + c(i)
             if(i.gt. n_fine_bins/4) exit
          enddo
          call coop_array_get_threshold(cc, 0.3, maxc)
          if(sum(c(1:i+n_cc))/(i+n_cc) .ge. maxc )then  
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
             if(i .lt. n_fine_bins*3/4) exit             
          enddo
          if(i.lt.n_fine_bins)then
             acc = acc-c(i)
             i = i + 1
          endif
          if( sum(c(i-n_cc:n_fine_bins))/(n_fine_bins+n_cc+1-i) .ge. maxc )then
             mc%right_is_tail(ip) = .false.
             mc%plotupper(ip) = mc%upper(ip)
          else
             mc%right_is_tail(ip) = .true.
             mc%plotupper(ip) = mc%upper(ip) - dx*(n_fine_bins - i+1)
          endif
          if(mc%upper(ip) - mc%plotupper(ip)  .lt. (mc%upper(ip)-mc%lower(ip))*0.05)then
             mc%plotupper(ip) = mc%upper(ip)
          endif
          if(mc%plotlower(ip) - mc%lower(ip)  .lt. (mc%upper(ip)-mc%lower(ip))*0.05)then
             mc%plotlower(ip) = mc%lower(ip)
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
    mc%nprim = count(mc%vary .and. .not. mc%derived)
    
    COOP_DEALLOC(mc%prim)
    COOP_DEALLOC(mc%used)
    COOP_DEALLOC(mc%corrmat)
    
    mc%np_pca = 0 !!default

    
    allocate(mc%used(mc%np_used), mc%corrmat(mc%np_used, mc%np_used), mc%prim(mc%nprim))
    i = 1
    j=1
    do ip=1, mc%np
       if(mc%vary(ip))then
          mc%used(i) = ip
          mc%map2used(ip) = i
          i = i +1                    
          if(.not. mc%derived(ip))then
             mc%prim(j) = ip  
             j = j +1
          endif
       else
          mc%map2used(ip) = 0
       endif
    enddo
    mc%covmat = 0.d0
    do i=1, mc%np_used
       do j= 1, i-1
          mc%covmat(mc%used(i), mc%used(j))  = sum((mc%params(:, mc%used(i))-mc%mean(mc%used(i)))*(mc%params(:, mc%used(j))-mc%mean(mc%used(j)))*mc%mult)/mc%totalmult
          mc%covmat(mc%used(j), mc%used(i)) = mc%covmat(mc%used(i), mc%used(j))
          mc%corrmat(i, j) =  mc%covmat(mc%used(i), mc%used(j))/ (mc%std(mc%used(i))*mc%std(mc%used(j)))
          mc%corrmat(j, i) = mc%corrmat(i, j)
       enddo
       mc%covmat(mc%used(i), mc%used(i)) = mc%std(mc%used(i))**2
       mc%corrmat(i, i) =  1.
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
    COOP_SINGLE::spec_ymin = 0.35, spec_ymax = 250., eps_ymax = 0.016
    COOP_INT,parameter::nk = 128
    COOP_INT,parameter::distlss = 13893.
    COOP_INT, parameter::num_1sigma_trajs = 50
    COOP_INT::num_samples_to_get_mean
    COOP_INT,parameter::lmin = coop_pp_lmin, lmax = coop_pp_lmax, num_cls_samples = 50
    COOP_REAL,parameter::low_ell_cut = 50
    COOP_REAL,parameter::low_k_cut = low_ell_cut/distlss
    COOP_STRING::allnames
    COOP_UNKNOWN_STRING output
    COOP_SHORT_STRING::rval=""
    type(coop_file)::fcl
    type(coop_asy) fp, fig_spec, fig_pot, fig_eps
    type(coop_asy_path) path    
    COOP_INT i , ip, j,  j2, k, ik, ik1,  ik2, ndof, l, junk, num_params,  pp_location, icontour, numpp, num_trajs, ind_lowk, isam, index_H, num_convex, ind_pivot
    COOP_SINGLE total_mult, cltraj_weight, x(mc%nb), ytop
    logical first_1sigma, inflation_consistency, do_dcl
    COOP_REAL  norm, lnkmin, lnkmax, cltt, errup, errdown, mean_lnAs, hubble, dns_trial, mean_ns
    COOP_REAL  lnk(nk),  standard_lnps(nk), kMpc(nk), ps(nk), pt(nk), lnpsmean(nk), lnptmean(nk), lnpscov(nk, nk), lnps(nk), lnpt(nk),  mineig, clnps(mypp_n), clnpt(mypp_n), lnps_bounds(-2:2,nk), lnpt_bounds(0:2,nk), eps_bounds(0:2,nk), lnV_bounds(-2:2, nk), phi_rs(nk)
    COOP_REAL, dimension(:,:),allocatable::pcamat, eps_samples, phi_samples, lnV_samples, lnps_samples, lnpt_samples    
    COOP_REAL, dimension(:),allocatable::eig, ipca
    COOP_REAL:: ps_trajs(nk, num_1sigma_trajs), pt_trajs(nk, num_1sigma_trajs), eps_trajs(nk, num_1sigma_trajs), phi_trajs(nk, num_1sigma_trajs), lnV_trajs(nk, num_1sigma_trajs)
!!knots statistics    
    COOP_REAL,dimension(:),allocatable::lnk_knots, lnps_knots, k_knots, lnps_mean_knots, lnps_standard_knots, mult_samples
    COOP_REAL,dimension(:,:),allocatable::cov_knots, cov_lowk, cov_highk, cov_all
    logical,dimension(:),allocatable::used
    COOP_STRING:: inline
    COOP_REAL::best_ns, best_ns_chisq, this_ns_chisq, best_ns_lowk_chisq, best_ns_highk_chisq

    num_samples_to_get_mean = min(10000, mc%n/3)
    allocate(eps_samples(num_samples_to_get_mean, nk), phi_samples(num_samples_to_get_mean, nk), lnV_samples(num_samples_to_get_mean, nk), lnps_samples(num_samples_to_get_mean, nk), lnpt_samples(num_samples_to_get_mean, nk), mult_samples(num_samples_to_get_mean))
    allocate(used(mc%n))
    used = .false.
    mc%output = trim(adjustl(output))//trim(coop_file_name_of(mc%prefix))
    !! =================================================================!!



    !!------------------------------------------------------------------!!
    call fp%open(trim(mc%output)//"_1sig.samples", "w")
    write(fp%unit, "("//trim(coop_num2str(mc%np))//"G14.5)") mc%params(mc%ibest, :)

    select case(mc%extmode)
    case(PK_SPLINE, PK_WAVELET)
       if(mypp_nknots .eq. 0) stop "Error: cannot find number of knots for power spectra extensions"      
       if(mc%index_of("ns") .ne. 0)then
          call coop_dictionary_lookup(mc%settings, "m2phi2_ns", mean_ns, dble(mc%mean(mc%index_of("ns"))))
       else
          mean_ns = mc%default_ns
       endif
       call coop_dictionary_lookup(mc%settings, "m2phi2_logA", mean_lnAs, dble(mc%mean(mc%index_of("logA"))))
       call coop_feedback("Generating primordial power spectra trajectories")
       write(*,*) "fiducial trajectory logA = ", mean_lnAs, " ns =", mean_ns       
       call fig_spec%open(trim(mc%output)//"_power_trajs.txt", "w")
       call fig_pot%open(trim(mc%output)//"_potential_trajs.txt", "w")
       call fig_eps%open(trim(mc%output)//"_eps_trajs.txt", "w")
       lnpsmean = 0
       lnptmean = 0
       clnps = 0
       clnpt = 0
       lnpscov = 0
       if(mc%extmode .eq. 1)then
          if(mc%index_r.eq.0)then
             call mypp_setup_pp(As = exp(dble(mc%Params(1, mc%index_logA)))*1.d-10, ns = mean_ns, nknots = mypp_nknots, dlnps = dble(mc%Params(1, mc%index_pp: mc%index_pp+mypp_nknots-1)), r = mc%default_r )
          else
             call mypp_setup_pp(As = exp(dble(mc%Params(1, mc%index_logA)))*1.d-10, ns = mean_ns, nknots = mypp_nknots, dlnps = dble(mc%Params(1, mc%index_pp: mc%index_pp+mypp_nknots-1)), r = dble(mc%params(1, mc%index_r)) )
          endif
       else
          if(mc%index_r.eq.0)then
             call mypp_setup_pp(As = exp(dble(mc%Params(1, mc%index_logA)))*1.d-10, ns = dble(mc%Params(1, mc%index_ns)), nknots = mypp_nknots, dlnps = dble(mc%Params(1, mc%index_pp: mc%index_pp+mypp_nknots-1)), r = mc%default_r )
          else
             call mypp_setup_pp(As = exp(dble(mc%Params(1, mc%index_logA)))*1.d-10, ns = dble(mc%Params(1, mc%index_ns)), nknots = mypp_nknots, dlnps = dble(mc%Params(1, mc%index_pp: mc%index_pp+mypp_nknots-1)), r = dble(mc%params(1, mc%index_r)) )
          endif
       endif
       if(mypp_nknots .gt. 4 .and. mc%extmode .eq. 1)then
          allocate(lnk_knots(0:mypp_nknots), lnps_knots(0:mypp_nknots), cov_knots(0:mypp_nknots, 0:mypp_nknots), k_knots(0:mypp_nknots), lnps_mean_knots(0:mypp_nknots), lnps_standard_knots(0:mypp_nknots))
          call coop_set_uniform(mypp_nknots+1, lnk_knots(0:mypp_nknots), mypp_lnkmin, mypp_lnkmax)
          k_knots = exp(lnk_knots)
          cov_knots = 0.d0
          lnps_mean_knots = 0.d0
          lnps_standard_knots = mean_lnAs+(mean_ns-1.)*(lnk_knots-mypp_lnkpiv)-log(1.d10)          
       endif
       lnkmin = mypp_lnkmin
       lnkmax = mypp_lnkmax
       call coop_set_uniform(nk, lnk, lnkmin, lnkmax)
       kMpc = exp(lnk)
       standard_lnps = mean_lnAs+(mean_ns -1.)*(lnk-mypp_lnkpiv)
       if(mc%index_r .eq. 0)then
          spec_ymin = 10.
          spec_ymax = 100.
       endif
       call fig_spec%init(xlabel="$ k [{\rm Mpc}^{-1}]$", ylabel = "$10^{10}\mathcal{P}_{{\cal R},\mathrm{t}}$", xlog=.true., ylog = .true., xmin = real(exp(mypp_lnkmin-0.08)), xmax = real(exp(mypp_lnkmax + 0.08)), ymin = spec_ymin, ymax =spec_ymax, doclip = .true.)          

       call coop_asy_topaxis(fig_spec, xmin = real(exp(mypp_lnkmin-0.08))*distlss,  xmax = real(exp(mypp_lnkmax + 0.08))*distlss, islog = .true. , label = "$\ell_k\equiv  k D_{\rm rec}$")

       if(mc%extmode .eq. 1)then
          call fig_pot%init(xlabel="$(\phi - \phi_{\rm pivot})/M_p$", ylabel = "$\ln (V/V_{\rm pivot})$", xmin = -0.8, xmax = 0.5, ymin = -0.06, ymax = 0.08, doclip = .true.)
          call fig_eps%init(xlabel = "$ k [{\rm Mpc}^{-1}]$", ylabel = "$\epsilon\equiv -\dot H/H^2$", xlog = .true. ,  xmin = real(exp(mypp_lnkmin-0.08)), xmax = real(exp(mypp_lnkmax + 0.08)), ymin = -eps_ymax/10, ymax = eps_ymax, doclip = .true.)
          call coop_asy_topaxis(fig_eps, xmin = real(exp(mypp_lnkmin-0.08))*distlss,  xmax = real(exp(mypp_lnkmax + 0.08))*distlss, islog = .true. , label = "$\ell_k\equiv  k D_{\rm rec}$")             
       endif
       num_trajs = 0
       first_1sigma = .true.
       index_H = mc%usedparams%index("H0*")
       isam = 0
       do while(isam .lt. num_samples_to_get_mean)
          isam = isam + 1
          if(isam .eq. 1)then
             j = mc%ibest
             used(j) = .true.
          else
             do while(used(j))
                j = coop_random_index(mc%n)
             enddo
             used(j) = .true.
          endif
          if(mc%extmode .eq. 1)then
             if(mc%index_r .eq. 0)then
                call mypp_setup_pp(As = exp(dble(mc%Params(j, mc%index_logA)))*1.d-10, ns = mc%default_ns, nknots = mypp_nknots, dlnps = dble(mc%Params(j, mc%index_pp: mc%index_pp+mypp_nknots-1)), r = mc%default_r )
             else
                call mypp_setup_pp(As = exp(dble(mc%Params(j, mc%index_logA)))*1.d-10, ns = mc%default_ns, nknots = mypp_nknots, dlnps = dble(mc%Params(j, mc%index_pp: mc%index_pp+mypp_nknots-1)), r = dble(mc%params(j, mc%index_r)))
             endif
          else
             if(mc%index_r.eq.0)then
                call mypp_setup_pp(As = exp(dble(mc%Params(j, mc%index_logA)))*1.d-10, ns = dble(mc%Params(j, mc%index_ns)), nknots = mypp_nknots, dlnps = dble(mc%Params(j, mc%index_pp: mc%index_pp+mypp_nknots-1)), r = mc%default_r )
             else
                call mypp_setup_pp(As = exp(dble(mc%Params(j, mc%index_logA)))*1.d-10, ns = dble(mc%Params(j, mc%index_ns)), nknots = mypp_nknots, dlnps = dble(mc%Params(j, mc%index_pp: mc%index_pp+mypp_nknots-1)), r = dble(mc%params(1, mc%index_r)) )
             endif
          endif
          if(mc%extmode .eq. 1) call mypp_get_potential()
          
          !$omp parallel do 
          do ik = 1, nk
             lnps_samples(isam, ik) = mypp_primordial_lnps(kMpc(ik))
             lnpt_samples(isam, ik) = mypp_primordial_lnpt(kMpc(ik))
             if(mc%extmode .eq. 1)call mypp_get_eps_phi_lnV(kMpc(ik), eps_samples(isam, ik), phi_samples(isam, ik), lnV_samples(isam, ik))
          enddo
          !$omp end parallel do
          if(mc%extmode .eq. 1)lnV_samples(isam,:) = lnV_samples(isam,:) - mypp_lnV(mypp_ipivot)
          
          mult_samples(isam) = mc%mult(j)
          clnps = clnps + mypp_lnps*mc%mult(j)
          clnpt = clnpt + mypp_lnpt*mc%mult(j)
          if(num_trajs .lt. num_1sigma_trajs .and. mc%like(j) .lt. mc%likecut(1))then
             write(fp%unit, "("//trim(coop_num2str(mc%np))//"G14.5)") mc%params(j, :)
             num_trajs = num_trajs+1
             ps = 1.e10*exp(lnps_samples(isam, :))
             pt = 1.e10*exp(lnpt_samples(isam, :))
             ps_trajs(:, num_trajs) = ps
             pt_trajs(:, num_trajs) = pt
             if(mc%extmode .eq. 1)then             
                eps_trajs(:, num_trajs) = eps_samples(isam, :)
                phi_trajs(:, num_trajs) = phi_samples(isam, :)
                lnv_trajs(:, num_trajs) = lnV_samples(isam, :)
             endif
          endif          
       enddo
       call fp%close()

       call coop_set_uniform(nk, phi_rs, minval(phi_samples(:, 1)), maxval(phi_samples(:, nk)))
       write(*,*) "Extrapolating V(phi) for phi from "//COOP_STR_OF(phi_rs(1))//" to "//COOP_STR_OF(phi_rs(nk))
       ind_pivot = coop_minloc(abs(phi_rs))
       num_convex = 0
       if(mc%extmode .eq. 1)then       
          do isam = 1, num_samples_to_get_mean
             call coop_cheb_resample(nk, phi_samples(isam, :), lnV_samples(isam, :), phi_rs, 5)
             if(sum(lnV_samples(isam, ind_pivot-2:ind_pivot+2))/5.d0 .gt.  lnV_samples(isam, ind_pivot)) num_convex = num_convex + 1
          enddo
       endif
       !!now plot the mean
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
       mypp_lnps = clnps
       mypp_lnpt = clnpt

       
       call mypp_get_potential()
       do ik = 1, nk
          call coop_get_bounds(lnps_samples(:, ik), (/ 0.023d0, 0.1585d0, 0.5d0, 0.8415d0, 0.977d0 /), lnps_bounds(-2:2, ik), mult_samples)
          call coop_get_bounds(lnpt_samples(:, ik), (/ 0.683d0, 0.954d0 /), lnpt_bounds(1:2, ik), mult_samples)
          if(mc%extmode .eq. 1)call coop_get_bounds(eps_samples(:, ik), (/ 0.683d0, 0.955d0 /), eps_bounds(1:2, ik), mult_samples)
          if(mc%extmode .eq. 1)call coop_get_bounds(lnV_samples(:, ik), (/ 0.023d0, 0.1585d0, 0.5d0, 0.8415d0, 0.977d0 /), lnV_bounds(-2:2, ik), mult_samples)
       enddo
       lnpt_bounds(0,:) = spec_ymin
       eps_bounds(0, :) = 0.
       call fig_spec%band(kmpc, 1.d10*exp(lnps_bounds(-2,:)), 1.d10*exp(lnps_bounds(2,:)), colorfill = trim(coop_asy_gray_color(0.67)), linecolor="invisible")
       call fig_spec%band(kmpc, 1.d10*exp(lnps_bounds(-1,:)), 1.d10*exp(lnps_bounds(1,:)), colorfill = trim(coop_asy_gray_color(0.42)), linecolor="invisible")
       call fig_spec%curve(kmpc, ps_trajs(:,1), color="HEX:006FED", linetype="dashed", linewidth=0.5,legend="$\mathcal{P}_{\cal R}$ samples")       
       if(mc%index_r .ne. 0)then
          call fig_spec%band(kmpc, lnpt_bounds(0,:), 1.d10*exp(lnpt_bounds(2,:)), colorfill = trim(coop_asy_gray_color(0.67)), linecolor="invisible")
          call fig_spec%band(kmpc, lnpt_bounds(0,:), 1.d10*exp(lnpt_bounds(1,:)), colorfill = trim(coop_asy_gray_color(0.42)), linecolor="invisible")
          call fig_spec%curve(kmpc, pt_trajs(:, 1), color="HEX:8CD3F5", linetype="dotted", linewidth=0.8, legend="$\mathcal{P}_{\mathrm{t}}$ samples")          
       endif

       if(mc%extmode .eq. 1)then       
          call fig_pot%band(phi_rs, lnV_bounds(-2,:), lnV_bounds(2,:), colorfill = trim(coop_asy_gray_color(0.67)), linecolor="invisible")
          call fig_pot%band(phi_rs, lnV_bounds(-1,:), lnV_bounds(1,:), colorfill = trim(coop_asy_gray_color(0.42)), linecolor="invisible")


          call fig_eps%band(kmpc, eps_bounds(0,:),  eps_bounds(2, :), colorfill = trim(coop_asy_gray_color(0.67)), linecolor="invisible")
          call fig_eps%band(kmpc, eps_bounds(0,:), eps_bounds(1, :),  colorfill = trim(coop_asy_gray_color(0.42)), linecolor="invisible")       
       
          call fig_eps%curve(kmpc, eps_trajs(:,1), color="HEX:006FED", linetype="dotted", linewidth=0.8,legend="inflation $\epsilon$ samples")
          call fig_pot%curve(phi_trajs(:,1), lnV_trajs(:, 1), color="HEX:006FED", linetype="dotted", linewidth=1.2, legend="inflation potential samples")
       endif       
       do j=2, num_trajs
          call fig_spec%curve(kmpc, ps_trajs(:,j), color="HEX:006FED", linetype="dashed", linewidth=0.5)
          if(mc%index_r .ne. 0) call fig_spec%curve(kmpc, pt_trajs(:, j), color="HEX:8CD3F5", linetype="dotted", linewidth=0.8)
          if(mc%extmode .eq. 1)then                 
             call fig_eps%curve(kmpc, eps_trajs(:,j), color="HEX:006FED", linetype="dotted", linewidth=0.8)
             call fig_pot%curve(phi_trajs(:,j), lnV_trajs(:, j), color="HEX:006FED", linetype="dotted", linewidth=0.8)
          endif
       enddo
       call fig_spec%curve(kmpc, ps, color = "red", linetype = "solid", linewidth = 1.5, legend="mean $\mathcal{P}_{\cal R}$")
       if(mc%index_r .ne. 0) call fig_spec%curve(kmpc, pt, color = "violet", linetype = "solid", linewidth = 1.2, legend="mean $\mathcal{P}_{\mathrm{t}}$")

       if(mc%extmode .eq. 1)then              
          call fig_pot%interpolate_curve(xraw = mypp_phi, yraw = mypp_lnV-mypp_lnV(mypp_ipivot), interpolate="LinearLinear", color = "red", linetype = "solid", linewidth = 1.5, legend="inflation potential mean")
          call fig_eps%interpolate_curve(xraw = exp(mypp_lnk), yraw = exp(mypp_lneps), interpolate = "LogLinear", color = "red", linetype = "solid", linewidth = 1.5, legend="inflation $\epsilon$ mean")
       endif


       call fig_spec%curve(kMpc, exp(standard_lnps), color = "black", linewidth=1.2, legend="fiducial model $\mathcal{P}_{\cal R}$")
       if(mc%index_r .ne. 0) call fig_spec%curve(kMpc, exp(mean_lnAs - 0.01625*lnk)*0.13, color = "cyan", linewidth=1.2, legend="fiducial model $\mathcal{P}_{\mathrm{t}}$")
       if(mypp_nknots .gt. 4 .and. mypp_nknots .lt. 20 .and. mc%extmode  .eq. 1)then
          ps(1:mypp_nknots+1) = spec_ymin*1.4
          call coop_asy_dots(fig_spec, k_knots, ps(1:mypp_nknots+1), "black", "$\Delta$")
          if(mc%extmode .eq. 1)then                 
             ps(1:mypp_nknots+1) = -eps_ymax/15.
             call coop_asy_dots(fig_eps, k_knots, ps(1:mypp_nknots+1), "black", "$\Delta$")
          endif
       endif
       if(mc%extmode .eq. 1) then
          call fig_pot%label( COOP_STR_OF(mypp_nknots+1)//" knots; p(convex) = "//trim(coop_num2str(real(num_convex)/num_samples_to_get_mean,"(F10.2)")),  0.06, 0.2)

          if(trim(mc%datasets) .eq. "")then
             if(mc%index_of("r") .ne. 0)then
                call coop_asy_label(fig_spec,  "BK15 + lowl + simall + plik TTTEEE + lensing", 0.013, 6., "black")
             else
                rval = trim(mc%inputparams%value("param[r]"))
                if(trim(rval) .ne. "")then
                   call coop_asy_label(fig_spec, "fixed $r="//COOP_STR_OF(coop_str2real(rval))//"$", 0.012, 8., "black")
                endif
             endif
          else
             call coop_asy_label(fig_spec, trim(mc%datasets), 0.013, 6., "black")
             call fig_pot%label(trim(mc%datasets), 0.06, 0.1)
             call fig_eps%label(trim(mc%datasets), 0.06, 0.93)         
          endif
       endif       
       call coop_asy_legend_advance(fig_spec, real(exp(lnkmin + 0.6)), spec_ymax*0.9, "invisible", 0., 0., 0.8, 0.9, 0.9, 2)
       call fig_spec%close()       
       if(mc%extmode .eq. 1)then              
          call coop_asy_legend_advance(fig_pot, -0.4, 0.072, "invisible", 0., 0., 0.8, 0.9, 0.9, 1)
          call coop_asy_legend_advance(fig_eps,  real(exp(mypp_lnkmin +2.)), eps_ymax*0.84, "invisible", 0., 0., 0.8, 0.9, 0.9, 1)       
          call fig_eps%close()
          call fig_pot%close()
       endif

       if(mc%extmode .eq. 1)then
          !!now do eigen modes    
          do ik=1, nk
             lnpscov(ik,ik) = lnpscov(ik,ik) + 1.d-6
          enddo
          call fp%open(trim(mc%output)//"_pwtraj_eig.txt","w")
          call fp%init(xlabel="$ k [{\rm Mpc}^{-1}]$", ylabel = "$\delta \ln \Delta_S^2$", xlog = .true. , xmin = real(exp(lnkmin - 0.01)), xmax = real(exp(lnkmax + 0.01)), width = 7.2, height = 6.)
          call coop_matsym_diagonalize(nk, nk, lnpscov, lnps)
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
       endif

       if(allocated(lnk_knots))deallocate(lnk_knots, k_knots, cov_knots, lnps_knots, lnps_mean_knots, lnps_standard_knots)
    case default
       write(*,*) "pp mode ", mc%extmode, " has not been implemented"
    end select
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


!!_good.ini    
    call fp%open(trim(mc%output)//"_good.ini", "w")
    do i = 1, mc%np
       if(.not. mc%derived(i) .and. mc%vary(i))then
          write(fp%unit, "(A, 5G16.7)") "param["//trim(mc%name(i))//"] = ",  mc%mean(i), max(mc%mean(i)- 5.*mc%std(i), mc%lower(i)), min(mc%mean(i)+5.*mc%std(i), mc%upper(i)), mc%std(i)*0.8, mc%std(i)*0.8
       endif
    enddo
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

    mc%nwant = count(mc%want_1d_output)
    COOP_DEALLOC(mc%want)
    allocate(mc%want(mc%nwant))
    i=0
    do ip = 1, mc%np_used
       if(mc%want_1d_output(ip))then
          i = i+1
          mc%want(i) = ip
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
    allnames = "# "
    do ip=1, mc%nprim
       allnames = trim(allnames)//" "//trim(mc%name(mc%prim(ip)))
    enddo
    call fp%open(trim(mc%output)//".covmat", "w")    
    write(fp%unit, "(A)") trim(allnames)
    call coop_write_matrix(fp%unit, mc%covmat(mc%prim, mc%prim), mc%nprim, mc%nprim)
    call fp%close()

    allnames = "# "
    do ip=1, mc%np_used
       allnames = trim(allnames)//" "//trim(mc%name(mc%used(ip)))
    enddo
    call fp%open(trim(mc%output)//".usedcov", "w")    
    write(fp%unit, "(A)") trim(allnames)
    call coop_write_matrix(fp%unit, mc%covmat(mc%used, mc%used), mc%np_used, mc%np_used)
    call fp%close()




    call fp%open(trim(mc%output)//"_corr.txt")
    call fp%init(xlabel="", ylabel="", width=10., height=8., xmin=0., xmax=real(mc%nwant), ymin=0., ymax=real(mc%nwant), nxticks=-1, nyticks=-1)
    do i=1, mc%np_used
       if(mc%want_1d_output(i))then
          call fp%text(mc%label(mc%used(i)), i-0.5d0, -0.4d0)
          call fp%text(mc%label(mc%used(i)), -0.1d0, i-0.6d0, alignment="left")
       endif
    enddo
    call fp%density(dble(mc%corrmat(mc%want, mc%want)), 0.d0, dble(mc%nwant), 0.d0, dble(mc%nwant), "correlation", 0.d0, 1.d0)
    do i=1, mc%nwant
       do j=1, mc%nwant
          call fp%text("{\bf "//trim(coop_num2str(mc%corrmat(mc%want(i), mc%want(j)), "(F10.2)"))//"}", i-0.5d0, j-0.6d0)
       enddo
    enddo
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

    if(mc%np_pca .gt. 0)then
       call coop_feedback( "Doing PCA")
       allocate(pcamat(mc%np_pca, mc%np_pca),eig(mc%np_pca), ipca(mc%np_pca))
       pcamat = mc%covmat(mc%pca, mc%pca)
       call coop_set_uniform(mc%np_pca, ipca, 1.d0, 1.d0*mc%np_pca)
       call coop_matsym_diagonalize(mc%np_pca, mc%np_pca, pcamat, eig) 
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
       call coop_asy_curve(fp, ipca, pcamat(:,1),  color = "red", linetype = "solid", linewidth = 2., legend="$\sigma_1="//trim(coop_num2str(eig(1),"(G11.2)"))//"$")
       if(mc%np_pca .ge. 2)then
          call coop_asy_curve(fp, ipca, pcamat(:,2),  color = "blue", linetype = "dashed", linewidth = 1.5, legend = "$\sigma_2="//trim(coop_num2str(eig(2), "(G11.2)"))//"$")
       endif
       if(mc%np_pca .ge. 3)then
          call coop_asy_curve(fp, ipca, pcamat(:,3),  color = "black", linetype = "dotted", linewidth = 1., legend = "$\sigma_3="//trim(coop_num2str(eig(3), "(G11.2)"))//"$")
       endif
       if(mc%np_pca .ge. 4)then
          call coop_asy_curve(fp, ipca, pcamat(:,4),  color = "green", linetype = "dashdotted", linewidth = 0.8, legend =  "$\sigma_4="//trim(coop_num2str(eig(4), "(G11.2)"))//"$")
       endif
       if(mc%np_pca .ge. 5)then
          call coop_asy_curve(fp, ipca, pcamat(:,5),  color = "gray", linetype = "longdashdotted", linewidth = 0.8, legend =  "$\sigma_5="//trim(coop_num2str(eig(5), "(G11.2)"))//"$")
       endif
       call coop_asy_legend(fp)
       call fp%close()
       deallocate(pcamat , eig , ipca)
    endif
    do ip = 1, mc%np_used
       if(mc%want_1d_output(ip))then
          call fp%open(trim(mc%output)//"_"//trim(mc%simplename(mc%used(ip)))//"_1D.txt", "w")
          call fp%init( xlabel = trim(mc%label(mc%used(ip))), ylabel = "$P/P_{\max}$", width=6.4, height=4.8, ymin=0., ymax=1.08)          
          if(abs(mc%skewness(mc%used(ip))) .lt. mc%fit_skewness_threshold .and. abs(mc%kurtosis(mc%used(ip))) .lt. mc%fit_kurtosis_threshold .and. mc%right_is_tail(mc%used(ip)) .and. mc%left_is_tail(mc%used(ip)) )then
             call coop_asy_plot_Gaussianfit(fp, middle = mc%base(mc%used(ip)), mean =  mc%mean(mc%used(ip)), rms = mc%std(mc%used(ip)), skewness = mc%skewness(mc%used(ip)), lbd=mc%lower(mc%used(ip)), ubd = mc%upper(mc%used(ip)))
          else
             write(*,*) trim(mc%name(mc%used(ip))), " skewness = ", mc%skewness(mc%used(ip)), "; kurtosis = ", mc%kurtosis(mc%used(ip))             
             do i = 1, mc%nb
                x(i) = mc%plotlower(mc%used(ip)) + mc%dx(mc%used(ip))*(i-1)
             enddo
             call coop_asy_plot_likelihood(fp, x, mc%c1d(:, ip)/maxval(mc%c1d(:, ip)), left_tail = .true., right_tail = .true., linewidth=1.5)
          endif
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
                   call path%from_array_gaussianfit(mc%c2d(:, :, k), mc%plotlower(mc%used(j)), mc%plotupper(mc%used(j)), mc%plotlower(mc%used(j2)), mc%plotupper(mc%used(j2)), mc%cut2d(icontour, k))

                   call coop_asy_contour(fp, path, colorfill = trim(mc%color2d(icontour)),  linecolor = "black", linetype = "solid")
                enddo
                call fp%close()
             endif
             if(j2 .ne. j .and. mc%want_2d_output(j2, j))then
                call fp%open(trim(mc%output)//"_"//trim(mc%simplename(mc%used(j2)))//"_"//trim(mc%simplename(mc%used(j)))//"_2D.txt", "w")
                call fp%init( xlabel = trim(mc%label(mc%used(j2))), ylabel = trim(mc%label(mc%used(j))), xmin=mc%plotlower(mc%used(j2)), xmax = mc%plotupper(mc%used(j2)), ymin=mc%plotlower(mc%used(j)), ymax = mc%plotupper(mc%used(j)) )
                do icontour = coop_postprocess_num_contours, 1, -1
                   call path%from_array_gaussianfit(transpose(mc%c2d(:, :, k)), mc%plotlower(mc%used(j2)), mc%plotupper(mc%used(j2)),  mc%plotlower(mc%used(j)), mc%plotupper(mc%used(j)), mc%cut2d(icontour, k))

                   call coop_asy_contour(fp, path, colorfill = trim(mc%color2d(icontour)),  linecolor = "black", linetype = "solid")
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
    COOP_INT nparams    
    COOP_INT i, ind
    nparams = size(Params)
    do i = 1, nparams
       params(i) = this%name2value(ind,  trim(this%allparams%key(i)))
    enddo
  end subroutine coop_mcmc_chain_getCosmoMCParams


  function coop_mcmc_chain_all_index_of(mc, name) result(ind)
    class(coop_mcmc_chain) mc
    COOP_INT ind
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
    COOP_INT isample
    COOP_UNKNOWN_STRING name
    COOP_SINGLE val
    COOP_INT ind
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
    call coop_feedback("preprocessing...")
    if(this%index_H0 .ne. 0)then
       do i=1, this%n
          if(this%params(i, this%index_H0) <55 .or. this%params(i, this%index_H0) >80)then
             this%params(i, this%index_H0) = 70.
             this%mult(i) = 0.
          endif
       enddo
    endif
       
  end subroutine coop_mcmc_chain_preprocess





End Module Coop_statchains_mod
