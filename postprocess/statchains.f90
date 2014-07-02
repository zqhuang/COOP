module statchains
  use coop_wrapper_utils
  use cosmolibwrap

  implicit none

#include "constants.h"

  integer,parameter::mcmc_stat_num_cls = 3
  real,dimension(mcmc_stat_num_cls)::mcmc_stat_cls = (/ 0.683, 0.954, 0.997 /)
  COOP_STRING :: measured_cltt_file = ""
  COOP_STRING :: measured_clee_file = ""
  COOP_STRING :: measured_clbb_file = ""
  COOP_STRING :: measured_clte_file = ""
  COOP_STRING :: bestfit_cl_file = ""


  type MCMC_chain
     COOP_STRING prefix, output
     real bestlike, worstlike, totalmult, likecut(mcmc_stat_num_cls)
     integer ibest, iworst
     integer np
     COOP_SHORT_STRING,dimension(:),allocatable::label, name, simplename
     real,dimension(:),allocatable::lower, upper, mean, std, plotlower, plotupper, dx, base !!plotlower and plotupper are used to make plots, with ~0 long tail truncated
     real,dimension(:,:),allocatable::lowsig, upsig
     logical,dimension(:),allocatable::vary, left_is_tail, right_is_tail
     integer,dimension(:),allocatable::map2used
     real,dimension(:,:),allocatable::covmat, corrmat, cov_used

     integer np_used
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
     COOP_SHORT_STRING Color2D_light, color2d_dark
     type(coop_dictionary) inputparams
  end type MCMC_chain

contains

  subroutine load_chain(mc, prefix, ignore_percent)
    type(MCMC_chain) mc
    COOP_UNKNOWN_STRING prefix
    integer,optional::ignore_percent
    integer nfiles, i, j, lens, k
    integer,dimension(:),allocatable::nskip, nlines
    COOP_LONG_STRING fname
    COOP_STRING inline
    integer ispace
    type(coop_file) fp
    mc%prefix = trim(prefix)
    nfiles = 0
    do 
       fname = trim(prefix)//"_"//trim(coop_num2str(nfiles+1))//".txt"
       if(coop_file_exists(fname))then
          nfiles = nfiles+1
          if(nfiles.eq.1)then
             mc%np = coop_file_numcolumns(fname) -2 
             if(allocated(mc%label))deallocate(mc%name, mc%simplename, mc%label, mc%std, mc%mean, mc%lower, mc%upper, mc%vary, mc%covmat,  mc%plotlower, mc%plotupper, mc%dx, mc%map2used, mc%lowsig, mc%upsig, mc%left_is_tail, mc%right_is_tail, mc%base)
             allocate(mc%label(mc%np), mc%simplename(mc%np), mc%std(mc%np), mc%mean(mc%np), mc%lower(mc%np), mc%upper(mc%np),  mc%plotlower(mc%np), mc%plotupper(mc%np), mc%vary(mc%np), mc%covmat(mc%np, mc%np), mc%dx(mc%np), mc%map2used(mc%np), mc%name(mc%np), mc%lowsig(mcmc_stat_num_cls, mc%np), mc%upsig(mcmc_stat_num_cls, mc%np), mc%left_is_tail(mc%np), mc%right_is_tail(mc%np), mc%base(mc%np))
          endif
       else
          exit
       endif
    enddo
    if(nfiles.eq.0)then
       write(*,*) trim(prefix)//"_1.txt is not found on the disk"
       stop
    else
       write(*,*) "found "//trim(coop_num2str(nfiles))//" chain files on the disck"
    endif
    fname = trim(prefix)//".inputparams"
    if(coop_file_exists(fname))then
       call coop_load_dictionary(trim(fname), mc%inputparams)
       pp_mode = mc%inputparams%value("primordial_power_mode")
       if(trim(pp_mode).eq."elldip")then
          call coop_dictionary_lookup(mc%inputparams, "pp_dip_lmin", pp_dip_lmin, 20)
          call coop_dictionary_lookup(mc%inputparams, "pp_dip_lmax", pp_dip_lmax, 30)
       endif
    else
       write(*,*) "Warning: inputparams file not found; cannot set default values"
       pp_mode = "standard"
    endif
    fname = trim(prefix)//".paramnames"
    if(coop_file_exists(fname))then
       call fp%open(fname, "r")
       do i =1, mc%np
          if(Coop_file_readOneLine(fp, inline))then
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
       call close_file(fp)
    else
       do i = 1, mc%np
          mc%label(i) = "param"//trim(coop_num2str(i))
          mc%name(i) = mc%label(i)
       enddo
    endif
    do i = 1, mc%np
       mc%simplename(i) = trim(coop_str_numalpha(mc%name(i)))
    enddo
    allocate(nskip(nfiles),nlines(nfiles))
    do i = 1, nfiles
       fname = trim(prefix)//"_"//trim(coop_num2str(i))//".txt"
       nlines(i)= coop_file_numlines(fname)          
       write(*,*) "found "//trim(coop_num2str(nlines(i)))//" lines in "//trim(fname)
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
          call coop_file_skiplines(fp, nskip(i))
          do j=nskip(i)+1, nlines(i)
             k = k + 1
             read(fp%unit, *) mc%mult(k), mc%like(k), mc%params(k,1:mc%np)
          enddo
          call close_file(fp)
       endif
    enddo
    if( k .eq. mc%n)then
       write(*,*) "Totally "//trim(coop_num2str(k))//" samples are used."
    else
       stop "Error: counting failed in load_chain"
    endif
    deallocate(nskip, nlines)
    write(*, *) "Analysing the chain ... "
    call analyze_chain(mc)
    write(*, *) "Chain analysed."
  end subroutine load_chain

  subroutine analyze_chain(mc)
    integer,parameter::n_fine_bins = 2048
    real c(n_fine_bins), dx,  multcut, acc
    real,dimension(:),allocatable::c2dlist
    type(mcmc_chain) mc
    integer ip, i, loc, j, j2, ip2, k, loc2, icl
    mc%totalmult = sum(mc%mult)
    if(mc%totalmult .le. 0 .or. mc%n.eq.0) stop "analyze_chains: found no samples"

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
       
       if (mc%upper(ip) .eq. mc%lower(ip)  )then
          mc%mean(ip) = mc%upper(ip)
          mc%std(ip) = 0.
          mc%vary(ip) = .false.
          mc%plotupper(ip) = mc%upper(ip)
          mc%plotlower(ip) = mc%lower(ip)
          mc%dx(ip) = 0.
          mc%lowsig(:, ip) = mc%upper(ip)
          mc%upsig(:, ip) = mc%upper(ip)
          mc%left_is_tail(ip) = .false.
          mc%right_is_tail(ip) = .false.
       else
          mc%mean(ip) = sum(mc%params(:,ip)*mc%mult)/mc%totalmult
          mc%std(ip) = sqrt(sum((mc%params(:,ip)-mc%mean(ip))**2*mc%mult)/mc%totalmult)
          mc%vary(ip) = .true.
          dx = (mc%upper(ip) - mc%lower(ip))/n_fine_bins*1.e-6
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
          multcut = mc%totalmult*max(min(sqrt(0.03/mc%n), 0.005), 0.001)
          do while(acc.lt.multcut)
             i = i + 1
             acc = acc + c(i)
          enddo
          if(i.gt.1)then
             acc = acc - c(i)
             i = i - 1
          endif
          if( acc/i  .ge. mc%totalmult/n_fine_bins )then  
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
          if(acc/(n_fine_bins-i+1) .ge. mc%totalmult/n_fine_bins)then
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
    mc%nb = min(max(13, ceiling(sqrt(mc%n/100.))), 30)
    if(allocated(mc%c1d))deallocate(mc%c1d, mc%c2d, mc%cut2d, mc%want_2d_output)
    allocate(mc%c1d(mc%nb, mc%np_used), mc%c2d(mc%nb, mc%nb, mc%np_used*(mc%np_used+1)/2), mc%cut2d(mcmc_stat_num_cls, mc%np_used*(mc%np_used+1)/2), mc%want_2d_output(mc%np_used, mc%np_used) )
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
          call quicksort(c2dlist)
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
  end subroutine analyze_chain


  subroutine export_stats(mc, output)
    integer::num_1sigma_trajs
    integer::num_samples_to_get_mean

    integer,parameter::nk = 71
    COOP_REAL  norm
    logical,parameter::do_traj_cov = .true.
    integer,parameter::distlss = 13893.
    integer,parameter::index_TT = 1, index_EE = 2, index_BB = 3, index_TE = 4
    integer,parameter:: lmax = 2500
    integer, parameter::dcl_lmax = 300
    real, parameter::dcl_legend_loc = sqrt(dcl_lmax*2.)
    COOP_REAL  Cls(4, 2:lmax), rl(2:lmax), cls_mean(4, 2:lmax),  maxcls(4), mincls(4), bestCls(4, 2:dcl_lmax), maxdcls(4), mindcls(4)
    
    COOP_REAL  lnk(nk), ps(nk), pt(nk), lnpsmean(nk), lnptmean(nk), lnpscov(nk, nk), lnps(nk), lnpt(nk), lnps_shift, lnpt_shift, kmpc(nk),  mineig !, lnptcov(nk, nk)
    type(mcmc_chain) mc
    COOP_UNKNOWN_STRING output
    type(coop_asy) fp, fp2, fpv, fpeps
    type(coop_asy) fpclTT, fpclEE, fpclBB, fpclTE
    type(coop_asy) fpdclTT, fpbest, fpdclEE, fpdclBB, fpdclTE
    integer i , ip, j,  j2, k, ik, numpiv, ltmp, ik2, ndof
    real mult
    real x(mc%nb), lnkmin, lnkmax
    logical do_pp, do_cl
    type(coop_asy_path) path
    logical first_1sigma
    COOP_REAL ,dimension(:,:),allocatable::pcamat
    COOP_REAL ,dimension(:),allocatable::eig, ipca, initpower
    integer num_initpower, numpp
    logical inflation_consistency
    real kpiv, ytop
    COOP_REAL  :: cltt, errup, errdown
    integer junk, l
    real cltraj_weight
    COOP_STRING inline

    do i = 2, lmax
       rl(i) = i
    enddo
    mc%output = trim(adjustl(output))//trim(coop_file_name_of(mc%prefix))
    call fp%open(trim(mc%output)//".likes", "w")
    write(fp%unit, "(A, G14.5)") "Best -lnlike = ", mc%bestlike
    write(fp%unit, "(A, G14.5)") "Worst -lnlike = ", mc%worstlike
    do j = 1, mcmc_stat_num_cls
       write(fp%unit, "(A, F14.3, A, G14.5)") "prob = ", mcmc_stat_cls(j)," truncation like = ", mc%likecut(j)
    enddo
    call close_file(fp)
    do_pp  = (trim(pp_mode) .ne. "")
    !! =================================================================!!
    call coop_dictionary_lookup(mc%inputparams, "num_init_power", num_initpower, 6)
    call coop_dictionary_lookup(mc%inputparams, "inflation_consistency", inflation_consistency, .true.)
    call coop_dictionary_lookup(mc%inputparams, "pivot_k", kpiv, 0.05)

    allocate(initpower(num_initpower))
    !!------------------------------------------------------------------!!
    call fp%open(trim(mc%output)//"_1sig.samples", "w")
    write(fp%unit, "("//trim(coop_num2str(mc%np))//"G14.5)") mc%params(mc%ibest, :)
    do_cl = (trim(measured_cltt_file).ne."" .or.trim(measured_clee_file).ne."" .or. trim(measured_clbb_file).ne."" .or. trim(measured_clte_file).ne."")
    if(do_pp)then
       write(*,*) "Generating primordial power spectra trajectories"
       select case(trim(pp_mode))
       case("bump")
          num_1sigma_trajs  = 70
          num_samples_to_get_mean = 5
       case default
          num_1sigma_trajs  = 80
          num_samples_to_get_mean = 5000
       end select
       if(do_cl) num_1sigma_trajs = num_1sigma_trajs  * 4/5

       if(do_cl)then
          if(trim(bestfit_cl_file) .ne. "")then
             call fpbest%open(trim(bestfit_cl_file), "r")
             do l=2, dcl_lmax
                read(fpbest%unit, *) ltmp, bestCls(:, l)
             enddo
             call close_file(fpbest)
          else
             bestCls = 0.d0
          endif
          write(*,*) "Generating Cl trajectories"
          call fpdclTT%open(trim(mc%output)//"_dclTT_trajs.txt", "w")
          call fpdclEE%open(trim(mc%output)//"_dclEE_trajs.txt", "w")
          call fpdclBB%open(trim(mc%output)//"_dclBB_trajs.txt", "w")
          call fpdclTE%open(trim(mc%output)//"_dclTE_trajs.txt", "w")
          call fpclTT%open(trim(mc%output)//"_clTT_trajs.txt", "w")
          call fpclEE%open(trim(mc%output)//"_clEE_trajs.txt", "w")
          call fpclBB%open(trim(mc%output)//"_clBB_trajs.txt", "w")
          call fpclTE%open(trim(mc%output)//"_clTE_trajs.txt", "w")
       endif
       call fp2%open(trim(mc%output)//"_power_trajs.txt", "w")
       call fpv%open(trim(mc%output)//"_potential_trajs.txt", "w")
       call fpeps%open(trim(mc%output)//"_epsilon_trajs.txt", "w")
       call get_initpower(mc, mc%ibest, initpower, num_initpower)
       call  pp_cosmomc_init(initpower, inflation_consistency = inflation_consistency, k_pivot = kpiv)
       lnkmin = pp_lnk_min
       lnkmax = pp_lnk_max
       call coop_set_uniform(nk, lnk, dble(lnkmin), dble(lnkmax))
       kmpc = exp(lnk)*kpiv
       if(do_cl)then
          call fpdclTT%init( xlabel="$\ell$", ylabel = "$\frac{\ell (\ell + 1)}{2\pi}\Delta C_\ell^{TT}$", xlog=.true., xmin = 1.8, xmax = real(dcl_lmax+1))
          call fpdclEE%init(xlabel="$\ell$", ylabel = "$\frac{\ell (\ell + 1)}{2\pi}\Delta C_\ell^{EE}$", xlog=.true., xmin = 1.8, xmax = real(dcl_lmax+1))
          call fpdclTE%init( xlabel="$\ell$", ylabel = "$\frac{\ell (\ell + 1)}{2\pi}\Delta C_\ell^{TE}$", xlog=.true., xmin = 1.8, xmax = real(dcl_lmax+1))
          call fpdclBB%init( xlabel="$\ell$", ylabel = "$\frac{\ell (\ell + 1)}{2\pi}\Delta C_\ell^{BB}$", xlog=.true., xmin = 1.8, xmax = real(dcl_lmax+1))
          call fpclTT%init( xlabel="$\ell$", ylabel = "$\frac{\ell (\ell + 1)}{2\pi}C_\ell^{TT}$", xlog=.true., xmin = 1.9, xmax = real(lmax))
          call fpclEE%init( xlabel="$\ell$", ylabel = "$\frac{\ell (\ell + 1)}{2\pi}C_\ell^{EE}$", xlog=.true., xmin = 1.9, xmax = real(lmax))
          call fpclBB%init( xlabel="$\ell$", ylabel = "$\frac{\ell (\ell + 1)}{2\pi}C_\ell^{BB}$", xlog=.true., xmin = 1.9, xmax = real(lmax))
          call fpclTE%init( xlabel="$\ell$", ylabel = "$\frac{\ell (\ell + 1)}{2\pi}C_\ell^{TE}$", xlog=.true., xmin = 1.9, xmax = real(lmax))
       endif

       call fp2%init(xlabel="$ k ({\rm Mpc}^{-1})$", ylabel = "$10^{10}\Delta_{S,T}^2$", xlog=.true., ylog = .true., xmin = exp(lnkmin - 0.01)*kpiv, xmax = kpiv*exp(lnkmax + 0.01), ymin = 1., ymax = 200., doclip = .true.)
       call coop_asy_topaxis(fp2, xmin =  exp(lnkmin - 0.01)*kpiv*distlss, xmax = kpiv*exp(lnkmax + 0.01)*distlss, islog = .true., label = "$\ell\equiv  k D_{\rm rec}$")
       call fpv%init(xlabel = "$ (\phi - \phi_{\rm pivot})/M_p$", ylabel = "$ \ln(V/V_{\rm pivot})$", xmin = -1.5, xmax = 0.5, ymin = -0.2, ymax=0.6, doclip = .true.)
       call fpeps%init(xlabel="$ k ({\rm Mpc}^{-1})$", ylabel = "$\epsilon$", xlog = .true. , xmin = exp(lnkmin - 0.01)*kpiv, xmax = kpiv*exp(lnkmax + 0.01), ymin = 0., ymax = 0.15, doclip = .true.)
       call coop_asy_topaxis(fpeps, xmin =  exp(lnkmin - 0.01)*kpiv*distlss, xmax = kpiv*exp(lnkmax + 0.01)*distlss, islog = .true., label = "$\ell\equiv  k D_{\rm rec}$")



       lnpsmean = 0
       lnptmean = 0
       if(do_traj_cov)then
          lnpscov = 0
!          lnptcov = 0
       endif
       mult = 0
      
       do j = 1, mc%n, max(mc%n/num_samples_to_get_mean, 1)
          call get_initpower(mc, j, initpower, num_initpower)
          call  pp_cosmomc_init(initpower, inflation_consistency = inflation_consistency, k_pivot = kpiv)
          !$omp parallel do 
          do ik = 1, nk
             lnps(ik) = log(pp_scalar_power(lnk(ik))) 
             lnpt(ik) = log(pp_tensor_power(lnk(ik))) 
          enddo
          !$omp end parallel do
          if(j.eq.1)then
             lnps_shift = - sum(lnps)/nk
             lnpt_shift = - sum(lnpt)/nk
          endif
          lnps = lnps + lnps_shift
          lnpt = lnpt + lnpt_shift
          lnpsmean = lnpsmean + lnps * mc%mult(j)
          lnptmean = lnptmean + lnpt * mc%mult(j)
          if(do_traj_cov)then
             do ik=1, nk
                do ik2 = ik, nk
                   lnpscov(ik, ik2) = lnpscov(ik, ik2) + lnps(ik) * lnps(ik2) * mc%mult(j)
!                   lnptcov(ik, ik2) = lnptcov(ik, ik2) + lnpt(ik) * lnpt(ik2) * mc%mult(j)
                enddo
             enddo
          endif
          mult = mult + mc%mult(j)
       enddo
    endif
    if(do_cl)then
       maxcls = -1.e30
       maxdcls = -1.e30
       mincls = 1.e30
       mindcls = 1.e30
       cls_mean = 0.
       cltraj_weight = 0.
    endif
    first_1sigma = .true.
    do i = 1, num_1sigma_trajs 
       j = coop_random_index(mc%n)
       if(do_pp)then
          call get_initpower(mc, j, initpower, num_initpower)
          call  pp_cosmomc_init(initpower, inflation_consistency = inflation_consistency, k_pivot = kpiv)
          !$omp parallel do 
          do ik = 1, nk
             lnps(ik) = log(pp_scalar_power(lnk(ik))) * mc%mult(j)
             lnpt(ik) = log(pp_tensor_power(lnk(ik))) * mc%mult(j)
          enddo
          !$omp end parallel do
          lnps = lnps + lnps_shift
          lnpt = lnpt + lnpt_shift
          lnpsmean = lnpsmean + lnps * mc%mult(j)
          lnptmean = lnptmean + lnpt * mc%mult(j)
          if(do_traj_cov)then
             do ik=1, nk
                do ik2 = ik, nk
                   lnpscov(ik, ik2) = lnpscov(ik, ik2) + lnps(ik) * lnps(ik2) * mc%mult(j)
!                   lnptcov(ik, ik2) = lnptcov(ik, ik2) + lnpt(ik) * lnpt(ik2) * mc%mult(j)
                enddo
             enddo
          endif
          mult = mult + mc%mult(j)
          if(do_cl)then
             call get_camb_cls(mc, j, lmax, cls)
             cls_mean = cls_mean + cls*mc%mult(j)
             cltraj_weight = cltraj_weight + mc%mult(j)
          endif
          if(trim(pp_mode).eq."bump" .and. mod(i, 5).eq. 0) write(*,*) i, " trajs tried"
       endif
       if(mc%like(j) .lt. mc%likecut(1))then
          write(fp%unit, "("//trim(coop_num2str(mc%np))//"G14.5)") mc%params(j, :)
          if(do_pp)then
             if(do_cl)then
                do k=1,4
                   maxcls(k) = max(maxval(cls(k,:)), maxcls(k))
                   maxdcls(k) = max(maxval(cls(k,2:dcl_lmax) - bestcls(k, 2:dcl_lmax)), maxdcls(k))
                   mincls(k) = min(minval(cls(k,:)), mincls(k))
                   mindcls(k) = min(minval(cls(k,2:dcl_lmax)-bestcls(k, 2:dcl_lmax)), mindcls(k))
                enddo
             endif
             do ik = 1, nk
                ps(ik) = 1.e10*pp_scalar_power(lnk(ik))
                pt(ik) = 1.e10*pp_tensor_power(lnk(ik))
             enddo
             call pp_get_potential()
             if(first_1sigma)then
                first_1sigma = .false.
                if(do_cl)then
                   call coop_asy_interpolate_curve(fpdclTT, rl(2:dcl_lmax), cls(index_TT, 2:dcl_lmax)-bestCls(index_TT, 2:dcl_lmax), interpolate="LogLinear", color = "blue", linetype = "dotted", legend="1-$\sigma$ trajs.")
                   call coop_asy_interpolate_curve(fpdclEE, rl(2:dcl_lmax), cls(index_EE, 2:dcl_lmax)-bestCls(index_EE, 2:dcl_lmax), interpolate="LogLinear", color = "blue", linetype = "dotted", legend="1-$\sigma$ trajs.")
                   call coop_asy_interpolate_curve(fpdclBB, rl(2:dcl_lmax), cls(index_BB, 2:dcl_lmax)-bestCls(index_BB, 2:dcl_lmax), interpolate="LogLinear", color = "blue", linetype = "dotted", legend="1-$\sigma$ trajs.")
                   call coop_asy_interpolate_curve(fpdclTE, rl(2:dcl_lmax), cls(index_TE, 2:dcl_lmax)-bestCls(index_TE, 2:dcl_lmax), interpolate="LogLinear", color = "blue", linetype = "dotted", legend="1-$\sigma$ trajs.")
                   call coop_asy_interpolate_curve(fpclTT, rl, cls(index_TT, :), interpolate="LogLinear", color = "blue", linetype = "dotted", legend="1-$\sigma$ trajs.")
                   call coop_asy_interpolate_curve(fpclEE, rl, cls(index_EE, :), interpolate="LogLinear", color = "blue", linetype = "dotted", legend="1-$\sigma$ trajs.")
                   call coop_asy_interpolate_curve(fpclBB, rl, cls(index_BB, :), interpolate="LogLinear", color = "blue", linetype = "dotted", legend="1-$\sigma$ trajs.")
                   call coop_asy_interpolate_curve(fpclTE, rl, cls(index_TE, :), interpolate="LogLinear", color = "blue", linetype = "dotted", legend="1-$\sigma$ trajs.")
                endif
                call coop_asy_curve(fp2, kmpc, ps, smooth = .false., color = "blue", linetype = "dashed", legend="1-$\sigma$ trajs.")
                call coop_asy_curve(fp2, kmpc, pt, smooth = .false., color = "skyblue", linetype = "dotted")
                call coop_asy_interpolate_curve(fpv, pp_phi, pp_lnV, interpolate = "LinearLinear", color = "blue", linetype = "dashed", legend="1-$\sigma$ trajs.")
                call coop_asy_interpolate_curve(fpeps, kpiv*exp(pp_lnk), pp_epsilon, interpolate="LogLinear", color = "blue", linetype = "dashed", legend="1-$\sigma$ trajs.")
             else
                if(do_cl)then
                   call coop_asy_interpolate_curve(fpdclTT, rl(2:dcl_lmax), cls(index_TT, 2:dcl_lmax)-bestCls(index_TT, 2:dcl_lmax), interpolate="LogLinear", color = "blue", linetype = "dotted")
                   call coop_asy_interpolate_curve(fpdclEE, rl(2:dcl_lmax), cls(index_EE, 2:dcl_lmax)-bestCls(index_EE, 2:dcl_lmax), interpolate="LogLinear", color = "blue", linetype = "dotted")
                   call coop_asy_interpolate_curve(fpdclBB, rl(2:dcl_lmax), cls(index_BB, 2:dcl_lmax)-bestCls(index_BB, 2:dcl_lmax), interpolate="LogLinear", color = "blue", linetype = "dotted")
                   call coop_asy_interpolate_curve(fpdclTE, rl(2:dcl_lmax), cls(index_TE, 2:dcl_lmax)-bestCls(index_TE, 2:dcl_lmax), interpolate="LogLinear", color = "blue", linetype = "dotted")
                   call coop_asy_interpolate_curve(fpclTT, rl, cls(index_TT, :), interpolate="LogLinear", color = "blue", linetype = "dotted")
                   call coop_asy_interpolate_curve(fpclEE, rl, cls(index_EE, :), interpolate="LogLinear", color = "blue", linetype = "dotted")
                   call coop_asy_interpolate_curve(fpclBB, rl, cls(index_BB, :), interpolate="LogLinear", color = "blue", linetype = "dotted")
                   call coop_asy_interpolate_curve(fpclTE, rl, cls(index_TE, :), interpolate="LogLinear", color = "blue", linetype = "dotted")
                endif
                call coop_asy_curve(fp2, kmpc, ps, smooth = .false., color = "blue", linetype = "dashed")
                call coop_asy_curve(fp2, kmpc, pt, smooth = .false., color = "skyblue", linetype = "dotted")
                call coop_asy_interpolate_curve(fpv, pp_phi, pp_lnV, interpolate = "LinearLinear", color = "blue", linetype = "dashed")
                call coop_asy_interpolate_curve(fpeps, kpiv*exp(pp_lnk), pp_epsilon, interpolate="LogLinear", color = "blue", linetype = "dashed")
             endif
          endif
       endif
    enddo
    call close_file(fp)

    !!now plot the mean
    if(do_pp)then
       lnpsmean = lnpsmean / mult 
       lnptmean = lnptmean / mult
       lnpscov = lnpscov/mult
!       lnptcov = lnptcov/mult
       if(do_traj_cov)then
          do ik=1, nk
             do ik2 = ik, nk
                lnpscov(ik, ik2) = lnpscov(ik, ik2) - lnpsmean(ik)*lnpsmean(ik2)
!                lnptcov(ik, ik2) = lnptcov(ik, ik2) - lnptmean(ik)*lnptmean(ik2)
             enddo
          enddo
          do ik=1,nk
             do ik2 = 1, ik-1
                lnpscov(ik, ik2) = lnpscov(ik2, ik) 
!                lnptcov(ik, ik2) = lnptcov(ik2, ik)
             enddo
          enddo
       endif
       lnpsmean = lnpsmean - lnps_shift
       lnptmean = lnptmean - lnpt_shift
       if(do_cl)then
          cls_mean = cls_mean / cltraj_weight   
          call coop_asy_interpolate_curve(fpdclTT, rl(2:dcl_lmax), cls_mean(index_TT, 2:dcl_lmax)-bestCls(index_TT, 2:dcl_lmax), interpolate="LogLinear", color = "red", linetype = "solid", linewidth=1.2, legend="mean traj")
          call coop_asy_interpolate_curve(fpdclEE, rl(2:dcl_lmax), cls_mean(index_EE, 2:dcl_lmax)-bestCls(index_EE, 2:dcl_lmax), interpolate="LogLinear", color = "red", linetype = "solid", linewidth=1.2, legend="mean traj")
          call coop_asy_interpolate_curve(fpdclBB, rl(2:dcl_lmax), cls_mean(index_BB, 2:dcl_lmax)-bestCls(index_BB, 2:dcl_lmax), interpolate="LogLinear", color = "red", linetype = "solid", linewidth=1.2, legend="mean traj")
          call coop_asy_interpolate_curve(fpdclTE, rl(2:dcl_lmax), cls_mean(index_TE, 2:dcl_lmax)-bestCls(index_TE, 2:dcl_lmax), interpolate="LogLinear", color = "red", linetype = "solid", linewidth=1.2, legend="mean traj")

          call coop_asy_interpolate_curve(fpclTT, rl, cls_mean(index_TT, :), interpolate="LogLinear", color = "red", linetype = "solid", linewidth=1.2, legend="mean traj")
          call coop_asy_interpolate_curve(fpclEE, rl, cls_mean(index_EE, :), interpolate="LogLinear", color = "red", linetype = "solid", linewidth=1.2, legend="mean traj")
          call coop_asy_interpolate_curve(fpclBB, rl, cls_mean(index_BB, :), interpolate="LogLinear", color = "red", linetype = "solid", linewidth=1.2, legend="mean traj")
          call coop_asy_interpolate_curve(fpclTE, rl, cls_mean(index_TE, :), interpolate="LogLinear", color = "red", linetype = "solid", linewidth=1.2, legend="mean traj")

          if(trim(bestfit_cl_file) .ne. "")then
             call coop_asy_curve_from_file(fpclTT, bestfit_cl_file, xcol = 1, ycol = index_TT + 1, interpolate = "LogLinear", color = "green", linewidth= 0.5, legend = "$\Lambda$CDM best-fit")
             call coop_asy_curve_from_file(fpclTE, bestfit_cl_file, xcol = 1, ycol = index_TE + 1, interpolate = "LogLinear", color = "green", linewidth= 0.5, legend = "$\Lambda$CDM best-fit")
             call coop_asy_curve_from_file(fpclEE, bestfit_cl_file, xcol = 1, ycol = index_EE + 1, interpolate = "LogLinear", color = "green", linewidth= 0.5, legend = "$\Lambda$CDM best-fit")
             call coop_asy_curve_from_file(fpclBB, bestfit_cl_file, xcol = 1, ycol = index_BB + 1, interpolate = "LogLinear", color = "green", linewidth= 0.5, legend = "$\Lambda$CDM best-fit")
          endif
          if(trim(measured_cltt_file).ne."")then
             call fp%open(measured_cltt_file, "r")
             do 
                if( coop_file_readoneline_string(fp, inline) )then
                   read(inline, *) l, junk, junk, cltt, errup, errdown
                else
                   exit
                endif
                maxcls(index_TT) = max(maxcls(index_TT), cltt+errup)
                mincls(index_TT) = min(mincls(index_TT), cltt-errdown)
                if(l.le. dcl_lmax)then
                   maxdcls(index_TT) = max(maxdcls(index_TT), cltt+errup-bestcls(index_TT, l))
                   mindcls(index_TT) = min(mindcls(index_TT), cltt-errdown-bestcls(index_TT, l))
                   call coop_asy_error_bar(fpdclTT, dble(l), cltt-bestcls(index_TT,l), dy_plus = errup, dy_minus = errdown)
                endif
                call coop_asy_error_bar(fpclTT, dble(l), cltt, dy_plus = errup, dy_minus = errdown)
             enddo
             call close_file(fp)
          endif
          if(trim(measured_clte_file).ne."")then
             call fp%open(measured_clte_file, "r")
             do 
                if( coop_file_readoneline_string(fp, inline) )then
                   read(inline, *) l, junk, junk, cltt, errup, errdown
                else
                   exit
                endif
                if(l.le. dcl_lmax)then
                   maxdcls(index_TE) = max(maxdcls(index_TE), cltt+errup-bestcls(index_TE, l))
                   mindcls(index_TE) = min(mindcls(index_TE), cltt-errdown-bestcls(index_TE, l))
                   call coop_asy_error_bar(fpdclTE, dble(l), cltt-bestcls(index_TE,l), dy_plus = errup, dy_minus = errdown)
                endif

                maxcls(index_TE) = max(maxcls(index_TE), cltt+errup)
                mincls(index_TE) = min(mincls(index_TE), cltt-errdown)
                call coop_asy_error_bar(fpclTE, dble(l), cltt, dy_plus = errup, dy_minus = errdown)
             enddo
             call close_file(fp)
          endif
          if(trim(measured_clEE_file).ne."")then
             call fp%open(measured_clEE_file, "r")
             do 
                if( coop_file_readoneline_string(fp, inline) )then
                   read(inline, *) l, junk, junk, cltt, errup, errdown
                else
                   exit
                endif
                if(l.le. dcl_lmax)then
                   maxdcls(index_EE) = max(maxdcls(index_EE), cltt+errup-bestcls(index_EE, l))
                   mindcls(index_EE) = min(mindcls(index_EE), cltt-errdown-bestcls(index_EE, l))
                   call coop_asy_error_bar(fpdclEE, dble(l), cltt-bestcls(index_EE,l), dy_plus = errup, dy_minus = errdown)
                endif

                maxcls(index_EE) = max(maxcls(index_EE), cltt+errup)
                mincls(index_EE) = min(mincls(index_EE), cltt-errdown)
                call coop_asy_error_bar(fpclEE, dble(l), cltt, dy_plus = errup, dy_minus = errdown)
             enddo
             call close_file(fp)
          endif
          if(trim(measured_clBB_file).ne."")then
             call fp%open(measured_clBB_file, "r")
             do 
                if( coop_file_readoneline_string(fp, inline) )then
                   read(inline, *) l, junk, junk, cltt, errup, errdown
                else
                   exit
                endif
                if(l.le. dcl_lmax)then
                   maxdcls(index_BB) = max(maxdcls(index_BB), cltt+errup-bestcls(index_BB, l))
                   mindcls(index_BB) = min(mindcls(index_BB), cltt-errdown-bestcls(index_BB, l))
                   call coop_asy_error_bar(fpdclBB, dble(l), cltt-bestcls(index_BB,l), dy_plus = errup, dy_minus = errdown)
                endif

                maxcls(index_BB) = max(maxcls(index_BB), cltt+errup)
                mincls(index_BB) = min(mincls(index_BB), cltt-errdown)
                call coop_asy_error_bar(fpclBB, dble(l), cltt, dy_plus = errup, dy_minus = errdown)
             enddo
             call close_file(fp)
          endif
       endif

       call pp_init(As = 1.d0, ns = 0.d0, nrun = 0.d0, r = 1.d0, nt = 0.d0, lnkmin = dble(lnkmin),  lnkmax = dble(lnkmax), psparams = lnpsmean, ptparams = lnptmean, tmp_mode = "spline0")
       call pp_get_potential()
       ps = 1.e10 * exp(lnpsmean)
       pt = 1.e10 * exp(lnptmean)

       call coop_asy_curve(fp2, kmpc, ps, smooth = .false., color = "red", linetype = "solid", linewidth = 2., legend="mean scalar")
       call coop_asy_curve(fp2, kmpc, pt, smooth = .false., color = "violet", linetype = "solid", linewidth = 1.2, legend="mean tensor")
       call coop_asy_curve(fp2, (/ kpiv*exp(lnkmin), kpiv*exp(lnkmax) /), (/ exp(3.091+(0.967-1.)*lnkmin),  exp(3.091+(0.967-1.)*lnkmax) /), smooth = .false., color = "black", linewidth=2., legend="$m^2\phi^2$ scalar")
       call coop_asy_curve(fp2, (/ kpiv*exp(lnkmin), kpiv*exp(lnkmax) /), (/ exp(3.091 - 0.01625*lnkmin)*0.13,  exp(3.091-0.01625*lnkmax)*0.13 /), smooth = .false., color = "cyan", linewidth=1.2, legend="$m^2\phi^2$ tensor")
       if(index(mc%prefix, "nobicep").ne.0)then
          ps = exp(3.114+(0.96-1. + (-0.013)/2.*lnk)*lnk)
       else
          ps = exp(3.114+(0.9593-1. + (-0.0285)/2.*lnk)*lnk)
       endif
       call coop_asy_curve(fp2, kmpc, ps, smooth = .false., color = "green", linetype = "solid", linewidth = 1.2, legend="$n_{\rm run}$ best-fit")
       select case(trim(pp_mode))
       case("spline0", "linear0")
          numpp = num_initpower - 5
       case("spline1", "linear1")
          numpp = num_initpower - 4
       case("spline2", "linear2")
          numpp = num_initpower - 3
       case default
          numpp = 0
       end select
       if(numpp .gt. 1)then
          call coop_set_uniform(numpp, lnk(1:numpp), dble(lnkmin), dble(lnkmax))
          ps(1:numpp) = 1.3
          call coop_asy_dots(fp2, kpiv*exp(lnk(1:numpp)), ps(1:numpp), "black", "$\Delta$")
       endif

       call coop_asy_interpolate_curve(fpv, pp_phi, pp_lnV, interpolate="LinearLinear", color = "red", linetype = "solid", linewidth = 2., legend = "mean traj.")
       call coop_asy_interpolate_curve(fpeps, kpiv*exp(pp_lnk), pp_epsilon, interpolate="LogLinear", color = "red", linetype = "solid", linewidth=2., legend="mean traj.")

       call coop_asy_legend(fp2, kpiv*exp(lnkmin + 1.), 116., 2)
       call coop_asy_legend(fpv, -0.2, 0.35, 1)
       call coop_asy_legend(fpeps, kpiv*exp(lnkmin+4.), 0.12, 1)
       call close_file(fp2)
       call close_file(fpv)
       call close_file(fpeps)
       if(do_cl)then
          if(any(bestcls(:, 2:dcl_lmax) .ne. 0.d0))then
             call coop_asy_line(fpdclTT, 1.8, 0., dcl_lmax+1., 0., linewidth = 0.5)
             call coop_asy_label(fpdclTT, "{\small ref: $\Lambda$CDM best-fit}", dcl_legend_loc,  real((maxdcls(index_TT) - mindcls(index_TT))*0.12 + mindcls(index_TT)))
             call coop_asy_line(fpdclEE, 1.8, 0., dcl_lmax+1., 0., linewidth = 0.5)
             call coop_asy_label(fpdclEE, "{\small ref: $\Lambda$CDM best-fit}", dcl_legend_loc,  real((maxdcls(index_EE) - mindcls(index_EE))*0.12 + mindcls(index_EE)))
             call coop_asy_line(fpdclBB, 1.8, 0., dcl_lmax+1., 0., linewidth = 0.5)
             call coop_asy_label(fpdclBB, "{\small ref: $\Lambda$CDM best-fit}", dcl_legend_loc,  real((maxdcls(index_BB) - mindcls(index_BB))*0.15 + mindcls(index_BB)))

             call coop_asy_line(fpdclTE, 1.8, 0., dcl_lmax+1., 0., linewidth = 0.5)
             call coop_asy_label(fpdclTE, "{\small ref: $\Lambda$CDM best-fit}", dcl_legend_loc,  real((maxdcls(index_TE) - mindcls(index_TE))*0.12 + mindcls(index_TE)))
          endif
          call coop_asy_legend(fpdclTT, dcl_legend_loc, real((maxdcls(index_TT) - mindcls(index_TT))*0.92 + mindcls(index_TT)), 1)
          call coop_asy_legend(fpdclEE, dcl_legend_loc, real((maxdcls(index_EE) - mindcls(index_EE))*0.92 + mindcls(index_EE)), 1)
          call coop_asy_legend(fpdclBB, dcl_legend_loc, real((maxdcls(index_BB) - mindcls(index_BB))*0.92 + mindcls(index_BB)), 1)
          call coop_asy_legend(fpdclTE, dcl_legend_loc, real((maxdcls(index_TE) - mindcls(index_TE))*0.92 + mindcls(index_TE)), 1)

          call coop_asy_legend(fpclTT, min(3.5, lmax/2.), real((maxcls(index_TT) - mincls(index_TT))*0.92 + mincls(index_TT)), 1)
          call coop_asy_legend(fpclEE, min(3.5, lmax/2.), real((maxcls(index_EE) - mincls(index_EE))*0.92 + mincls(index_EE)), 1)
          call coop_asy_legend(fpclBB, min(3.5, lmax/2.), real((maxcls(index_BB) - mincls(index_BB))*0.92 + mincls(index_BB)), 1)
          call coop_asy_legend(fpclTE, min(3.5, lmax/2.), real((maxcls(index_TE) - mincls(index_TE))*0.92 + mincls(index_TE)), 1)
          call close_file(fpdclTT)
          call close_file(fpdclEE)
          call close_file(fpdclBB)
          call close_file(fpdclTE)
          call close_file(fpclTT)
          call close_file(fpclEE)
          call close_file(fpclBB)
          call close_file(fpclTE)
       endif
       if(do_traj_cov)then
          lnps_shift = sum(lnpscov**2)/nk/nk*1.d-18
          do ik=1, nk
             lnpscov(ik,ik) = lnpscov(ik,ik) + lnps_shift
          enddo
          call fp%open(trim(mc%output)//"_pwtraj_eig.txt","w")
          call fp%init(xlabel="$ k ({\rm Mpc}^{-1})$", ylabel = "$\delta \ln \Delta_S^2$", xlog = .true. , xmin = exp(lnkmin - 0.01)*kpiv, xmax = kpiv*exp(lnkmax + 0.01), width = 7.2, height = 6.)
          call matsym_diagonalize(lnpscov, lnps, sort = .true.)
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
          write(*,*) "found "//trim(coop_num2str(ndof))//" degrees of freedom in scalar power spectrum"
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


          call coop_asy_legend(fp, 0.0001, ytop*0.9)
          call close_file(fp)
       endif
    endif
    
    deallocate(initpower)

    call fp%open(trim(mc%output)//".margstats", "w")
    select case(mcmc_stat_num_cls)
    case(3:)
       write(fp%unit, "(2A32, 10A16)") "name", "label", "best", "mean", "std", "base", "-1sigma", "+1sigma", "-2sigma", "+2sigma", "-3sigma", "+3sigma"
       do i=1, mc%np
          write(fp%unit, "(2A32, 10G16.7)") trim(mc%name(i)), trim(mc%label(i)), mc%params(mc%ibest, i), mc%mean(i), mc%std(i), mc%base(i), mc%lowsig(1, i), mc%upsig(1, i), mc%lowsig(2, i), mc%upsig(2, i),mc%lowsig(3, i), mc%upsig(3, i)
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
    call close_file(fp)
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
    do ip = 1, mc%np
       if(mc%vary(ip))then
          write(fp%unit, "(A)") trim(mc%label(ip))//" & "//trim(latex_range(mc%base(ip), mc%lowsig(:, ip), mc%upsig(:, ip)))//" \\ "
          Write(fp%unit,"(a)")"\hline"
       endif
    end do
    Write(fp%unit,"(a)")"\end{tabular}"
    Write(fp%unit,"(a)")"\end{center}"
    Write(fp%unit,"(a)")"\end{table}"
    Write(fp%unit,"(a)")"\end{document}"
    call close_file(fp)

    call fp%open(trim(mc%output)//".covmat", "w")
    call write_matrix(fp%unit, mc%covmat, mc%np, mc%np)
    call close_file(fp)
    call fp%open(trim(mc%output)//".corr", "w")
    call write_matrix(fp%unit, mc%corrmat, mc%np_used, mc%np_used)
    call close_file(fp)
    call fp%open(trim(mc%output)//".covused", "w")
    call write_matrix(fp%unit, mc%cov_used, mc%np_used, mc%np_used)
    call close_file(fp)
    if(mc%np_pca .gt. 0)then
       allocate(pcamat(mc%np_pca, mc%np_pca),eig(mc%np_pca), ipca(mc%np_pca))
       pcamat = mc%covmat(mc%pca, mc%pca)
       call coop_set_uniform(mc%np_pca, ipca, 1.d0, 1.d0*mc%np_pca)
       call matsym_diagonalize(pcamat, eig, sort=.true.) !!sort eigen values
       eig = sqrt(eig)
       call fp%open(trim(mc%output)//".pcamat", "w")
       write(fp%unit, "(A)") "# format is  i, sigma_i (newline) eigen vector (i = 2, 3, ...)"
       do i = 1, mc%np_pca
          write(fp%unit, "(I8, G14.5)") i, eig(i)
          write(fp%unit, "("//trim(coop_num2str(mc%np_pca))//"G14.5)") pcamat(:, i)
       enddo
       call close_file(fp)
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
       call close_file(fp)
       deallocate(pcamat , eig , ipca)
    endif
    do ip = 1, mc%np_used
       call fp%open(trim(mc%output)//"_"//trim(mc%simplename(mc%used(ip)))//"_1D.txt", "w")
       do i = 1, mc%nb
          x(i) = mc%plotlower(mc%used(ip)) + mc%dx(mc%used(ip))*(i-1)
       enddo
       call fp%init( xlabel = trim(mc%label(mc%used(ip))), ylabel = "P")
       call coop_asy_curve(fp, x, mc%c1d(:, ip)/maxval(mc%c1d(:, ip)))
       call close_file(fp)
    enddo

    do j = 1, mc%np_used
       do j2 = 1, j
          k = j*(j-1)/2 + j2
          if(mc%want_2d_output(j, j2))then
             call fp%open(trim(mc%output)//"_"//trim(mc%simplename(mc%used(j)))//"_"//trim(mc%simplename(mc%used(j2)))//"_2D.txt", "w")
             call fp%init( xlabel = trim(mc%label(mc%used(j))), ylabel = trim(mc%label(mc%used(j2))), xmin=mc%plotlower(mc%used(j)), xmax = mc%plotupper(mc%used(j)), ymin=mc%plotlower(mc%used(j2)), ymax = mc%plotupper(mc%used(j2)) )
             call coop_asy_path_from_array(path, mc%c2d(:, :, k), mc%plotlower(mc%used(j)), mc%plotupper(mc%used(j)), mc%plotlower(mc%used(j2)), mc%plotupper(mc%used(j2)), mc%cut2d(2, k))
             call coop_asy_contour_path(fp, path, colorfill = trim(mc%color2d_light), smooth = .false., linecolor = "black", linetype = "solid")
             call coop_asy_path_from_array(path, mc%c2d(:, :, k), mc%plotlower(mc%used(j)), mc%plotupper(mc%used(j)), mc%plotlower(mc%used(j2)), mc%plotupper(mc%used(j2)), mc%cut2d(1, k))
             call coop_asy_contour_path(fp, path, colorfill = trim(mc%color2d_dark), smooth = .false., linecolor = "black", linetype = "solid")
             call close_file(fp)
          endif
          if(j2 .ne. j .and. mc%want_2d_output(j2, j))then
             call fp%open(trim(mc%output)//"_"//trim(mc%simplename(mc%used(j2)))//"_"//trim(mc%simplename(mc%used(j)))//"_2D.txt", "w")
             call fp%init( xlabel = trim(mc%label(mc%used(j2))), ylabel = trim(mc%label(mc%used(j))), xmin=mc%plotlower(mc%used(j2)), xmax = mc%plotupper(mc%used(j2)), ymin=mc%plotlower(mc%used(j)), ymax = mc%plotupper(mc%used(j)) )
             call coop_asy_path_from_array(path, transpose(mc%c2d(:, :, k)),  mc%plotlower(mc%used(j2)), mc%plotupper(mc%used(j2)), mc%plotlower(mc%used(j)), mc%plotupper(mc%used(j)), mc%cut2d(2, k))
             call coop_asy_contour_path(fp, path, colorfill = trim(mc%color2d_light), smooth = .false., linecolor = "black", linetype = "solid")
             call coop_asy_path_from_array(path, transpose(mc%c2d(:, :, k)), mc%plotlower(mc%used(j2)), mc%plotupper(mc%used(j2)),  mc%plotlower(mc%used(j)), mc%plotupper(mc%used(j)), mc%cut2d(1, k))
             call coop_asy_contour_path(fp, path, colorfill = trim(mc%color2d_dark), smooth = .false., linecolor = "black", linetype = "solid")
             call close_file(fp)

          endif
       enddo
    enddo    
  end subroutine export_stats


  subroutine get_initpower(mc, ind, initp, np)
    integer np, ind, i
    type(mcmc_chain) mc
    COOP_REAL  initp(np)
    initp(1) = name2value(mc, ind, "ns")
    initp(2) = name2value(mc, ind, "nt")
    initp(3) = name2value(mc, ind, "nrun")
    initp(4) = name2value(mc, ind, "logA")
    initp(5) = name2value(mc, ind, "r")
    initp(6) = name2value(mc, ind, "Aphiphi")
    select case(trim(pp_mode))
    case("standard")
       return
    case("elldip")
       pp_dip_Amp = name2value(mc, ind, "dipamp")
       return
    case("bump")
       initp(7) = name2value(mc, ind, "mphi")
       initp(8) = name2value(mc, ind, "nefolds")
       initp(9) = name2value(mc, ind, "bumpdn")
       initp(10) = name2value(mc, ind, "bumpamp")
       initp(11) = name2value(mc, ind, "bumpwidth")
    case("spline0", "spline1", "spline2", "linear0", "linear1", "linear2")
       do i = 1, np - 6
          initp(i+6) = name2value(mc, ind, "pp"//trim(coop_num2str(i)))
       enddo
    case("m2phi2")
       initp(7) = name2value(mc, ind, "invm2")
       initp(8) = name2value(mc, ind, "nefolds")
    case("lambdaphi4")
       initp(7) = name2value(mc, ind, "lambda")
       initp(8) = name2value(mc, ind, "nefolds")
    case("higgs")
       initp(7) = name2value(mc, ind, "lambda")
       initp(8) = name2value(mc, ind, "nefolds")
       initp(9) = name2value(mc, ind, "xih")
       initp(10) = name2value(mc, ind, "alphah")
    case default
       stop "Unknown pp_mode"
    end select
  end subroutine get_initpower


  function chain_index_of_name(mc, name) result(ind)
    type(mcmc_chain) mc
    COOP_UNKNOWN_STRING name
    integer ind, ip
    COOP_STRING sname
    sname = trim(coop_str_numalpha(name))
    if(trim(sname).eq."")then
       ind = 0
       return
    endif
    do ip = 1, mc%np
       if(trim(sname).eq. trim(mc%simplename(ip)))then
          ind = ip
          return
       endif
    enddo
    ind = 0
    return
  end function chain_index_of_name

  subroutine chain_preprocess(mc)
    type(mcmc_chain) mc
  end subroutine chain_preprocess


!!Because cosmomc changes every two months, it is not a good idea to link to it.
!!The idea is to 1. Produce an ini file for camb; 2. call camb to produce Cls. 3. Read camb output.
  subroutine get_camb_cls(mc, isample, lmax, cls)
    COOP_UNKNOWN_STRING,parameter::camb_path = "../../cosmomc/camb/"
    type(mcmc_chain) mc
    integer isample, lmax
    integer il, l
    COOP_REAL  cls(4, 2:lmax)  !!TT, EE, BB, TE
    type(coop_file) fp
    COOP_REAL ,dimension(:),allocatable:: initpower
    integer num_initpower
    character ans
    call fp%open("tmp.ini", "w")
    write(fp%unit, "(A)") "output_root = tmp"
    write(fp%unit, "(A)") "get_scalar_cls = T"
    write(fp%unit, "(A)") "get_vector_cls = F"
    write(fp%unit, "(A)") "get_tensor_cls = T"
    write(fp%unit, "(A)") "get_transfer = F"
    write(fp%unit, "(A)") "do_lensing = T"
    write(fp%unit, "(A)") "do_nonlinear = 0"
    write(fp%unit, "(A)") "l_max_scalar = "//trim(coop_num2str(max(lmax+200, 2000)))
    write(fp%unit, "(A)") "l_max_tensor = 1500"
    write(fp%unit, "(A)") "k_eta_max_tensor = 3000"
    write(fp%unit, "(A)") "use_physical = T"
    write(fp%unit, "(A)") "ombh2 = "//trim(coop_num2str(name2value(mc, isample, "omegabh2")))
    write(fp%unit, "(A)") "omch2 = "//trim(coop_num2str(name2value(mc, isample, "omegach2")))
    write(fp%unit, "(A)") "omnuh2 = "//trim(coop_num2str(name2value(mc, isample, "omeganuh2")))
    write(fp%unit, "(A)") "omk = "//trim(coop_num2str(name2value(mc, isample, "omegak")))
    write(fp%unit, "(A)") "hubble = "//trim(coop_num2str(name2value(mc, isample, "H0")))
    write(fp%unit, "(A)") "w = "//trim(coop_num2str(name2value(mc, isample, "w")))
    write(fp%unit, "(A)") "cs2_lam = 1"
    write(fp%unit, "(A)") "temp_cmb = 2.7255"
    write(fp%unit, "(A)") "helium_fraction = "//trim(coop_num2str(name2value(mc, isample, "yheused")))

    write(fp%unit, "(A)") "primordial_power_mode = "//trim(pp_mode)
    call coop_dictionary_lookup(mc%inputparams, "num_init_power", num_initpower , 6)
    write(fp%unit, "(A)")  "num_initpower = "//trim(coop_num2str(num_initpower))
    allocate(initpower(num_initpower))
    call get_initpower(mc, isample, initpower, num_initpower)
    write(fp%unit, "(A, "//trim(coop_num2str(num_initpower))//"G15.6)") "initpower =", initpower
    deallocate(initpower)
    write(fp%unit, "(A)") "inflation_consistency = "//trim(mc%inputparams%value("inflation_consistency"))
    write(fp%unit, "(A)") "k_pivot = "//trim(mc%inputparams%value( "pivot_k"))
    write(fp%unit, "(A)") "massless_neutrinos = 2.046"
    write(fp%unit, "(A)") "nu_mass_eigenstates = 1"
    write(fp%unit, "(A)") "massive_neutrinos  = 1"
    write(fp%unit, "(A)") "share_delta_neff = T"
    write(fp%unit, "(A)") "nu_mass_fractions = 1"
    write(fp%unit, "(A)") "nu_mass_degeneracies = "
    write(fp%unit, "(A)") "reionization         = T"
    write(fp%unit, "(A)") "re_use_optical_depth = T"
    write(fp%unit, "(A)") "re_optical_depth     = "//trim(coop_num2str(name2value(mc, isample, "tau")))
    write(fp%unit, "(A)") "re_delta_redshift    = 1.5"
    write(fp%unit, "(A)") "re_ionization_frac   = -1"
    write(fp%unit, "(A)") "RECFAST_fudge = 1.14"
    write(fp%unit, "(A)") "RECFAST_fudge_He = 0.86"
    write(fp%unit, "(A)") "RECFAST_Heswitch = 6"
    write(fp%unit, "(A)") "RECFAST_Hswitch  = T"
    write(fp%unit, "(A)") "initial_condition   = 1"
    write(fp%unit, "(A)") "initial_vector = -1 0 0 0 0"
    write(fp%unit, "(A)") "vector_mode = 0"
    write(fp%unit, "(A)") "COBE_normalize = F"
    write(fp%unit, "(A)") "CMB_outputscale = 7.42835025e12"
    write(fp%unit, "(A)") "scalar_output_file = scalCls.dat"
    write(fp%unit, "(A)") "vector_output_file = vecCls.dat"
    write(fp%unit, "(A)") "tensor_output_file = tensCls.dat"
    write(fp%unit, "(A)") "total_output_file  = totCls.dat"
    write(fp%unit, "(A)") "lensed_output_file = lensedCls.dat"
    write(fp%unit, "(A)") "lensed_total_output_file  =lensedtotCls.dat"
    write(fp%unit, "(A)") "lens_potential_output_file = lenspotentialCls.dat"
    write(fp%unit, "(A)") "FITS_filename      = scalCls.fits"
    write(fp%unit, "(A)") "do_lensing_bispectrum = F"
    write(fp%unit, "(A)") "do_primordial_bispectrum = F"
    write(fp%unit, "(A)") "feedback_level = 0"
    write(fp%unit, "(A)") "derived_parameters = T"
    write(fp%unit, "(A)") "lensing_method = 1"
    write(fp%unit, "(A)") "accurate_BB = F"
    write(fp%unit, "(A)") "massive_nu_approx = 1"
    write(fp%unit, "(A)") "accurate_polarization   = T"
    write(fp%unit, "(A)") "accurate_reionization   = T"
    write(fp%unit, "(A)") "do_tensor_neutrinos     = T"
    write(fp%unit, "(A)") "do_late_rad_truncation   = T"
    write(fp%unit, "(A)") "number_of_threads       = 0"
    write(fp%unit, "(A)") "high_accuracy_default=T"
    write(fp%unit, "(A)") "accuracy_boost          = 1"
    write(fp%unit, "(A)") "l_accuracy_boost        = 1"
    write(fp%unit, "(A)") "l_sample_boost          = 1"
    call close_file(fp)
    call system(camb_path//"camb tmp.ini")
    call fp%open("tmp_lensedtotCls.dat", "r")
    do l=2, lmax
       read(fp%unit, *) il, Cls(:, l)
       if(il.ne.l) stop "camb output broken"
    enddo
    call close_file(fp)
    call clean_camb_output()
  end subroutine get_camb_cls

  subroutine clean_camb_output()
    call delete_file("tmp.ini")
    call delete_file("tmp_*Cls.dat")
  end subroutine clean_camb_output

  function name2value(mc, isample, name) result(val)
    type(mcmc_chain) mc
    integer isample
    COOP_UNKNOWN_STRING name
    real val
    integer ind
    COOP_SHORT_STRING str
    if(isample .le. 0 .or. isample .gt. mc%n)  stop "name2value: index overflow"
    ind = chain_index_of_name(mc, name)
    if(ind.eq.0)then
       str = mc%inputparams%value("param["//trim(name)//"]")
       if(trim(str).eq."") stop "name2value: cannot find the name"
       read(str, *) val
    else
       val = mc%params( isample, ind )
    end if
  end function name2value



End Module Statchains
