program getdist
  use coop_statchains_mod
  use coop_wrapper_utils
  use coop_inifile_mod
  implicit none
#include "constants.h"
  type(coop_mcmc_chain) mc
  logical error
  COOP_STRING fini, prefix, outdir, p, name1, name2
  COOP_LONG_STRING inline
  integer discard_percent
  type(coop_list_string) sl
  type(coop_list_integer) pcal
  integer i, j, ip1,ip2, if1, if2
  call coop_random_init()
  fini = coop_InputArgs(1)
  if(fini(len_trim(fini)-3:len_trim(fini)).ne.".ini")then
     fini=trim(fini)//".ini"
  endif
  
  call coop_load_dictionary(fini, mc%settings)  
  call Ini_open(trim(fini), 7, error)
  if(error)then
     print*,"can not find ini file:"//trim(fini)
     stop
  endif
  prefix = ini_read_string("root", .false.)
  mc%do_preprocess = ini_read_logical("preprocess", .false.)
  mc%extmode = ini_read_int("extension_mode", 0)
  if(mc%extmode > 0)then
     write(*,*) "doing extensions with mode ", mc%extmode
     mypp_model = mc%extmode
  endif
  if(trim(prefix).eq."") stop "You need to specify the key 'root' in ini file"
  if(trim(Coop_InputArgs(2)).ne."")then
     prefix = trim(adjustl(coop_file_path_of(prefix)))//trim(adjustl(coop_inputArgs(2)))
  endif
  if(trim(coop_InputArgs(3)).ne."")then
     mc%datasets = trim(coop_InputArgs(3))
  else
     mc%datasets = ini_read_string("datasets")
  endif
  outdir = ini_read_string("output", .false.)
  coop_postprocess_nbins = ini_read_int("num_bins", 0)
  coop_postprocess_num_contours = ini_read_int("num_contours", 2)
  coop_postprocess_num_contours  = min(coop_postprocess_num_contours, mcmc_stat_num_cls)
  coop_postprocess_do_cls = ini_read_logical("do_cls", .false.)  
  if(trim(outdir) .eq. "") stop "You need to specify the key 'output' in ini file"
  discard_percent = ini_read_int("discard_percent", 30)
  mc%fit_skewness_threshold = ini_read_real("fit_skewness_threshold", 0.1)
  mc%fit_kurtosis_threshold = ini_read_real("fit_kurtosis_threshold", 0.02)
  if(discard_percent .lt. 100)then
     write(*,*) "discarding "//trim(coop_num2str(discard_percent))//"% samples at the beginning of chains"
  else
     write(*,*) "discarding "//trim(coop_num2str(discard_percent))//"lines at the beginning of chains"     
  endif

  if(mcmc_stat_num_cls .ge.6)then
     stop "too many contour levels, please adjust mcmc_stat_num_cls in statchains.f90"
  endif

  do i=1, mcmc_stat_num_cls
     mcmc_stat_cls(i) = ini_read_real("confidence_level_"//trim(coop_num2str(i))//"sigma", mcmc_stat_cls(i))
     mc%color2d(i) = ini_read_string("color2d_"//trim(coop_num2str(i))//"sigma", .false.)
     if(trim(mc%color2d(i)).eq."") mc%color2d(i) = trim(adjustl(coop_asy_rgb_color(0.8*(i-1)/mcmc_stat_num_cls, min(1.2*(i-1)/mcmc_stat_num_cls, 1.), 1.)))
     print*, "color2d_"//trim(coop_num2str(i))//"sigma: "//trim(mc%color2d(i))
  enddo
  measured_cltt_file = ini_read_string("measured_cltt_file", .false.)
  bestfit_cl_file = ini_read_string("bestfit_cl_file", .false.)
  bestfit_run_file = ini_read_string("bestfit_run_file", .false.)
  bestfit_varytau_file = ini_read_string("bestfit_varytau_file", .false.)
  call mc%load(prefix, discard_percent)

  inline = trim(adjustl(ini_read_string("want2d", .false.)))
  if(trim(inline).ne. "")then
     call coop_string_to_list(inline, sl, ";")
     do i=1, sl%n
        call coop_list_get_element(sl, i, p)
        j = scan(p, ",")
        if(j.eq.0)cycle
        name1 = trim(adjustl(p(1:j-1)))
        name2 = trim(adjustl(p(j+1:)))
        ip1 = mc%index_of(name1)
        ip2 = mc%index_of(name2)
        if(ip1 .gt. 0.and. ip2.gt.0)then
           ip1 = mc%map2used(ip1)
           ip2 = mc%map2used(ip2)
           if(ip1.ne.0 .and. ip2 .ne.0)then
              write(*,*) " 2d contours for "//trim(name1)//" "//trim(name2)
              mc%want_2d_output(ip1, ip2) = .true.
           endif
        endif
     enddo
  endif

  inline = trim(adjustl(ini_read_string("want1d", .false.)))
  if(trim(inline).ne. "")then
     mc%want_1d_output = .false.
     call coop_string_to_list(inline, sl, ",")
     do i=1, sl%n
        call coop_list_get_element(sl, i, p)
        name1 = trim(adjustl(p))
        ip1 = mc%index_of(name1)
        if(ip1 .gt. 0)then
           ip1 = mc%map2used(ip1)
           if(ip1 .ne. 0) mc%want_1d_output(ip1) = .true.
        endif
     enddo
  else
     mc%want_1d_output = .true.  !!simply plot all
  endif
  
  inline = trim(adjustl(ini_read_string("wantpca", .false.)))
  if(trim(inline).ne."")then
     call coop_string_to_list(inline, sl, ",")
     mc%np_pca = sl%n
     if(allocated(mc%pca))deallocate(mc%pca)
     allocate(mc%pca(mc%np_pca))
     do i = 1, sl%n
        mc%pca(i) = mc%index_of(sl%element(i))
        if(mc%pca(i).eq.0)then
           write(*,*) "PCA name not matched"
           mc%np_pca = 0
           deallocate(mc%pca)
           exit
        endif
        if(.not. mc%vary(mc%pca(i)))then
           write(*,*) "PCA name not matched"
           mc%np_pca = 0
           deallocate(mc%pca)
           exit
        endif
     enddo
     write(*,*) "Doing PCA for "//trim(inline)
  endif
  call mc%export_stats(outdir)

end program getdist


