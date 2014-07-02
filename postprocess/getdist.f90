program getdist
  use inifile_utils
  use wrap_utils
  use asy_utils
  use statchains
  use cosmolibwrap
  implicit none
#include "utils.h"
  type(mcmc_chain) mc
  logical error
  STRING fini, prefix, outdir, p, name1, name2
  LONG_STRING inline
  integer discard_percent
  type(list_string) sl
  type(list_integer) pcal
  integer i, j, ip1,ip2, if1, if2
  call global_initialization(nobessel = .true.)
  fini = GetIniParams(1)
  if(fini(len_trim(fini)-3:len_trim(fini)).ne.".ini")then
     fini=trim(fini)//".ini"
  endif  
  call Ini_open(trim(fini), tmp_file_unit, error)
  if(error)then
     print*,"can not find ini file:"//trim(fini)
     stop
  endif
  prefix = ini_read_string("root", .false.)
  if(trim(prefix).eq."") stop "You need to specify the key 'root' in ini file"
  if(trim(getparam(2)).ne."")then
     prefix = trim(adjustl(file_path_of(prefix)))//trim(adjustl(getparam(2)))
  endif
  outdir = ini_read_string("output", .false.)
  if(trim(outdir) .eq. "") stop "You need to specify the key 'output' in ini file"
  discard_percent = ini_read_int("discard_percent", 30)
  write(*,*) "discarding "//trim(num2str(discard_percent))//"% samples at the beginning of chains"
  mcmc_stat_cls(1) = ini_read_real("1sigma_threshold", mcmc_stat_cls(1))
  mcmc_stat_cls(2) = ini_read_real("2sigma_threshold", max(mcmc_stat_cls(2), mcmc_stat_cls(1)+0.001))
  mcmc_stat_cls(3) = ini_read_real("3sigma_threshold", max(mcmc_stat_cls(3), mcmc_stat_cls(2)+0.001))
  if(mcmc_stat_cls(3).gt.1.) stop "your sigma_threshold is too high"
  
  mc%color2d_light = ini_read_string("color2d_light", .false.)
  mc%color2d_dark = ini_read_string("color2d_dark", .false.)
  if(trim(mc%color2d_light).eq."") mc%color2d_light = "skyblue"
  if(trim(mc%color2d_dark).eq."") mc%color2d_dark = "blue"

  measured_cltt_file = ini_read_string("measured_cltt_file", .false.)
  measured_clee_file = ini_read_string("measured_clee_file", .false.)
  measured_clbb_file = ini_read_string("measured_clbb_file", .false.)
  measured_clte_file = ini_read_string("measured_clte_file", .false.)
  bestfit_cl_file = ini_read_string("bestfit_cl_file", .false.)
  call load_chain(mc, prefix, discard_percent)

  inline = trim(adjustl(ini_read_string("want2d", .false.)))
  if(trim(inline).ne. "")then
     call string_to_list_string(inline, sl, ";")
     do i=1, sl%n
        call list_get_element(sl, i, p)
        j = scan(p, ",")
        if(j.eq.0)cycle
        name1 = trim(adjustl(p(1:j-1)))
        name2 = trim(adjustl(p(j+1:)))
        ip1 = chain_index_of_name(mc, name1)
        ip2 = chain_index_of_name(mc, name2)
        if(ip1 .gt. 0.and. ip2.gt.0)then
           if(mc%vary(ip1) .and. mc%vary(ip2))then
              write(*,*) " 2d contours for "//trim(name1)//" "//trim(name2)
              mc%want_2d_output(ip1, ip2) = .true.
           endif
        endif
     enddo
  endif
  inline = trim(adjustl(ini_read_string("wantpca", .false.)))
  if(trim(inline).ne."")then
     call string_to_list_string(inline, sl, ",")
     mc%np_pca = sl%n
     if(allocated(mc%pca))deallocate(mc%pca)
     allocate(mc%pca(mc%np_pca))
     do i = 1, sl%n
        mc%pca(i) = chain_index_of_name(mc, list_element(sl, i))
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
  call export_stats(mc, outdir)

end program getdist


