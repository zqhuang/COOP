program stackth
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod  
  use coop_fitsio_mod
  use coop_healpix_mod
#ifdef HAS_HEALPIX  
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
#endif  
  implicit none
#include "constants.h"

#ifdef HAS_HEALPIX
  logical,parameter::flat = .false. !!use nonflat is actually faster  
  logical::want_noise = .false.
  logical::normalize_sigma0, hot
  COOP_REAL::sigma0_norm = 87.44  !#59.141857614604476d0
  
  COOP_REAL:: r_degree, width, dr, nu, fwhm_arcmin
  COOP_INT::lmax, hpauto_lowl, hpauto_highl, n, ind_auto, ind_cross, hpcross_lowl, hpcross_highl
  COOP_STRING::clfile, field, peak, orient, output_prefix, colortable
  type(coop_file)::fp, fpfr
  type(coop_asy)::figCr, fpI(0:2), fpQU(0:2), fpIdat(0:2), fpQUdat(0:2)
  integer i, j, il, l, iomega, m
  type(coop_arguments)::args
  type(coop_healpix_patch)::patchI, patchQU, patchQrUr
  COOP_REAL, dimension(:), allocatable::Pl0, Pl2, Pl4, r, ell
  COOP_REAL:: sigma, sigma2, sigma0, sigma1, cosbeta, j2, j4, j0, omega, weights(4), pomega, phi, romega, kr
  COOP_REAL,dimension(:,:),allocatable:: cls, l2cls,  frI, frQU
  COOP_REAL,dimension(:,:,:),allocatable:: cr
  call coop_MPI_Init()
  if(iargc() .lt. 2)then
     write(*,"(A)") "Syntax:"
     write(*,"(A)") "./GetTheo -cl CL_FILE -auto CL_COLUMN_AUTO_CORRELATION -cross CL_COLUMN_CROSS_CORRELATION -field [T|E|B|zeta] -peak [RANDOM|T|E|B|P_T|P|\zeta|P_\zeta] -out OUTPUT_PREFIX"
     write(*, "(A)") "Optional:"
     write(*, "(A)") "-orient [RANDOM|(Q_T,U_T)|(Q,U)]"
     write(*, "(A)") "-nu THRESHOLD"
     write(*, "(A)") "-fwhm FWHM_IN_ARCMIN"
     write(*, "(A)") "-hot [T/F]"
     write(*, "(A)") "-want_caption [T/F]"
     write(*, "(A)") "-want_label [T/F]"
     write(*, "(A)") "-want_arrow [T/F]"
     write(*, "(A)") "-lmax LMAX"     
     write(*, "(A)") "-radius RAIDUS_IN_DEGREE"
     write(*, "(A)") "-res RESOLUTION[10-200]"
     write(*, "(A)") "-normrms NORMALIZE_SIGMA_0[T/F]"
     write(*, "(A)") "-hpauto_lowl HIGHPASS_AUTO_LOWER"
     write(*, "(A)") "-hpauto_highl HIGHPASS_AUTO_UPPER"
     write(*, "(A)") "-hpcross_lowl HIGHPASS_CROSS_LOWER"
     write(*, "(A)") "-hpcross_highl HIGHPASS_CROSS_UPPER"
     write(*, "(A)") "-color_table [Rainbow/Planck]"     
     
     stop
  endif
  call coop_get_command_line_argument(key = 'cl', arg = clfile)
  call coop_get_command_line_argument(key = 'auto', arg = ind_auto)
  call coop_get_command_line_argument(key = 'cross', arg = ind_cross)  
  call coop_get_command_line_argument(key = 'field', arg = field)
  call coop_get_command_line_argument(key = 'peak', arg = peak)
  call coop_get_command_line_argument(key = 'out', arg = output_prefix)  
  call coop_get_command_line_argument(key = 'orient', arg = orient, default = 'RANDOM')
  call coop_get_command_line_argument(key = 'nu', arg = nu, default = 0.d0)
  call coop_get_command_line_argument(key = 'fwhm', arg = fwhm_arcmin, default = 5.d0)
    call coop_get_command_line_argument(key = 'want_noise', arg = want_noise, default = .false.)
  call coop_get_command_line_argument(key = 'want_caption', arg = coop_healpix_patch_default_want_caption, default = .true.)
  call coop_get_command_line_argument(key = 'want_label', arg = coop_healpix_patch_default_want_label, default  = .true.)
  call coop_get_command_line_argument(key = 'want_arrow', arg = coop_healpix_patch_default_want_arrow, default  = .true.)
  call coop_get_command_line_argument(key = 'lmax', arg = lmax, default = 2500)
  call coop_get_command_line_argument(key = 'hpauto_lowl', arg = hpauto_lowl, default =-1)
  call coop_get_command_line_argument(key = 'hpauto_highl', arg = hpauto_highl, default = 0)
  call coop_get_command_line_argument(key = 'hpcross_lowl', arg = hpcross_lowl, default = -1)
  call coop_get_command_line_argument(key = 'hpcross_highl', arg = hpcross_highl, default = 0)
  call coop_get_command_line_argument(key = 'radius', arg = r_degree, default = 2.d0)
  call coop_get_command_line_argument(key = "res", arg = n, default = 50)
  call coop_get_command_line_argument(key = 'normrms', arg = normalize_sigma0, default = .false.)
  call coop_get_command_line_argument(key = 'hot', arg = hot, default = .true.)
  call coop_get_command_line_argument(key = 'color_table', arg = colortable, default = "Rainbow")  
  call coop_random_init()
  
  allocate(Pl0(0:lmax), Pl2(0:lmax), Pl4(0:lmax), cls(4, 2:lmax), l2cls(4, 2:lmax), cr(0:1, 0:2, 0:n*3/2), frI(0:2, 0:n*3/2), frQU(0:2, 0:n*3/2) , ell(2:lmax), r(0:n*3/2))
  width = 2.d0*sin(r_degree*coop_SI_degree/2.d0)
  dr = width/n  
  sigma = fwhm_arcmin*coop_sigma_by_fwhm*coop_SI_arcmin
  if(hpauto_lowl .gt. 0 .or. hpcross_lowl.gt.0)then
     write(*,*) "doing highpass filtering", hpauto_lowl, hpauto_highl, hpcross_lowl, hpcross_highl
  endif
  call fp%open_skip_comments(clfile)
  do l=2, lmax
     read(fp%unit, *, ERR=100, END=100) il, l2cls(:, l)
     if(want_noise)then
        l2Cls(coop_index_ClTT, l)   =     l2Cls(coop_index_ClTT, l)  +  coop_Planck_TNoise(l)*(l*(l+1.d0))/coop_2pi*(COOP_DEFAULT_TCMB *1.d6)**2
        l2Cls(coop_index_ClEE, l)   =     l2Cls(coop_index_ClEE, l)  +  coop_Planck_ENoise(l)*(l*(l+1.d0))/coop_2pi*(COOP_DEFAULT_TCMB *1.d6)**2
     endif
     ell(l)  = l
     l2cls(:,l) = l2cls(:, l)*(coop_2pi*exp(-l*(l+1.d0)*sigma**2))
     l2cls(ind_cross, l) = l2cls(ind_cross, l)*coop_highpass_filter(hpcross_lowl, hpcross_highl, l)*coop_highpass_filter(hpauto_lowl, hpauto_highl, l)
     if(ind_auto .ne. ind_cross)then
        l2cls(ind_auto, l) = l2cls(ind_auto, l)*coop_highpass_filter(hpauto_lowl, hpauto_highl, l)**2
     endif
     cls(:,l) = l2cls(:,l)/(l*(l+1.d0))
     if(il.ne.l) stop "Error in Cl file"
  enddo
  goto 200
100 stop "Error in Cl file"
200 call fp%close()

  sigma0 = sqrt(sum(Cls(ind_auto,:)*(ell+0.5d0))/coop_2pi)
  write(*,*) "Expected rms = ", sigma0
  if(normalize_sigma0)then
     Cls = Cls*(sigma0_norm/sigma0)**2
     l2Cls = l2Cls*(sigma0_norm/sigma0)**2
     sigma0 = sigma0_norm
  endif
  sigma1 = sqrt(sum(l2Cls(ind_auto,:)*(ell+0.5d0))/coop_2pi)
  sigma2 = sqrt(sum(l2Cls(ind_auto,:)*(ell+0.5d0)*(ell*(ell+1.d0)))/coop_2pi)
  cosbeta = sigma1**2/sigma0/sigma2
  call coop_gaussian_npeak_set_args(args, 2, sigma0, sigma1, sigma2)

  select case(trim(peak))
  case("T", "E", "B", "zeta")
     select case(trim(orient))
     case("NULL", "RANDOM")
        call coop_gaussian_get_nonoriented_stacking_weights(nu, args, weights)
     case default
        call coop_gaussian_get_oriented_stacking_weights(nu, args, weights)
     end select
  case("P", "P_T", "P_\zeta")
     call coop_gaussian_get_pmax_stacking_weights(nu, args, weights)
  case("RANDOM")
     select case(trim(orient))
     case("NULL", "RANDOM")     
        call coop_gaussian_get_fieldpoints_nonoriented_stacking_weights(nu, args, weights)
     case default
        call coop_gaussian_get_fieldpoints_oriented_stacking_weights(nu, args, weights)
     end select
  case default
     write(*,*) trim(peak)
     stop "Unknown peak type"
  end select

  do m=0,2
     call fpI(m)%open(trim(output_prefix)//"_I_fwhm"//COOP_STR_OF(nint(fwhm_arcmin))//"_m"//COOP_STR_OF(m*2)//".txt")
     call fpIdat(m)%open(trim(output_prefix)//"_I_fwhm"//COOP_STR_OF(nint(fwhm_arcmin))//"_m"//COOP_STR_OF(m*2)//".dat")
     call fpQU(m)%open(trim(output_prefix)//"_QU_fwhm"//COOP_STR_OF(nint(fwhm_arcmin))//"_m"//COOP_STR_OF(m*2)//".txt")
     call fpQUdat(m)%open(trim(output_prefix)//"_QU_fwhm"//COOP_STR_OF(nint(fwhm_arcmin))//"_m"//COOP_STR_OF(m*2)//".dat")     
     call fpI(m)%init(xlabel="$\varpi$", ylabel="$T_"//COOP_STR_OF(m*2)//"$")
     call fpQU(m)%init(xlabel="$\varpi$", ylabel="$P_"//COOP_STR_OF(m*2)//"$")
     write(fpI(m)%unit, "(A)") "CURVE"
     write(fpQU(m)%unit, "(A)") "CURVE"
     write(fpI(m)%unit, "(I6)") n+1
     write(fpQU(m)%unit, "(I6)") n+1
     write(fpI(m)%unit, "(A)") "Theory"
     write(fpQU(m)%unit, "(A)") "Theory"
     write(fpI(m)%unit, "(A)") "blue_dashed"
     write(fpQU(m)%unit, "(A)") "blue_dashed"
     write(fpI(m)%unit, "(A)") "0"
     write(fpQU(m)%unit, "(A)") "0"
  enddo

  cr = 0.d0       
  do i = 0, n*3/2
     omega = dr*i
     r(i) = omega
     if(.not.flat)then  !!get Plms
        call coop_get_normalized_Plm_array(m=0, lmax=lmax, x = 1.d0-omega**2/2.d0, Plms = Pl0)
        call coop_get_normalized_Plm_array(m=2, lmax=lmax, x = 1.d0-omega**2/2.d0, Plms = Pl2)
        call coop_get_normalized_Plm_array(m=4, lmax=lmax, x = 1.d0-omega**2/2.d0, Plms = Pl4)
     endif
     do l = 0, lmax
        if(flat)then
           kr = sqrt(l*(l+1.d0))*omega
           j0 = coop_bessj(0, kr)
           j2 = coop_bessj(2, kr)
           j4 = coop_bessj(4, kr)
        else
           j0 = Pl0(l)
           j2 = Pl2(l)
           j4 = Pl4(l)
        endif
        cr(0, 0, i) = cr(0, 0, i) + (l+0.5d0)*j0*cls(ind_cross, l)
        cr(1, 0, i) = cr(1, 0, i) - (l+0.5d0)*j0*l2cls(ind_cross, l)
        cr(0, 1, i) = cr(0, 1, i) + (l+0.5d0)*j2*cls(ind_cross, l)
        cr(1, 1, i) = cr(1, 1, i) - (l+0.5d0)*j2*l2cls(ind_cross, l)
        cr(0, 2, i) = cr(0, 2, i) + (l+0.5d0)*j4*cls(ind_cross, l)
        cr(1, 2, i) = cr(1, 2, i) - (l+0.5d0)*j4*l2cls(ind_cross, l)        
     enddo
     cr(:,:,i) = cr(:,:,i)/coop_2pi

     call coop_gaussian_radial_modes_I(weights, cr(:,:,i), frI(:, i))
     call coop_gaussian_radial_modes_QU(weights, cr(:,:,i), frQU(:, i))
  enddo

  if(.not. hot)then
     frI(0,:) = -frI(0,:)
     frQU(1,:) = -frQU(1,:)
  endif
  call patchI%init(trim(field), n, dr)
  select case(trim(field))
  case("T", "I")
     call patchQU%init("QTUT", n, dr)
     call patchQrUr%init("QTrUTr", n, dr)
  case("E", "B")
     call patchQU%init("QU", n, dr)
     call patchQrUr%init("QrUr", n, dr)
  case("zeta")
     call patchQU%init("QZUZ", n, dr)
     call patchQrUr%init("QZrUZr", n, dr)
  case ("LT")
     call patchQU%init("QLTULT", n, dr)
     call patchQrUr%init("QLTrULTr", n, dr)     
  case default
     write(*,*) trim(field)
     stop "unknown field name."
  end select
     
  patchI%color_table = trim(colortable)
  patchQU%color_table = trim(colortable)
  patchQrUr%color_table = trim(colortable)    
  patchI%nstack = 1.d0
  patchQU%nstack = 1.d0
  patchQrUr%nstack = 1.d0
  patchI%nstack_raw = 1.d0
  patchQU%nstack_raw = 1.d0
  patchQrUr%nstack_raw = 1.d0

  do i=-n, n
     do j=-n, n
        if(i.eq.0.and.j.eq.0)then
           patchI%image(i,j,1) = frI(0,0)
           patchQU%image(i,j,1) = frQU(0,0)
           patchQU%image(i,j,2) = 0.d0
           patchQrUr%image(i,j,1) = frQU(0,0)
           patchQrUr%image(i,j,2) = 0.d0
        else
           phi = atan2(dble(j), dble(i))
           omega = sqrt(dble(i)**2+dble(j)**2)
           iomega = floor(omega)
           romega = omega - iomega
           patchI%image(i,j,1) = (frI(0, iomega)+ frI(1, iomega)*cos(2*phi) + frI(2, iomega)*cos(4*phi) ) *(1.d0-romega)+(frI(0, iomega+1)+ frI(1, iomega+1)*cos(2*phi) + frI(2, iomega+1)*cos(4*phi) ) * romega
           patchQU%image(i, j,1) = (frQU(0, iomega)+ frQU(1, iomega)*cos(2*phi) + frQU(2, iomega)*cos(4*phi) ) *(1.d0-romega)+(frQU(0, iomega+1)+ frQU(1, iomega+1)*cos(2*phi) + frQU(2, iomega+1)*cos(4*phi) ) * romega
           patchQU%image(i, j, 2) = (frQU(1, iomega)*sin(2*phi) + frQU(2, iomega)*sin(4*phi) ) *(1.d0-romega)+(frQU(0, iomega+1)+ frQU(1, iomega+1)*sin(2*phi) + frQU(2, iomega+1)*sin(4*phi) ) * romega
           patchQrUr%image(i, j, 1) = -patchQU%image(i, j,1)*cos(2*phi) - patchQU%image(i, j, 2)*sin(2*phi)  !!urgh, the stupid minus sign...
           patchQrUr%image(i, j, 2) = patchQU%image(i, j,1)*sin(2*phi) - patchQU%image(i, j, 2)*cos(2*phi)
        endif
     enddo
  enddo
  
  patchI%caption = "$\Lambda$CDM, $\nu="//COOP_STR_OF(nint(nu))//"$"
  patchQU%caption = "$\Lambda$CDM, $\nu="//COOP_STR_OF(nint(nu))//"$"
  patchQrUr%caption = "$\Lambda$CDM, $\nu="//COOP_STR_OF(nint(nu))//"$"
  call patchI%plot(1, trim(output_prefix)//"_I_stack.txt")
  call patchQU%plot(1, trim(output_prefix)//"_Q_stack.txt")
  call patchQU%plot(2, trim(output_prefix)//"_U_stack.txt")
  call patchQrUr%plot(1, trim(output_prefix)//"_Qr_stack.txt")
  call patchQrUr%plot(2, trim(output_prefix)//"_Ur_stack.txt")
  
  !!test the integrator
  call patchI%get_all_radial_profiles()
  call patchQU%get_all_radial_profiles()
  call fpfr%open(trim(output_prefix)//"frI.txt", "w")
  write(fpfr%unit, *) patchI%fr
  call fpfr%close()
  call fpfr%open(trim(output_prefix)//"frQU.txt", "w")
  write(fpfr%unit, *) patchQU%fr
  call fpfr%close()
  do i=0, n
     do m=0,2
        write(fpI(m)%unit, "(2E16.7)") dr*i, frI(m, i)
        write(fpQU(m)%unit, "(2E16.7)") dr*i, frQU(m, i)
        write(fpIdat(m)%unit, "(2E16.7)") dr*i, frI(m, i)
        write(fpQUdat(m)%unit, "(2E16.7)") dr*i, frQU(m, i)
     enddo
  enddo
  do m=0,2
     call fpI(m)%close()
     call fpIdat(m)%close()
     call fpQU(m)%close()
     call fpQUdat(m)%close()
  enddo
  call figCr%open(trim(output_prefix)//"cr.txt")
  call figCr%init(xlabel="$\varpi$", ylabel="$c_{m,n}(\varpi)$")
  call figCr%curve(x = r, y = cr(0,0,:)/sigma0, color="red", linewidth=1.8, linetype="solid", legend="$c_{0,0}(\varpi)/\sigma_0$")
  call figCr%curve(x = r, y = cr(0,1,:)/sigma0, color="orange", linewidth=1.5, linetype="dotted", legend="$c_{0,1}(\varpi)/\sigma_0$")
  call figCr%curve(x = r, y = cr(0,2,:)/sigma0, color="violet", linewidth=1.2, linetype="dashed", legend="$c_{0,2}(\varpi)/\sigma_0$")  
  call figCr%curve(x = r, y = -cr(1,0,:)/sigma2, color="blue", linewidth=1.8, linetype="solid", legend="$-c_{1,0}(\varpi)/\sigma_2$")  
  
  call figCr%curve(x = r, y = -cr(1,1,:)/sigma2, color="black", linewidth=1.5, linetype="dotted", legend="$-c_{1,1}(\varpi)/\sigma_2$")  

  call figCr%curve(x = r, y = -cr(1,2,:)/sigma2, color="skblue", linewidth=1.2, linetype="dashed", legend="$-c_{1,2}(\varpi)/\sigma_2$")  
  
  call figCr%legend(0.5, 0.9, 1)
  call figCr%close()

#else
  stop "you need to install healpix"
#endif  

end program stackth
