program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  type(coop_cosmology_firstorder)::cosmology
  COOP_INT,parameter::lmax = 3000, nlogl = 150, nchi = 400
  COOP_REAL,dimension(:,:,:),allocatable::trans
  type(coop_file)::fp
  type(coop_asy)::fig
  COOP_INT::nl, il, l, ichi
  COOP_REAL,dimension(:,:),allocatable::ClC
  COOP_REAL,dimension(:,:,:),allocatable::ClZ  
  COOP_INT,dimension(:),allocatable::ells
  COOP_REAL,dimension(:),allocatable::rl
  COOP_REAL::chi(nchi), logl(nlogl), info(nchi, nlogl), vis(nchi), dphidlna(nchi)
  COOP_INT::ik, idense
  logical::do_noise = .false.
  
#if DO_ZETA_TRANS
  call cosmology%Set_Planck_bestfit()
  call cosmology%compute_source(0)
  print*, "z= 0.5, chi =", (cosmology%tau0 - cosmology%tauofa(1.d0/(1.d0 + 0.5d0)))/cosmology%distlss
  print*, "z = 1, chi = ", (cosmology%tau0 - cosmology%tauofa(1.d0/(1.d0 + 1.d0)))/cosmology%distlss
  print*, "z = 10, chi = ", (cosmology%tau0 - cosmology%tauofa(1.d0/(1.d0 + 10.d0)))/cosmology%distlss  
  print*, "z = ", cosmology%zre, ", chi = ", (cosmology%tau0 - cosmology%tauofa(1.d0/(1.d0 + cosmology%zre)))/cosmology%distlss
  allocate(trans(cosmology%source(0)%nsrc, coop_k_dense_fac, cosmology%source(0)%nk))
  call coop_set_uniform(nchi, chi, cosmology%distlss*1.d-2, cosmology%distlss*1.2d0)
  do ichi = 1, nchi
     if(chi(ichi) .lt. cosmology%tau0*0.99)then
        vis(ichi) = cosmology%visofa(cosmology%aoftau(cosmology%tau0 - chi(ichi)))
        
     else
        vis(ichi) = 0.d0
     endif
  enddo
  vis = max(vis, 0.01d0)
  nl = coop_nl_to_lmax(lmax)
  allocate(ells(nl), rl(nl))
  call coop_set_ells(ells, lmax)
  rl = dble(ells)
  allocate(ClC(3, nl), ClZ(3, nl, nchi))

  if(coop_file_exists("zetacorr.dat"))then
     call fp%open("zetacorr.dat", "ru")
     read(fp%unit) nl
     read(fp%unit) ells
     read(fp%unit) il
     if(il.ne.nchi) stop "seems you have changed nchi"
     read(fp%unit) chi     
     read(fp%unit) ClC
     read(fp%unit) ClZ  
     call fp%close()
  else
     do il = 1, nl
        l = ells(il)
        write(*,*) il, "/", nl, ", ell = ", l
        call cosmology%source(0)%get_transfer(l, trans)
        !!TT, EE, TE
        ClC(1, il) = sum(cosmology%source(0)%ws_dense * trans(coop_index_source_T, :, :)**2)*coop_4pi  
        ClC(2, il) = (l+2.d0)*(l+1.d0)*l*(l-1.d0)*sum(cosmology%source(0)%ws_dense * trans(coop_index_source_E,:,:)**2)*coop_4pi     
        ClC(3, il) = sqrt((l+2.d0)*(l+1.d0)*l*(l-1.d0))*sum(cosmology%source(0)%ws_dense * trans(coop_index_source_T, :, :)*trans(coop_index_source_E, :, :))*coop_4pi

        do ichi = 1, nchi
           coop_zeta_single_slice_chi = chi(ichi)
           !$omp parallel do private(ik , idense)
           do ik = 1, cosmology%source(0)%nk
              do idense = 1, coop_k_dense_fac
                 trans(coop_index_source_zeta, idense, ik) = -coop_jl(l, cosmology%source(0)%k_dense(idense, ik)*coop_zeta_single_slice_chi)
              enddo
           enddo
           !$omp end parallel do
           !!TZ, EZ, ZZ
           ClZ(1, il, ichi) = sum(cosmology%source(0)%ws_dense * trans(coop_index_source_T, :, :)*trans(coop_index_source_zeta,:,:))*coop_4pi
           ClZ(2, il, ichi) = sqrt((l+2.d0)*(l+1.d0)*l*(l-1.d0))*sum(cosmology%source(0)%ws_dense * trans(coop_index_source_E, :, :) * trans(coop_index_source_zeta,:,:))*coop_4pi
           ClZ(3, il, ichi) = cosmology%clzetazeta_at_R(l, coop_zeta_single_slice_chi) !sum(cosmology%source(0)%ws_dense * trans(coop_index_source_zeta, :, :)**2)*coop_4pi
        enddo
     enddo
     call fp%open("zetacorr.dat", "u")
     write(fp%unit) nl
     write(fp%unit) ells
     write(fp%unit) nchi
     write(fp%unit) chi
     write(fp%unit) ClC
     write(fp%unit) ClZ  
     call fp%close()
  endif
  !!make figure
  call coop_set_uniform(nlogl, logl, log10(2.d0), log10(1.d0*lmax))

  call fig%open("ztinfo.txt")
  call fig%init(xlabel = "$\chi/\chi_{\rm rec}$", ylabel = "$\log_{10}\ell$")
  do il = 1, nlogl
     do ichi = 1, nchi
        info(ichi, il) =  TS2N(ichi, logl(il))
     enddo
  enddo
  call coop_asy_density(fig, z = info, xmin = chi(1)/cosmology%distlss, xmax = chi(nchi)/cosmology%distlss, ymin = logl(1), ymax = logl(nlogl), label = "$(S/N)_{T}$",  color_table = "BWRainbow")  
  call fig%close()

  call fig%open("zeinfo.txt")
  call fig%init(xlabel = "$\chi/\chi_{\rm rec}$", ylabel = "$\log_{10}\ell$")
  do il = 1, nlogl
     do ichi = 1, nchi
        info(ichi, il) =  ES2N(ichi, logl(il))
     enddo
  enddo
  call coop_asy_density(fig, z = info, xmin = chi(1)/cosmology%distlss, xmax = chi(nchi)/cosmology%distlss, ymin = logl(1), ymax = logl(nlogl), label = "$(S/N)_{EE}$",  color_table = "BWRainbow")
  call fig%close()


  call fig%open("zetinfo.txt")
  call fig%init(xlabel = "$\chi/\chi_{\rm rec}$", ylabel = "$log_{10}\ell$")
  do il = 1, nlogl
     do ichi = 1, nchi
        info(ichi, il) = TES2N(ichi, logl(il))
     enddo
  enddo
  call coop_asy_density(fig, z = info, xmin = chi(1)/cosmology%distlss, xmax = chi(nchi)/cosmology%distlss, ymin = logl(1), ymax = logl(nlogl), label = "$(S/N)_{TE}$", color_table = "BWRainbow")
  call fig%close()  

  call fig%open("vis.txt")
  call fig%init(xlabel = "$\chi/\chi_{\rm rec}$", ylabel = "$\dot\kappa e^{-\kappa}$", ylog = .true.)
  call fig%curve(chi/cosmology%distlss, vis)
  call fig%close()

  
#else
  stop "you need to switch zeta transfer on (set DO_ZETA_TRANS=COOP_YES in include/constants.h)"
#endif

contains

  function ES2N(ichi, logl)
    COOP_REAL::logl, thisl, ES2N, CTT, CTE, CEE, CTZ, CEZ, CZZ, delta, coef_T, coef_E, S2
    COOP_INT::il, ichi
    COOP_REAL::al,  bl
    thisl = 10.d0**(logl)
    il = max(min(coop_left_index(nl, rl, thisl), nl-1),1)
    bl = (thisl - rl(il))/(rl(il+1)-rl(il))
    al = 1.d0 - bl
    CTT = ClC(1, il)*al + ClC(1, il+1)*bl
    CEE = ClC(2, il)*al + ClC(2, il+1)*bl
    CTE = ClC(3, il)*al + ClC(3, il+1)*bl
    CTZ = ClZ(1, il, ichi)*al + ClZ(1, il+1, ichi)*bl
    CEZ = ClZ(2, il, ichi)*al + ClZ(2, il+1, ichi)*bl
    CZZ = ClZ(3, il, ichi)*al + ClZ(3, il+1, ichi)*bl
    if(do_noise)then
       CTT = CTT + coop_Planck_TNoise(thisl)
       CEE = CEE + coop_Planck_ENoise(thisl)
    endif
    
    coef_E = CEZ/CEE
    S2 = coef_E**2*CEE
    ES2N = sqrt(min(max(S2/(CZZ-S2), 1.d-4), 1.d4))
  end function ES2N

  function TS2N(ichi, logl)
    COOP_REAL::logl, thisl, TS2N, CTT, CTE, CEE, CTZ, CEZ, CZZ, delta, coef_T, coef_E, S2
    COOP_INT::il, ichi
    COOP_REAL::al,  bl
    thisl = 10.d0**(logl)
    il = max(min(coop_left_index(nl, rl, thisl), nl-1),1)
    bl = (thisl - rl(il))/(rl(il+1)-rl(il))
    al = 1.d0 - bl
    CTT = ClC(1, il)*al + ClC(1, il+1)*bl
    CEE = ClC(2, il)*al + ClC(2, il+1)*bl
    CTE = ClC(3, il)*al + ClC(3, il+1)*bl
    CTZ = ClZ(1, il, ichi)*al + ClZ(1, il+1, ichi)*bl
    CEZ = ClZ(2, il, ichi)*al + ClZ(2, il+1, ichi)*bl
    CZZ = ClZ(3, il, ichi)*al + ClZ(3, il+1, ichi)*bl
    if(do_noise)then
       CTT = CTT + coop_Planck_TNoise(thisl)
       CEE = CEE + coop_Planck_ENoise(thisl)
    endif
    
    coef_T = CTZ/CTT
    S2 = coef_T**2*CTT
    TS2N = sqrt(min(max(S2/(CZZ-S2), 1.d-4), 1.d4))
  end function TS2N


  function TES2N(ichi, logl)
    COOP_REAL::logl, thisl, TES2N, CTT, CTE, CEE, CTZ, CEZ, CZZ, delta, coef_T, coef_E, S2
    COOP_INT::il, ichi
    COOP_REAL::al,  bl
    thisl = 10.d0**(logl)
    il = max(min(coop_left_index(nl, rl, thisl), nl-1),1)
    bl = (thisl - rl(il))/(rl(il+1)-rl(il))
    al = 1.d0 - bl
    CTT = ClC(1, il)*al + ClC(1, il+1)*bl
    CEE = ClC(2, il)*al + ClC(2, il+1)*bl
    CTE = ClC(3, il)*al + ClC(3, il+1)*bl
    CTZ = ClZ(1, il, ichi)*al + ClZ(1, il+1, ichi)*bl
    CEZ = ClZ(2, il, ichi)*al + ClZ(2, il+1, ichi)*bl
    CZZ = ClZ(3, il, ichi)*al + ClZ(3, il+1, ichi)*bl
    
    if(do_noise)then
       CTT = CTT + coop_Planck_TNoise(thisl)
       CEE = CEE + coop_Planck_ENoise(thisl)
    endif
    delta = max(CTT*CEE -  CTE**2, 1.d-40)
    coef_T = (CTZ*CEE - CEZ*CTE)/delta
    coef_E = (-CTZ*CTE + CEZ*CTT)/delta
    S2 =  CTT*coef_T**2 +CEE*coef_E**2 +  2.d0*coef_E*coef_T*CTE
    TES2N = sqrt(min(max(S2/(CZZ-S2), 1.d-4), 1.d4))
  end function TES2N
  

  function coop_Planck_TNoise(l) result(Nl)
    COOP_REAL l
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
    Nl = exp(Nl)/(2.726) **2
  end function Coop_Planck_TNoise

  function coop_Planck_ENoise(l) result(Nl)
    COOP_REAL l
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
    Nl = exp(Nl)/(2.726)**2
  end function Coop_Planck_ENoise
  
end program test
