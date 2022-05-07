module coop_latex_mod
  use coop_wrapper_utils
  implicit none
  

#include "constants.h"

contains


  function latex_range(base, lowsig, upsig) result(str)
    real base, lowsig(:), upsig(:)
    COOP_STRING str
    if(lowsig(1) .gt. 0. .and. upsig(1) .gt. 0.)then
       str = getProperstr(dble(base), dble(upsig(1)), dble(lowsig(1)))
    elseif(lowsig(1) .gt. 0.)then
       str = getProperstr_Onetail(dble(base), dble(lowsig(1)), dble(lowsig(2)), no_left_tail = .false.)
    else
       str = getProperstr_Onetail(dble(base), dble(upsig(1)), dble(upsig(2)), no_left_tail = .true.)
    endif
  end function latex_range



  function ReadDigits(Number,Digit)
    COOP_REAL Number
    COOP_REAL  X
    integer Digit
    COOP_SHORT_STRING ReadDigits
    X=abs(Number)
    if(Digit.eq.0)THEN
       ReadDigits=Trim(coop_num2str(NINT(X)))
    else
       ReadDigits=Trim(coop_num2str(X,"(F15."//Trim(coop_num2str(Digit))//")"))
    endif
    if(ReadDigits(1:1).eq.".")then
       ReadDigits="0"//Trim(ReadDigits)
    endif
    if(number .lt.0.)then
       ReadDigits="-"//Trim(ReadDigits)
    end if
  end function ReadDigits



  function getProperstr(Mean,Psigma,Msigma)
    COOP_REAL  Mean,Psigma,Msigma
    COOP_REAL  X
    COOP_STRING getProperstr
    integer digitP,digitM,digit
    if(Psigma.lT.1.e-20 .OR. Psigma.gt.1.e20 .OR. Msigma.lT.1.e-20 .OR. Msigma.gt.1.e20)Then
       print*,"sigmas too large or too small",Psigma,Msigma
       getProperstr=""
       return
    endif
    X=Psigma
    digitP=0
    do While(X.lT.0.95)
       X=X*10.
       digitP=digitP+1
    enddo
    if(abs(X-nint(X)).ge. 0.05*X)digitP=digitP+1
    X=Msigma
    digitM=0
    do While(X.lT.0.95)
       X=X*10.
       digitM=digitM+1
    enddo
    if(abs(X-nint(X)).ge. 0.05*X)digitM=digitM+1
    digit=MaX(digitM,digitP)
    X=abs(Mean-nint(Mean*10.**digit)/10.**digit)
    if(X.gt.Min(Psigma,Msigma)*0.1)digit=digit+1
    getProperstr="$"//TRiM(Readdigits(Mean,digit))//  &
         "^{+"//TRiM(Readdigits(Psigma,digit))// &
         "}_{-"//TRiM(Readdigits(Msigma,digit))//"}$"
  end function getProperstr

  
  function getProperstr_OneTail(bound,sigma,sigma2,No_left_tail)
    COOP_REAL  bound,sigma,sigma2
    logical No_left_tail
    COOP_REAL  X
    COOP_STRING getProperstr_OneTail
    integer digit,digit2
    if(sigma2.lT.1.e-8 .OR. sigma2.gt.1.e8 .or. sigma.lt.1.e-8 .or. sigma .ge. 1.e8)Then
       print*,"sigma too large or too small",sigma,sigma2
       getProperstr_OneTail=""
       ReTURn
    endif
    X=sigma
    digit=0
    do While(X.lT.0.95)
       X=X*10.
       digit=digit+1
    enddo
    if(abs(X-nint(X)).ge. 0.05*X)digit=digit+1
    X=sigma2
    digit2=0
    do While(X.lT.0.95)
       X=X*10.
       digit2=digit2+1
    enddo
    if(abs(X-nint(X)).ge. 0.05*X)digit2=digit2+1
    digit=Max(digit,digit2)
    X=abs(bound-nint(bound*10.**digit)/10.**digit)
    if(X.gt.sigma*0.1)digit=digit+1
    if(No_left_tail)Then
       getProperstr_OneTail="$"//TRiM(Readdigits(bound,digit))//  &
            "^{+"//TRiM(Readdigits(sigma,digit))//"+"//Trim(Readdigits(sigma2,digit))//"}$"
    else
       getProperstr_OneTail="$"//TRiM(Readdigits(bound,digit))//  &
            "^{-"//TRiM(Readdigits(sigma,digit))//"-"//Trim(Readdigits(sigma2,digit))//"}$"
    endif

  end function getProperstr_OneTail


  
end module coop_latex_mod
