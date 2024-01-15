!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!Module to READ in name/value pairs from a file, with each line of the form line 'name = value'
!!Should correctly interpret FITS headers
!!by Antony Lewis (http://cosmologist.info/). Released to the public domain.
!!This version March 2005.
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE coop_InIFile_mod
  use coop_wrapper_typedef
  IMPLICIT NONE
  INTEGER, PARAMETER :: Ini_max_name_len = 128
  INTEGER, PARAMETER :: Ini_max_string_len = 1024
  LOGICAL :: Ini_fail_on_not_found = .false.
  LOGICAL :: Ini_Echo_READ = .false.

  TYPE TNameValue
     !no known way to make CHARACTER string POINTERs..
     CHARACTER(Ini_max_name_len)  :: Name
     CHARACTER(Ini_max_string_len):: Value
  END TYPE TNameValue

  TYPE TNameValue_POINTER
     TYPE(TNameValue), POINTER :: P
  END TYPE TNameValue_POINTER

  TYPE TNameValueList
     INTEGER Count
     INTEGER Delta
     INTEGER Capacity
     TYPE(TNameValue_POINTER), DIMENSION(:), POINTER :: Items
  END TYPE TNameValueList

  TYPE TInIFile
     LOGICAL SlashComments
     TYPE (TNameValueList) :: L, READValues
  END TYPE TInIFile

  TYPE(TInIFile) :: DefIni

CONTAINS

  SUBROUTINE TNameValueList_Init(L)
    TYPE (TNameValueList) :: L
    L%Count = 0
    L%Capacity = 0
    L%Delta = 128
    NULLIFY(L%Items)
  END SUBROUTINE TNameValueList_Init

  SUBROUTINE TNameValueList_Clear(L)
    TYPE (TNameValueList) :: L
    INTEGER i, status
    DO i=L%count,1,-1
       DEALLOCATE (L%Items(i)%P, stat = status)
    ENDDO
    DEALLOCATE (L%Items, stat = status)
    CALL TNameValueList_Init(L)
  END SUBROUTINE TNameValueList_Clear

  SUBROUTINE TNameValueList_ValueOf(L, AName, AValue)
    TYPE (TNameValueList) :: L
    CHARACTER(LEN=*), INTENT(in) :: AName
    CHARACTER(LEN=*) :: AValue
    INTEGER i
    DO i=1, L%Count
       IF (trim(L%Items(i)%P%Name) .eq. trim(AName)) THEN
          AValue = L%Items(i)%P%Value 
          RETURN
       ENDIF
    ENDDO
    AValue = ''
  END SUBROUTINE TNameValueList_ValueOf

  SUBROUTINE TNameValueList_Add(L, AName, AValue)
    TYPE (TNameValueList) :: L
    CHARACTER(LEN=*), INTENT(in) :: AName, AValue
    IF (L%Count .EQ. L%Capacity) CALL TNameValueList_SetCapacity(L, L%Capacity + L%Delta)
    L%Count = L%Count + 1
    ALLOCATE(L%Items(L%Count)%P)
    L%Items(L%Count)%P%Name = AName
    L%Items(L%Count)%P%Value = AValue
  END SUBROUTINE TNameValueList_Add

  SUBROUTINE TNameValueList_SetCapacity(L, C)
    TYPE (TNameValueList) :: L
    INTEGER C
    TYPE(TNameValue_POINTER), DIMENSION(:), POINTER :: TmpItems
    IF (L%Count .GT. 0) THEN
       IF (C .LT. L%Count) STOP 'TNameValueList_SetCapacity: smaller than Count'
       ALLOCATE(TmpItems(L%Count))
       TmpItems = L%Items(1:L%Count)
       DEALLOCATE(L%Items)
       ALLOCATE(L%Items(C))
       L%Items(1:L%Count) = TmpItems
       DEALLOCATE(TmpItems)
    ELSE
       ALLOCATE(L%Items(C))
    ENDIF
    L%Capacity = C
  END SUBROUTINE TNameValueList_SetCapacity

  SUBROUTINE TNameValueList_Delete(L, i)
    TYPE (TNameValueList) :: L
    INTEGER, INTENT(in) :: i
    DEALLOCATE(L%Items(i)%P)
    IF (L%Count .GT. 1) L%Items(i:L%Count-1) = L%Items(i+1:L%Count)
    L%Count = L%Count -1
  END SUBROUTINE TNameValueList_Delete

  SUBROUTINE Ini_NameValue_Add(Ini,AInLine)
    TYPE(TInIFile) :: Ini
    CHARACTER (LEN=*), INTENT(IN) :: AInLine
    INTEGER EqPos, slashpos, lastpos
    CHARACTER (LEN=len(AInLine)) :: AName, S, InLine
    InLine=trim(adjustl(AInLine))
    EqPos = scan(InLine,'=')
    IF (EqPos/=0 .and. InLine(1:1)/='#' .and. InLine(1:7) /= 'COMMENT' ) THEN
       AName = trim(InLine(1:EqPos-1))
       S = adjustl(InLine(EqPos+1:)) 
       IF (Ini%SlashComments) THEN
          slashpos=scan(S,'/')
          IF (slashpos /= 0) THEN
             S  = S(1:slashpos-1)
          ENDIF
       ENDIF
       lastpos=len_trim(S)
       IF (lastpos>1) THEN
          IF (S(1:1)=='''' .and. S(lastpos:lastpos)=='''') THEN
             S = S(2:lastpos-1)
          ENDIF
       ENDIF
       CALL TNameValueList_Add(Ini%L, AName, S)
    ENDIF
  END SUBROUTINE Ini_NameValue_Add


  SUBROUTINE Ini_OPEN(filename, unit_id,  error, slash_comments)
    CHARACTER (LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: unit_id
    LOGICAL, OPTIONAL, INTENT(OUT) :: error
    LOGICAL, OPTIONAL, INTENT(IN) :: slash_comments
    LOGICAL aerror
    CALL TNameValueList_Init(DefIni%L)
    CALL TNameValueList_Init(DefIni%READValues)
    IF (PRESENT(slash_comments)) THEN
       CALL Ini_OPEN_File(DefIni,filename,unit_id,aerror,slash_comments)
    ELSE
       CALL Ini_OPEN_File(DefIni,filename,unit_id,aerror)
    ENDIF
    IF (PRESENT(error)) THEN
       error = aerror
    ELSE
       IF (aerror) THEN
          WRITE (*,*) 'Ini_OPEN: Error OPENing file ' // trim(filename)
          STOP
       ENDIF
    ENDIF
  END SUBROUTINE Ini_OPEN


  SUBROUTINE Ini_OPEN_File(Ini, filename, unit_id,  error, slash_comments)
    TYPE(TInIFile) :: Ini
    CHARACTER (LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: unit_id
    LOGICAL, INTENT(OUT) :: error
    LOGICAL, OPTIONAL, INTENT(IN) :: slash_comments
    CHARACTER (LEN=256) :: InLine
    CALL TNameValueList_Init(Ini%L)
    CALL TNameValueList_Init(Ini%READValues)
    IF (PRESENT(slash_comments)) THEN
       Ini%SlashComments = slash_comments
    ELSE
       Ini%SlashComments = .false.
    ENDIF
    OPEN(unit=unit_id,file=filename,form='formatted',status='old', err=500)
    DO 
       READ (unit_id,'(a)',END=400) InLine
       InLine=trim(adjustl(InLine))
       IF (InLine == 'END') exit;
       IF (InLine /= '') CALL Ini_NameValue_Add(Ini,InLine) 
    ENDDO
400 close(unit_id)
    error=.false.
    RETURN
500 error=.true.
  END SUBROUTINE Ini_OPEN_File

  SUBROUTINE Ini_OPEN_Fromlines(Ini, Lines, NumLines, slash_comments)
    TYPE(TInIFile) :: Ini
    INTEGER, INTENT(IN) :: NumLines
    CHARACTER (LEN=*), DIMENSION(NumLines), INTENT(IN) :: Lines
    LOGICAL, INTENT(IN) :: slash_comments
    INTEGER i
    CALL TNameValueList_Init(Ini%L)
    Ini%SlashComments = slash_comments
    DO i=1,NumLines
       CALL Ini_NameValue_Add(Ini,Lines(i))
    ENDDO
  END  SUBROUTINE Ini_OPEN_Fromlines

  SUBROUTINE Ini_Close
    CALL Ini_close_File(DefIni)
  END SUBROUTINE Ini_Close


  SUBROUTINE Ini_Close_File(Ini)
    TYPE(TInIFile) :: Ini
    CALL TNameValueList_Clear(Ini%L)
    CALL TNameValueList_Clear(Ini%READValues)
  END  SUBROUTINE Ini_Close_File


  FUNCTION Ini_READ_String(Key, NotFoundFail) result(AValue)
    CHARACTER (LEN=*), INTENT(IN) :: Key
    LOGICAL, OPTIONAL, INTENT(IN) :: NotFoundFail
    CHARACTER(LEN=Ini_max_string_len) :: AValue
    IF (PRESENT(NotFoundFail)) THEN
       AValue = Ini_READ_String_File(DefIni, Key, NotFoundFail)
    ELSE
       AValue = Ini_READ_String_File(DefIni, Key)
    ENDIF
  END FUNCTION Ini_READ_String


  FUNCTION Ini_READ_String_File(Ini, Key, NotFoundFail) result(AValue)
    TYPE(TInIFile) :: Ini
    CHARACTER (LEN=*), INTENT(IN) :: Key
    LOGICAL, OPTIONAL, INTENT(IN) :: NotFoundFail
    CHARACTER(LEN=Ini_max_string_len) :: AValue
    CALL TNameValueList_ValueOf(Ini%L, Key, AValue)
    IF (AValue/='') THEN
       CALL  TNameValueList_Add(Ini%READValues, Key, AValue)
       IF (Ini_Echo_READ) WRITE (*,*) trim(Key)//' = ',trim(AValue)
       RETURN
    ENDIF
    IF (Ini_fail_on_not_found) THEN
       WRITE(*,*) 'key not found : '//Key
       STOP
    ENDIF
    IF (PRESENT(NotFoundFail)) THEN
       IF (NotFoundFail) THEN
          WRITE(*,*) 'key not found : '//Key
          STOP
       ENDIF
    ENDIF
  END FUNCTION Ini_READ_String_File


  FUNCTION Ini_READ_Int(Key, Default)
    INTEGER, OPTIONAL, INTENT(IN) :: Default
    CHARACTER (LEN=*), INTENT(IN) :: Key
    INTEGER Ini_READ_Int
    IF (PRESENT(Default)) THEN
       Ini_READ_Int = Ini_READ_Int_File(DefIni, Key, Default)
    ELSE
       Ini_READ_Int = Ini_READ_Int_File(DefIni, Key)
    ENDIF
  END FUNCTION Ini_READ_Int


  FUNCTION Ini_READ_Int_File(Ini, Key, Default)
    TYPE(TInIFile) :: Ini
    INTEGER Ini_READ_Int_File
    INTEGER, OPTIONAL, INTENT(IN) :: Default
    CHARACTER  (LEN=*), INTENT(IN) :: Key
    CHARACTER(LEN=Ini_max_string_len) :: S
    S = Ini_READ_String_File(Ini, Key,.not. PRESENT(Default))
    IF (S .EQ. '') THEN
       IF (.not. PRESENT(Default)) THEN
          WRITE(*,*) 'no value for key: '//Key
          STOP
       ENDIF
       Ini_READ_Int_File = Default
       WRITE (S,*) Default
       CALL  TNameValueList_Add(Ini%READValues, Key, S)
    ELSE
       IF (VERIFY(TRIM(S),'-+0123456789') /= 0) goto 10
       READ (S,*, err = 10) Ini_READ_Int_File
    ENDIF
    RETURN
10  WRITE (*,*) 'error READing INTEGER for key: '//Key
    STOP
  END FUNCTION Ini_READ_Int_File



  FUNCTION Ini_READ_Int_Arr(Key,n)
    Integer n
    CHARACTER (LEN=*), INTENT(IN) :: Key
    INTEGER Ini_READ_Int_Arr(n)
    Ini_READ_Int_Arr = Ini_READ_Int_Arr_File(DefIni, Key,n)
  END FUNCTION Ini_READ_Int_Arr


  FUNCTION Ini_READ_Int_Arr_File(Ini, Key, n)
    TYPE(TInIFile) :: Ini
    Integer n
    INTEGER Ini_READ_Int_Arr_File(n)
    CHARACTER  (LEN=*), INTENT(IN) :: Key
    CHARACTER(LEN=Ini_max_string_len) :: S
    S = Ini_READ_String_File(Ini, Key,.true.)
    IF (S .EQ. '') THEN
       WRITE(*,*) 'no value for key: '//Key
       STOP
    ENDIF
    IF (VERIFY(TRIM(S),'-+0123456789 ') /= 0) goto 10
    READ (S,*, err = 10) Ini_READ_Int_Arr_File(1:n)
    RETURN
10  WRITE (*,*) 'error READing INTEGER for key: '//Key
    STOP
  END FUNCTION Ini_READ_Int_Arr_File


  FUNCTION Ini_READ_DOUBLE(Key, Default)
    DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: Default
    CHARACTER (LEN=*), INTENT(IN) :: Key
    DOUBLE PRECISION Ini_READ_DOUBLE
    IF (PRESENT(Default)) THEN
       Ini_READ_DOUBLE = Ini_READ_DOUBLE_File(DefIni, Key, Default)
    ELSE
       Ini_READ_DOUBLE = Ini_READ_DOUBLE_File(DefIni, Key)
    ENDIF
  END FUNCTION Ini_READ_DOUBLE


  FUNCTION Ini_READ_DOUBLE_File(Ini,Key, Default)
    TYPE(TInIFile) :: Ini
    DOUBLE PRECISION Ini_READ_DOUBLE_File 
    DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: Default
    CHARACTER (LEN=*), INTENT(IN) :: Key
    CHARACTER(LEN=Ini_max_string_len) :: S
    S = Ini_READ_String_File(Ini,Key,.not. PRESENT(Default))
    IF (S == '') THEN
       IF (.not. PRESENT(Default)) THEN
          WRITE(*,*) 'no value for key: '//Key
          STOP
       ENDIF
       Ini_READ_DOUBLE_File = Default
       WRITE (S,*) Default
       CALL  TNameValueList_Add(Ini%READValues, Key, S)
    ELSE
       READ (S,*, err=10) Ini_READ_DOUBLE_File
    ENDIF
    RETURN
10  WRITE (*,*) 'error READing DOUBLE for key: '//Key
    STOP
  END FUNCTION Ini_READ_DOUBLE_File



  FUNCTION Ini_READ_DOUBLE_Arr(Key,n)
    Integer n
    CHARACTER (LEN=*), INTENT(IN) :: Key
    DOUBLE PRECISION Ini_READ_DOUBLE_arr(n)
    Ini_READ_DOUBLE_arr = Ini_READ_DOUBLE_Arr_File(DefIni, Key,n)
  END FUNCTION Ini_READ_DOUBLE_Arr


  FUNCTION Ini_READ_DOUBLE_Arr_File(Ini,Key,n)
    Integer n
    TYPE(TInIFile) :: Ini
    DOUBLE PRECISION Ini_READ_DOUBLE_Arr_File(n)
    CHARACTER (LEN=*), INTENT(IN) :: Key
    CHARACTER(LEN=Ini_max_string_len) :: S
    S = Ini_READ_String_File(Ini,Key,.true.)
    IF (S == '') THEN
       WRITE(*,*) 'no value for key: '//Key
       STOP
    ENDIF
    READ (S,*, err=10) Ini_READ_DOUBLE_Arr_File(1:n)
    RETURN
10  WRITE (*,*) 'error READing DOUBLE for key: '//Key
    STOP
  END FUNCTION Ini_READ_DOUBLE_Arr_File



  FUNCTION Ini_READ_REAL(Key, Default)
    REAL, OPTIONAL, INTENT(IN) :: Default
    CHARACTER (LEN=*), INTENT(IN) :: Key
    REAL Ini_READ_REAL
    IF (PRESENT(Default)) THEN
       Ini_READ_REAL = Ini_READ_REAL_File(DefIni, Key, Default)
    ELSE
       Ini_READ_REAL = Ini_READ_REAL_File(DefIni, Key)
    ENDIF
  END FUNCTION Ini_READ_REAL

  FUNCTION Ini_READ_REAL_File(Ini,Key, Default)
    TYPE(TInIFile) :: Ini
    REAL Ini_READ_REAL_File 
    REAL, OPTIONAL, INTENT(IN) :: Default
    CHARACTER (LEN=*), INTENT(IN) :: Key
    CHARACTER(LEN=Ini_max_string_len) :: S
    S = Ini_READ_String_File(Ini,Key,.not. PRESENT(Default))
    IF (S == '') THEN
       IF (.not. PRESENT(Default)) THEN
          WRITE(*,*) 'no value for key: '//Key
          STOP
       ENDIF
       Ini_READ_REAL_File = Default
       WRITE (S,*) Default
       CALL  TNameValueList_Add(Ini%READValues, Key, S)
    ELSE
       READ (S,*, err=10) Ini_READ_REAL_File
    ENDIF
    RETURN
10  WRITE (*,*) 'error READing DOUBLE for key: '//Key
    STOP
  END FUNCTION Ini_READ_REAL_File


  FUNCTION Ini_READ_LOGICAL(Key, Default)
    LOGICAL, OPTIONAL, INTENT(IN) :: Default
    CHARACTER (LEN=*), INTENT(IN) :: Key
    LOGICAL Ini_READ_LOGICAL
    IF (PRESENT(Default)) THEN
       Ini_READ_LOGICAL = Ini_READ_LOGICAL_File(DefIni, Key, Default)
    ELSE
       Ini_READ_LOGICAL = Ini_READ_LOGICAL_File(DefIni, Key)
    ENDIF
  END FUNCTION Ini_READ_LOGICAL

  FUNCTION Ini_READ_LOGICAL_File(Ini, Key, Default)
    TYPE(TInIFile) :: Ini
    LOGICAL Ini_READ_LOGICAL_File
    LOGICAL, OPTIONAL, INTENT(IN) :: Default
    CHARACTER  (LEN=*), INTENT(IN) :: Key
    CHARACTER(LEN=Ini_max_string_len) :: S
    S = Ini_READ_String_File(Ini,Key,.not. PRESENT(Default))
    IF (S == '') THEN
       IF (.not. PRESENT(Default)) THEN
          WRITE(*,*) 'no value for key: '//Key
          STOP
       ENDIF
       Ini_READ_LOGICAL_File = Default
       WRITE (S,*) Default
       CALL  TNameValueList_Add(Ini%READValues, Key, S)
    ELSE
       IF (verIFy(trim(S),'10TF') /= 0) goto 10  
       READ (S,*, err = 10) Ini_READ_LOGICAL_File
    ENDIF
    RETURN
10  WRITE (*,*) 'error READing LOGICAL for key: '//Key
    STOP
  END FUNCTION Ini_READ_LOGICAL_File


  SUBROUTINE Ini_SaveREADValues(afile,unit_id)
    CHARACTER(LEN=*)  :: afile
    INTEGER, INTENT(in) :: unit_id
    CALL Ini_SaveREADValues_File(DefIni, afile, unit_id)
  END SUBROUTINE Ini_SaveREADValues


  SUBROUTINE Ini_SaveREADValues_File(Ini, afile, unit_id)
    TYPE(TInIFile) :: Ini
    CHARACTER(LEN=*), INTENT(in) :: afile
    INTEGER, INTENT(in) :: unit_id
    INTEGER i
    OPEN(unit=unit_id,file=afile,form='formatted',status='replace', err=500)
    DO i=1, Ini%READValues%Count
       WRITE (unit_id,'(a)') trim(Ini%READValues%Items(i)%P%Name) // ' = ' &
            //trim(Ini%READValues%Items(i)%P%Value)
    ENDDO
    close(unit_id)
    RETURN
500 WRITE(*,*) 'Ini_SaveREADValues_File: Error creating '//trim(afile)
  END SUBROUTINE Ini_SaveREADValues_File

  FUNCTION GetIniParams(N)
    CHARACTER(LEN=120) GetIniParams
    !INTEGER IARGC
    INTEGER N
    IF(IARGC().ge. N)THEN
       call GetArg(N,GetIniParams)
    ELSE
       Write(*,*) "Enter the initial parameters:"
       READ(*,'(A)')GetIniParams
    ENDIF
  END FUNCTION GetIniParams

END MODULE Coop_InIFile_mod
