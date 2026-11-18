!================================================
!   The case file: one namelist file, named by the environment variable
!   AMR_CASE, that holds every run-time choice of a run (see cases/):
!     &run          the case and the model (Registry_Mod)
!     &boundaries   the kind of each side of the domain (Boundaries_Mod)
!     &<case>       the parameters of the case, e.g. &rising_bubble
!     &<model>      options of the model, e.g. &pre2017
!   Each module reads its own group: it opens the file with Open_Case_File,
!   reads the group and hands the status to Close_Case_File, which stops on
!   an error. A missing group is an error only for a required one (the
!   case's own group). Without AMR_CASE every group keeps its defaults.
!================================================
MODULE Case_File_Mod
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Case_File_Name, Open_Case_File, Close_Case_File, To_Lower

CONTAINS

!== the case file named by AMR_CASE, '' without it
FUNCTION Case_File_Name() RESULT( file )
    CHARACTER(len=:), ALLOCATABLE :: file

    INTEGER :: length, stat

    CALL GET_ENVIRONMENT_VARIABLE( 'AMR_CASE', LENGTH=length, STATUS=stat )
    IF( stat /= 0 .OR. length == 0 ) THEN
        file = ''
        RETURN
    END IF
    ALLOCATE( CHARACTER(len=length) :: file )
    CALL GET_ENVIRONMENT_VARIABLE( 'AMR_CASE', file )

END FUNCTION Case_File_Name

!******************************************************************************

!== open the case file at its start; found = F without AMR_CASE
SUBROUTINE Open_Case_File( unit, found )
    INTEGER, INTENT(OUT) :: unit
    LOGICAL, INTENT(OUT) :: found

    INTEGER :: stat

    found = ( LEN(Case_File_Name()) > 0 )
    unit  = 0
    IF( .NOT. found ) RETURN
    OPEN( NEWUNIT=unit, FILE=Case_File_Name(), STATUS='OLD', ACTION='READ', IOSTAT=stat )
    IF( stat /= 0 ) THEN
        PRINT '(2A)', ' AMR_CASE: cannot open ', Case_File_Name()
        ERROR STOP 'AMR_CASE: cannot open the case file'
    END IF

END SUBROUTINE Open_Case_File

!******************************************************************************

!== close the case file after reading the namelist group `group` with
!   IOSTAT=stat, IOMSG=msg; stop on a read error, and on a missing group
!   if it is required
SUBROUTINE Close_Case_File( unit, stat, msg, group, required )
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: IOSTAT_END
    INTEGER,          INTENT(IN) :: unit, stat
    CHARACTER(len=*), INTENT(IN) :: msg, group
    LOGICAL,          INTENT(IN) :: required

    CLOSE( unit )
    IF( stat == 0 ) RETURN
    IF( stat == IOSTAT_END ) THEN
        IF( .NOT. required ) RETURN
        PRINT '(4A)', ' AMR_CASE: no namelist &', group, ' in ', Case_File_Name()
        ERROR STOP 'AMR_CASE: missing namelist group'
    END IF
    PRINT '(4A)', ' AMR_CASE: error in namelist &', group, ': ', TRIM(msg)
    ERROR STOP 'AMR_CASE: error in a namelist group'

END SUBROUTINE Close_Case_File

!******************************************************************************

!== str in lower case (case-file names and kinds are case-insensitive)
PURE FUNCTION To_Lower( str ) RESULT( low )
    CHARACTER(len=*), INTENT(IN) :: str
    CHARACTER(len=LEN(str)) :: low
    INTEGER :: i
    low = str
    DO i = 1, LEN(str)
        IF( str(i:i) >= 'A' .AND. str(i:i) <= 'Z' ) low(i:i) = ACHAR( IACHAR(str(i:i)) + 32 )
    END DO
END FUNCTION To_Lower

END MODULE Case_File_Mod
