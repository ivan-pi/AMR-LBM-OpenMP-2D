!================================================
!   The models and test cases known to the driver, by name
!
!   The run is chosen at run time by the namelist &run of the case file
!   (Case_File_Mod):
!       &run  case = 'rayleigh_taylor', model = 'pre2017' /
!   and the environment variable AMR_MODEL, if set, overrides the model, so
!   one case file can be run with every model. Without a case file the run
!   is the rising bubble with the JCP 2016 model.
!
!   A new model or case is a module with a constructor function returning
!   an LBM_Model (Model_Mod) or an AMR_Case (Case_Mod), plus one entry in
!   the lists and SELECT CASE blocks below.
!================================================
MODULE Registry_Mod
    USE Model_Mod, ONLY: LBM_Model
    USE Case_Mod,  ONLY: AMR_Case
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Read_Run, Select_Model, Select_Case

    CHARACTER(len=*), PARAMETER :: model_names = 'jcp2016, pre2017'
    CHARACTER(len=*), PARAMETER :: case_names  = 'rising_bubble, rayleigh_taylor'

CONTAINS

!== the names of the case and the model of this run (&run, AMR_MODEL)
SUBROUTINE Read_Run( case_name, model_name )
    USE Case_File_Mod, ONLY: Case_File_Name, Open_Case_File, Close_Case_File, To_Lower
    CHARACTER(len=:), ALLOCATABLE, INTENT(OUT) :: case_name, model_name

    CHARACTER(len=64)  :: case, model, model_env
    CHARACTER(len=256) :: msg
    INTEGER :: u, stat, length
    LOGICAL :: found

    NAMELIST /run/ case, model

    case  = 'rising_bubble'
    model = 'jcp2016'

    CALL Open_Case_File( u, found )
    IF( found ) THEN
        PRINT '(/2A)', 'Case file: ', Case_File_Name()
        msg = ''
        READ( u, NML=run, IOSTAT=stat, IOMSG=msg )
        CALL Close_Case_File( u, stat, msg, 'run', required=.FALSE. )
    END IF

    CALL GET_ENVIRONMENT_VARIABLE( 'AMR_MODEL', model_env, length, stat )
    IF( stat == 0 .AND. length > 0 ) model = model_env

    case_name  = To_Lower( TRIM(ADJUSTL(case)) )
    model_name = To_Lower( TRIM(ADJUSTL(model)) )

END SUBROUTINE Read_Run

!******************************************************************************

FUNCTION Select_Model( name ) RESULT( model )
    USE Model_JCP2016_Mod, ONLY: JCP2016_Model
    USE Model_PRE2017_Mod, ONLY: PRE2017_Model
    CHARACTER(len=*), INTENT(IN) :: name
    TYPE(LBM_Model) :: model

    SELECT CASE( name )
    CASE( 'jcp2016' )
        model = JCP2016_Model()
    CASE( 'pre2017' )
        model = PRE2017_Model()
    CASE DEFAULT
        PRINT '(4A)', ' unknown model ''', name, '''; known: ', model_names
        ERROR STOP 'Select_Model: unknown model'
    END SELECT

END FUNCTION Select_Model

!******************************************************************************

FUNCTION Select_Case( name ) RESULT( problem )
    USE Rising_Bubble_Mod,   ONLY: Rising_Bubble_Case
    USE Rayleigh_Taylor_Mod, ONLY: Rayleigh_Taylor_Case
    CHARACTER(len=*), INTENT(IN) :: name
    TYPE(AMR_Case) :: problem

    SELECT CASE( name )
    CASE( 'rising_bubble' )
        problem = Rising_Bubble_Case()
    CASE( 'rayleigh_taylor' )
        problem = Rayleigh_Taylor_Case()
    CASE DEFAULT
        PRINT '(4A)', ' unknown case ''', name, '''; known: ', case_names
        ERROR STOP 'Select_Case: unknown case'
    END SELECT

END FUNCTION Select_Case

END MODULE Registry_Mod
