!================================================
!   2D two-phase AMR-LBM solver
!
!   One driver for every model and test case: the case file (environment
!   variable AMR_CASE) names both in its namelist &run, and AMR_MODEL can
!   override the model (Registry_Mod). Without a case file it runs the
!   rising bubble with the mass-conserving phase-field model of Fakhari,
!   Geier & Lee (JCP 2016).
!
!   usage: [AMR_CASE=file] [AMR_MODEL=name] amr_lbm [number_of_time_steps [output_interval_in_steps]]
!================================================
PROGRAM AMR_LBM_2D
    USE OMP_LIB, ONLY: OMP_GET_MAX_THREADS, OMP_GET_NUM_PROCS
    USE Solver_Mod, ONLY: Run_AMR_LBM
    USE Registry_Mod, ONLY: Read_Run, Select_Model, Select_Case
    IMPLICIT NONE

    CHARACTER(len=:), ALLOCATABLE :: case_name, model_name
    INTEGER(8) :: t1, t2, Rate

    ! The number of threads is controlled with the OMP_NUM_THREADS
    ! environment variable (default: all available hardware threads).
    PRINT '(/A,I5)', 'OpenMP number of threads used in this run:', OMP_GET_MAX_THREADS()
    PRINT '( A,I5)', 'Number of available processors:           ', OMP_GET_NUM_PROCS()

    CALL SYSTEM_CLOCK( t1, Rate )

    CALL Read_Run( case_name, model_name )
    CALL Run_AMR_LBM( Select_Model( model_name ), Select_Case( case_name ) )

    CALL SYSTEM_CLOCK( t2 )
    PRINT '(/A)',' *****************************************'
    PRINT '(A,F12.1)', ' Time Elapsed:', DBLE(t2-t1)/Rate
    PRINT '( A)',' *****************************************'

END PROGRAM AMR_LBM_2D
