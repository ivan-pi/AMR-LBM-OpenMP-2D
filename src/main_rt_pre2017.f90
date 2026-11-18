!================================================
!   2D two-phase AMR-LBM solver: Rayleigh-Taylor instability with the
!   improved-locality phase-field model of
!   Fakhari, Mitchell, Leonardi & Bolster (PRE 2017)
!================================================
PROGRAM Rayleigh_Taylor_AMR_2D_PRE2017
    USE OMP_LIB, ONLY: OMP_GET_MAX_THREADS, OMP_GET_NUM_PROCS
    USE Solver_Mod, ONLY: Run_AMR_LBM, Read_Command_Line
    USE Rayleigh_Taylor_Mod, ONLY: Rayleigh_Taylor_Case, Light_Phase_Refinement
    USE Model_PRE2017_Mod, ONLY: PRE2017_Model
    IMPLICIT NONE

    INTEGER(8) :: t1, t2, Rate

    CALL Read_Command_Line( 'Rayleigh-Taylor instability, model of Fakhari et al. (PRE 2017)' )

    ! The number of threads is controlled with the OMP_NUM_THREADS
    ! environment variable (default: all available hardware threads).
    PRINT '(/A,I5)', 'OpenMP number of threads used in this run:', OMP_GET_MAX_THREADS()
    PRINT '( A,I5)', 'Number of available processors:           ', OMP_GET_NUM_PROCS()

    CALL SYSTEM_CLOCK( t1, Rate )

    ! the refinement criterion (interface only, or also the light phase)
    ! is chosen by refine_light_phase in the case file
    CALL Run_AMR_LBM( PRE2017_Model( Light_Phase_Refinement() ), Rayleigh_Taylor_Case() )

    CALL SYSTEM_CLOCK( t2 )
    PRINT '(/A)',' *****************************************'
    PRINT '(A,F12.1)', ' Time Elapsed:', DBLE(t2-t1)/Rate
    PRINT '( A)',' *****************************************'

END PROGRAM Rayleigh_Taylor_AMR_2D_PRE2017
