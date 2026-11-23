!================================================
!   2D two-phase AMR-LBM solver: Rayleigh-Taylor instability with the
!   mass-conserving phase-field model of Fakhari, Geier & Lee (JCP 2016)
!================================================
PROGRAM Rayleigh_Taylor_AMR_2D_JCP2016
    USE OMP_LIB, ONLY: OMP_GET_MAX_THREADS, OMP_GET_NUM_PROCS
    USE Solver_Mod, ONLY: Run_AMR_LBM
    USE Rayleigh_Taylor_Mod, ONLY: Rayleigh_Taylor_Case
    USE Model_JCP2016_Mod, ONLY: JCP2016_Model
    USE Kinds_Mod, ONLY: i64
    IMPLICIT NONE

    INTEGER(i64) :: t1, t2, Rate

    ! The number of threads is controlled with the OMP_NUM_THREADS
    ! environment variable (default: all available hardware threads).
    PRINT '(/A,I5)', 'OpenMP number of threads used in this run:', OMP_GET_MAX_THREADS()
    PRINT '( A,I5)', 'Number of available processors:           ', OMP_GET_NUM_PROCS()

    CALL SYSTEM_CLOCK( t1, Rate )

    CALL Run_AMR_LBM( JCP2016_Model(), Rayleigh_Taylor_Case() )

    CALL SYSTEM_CLOCK( t2 )
    PRINT '(/A)',' *****************************************'
    PRINT '(A,F12.1)', ' Time Elapsed:', DBLE(t2-t1)/Rate
    PRINT '( A)',' *****************************************'

END PROGRAM Rayleigh_Taylor_AMR_2D_JCP2016
