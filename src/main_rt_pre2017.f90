!================================================
!   2D two-phase AMR-LBM solver: Rayleigh-Taylor instability with the
!   improved-locality phase-field model of
!   Fakhari, Mitchell, Leonardi & Bolster (PRE 2017)
!================================================
PROGRAM Rayleigh_Taylor_AMR_2D_PRE2017
    USE OMP_LIB, ONLY: OMP_GET_MAX_THREADS, OMP_GET_NUM_PROCS
    USE Solver_Mod, ONLY: Run_AMR_LBM
    USE Rayleigh_Taylor_Mod, ONLY: Rayleigh_Taylor_Case
    USE Model_PRE2017_Mod, ONLY: PRE2017_Model
    IMPLICIT NONE

    INTEGER(8) :: t1, t2, Rate

    ! The number of threads is controlled with the OMP_NUM_THREADS
    ! environment variable (default: all available hardware threads).
    PRINT '(/A,I5)', 'OpenMP number of threads used in this run:', OMP_GET_MAX_THREADS()
    PRINT '( A,I5)', 'Number of available processors:           ', OMP_GET_NUM_PROCS()

    CALL SYSTEM_CLOCK( t1, Rate )

    CALL Run_AMR_LBM( PRE2017_Model(), Rayleigh_Taylor_Case() )

    CALL SYSTEM_CLOCK( t2 )
    PRINT '(/A)',' *****************************************'
    PRINT '(A,F12.1)', ' Time Elapsed:', DBLE(t2-t1)/Rate
    PRINT '( A)',' *****************************************'

END PROGRAM Rayleigh_Taylor_AMR_2D_PRE2017
