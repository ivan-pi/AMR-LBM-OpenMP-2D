!================================================
!   2D two-species AMR-LBM solver: droplets and spinodal decomposition with
!   the chemical-potential diffusion model of Parker & Jovanovic (2009)
!================================================
PROGRAM Multicomponent_AMR_2D_Parker2009
    USE OMP_LIB, ONLY: OMP_GET_MAX_THREADS, OMP_GET_NUM_PROCS
    USE Solver_Mod, ONLY: Run_AMR_LBM, Read_Command_Line
    USE Parker_Cases_Mod, ONLY: Parker_Case
    USE Model_Parker2009_Mod, ONLY: Parker2009_Model
    USE Kinds_Mod, ONLY: i64
    IMPLICIT NONE

    INTEGER(i64) :: t1, t2, Rate

    CALL Read_Command_Line( 'droplets and spinodal decomposition, model of Parker & Jovanovic (2009)' )

    ! The number of threads is controlled with the OMP_NUM_THREADS
    ! environment variable (default: all available hardware threads).
    PRINT '(/A,I5)', 'OpenMP number of threads used in this run:', OMP_GET_MAX_THREADS()
    PRINT '( A,I5)', 'Number of available processors:           ', OMP_GET_NUM_PROCS()

    CALL SYSTEM_CLOCK( t1, Rate )

    CALL Run_AMR_LBM( Parker2009_Model(), Parker_Case() )

    CALL SYSTEM_CLOCK( t2 )
    PRINT '(/A)',' *****************************************'
    PRINT '(A,F12.1)', ' Time Elapsed:', DBLE(t2-t1)/Rate
    PRINT '( A)',' *****************************************'

END PROGRAM Multicomponent_AMR_2D_Parker2009
