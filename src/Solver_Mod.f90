!================================================
!   The AMR-LBM solver: initialization and time loop for a given LBM model
!   (Model_Mod) and test case (Case_Mod); the drivers main_*.f90 select both.
!   The domain is set in Config_Mod: periodic sides, walls at the bottom
!   and top.
!================================================
MODULE Solver_Mod
    USE Kinds_Mod, ONLY: wp, dp, i64
    USE Mesh_Mod, ONLY: Create_Domain
    USE Refine_Mod, ONLY: Refine_Derefine_Mesh
    USE LBM_Common_Mod, ONLY: Perfect_Shift_AMR, Streaming_LaxWendroff_AMR, Macroscopic_Values_h_AMR
    USE Model_Mod, ONLY: LBM_Model
    USE Case_Mod, ONLY: AMR_Case
    USE Ghost_C_Mod, ONLY: Fill_Ghost_Cells_C
    USE Ghost_hg_Mod, ONLY: Fill_Ghost_Cells_hg
    USE Output_Mod, ONLY: Select_Writer, Write_Output, Close_Output
    USE Viewer_Mod, ONLY: Viewer_Open, Viewer_Due, Viewer_Update, Viewer_Close
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Run_AMR_LBM

CONTAINS

!== run a test case with an LBM model; the case may adjust the model's mesh
!   policy (tolerances, island refinement) after its setup
SUBROUTINE Run_AMR_LBM( model_in, problem )
    USE Config_Mod, ONLY: minLevels, maxLevels
    USE Blocks_Mod, ONLY: nScalars
    IMPLICIT NONE
    TYPE(LBM_Model), INTENT(IN) :: model_in
    TYPE(AMR_Case),  INTENT(IN) :: problem

    TYPE(LBM_Model) :: model
    INTEGER    :: L
    INTEGER(i64) :: t1, t2, Rate
    REAL(wp)   :: t_end, time_per_step

    model = model_in
    IF( .NOT. ( ASSOCIATED(model%Initialize) .AND. ASSOCIATED(model%Collide) .AND. &
                ASSOCIATED(model%Moments_g)  .AND. ASSOCIATED(model%Moments_and_Collide) .AND. &
                ASSOCIATED(model%Refinement_Error) ) ) ERROR STOP 'Run_AMR_LBM: incomplete model'
    IF( .NOT. ( ASSOCIATED(problem%Setup) .AND. ASSOCIATED(problem%Initial_State) .AND. &
                ASSOCIATED(problem%Diagnostics) ) ) ERROR STOP 'Run_AMR_LBM: incomplete case'

    CALL SYSTEM_CLOCK( t1, Rate )

    PRINT '(/2A)', 'Case: ', TRIM(problem%name)
    CALL problem%Setup( t_end, time_per_step )
    IF( ASSOCIATED(problem%Adjust_Model) ) CALL problem%Adjust_Model( model )

    nScalars = model%n_scalars      ! before the block pool is allocated
    CALL Create_Domain

    CALL Initialize_State( model, problem )

    !== refine level by level, from the initial state at each level
    DO L = minLevels, maxLevels-1
        CALL Refine_Derefine_Mesh( model )
        CALL Initialize_State( model, problem )
    END DO

    CALL SYSTEM_CLOCK( t2 )

    CALL Main_Loop( model, problem, t_end, time_per_step )

    PRINT '(A,F12.3)', ' t_initialization:', DBLE(t2-t1)/Rate

END SUBROUTINE Run_AMR_LBM

!******************************************************************************

!== initial state of all leaf blocks: C, P, U from the case, h and g from the model
SUBROUTINE Initialize_State( model, problem )
    USE Blocks_Mod, ONLY: nLeaf, LeafByLevel
    IMPLICIT NONE
    TYPE(LBM_Model), INTENT(IN) :: model
    TYPE(AMR_Case),  INTENT(IN) :: problem

    INTEGER :: m, n

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED( LeafByLevel, nLeaf, model, problem )
    DO m = 1, nLeaf
        n = LeafByLevel(m)
        CALL problem%Initial_State( n )
        CALL model%Initialize( n )
    END DO
!$OMP END PARALLEL DO

END SUBROUTINE Initialize_State

!******************************************************************************

SUBROUTINE Main_Loop( model, problem, t_end, time_per_step )
    USE Config_Mod, ONLY: nx, ny, maxLevels, inv_dx_level
    USE Blocks_Mod, ONLY: nLeaf, Blk_h, Blk_g, Blk_C, Blk_C2, LeafByLevel, LevelStart, LevelMin, LevelMax
    USE Physics_Mod, ONLY: W, M0, taul, tauh, wc, Sigma
    USE OMP_LIB, ONLY: OMP_GET_MAX_THREADS
    USE ISO_FORTRAN_ENV, ONLY: OUTPUT_UNIT
    USE Mesh_Mod, ONLY: t_nbr_table, n_nbr_rebuilds
    IMPLICIT NONE
    TYPE(LBM_Model), INTENT(IN) :: model
    TYPE(AMR_Case),  INTENT(IN) :: problem
    REAL(wp),        INTENT(IN) :: t_end            ! end time, in the case's time unit
    REAL(wp),        INTENT(IN) :: time_per_step    ! the case's time unit per time step

    INTEGER :: T, n, m, L, Tf, nOutStep, nargs, ios
    CHARACTER(len=32) :: arg
    REAL(wp)    :: cfl
    LOGICAL    :: changed, moments_done, is_output, is_diagnostic
    INTEGER    :: nDiagStep
    INTEGER(i64) :: t1, Rate
    INTEGER(i64) :: t_collision, t_advection, t_Output, t_Ghost_C, t_Ghost_hg, t_refine
    INTEGER    :: nProgStep
    INTEGER(i64) :: tp0, tpLast, cellUpdates

    !==================================
    Tf   = NINT(t_end / time_per_step /20)*20
    !== optional command-line override of the number of time steps
    nargs = COMMAND_ARGUMENT_COUNT()
    IF( nargs >= 1 ) THEN
        CALL GET_COMMAND_ARGUMENT( 1, arg )
        READ( arg, *, IOSTAT=ios ) Tf
        IF( ios /= 0 .OR. Tf < 1 ) STOP 'usage: amr_lbm [number_of_time_steps [output_interval_in_steps]]'
    END IF
    !== optional output interval in time steps (default: 20 outputs over the run);
    !   the initial and the final state are always written
    nOutStep = MAX( Tf/20, 1 )
    IF( nargs >= 2 ) THEN
        CALL GET_COMMAND_ARGUMENT( 2, arg )
        READ( arg, *, IOSTAT=ios ) nOutStep
        IF( ios /= 0 .OR. nOutStep < 1 ) STOP 'usage: amr_lbm [number_of_time_steps [output_interval_in_steps]]'
    END IF
    CALL Select_Writer
    CALL Viewer_Open
    !== the case's diagnostics are recorded 120 times per run and at every output
    nDiagStep = MAX( Tf/120, 1 )
    !== a progress line (speed, estimated time left) 20 times per run
    nProgStep = MAX( Tf/20, 1 )
    !==================================

    t_Output    = 0
    t_collision = 0
    t_advection = 0
    t_Ghost_C   = 0
    t_Ghost_hg  = 0
    t_refine    = 0

    !==========================================================
    PRINT '(/2A)', 'Model: ', TRIM(model%name)
    PRINT '(/A/)', '    Tf     Sigma    W      M0    tauc   taul   tauh     Ref_Tols'
    PRINT '(I8,F8.4,F6.2,4F7.3,2F8.4)', Tf, Sigma, W, M0, 1/wc, taul, tauh, &
        model%refine_tol, model%deref_tol
    PRINT '(/A/)','        T   Nblks  nLeaf  Unused  Rho_min    Rho_max   U_max      Mass     Levels'
    !==========================================================

    CALL SYSTEM_CLOCK( t1, Rate )
    tp0 = t1
    tpLast = t1
    cellUpdates = 0

    !== One time step, starting from h, g and C (interior) after streaming:
    !   1) ghost layers of C, as many as the model's collision reads; they also
    !      serve the refinement criterion and the moments of g
    !   2) mesh adaptation (from the second step on: the initialization adapted
    !      the mesh); if it changed, the ghost layers are filled again
    !   3) output of the current state and the rise of the bubble (pressure and
    !      velocity computed first)
    !   4) one parallel region: moments of g and collision in one block loop,
    !      ghost cells of h and g, streaming together with C = sum(h).
    !      The work of each stage is shared with orphaned !$OMP DO constructs
    !      (also inside the Fill_Ghost_Cells_* routines); their implicit
    !      barriers separate the stages. Timers are updated by the primary
    !      thread right after each barrier.
    !   The moments of g are not recomputed before the first collision (the
    !   initial P, U are used) nor twice on output steps.
    DO T = 0, Tf

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED( model )
        CALL Fill_Ghost_Cells_C( model%ghost_C )
!$OMP END PARALLEL
        CALL Lap( t_Ghost_C )

        IF( T > 0 ) THEN
            CALL Refine_Derefine_Mesh( model, changed )
            IF( changed ) THEN
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED( model )
                CALL Fill_Ghost_Cells_C( model%ghost_C )
!$OMP END PARALLEL
            END IF
            CALL Lap( t_refine )
        END IF

        !== real-time view (optional, AMR_VIEW): C of this step, P and U of the last moments
        IF( Viewer_Due( T, Tf ) ) THEN
            CALL Viewer_Update( T, T*time_per_step )
            CALL Lap( t_Output )
        END IF

        moments_done  = ( T == 0 )
        is_output     = ( MOD(T,nOutStep) == 0 .OR. T == Tf )
        is_diagnostic = ( is_output .OR. MOD(T,nDiagStep) == 0 )
        IF( is_diagnostic ) THEN
            IF( .NOT. moments_done ) THEN
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED( LeafByLevel, nLeaf, model )
                DO m = 1, nLeaf
                    CALL model%Moments_g( LeafByLevel(m) )
                END DO
!$OMP END PARALLEL DO
                moments_done = .True.
            END IF
            IF( is_output ) THEN
                CALL Write_Output( T*time_per_step )
                CALL Print_Diagnostics( T*time_per_step )
            END IF
            CALL problem%Diagnostics( T*time_per_step )
            CALL Lap( t_Output )
            IF( T==Tf ) THEN
                CALL Close_Output
                IF( ASSOCIATED(problem%Finish) ) CALL problem%Finish
                EXIT
            END IF
        END IF

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED( Blk_h, Blk_g, Blk_C, Blk_C2, LeafByLevel, LevelStart, LevelMin, LevelMax, &
!$OMP   model, moments_done, &
!$OMP   t_collision, t_Ghost_hg, t_advection )

        !== pressure and velocity from g, then the collision, block by block.
        !   This loop and the streaming loop run level by level over LeafByLevel
        !   with the same static schedule as Fill_Ghost_Cells_hg, so each block
        !   is collided, ghost-filled and streamed by the same thread (OpenMP
        !   assigns iterations of static loops of equal length identically);
        !   NOWAIT: the levels are independent, one barrier follows the last.
        DO L = LevelMin, LevelMax
!$OMP DO SCHEDULE(STATIC)
        DO m = LevelStart(L), LevelStart(L+1)-1
            n = LeafByLevel(m)
            IF( moments_done ) THEN
                CALL model%Collide( n )
            ELSE
                CALL model%Moments_and_Collide( n )
            END IF
        END DO
!$OMP END DO NOWAIT
        END DO
!$OMP BARRIER
!$OMP MASKED
        CALL Lap( t_collision )
!$OMP END MASKED

        CALL Fill_Ghost_Cells_hg
!$OMP MASKED
        CALL Lap( t_Ghost_hg )
!$OMP END MASKED

        !== streaming, and the phase field C = sum(h) while h of the block is in cache:
        !   exact shift on the finest level, Lax-Wendroff with CFL number
        !   1/dx (a time step of 1) on the coarser ones
        DO L = LevelMin, LevelMax
        cfl = inv_dx_level( L )
!$OMP DO SCHEDULE(STATIC)
        DO m = LevelStart(L), LevelStart(L+1)-1
            n = LeafByLevel(m)
            IF( L == maxLevels ) THEN
                CALL Perfect_Shift_AMR( Blk_h(:,:,:,n) )
                CALL Perfect_Shift_AMR( Blk_g(:,:,:,n) )
            ELSE
                CALL Streaming_LaxWendroff_AMR( cfl, Blk_h(:,:,:,n) )
                CALL Streaming_LaxWendroff_AMR( cfl, Blk_g(:,:,:,n) )
            END IF
            CALL Macroscopic_Values_h_AMR( Blk_h(:,:,:,n), Blk_C(:,:,n) )
            IF( model%n_scalars == 2 ) CALL Macroscopic_Values_h_AMR( Blk_g(:,:,:,n), Blk_C2(:,:,n) )
        END DO
!$OMP END DO NOWAIT
        END DO
!$OMP BARRIER
!$OMP MASKED
        CALL Lap( t_advection )
!$OMP END MASKED

!$OMP END PARALLEL

        cellUpdates = cellUpdates + INT(nLeaf,i64)*nx*ny
        IF( MOD(T+1,nProgStep) == 0 .AND. T+1 < Tf ) CALL Progress( T+1 )

    END DO

    PRINT '(/A)', '------------------------------'
    PRINT '(A,F12.3)',   ' t_Output:        ', DBLE(t_Output   )/Rate
    PRINT '(A,F12.3,A)', ' t_collision:     ', DBLE(t_collision)/Rate, '   (moments of g and collision)'
    PRINT '(A,F12.3,A)', ' t_advection:     ', DBLE(t_advection)/Rate, '   (streaming and C)'
    PRINT '(A,F12.3)',   ' t_Ghost_C :      ', DBLE(t_Ghost_C  )/Rate
    PRINT '(A,F12.3)',   ' t_Ghost_hg:      ', DBLE(t_Ghost_hg )/Rate
    PRINT '(A,F12.3)',   ' t_refine:        ', DBLE(t_refine   )/Rate
    PRINT '(A,F12.3,A,I0,A,F0.2,A)', ' t_nbr_table:     ', DBLE(t_nbr_table)/Rate, &
        '   (in t_refine; ', n_nbr_rebuilds, ' rebuilds, ', &
        1.e6_dp*DBLE(t_nbr_table)/Rate/MAX(n_nbr_rebuilds,1), ' us each)'

    CALL Viewer_Close

CONTAINS

    !== progress after step Tdone: the speed over the last interval in million
    !   lattice-cell updates per second (MLUPS; every leaf cell counts once per
    !   step whatever its level) and the time left estimated from the average
    !   time per step so far (it changes with the number of refined blocks)
    SUBROUTINE Progress( Tdone )
        INTEGER, INTENT(IN) :: Tdone
        INTEGER(i64) :: tnow
        REAL(dp)   :: elapsed, mlups

        CALL SYSTEM_CLOCK( tnow )
        elapsed = REAL( tnow - tp0, dp ) / Rate
        mlups   = REAL( cellUpdates, dp ) / MAX( REAL(tnow-tpLast,dp)/Rate, 1.e-9_dp ) / 1.e6_dp
        PRINT '(A,I9,A,I9,A,F5.1,A,F8.3,A,I9,A,F8.2,A,F7.2,4A)', ' step', Tdone, ' /', Tf, &
            ' (', 100._dp*Tdone/Tf, '%)  t* =', Tdone*time_per_step, '  leaf cells', nLeaf*nx*ny, &
            '  MLUPS', mlups, ' (', mlups/OMP_GET_MAX_THREADS(), ' per thread)  elapsed ', &
            Hms(elapsed), '  left ~', Hms( elapsed/Tdone*(Tf-Tdone) )
        FLUSH( OUTPUT_UNIT )
        tpLast = tnow
        cellUpdates = 0
    END SUBROUTINE Progress

    !== seconds as h:mm:ss
    FUNCTION Hms( sec ) RESULT( str )
        REAL(dp), INTENT(IN) :: sec
        CHARACTER(len=10) :: str
        INTEGER :: isec
        isec = NINT( sec )
        WRITE( str, '(I4,":",I2.2,":",I2.2)' ) isec/3600, MOD(isec/60,60), MOD(isec,60)
    END FUNCTION Hms

    !== add the time elapsed since the last call to the accumulator tacc
    SUBROUTINE Lap( tacc )
        INTEGER(i64), INTENT(INOUT) :: tacc
        INTEGER(i64) :: tnow
        CALL SYSTEM_CLOCK( tnow )
        tacc = tacc + (tnow - t1)
        t1 = tnow
    END SUBROUTINE Lap

END SUBROUTINE Main_Loop

!******************************************************************************

!== one line of diagnostics: density range, maximum vertical velocity, total mass
SUBROUTINE Print_Diagnostics( Time )
    USE Config_Mod, ONLY: nx, ny, dx_level
    USE Blocks_Mod, ONLY: Blk, Nblks, LevelMin, LevelMax, nLeaf, LeafByLevel, nUnused, Blk_C, Blk_Uy
    USE Physics_Mod, ONLY: Rhol, Rhoh
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: Time    ! time of the current step (the case's unit)

    INTEGER :: m, n, dx, i, j
    REAL(wp) :: Rho, Rho_min, Rho_max, U_max
    REAL(dp) :: Mass    ! accumulated in double precision independent of wp

    Rho_min =  10
    Rho_max = -10
    Mass  = 0
    U_max = 0
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED( Blk, Blk_C, Blk_Uy, LeafByLevel, nLeaf, Rhol, Rhoh ) &
!$OMP   REDUCTION(min:Rho_min) REDUCTION(max:Rho_max,U_max) REDUCTION(+:Mass)
    DO m = 1, nLeaf
        n = LeafByLevel(m)

        dx = dx_level( Blk(n)%Level )
        DO j = 1, ny
        DO i = 1, nx
            Rho = Rhol + Blk_C(i,j, n) * (Rhoh - Rhol)
            Rho_min = MIN( Rho_min, Rho )
            Rho_max = MAX( Rho_max, Rho )
            Mass    = Mass + REAL(Rho,dp) * dx**2
        END DO
        END DO

        U_max = MAX( U_max, MAXVAL(ABS(Blk_Uy(:,:,n))) )
    END DO
!$OMP END PARALLEL DO
    WRITE(*,10) Time, Nblks, nLeaf, nUnused(), Rho_min, Rho_max, U_max, Mass, LevelMin, LevelMax
10  FORMAT(F10.2,2I7,I6,2F11.5,F9.5,F12.1,2I4)

END SUBROUTINE Print_Diagnostics

END MODULE Solver_Mod
