!================================================
!   Real-time view of the mesh and a field (optional)
!
!   Built with -DAMR_WITH_VIEWER=ON, the solver can show the leaf blocks
!   and one field in an OpenGL window while it runs (C library amrviz in
!   viewer/, SDL2). The command-line option --view[=N] turns it on and sets
!   the update interval in time steps, e.g. ./build/amr_lbm --view=20.
!   Without it, or in a build without the viewer, Viewer_Due is always
!   .FALSE. and nothing else is called. How the window is shown is the
!   viewer library's business, configured by its environment variables
!   (AMRVIZ_THREAD, see viewer/amrviz.h).
!
!   Each update copies C, P, Ux, Uy and the speed |U| of every leaf block
!   (as float) into the viewer's buffer; drawing runs on the viewer's own
!   thread and never holds up the solver, except while the user has paused
!   the run with the space bar. On macOS the window has to live on the main
!   thread, the solver's: each update then also handles the window's input
!   and draws, so the window reacts only at updates (keep N small).
!   Pressure and velocity are those of the last moments of g (one step
!   old, except on output and diagnostic steps).
!   Closing the window lets the run continue without the viewer; at the end
!   of the run the window stays open until it is closed.
!================================================
MODULE Viewer_Mod
    USE Kinds_Mod, ONLY: wp
#ifdef USE_AMRVIZ
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_LONG, C_DOUBLE, C_NULL_CHAR
    USE amrviz, ONLY: amrviz_open, amrviz_set_field_name, amrviz_is_open, amrviz_begin_snapshot, &
                      amrviz_add_block, amrviz_block_field, amrviz_end_snapshot, amrviz_wait_if_paused, &
                      amrviz_close
#endif
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Viewer_Open, Viewer_Due, Viewer_Update, Viewer_Close

    !== update interval in time steps; 0: viewer off
    INTEGER :: nViewStep = 0

#ifdef USE_AMRVIZ
    INTEGER, PARAMETER :: nFields = 5
    CHARACTER(len=*), PARAMETER :: FieldNames(nFields) = [character(len=5) :: 'C', 'P', 'Ux', 'Uy', '|U|']
#endif

CONTAINS

!******************************************************************************

!== open the window if the update interval nStep (time steps) is positive
SUBROUTINE Viewer_Open( nStep )
#ifdef USE_AMRVIZ
    USE Config_Mod, ONLY: nx, ny, L0, nBaseX, nBaseY
#endif
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nStep
#ifdef USE_AMRVIZ
    INTEGER :: f, win_h, win_w
#endif

    nViewStep = MAX( nStep, 0 )
    IF( nViewStep == 0 ) RETURN

#ifdef USE_AMRVIZ
    !== a window of the domain's aspect ratio, 900 pixels high
    win_h = 900
    win_w = MAX( 480, NINT( REAL(win_h) * nBaseX / nBaseY ) )
    IF( amrviz_open( win_w, win_h, nx, ny, REAL(nBaseX*L0, C_DOUBLE), REAL(nBaseY*L0, C_DOUBLE), &
                     nFields ) /= 0 ) THEN
        PRINT '(/A)', 'Viewer: no window (see the message above); running without it'
        nViewStep = 0
        RETURN
    END IF
    DO f = 1, nFields
        CALL amrviz_set_field_name( f-1, TRIM(FieldNames(f)) // C_NULL_CHAR )
    END DO
    PRINT '(/A,I0,A)', 'Viewer: update every ', nViewStep, ' steps (H in the window: controls)'
#else
    PRINT '(/A)', 'Viewer: --view ignored, built without the viewer (configure with -DAMR_WITH_VIEWER=ON)'
    nViewStep = 0
#endif

END SUBROUTINE Viewer_Open

!******************************************************************************

!== .TRUE. on the steps that update the viewer: every nViewStep steps and the last
LOGICAL FUNCTION Viewer_Due( T, Tf )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: T, Tf

    Viewer_Due = .FALSE.
    IF( nViewStep > 0 ) Viewer_Due = ( MOD(T,nViewStep) == 0 .OR. T == Tf )

END FUNCTION Viewer_Due

!******************************************************************************

!== hand the leaf blocks and their fields to the viewer; waits while the run is paused
SUBROUTINE Viewer_Update( T, Time )
#ifdef USE_AMRVIZ
    USE Config_Mod, ONLY: nx, ny, dx_level
    USE Blocks_Mod, ONLY: Blk, nLeaf, LeafByLevel, Blk_C, Blk_P, Blk_Ux, Blk_Uy
#endif
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: T       ! time step
    REAL(wp), INTENT(IN) :: Time    ! time in the case's unit
#ifdef USE_AMRVIZ
    INTEGER(C_INT), PARAMETER :: ldC = nx+4, offC = 2 + 2*(nx+4)   ! interior of C(-1:nx+2,-1:ny+2)
    INTEGER  :: m, b, dx
    INTEGER(C_INT) :: slot
    REAL(wp) :: speed(nx,ny)

    IF( amrviz_is_open() == 0 ) THEN    ! closed by the user
        PRINT '(A,I0,A)', ' Viewer closed at step ', T, '; the run continues'
        CALL amrviz_close( 0_C_INT )
        nViewStep = 0
        RETURN
    END IF

    CALL amrviz_begin_snapshot( REAL(Time, C_DOUBLE), INT(T, C_LONG) )
    DO m = 1, nLeaf
        b  = LeafByLevel(m)
        dx = dx_level( Blk(b)%Level )
        slot = amrviz_add_block( Blk(b)%Level, REAL(Blk(b)%I*dx*nx, C_DOUBLE), &
                                 REAL(Blk(b)%J*dx*ny, C_DOUBLE), REAL(dx, C_DOUBLE) )
        IF( slot < 0 ) EXIT
        speed = SQRT( Blk_Ux(:,:,b)**2 + Blk_Uy(:,:,b)**2 )
        CALL amrviz_block_field( slot, 0_C_INT, Blk_C(:,:,b),  ldC, offC )
        CALL amrviz_block_field( slot, 1_C_INT, Blk_P(:,:,b),  nx, 0_C_INT )
        CALL amrviz_block_field( slot, 2_C_INT, Blk_Ux(:,:,b), nx, 0_C_INT )
        CALL amrviz_block_field( slot, 3_C_INT, Blk_Uy(:,:,b), nx, 0_C_INT )
        CALL amrviz_block_field( slot, 4_C_INT, speed,         nx, 0_C_INT )
    END DO
    CALL amrviz_end_snapshot
    CALL amrviz_wait_if_paused
#endif

END SUBROUTINE Viewer_Update

!******************************************************************************

!== end of the run: keep the window open until the user closes it
SUBROUTINE Viewer_Close
    IMPLICIT NONE

#ifdef USE_AMRVIZ
    IF( amrviz_is_open() /= 0 ) THEN
        PRINT '(/A)', 'Viewer: run finished; close the window to exit'
        CALL amrviz_close( 1_C_INT )
    ELSE
        CALL amrviz_close( 0_C_INT )
    END IF
#endif
    nViewStep = 0

END SUBROUTINE Viewer_Close

END MODULE Viewer_Mod
