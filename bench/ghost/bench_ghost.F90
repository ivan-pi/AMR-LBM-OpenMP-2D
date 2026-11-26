!================================================
!   Benchmark of the ghost-cell fill of h and g (src/Ghost_hg_Mod.f90) and
!   its variants (ghost_template.inc) on the solver's own mesh: the rising
!   bubble of the build (-DAMR_L0, -DAMR_MIN_LEVEL), refined level by level
!   as the solver does before its first step.
!
!     bench_ghost [sweeps [trials]]        BENCH_ONLY=name1,name2
!
!   A sample is `sweeps` fills of all leaf blocks. Before every fill all
!   blocks are collided (not timed), so that the fill finds the caches as in
!   a time step. MLUPS = leaf cells per second of fill; the last column counts
!   the collision too (variants fused_*: every block pushes its edges into the
!   ghost cells of its same-level neighbours right after its collision, and
!   the fill does the rest).
!
!   Check: from the same saved state, fill, then stream every leaf block
!   (shift on the finest level, Lax-Wendroff below); err is the largest
!   difference of the streamed interior of h and g against the variant
!   `ref`, the solver's module before the study (src: the module in src/ now,
!   a fused variant: Push_Same_Level after the collision, then the fill).
!   Ghost values that no streaming reads may differ between variants.
!================================================
PROGRAM Bench_Ghost
    USE OMP_LIB, ONLY: OMP_GET_MAX_THREADS
    USE Kinds_Mod, ONLY: wp, dp, i64
    USE Config_Mod, ONLY: nx, ny, minLevels, maxLevels, inv_dx_level
    USE Blocks_Mod, ONLY: Blk, Blk_h, Blk_g, Nblks, nLeaf, LeafByLevel, LevelStart, LevelMin, LevelMax
    USE Mesh_Mod, ONLY: Create_Domain
    USE Refine_Mod, ONLY: Refine_Derefine_Mesh
    USE Model_Mod, ONLY: LBM_Model, Block_Kernel
    USE Case_Mod, ONLY: AMR_Case
    USE Model_JCP2016_Mod, ONLY: JCP2016_Model
    USE Rising_Bubble_Mod, ONLY: Rising_Bubble_Case
    USE Ghost_C_Mod, ONLY: Fill_Ghost_Cells_C
    USE LBM_Common_Mod, ONLY: Perfect_Shift_AMR, Streaming_LaxWendroff_AMR
#include "ghost_use.inc"
    IMPLICIT NONE

    ABSTRACT INTERFACE
        SUBROUTINE Fill_Kernel
        END SUBROUTINE Fill_Kernel
    END INTERFACE
    TYPE :: Variant_t
        CHARACTER(len=20) :: name = ''
        LOGICAL :: active = .TRUE.
        PROCEDURE(Fill_Kernel), POINTER, NOPASS :: run => NULL()
        !== fused variants: called for every block right after its collision
        PROCEDURE(Block_Kernel), POINTER, NOPASS :: post => NULL()
        REAL(dp), ALLOCATABLE :: total(:)       ! MLUPS of collision + fill
        REAL(dp) :: err = 0
        REAL(dp), ALLOCATABLE :: mlups(:)
    END TYPE Variant_t
#include "ghost_table.inc"

    TYPE(LBM_Model) :: model
    TYPE(AMR_Case)  :: problem
    REAL(wp), ALLOCATABLE :: h0(:,:,:,:), g0(:,:,:,:), href(:,:,:,:), gref(:,:,:,:)
    REAL(wp) :: t_end, time_per_step
    INTEGER  :: sweeps, trials, v, t, s, L, m, n, length, stat, cnt(3), k, nb
    INTEGER(i64) :: t0, t1, t2, tsum, tcol, Rate
    PROCEDURE(Block_Kernel), POINTER :: nopost => NULL()
    CHARACTER(len=512) :: arg, only

    sweeps = 20
    trials = 11
    IF( COMMAND_ARGUMENT_COUNT() >= 1 ) THEN
        CALL GET_COMMAND_ARGUMENT( 1, arg ); READ( arg, * ) sweeps
    END IF
    IF( COMMAND_ARGUMENT_COUNT() >= 2 ) THEN
        CALL GET_COMMAND_ARGUMENT( 2, arg ); READ( arg, * ) trials
    END IF
#include "ghost_set.inc"
    CALL GET_ENVIRONMENT_VARIABLE( 'BENCH_ONLY', only, length, stat )
    IF( stat == 0 .AND. length > 0 ) THEN
        DO v = 1, nVar
            var(v)%active = INDEX( ','//TRIM(only)//',', ','//TRIM(var(v)%name)//',' ) > 0 .OR. var(v)%name == 'ref'
        END DO
    END IF

    !== the solver's initialization (Run_AMR_LBM)
    model   = JCP2016_Model()
    problem = Rising_Bubble_Case()
    CALL problem%Setup( t_end, time_per_step )
    CALL Create_Domain
    CALL Initialize_State
    DO L = minLevels, maxLevels-1
        CALL Refine_Derefine_Mesh( model )
        CALL Initialize_State
    END DO
!$OMP PARALLEL
    CALL Fill_Ghost_Cells_C( model%ghost_C )
!$OMP END PARALLEL
    CALL Collide_All( nopost )      ! off equilibrium

    !== what the fill has to do on this mesh: same-level copies, coarse/fine fills, walls
    cnt = 0
    DO m = 1, nLeaf
        n = LeafByLevel(m)
        DO k = 1, 8
            nb = Blk(n)%nbr(k)
            IF( nb < 0 ) THEN
                cnt(3) = cnt(3) + 1
            ELSEIF( nb > 0 ) THEN
                IF( Blk(nb)%IsLeaf .AND. Blk(nb)%Level == Blk(n)%Level ) THEN
                    cnt(1) = cnt(1) + 1
                ELSE
                    cnt(2) = cnt(2) + 1
                END IF
            END IF
        END DO
    END DO
    PRINT '(/A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)', 'bench_ghost: block ', nx, 'x', ny, ', ', nLeaf, ' leaf blocks on levels ', &
        LevelMin, '-', LevelMax, ', threads ', OMP_GET_MAX_THREADS()
    PRINT '(A,10I6)', ' leaves per level:', ( LevelStart(L+1)-LevelStart(L), L = LevelMin, LevelMax )
    PRINT '(A,I0,A,I0,A,I0)', ' neighbour directions: same-level copies ', cnt(1), ', coarse/fine fills ', cnt(2), ', wall ', cnt(3)

    ALLOCATE( h0(0:nx+1,0:ny+1,9,Nblks), g0(0:nx+1,0:ny+1,9,Nblks), href(nx,ny,9,Nblks), gref(nx,ny,9,Nblks) )
    h0 = Blk_h(:,:,:,1:Nblks)
    g0 = Blk_g(:,:,:,1:Nblks)

    DO v = 1, nVar
        IF( var(v)%name /= 'ref' ) CYCLE
        CALL Fill_and_Stream( var(v)%run, nopost )
        href = Blk_h(1:nx,1:ny,:,1:Nblks)
        gref = Blk_g(1:nx,1:ny,:,1:Nblks)
    END DO
    DO v = 1, nVar
        IF( .NOT. var(v)%active ) CYCLE
        ALLOCATE( var(v)%mlups(trials), var(v)%total(trials) )
        CALL Fill_and_Stream( var(v)%run, var(v)%post )
        var(v)%err = MAX( MAXVAL( ABS( Blk_h(1:nx,1:ny,:,1:Nblks) - href ) )/MAXVAL( ABS(href) ), &
                          MAXVAL( ABS( Blk_g(1:nx,1:ny,:,1:Nblks) - gref ) )/MAXVAL( ABS(gref) ) )
    END DO

    CALL SYSTEM_CLOCK( t1, Rate )
    DO t = 1, trials
        DO v = 1, nVar
            IF( .NOT. var(v)%active ) CYCLE
            Blk_h(:,:,:,1:Nblks) = h0
            Blk_g(:,:,:,1:Nblks) = g0
            tsum = 0
            tcol = 0
            DO s = 1, sweeps
                CALL SYSTEM_CLOCK( t0 )
                CALL Collide_All( var(v)%post )
                CALL SYSTEM_CLOCK( t1 )
                tcol = tcol + (t1 - t0)
!$OMP PARALLEL
                CALL var(v)%run
!$OMP END PARALLEL
                CALL SYSTEM_CLOCK( t2 )
                tsum = tsum + (t2 - t1)
            END DO
            var(v)%mlups(t) = REAL(nLeaf,dp)*nx*ny*sweeps / ( REAL(tsum,dp)/Rate ) / 1.e6_dp
            var(v)%total(t) = REAL(nLeaf,dp)*nx*ny*sweeps / ( REAL(tsum+tcol,dp)/Rate ) / 1.e6_dp
        END DO
    END DO

    PRINT '(/A)', '# variant              fill MLUPS(best)  (median)  collision+fill(best)   err'
    DO v = 1, nVar
        IF( .NOT. var(v)%active ) CYCLE
        CALL Sort( var(v)%mlups )
        CALL Sort( var(v)%total )
        PRINT '(A,A22,3F14.2,ES12.2)', 'R ', var(v)%name, var(v)%mlups(trials), var(v)%mlups((trials+1)/2), &
            var(v)%total(trials), var(v)%err
    END DO

CONTAINS

    SUBROUTINE Initialize_State
        INTEGER :: m, n
!$OMP PARALLEL DO PRIVATE( n )
        DO m = 1, nLeaf
            n = LeafByLevel(m)
            CALL problem%Initial_State( n )
            CALL model%Initialize( n )
        END DO
!$OMP END PARALLEL DO
    END SUBROUTINE Initialize_State

    SUBROUTINE Collide_All( post )
        PROCEDURE(Block_Kernel), POINTER, INTENT(IN) :: post
        INTEGER :: L, m
!$OMP PARALLEL PRIVATE( L, m )
        DO L = LevelMin, LevelMax
!$OMP DO SCHEDULE(STATIC)
        DO m = LevelStart(L), LevelStart(L+1)-1
            CALL model%Collide( LeafByLevel(m) )
            IF( ASSOCIATED(post) ) CALL post( LeafByLevel(m) )
        END DO
!$OMP END DO NOWAIT
        END DO
!$OMP END PARALLEL
    END SUBROUTINE Collide_All

    !== from the saved state: fill, then stream all leaf blocks
    SUBROUTINE Fill_and_Stream( fill, post )
        PROCEDURE(Fill_Kernel) :: fill
        PROCEDURE(Block_Kernel), POINTER, INTENT(IN) :: post
        INTEGER :: L, m, n
        Blk_h(:,:,:,1:Nblks) = h0
        Blk_g(:,:,:,1:Nblks) = g0
        IF( ASSOCIATED(post) ) THEN
            DO m = 1, nLeaf
                CALL post( LeafByLevel(m) )
            END DO
        END IF
!$OMP PARALLEL PRIVATE( L, m, n )
        CALL fill
        DO L = LevelMin, LevelMax
!$OMP DO SCHEDULE(STATIC)
        DO m = LevelStart(L), LevelStart(L+1)-1
            n = LeafByLevel(m)
            IF( L == maxLevels ) THEN
                CALL Perfect_Shift_AMR( Blk_h(:,:,:,n) )
                CALL Perfect_Shift_AMR( Blk_g(:,:,:,n) )
            ELSE
                CALL Streaming_LaxWendroff_AMR( inv_dx_level(L), Blk_h(:,:,:,n) )
                CALL Streaming_LaxWendroff_AMR( inv_dx_level(L), Blk_g(:,:,:,n) )
            END IF
        END DO
!$OMP END DO NOWAIT
        END DO
!$OMP END PARALLEL
    END SUBROUTINE Fill_and_Stream

    SUBROUTINE Sort( a )
        REAL(dp), INTENT(INOUT) :: a(:)
        INTEGER :: i, j
        REAL(dp) :: x
        DO i = 2, SIZE(a)
            x = a(i)
            j = i - 1
            DO WHILE( j >= 1 )
                IF( a(j) <= x ) EXIT
                a(j+1) = a(j)
                j = j - 1
            END DO
            a(j+1) = x
        END DO
    END SUBROUTINE Sort

END PROGRAM Bench_Ghost
