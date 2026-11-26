!================================================
!   Benchmark of the ghost-cell fill of C (src/Ghost_C_Mod.f90) and its
!   variants (ghost_c_template.inc) on the solver's own mesh, set up as in
!   bench_ghost, with the ghost layers of the model (JCP 2016: two).
!
!     bench_ghost_c [sweeps [trials]]        BENCH_ONLY=name1,name2
!
!   Before every fill C = sum(h) is computed for all blocks, as the streaming
!   loop of a time step does (timed apart; the last column counts it too).
!   Variants with a `post` kernel call it for every block right after its
!   C = sum(h) (fused_push); a `pre` kernel runs for all blocks in the timed
!   fill (push). MLUPS = leaf cells per second.
!
!   Check: C of all leaf blocks with all ghost layers after the fill against
!   the variant `ref`, the solver's module before the study. The saved state
!   has invalid ghost layers, so a variant must fill all that `ref` fills.
!================================================
PROGRAM Bench_Ghost_C
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
    USE Blocks_Mod, ONLY: Blk_C
    USE LBM_Common_Mod, ONLY: Macroscopic_Values_h_AMR
#include "ghost_c_use.inc"
    IMPLICIT NONE

    ABSTRACT INTERFACE
        SUBROUTINE Fill_Kernel( ng )
            INTEGER, INTENT(IN) :: ng
        END SUBROUTINE Fill_Kernel
    END INTERFACE
    TYPE :: Variant_t
        CHARACTER(len=20) :: name = ''
        LOGICAL :: active = .TRUE.
        PROCEDURE(Fill_Kernel), POINTER, NOPASS :: run => NULL()
        !== called for every block right after its C = sum(h) / for all blocks in the timed fill
        PROCEDURE(Block_Kernel), POINTER, NOPASS :: post => NULL(), pre => NULL()
        REAL(dp), ALLOCATABLE :: total(:)       ! MLUPS of collision + fill
        REAL(dp) :: err = 0
        REAL(dp), ALLOCATABLE :: mlups(:)
    END TYPE Variant_t
#include "ghost_c_table.inc"

    TYPE(LBM_Model) :: model
    TYPE(AMR_Case)  :: problem
    REAL(wp), ALLOCATABLE :: C0(:,:,:), Cref(:,:,:)
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
#include "ghost_c_set.inc"
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
    PRINT '(/A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)', 'bench_ghost_c: block ', nx, 'x', ny, ', ', nLeaf, ' leaf blocks on levels ', &
        LevelMin, '-', LevelMax, ', threads ', OMP_GET_MAX_THREADS()
    PRINT '(A,10I6)', ' leaves per level:', ( LevelStart(L+1)-LevelStart(L), L = LevelMin, LevelMax )
    PRINT '(A,I0,A,I0,A,I0)', ' neighbour directions: same-level copies ', cnt(1), ', coarse/fine fills ', cnt(2), ', wall ', cnt(3)

    ALLOCATE( C0(-1:nx+2,-1:ny+2,Nblks), Cref(-1:nx+2,-1:ny+2,Nblks) )
    !== the saved state has invalid ghost layers: every fill has to produce them all
    CALL Sum_h_All( nopost )
    C0 = -7
    C0(1:nx,1:ny,:) = Blk_C(1:nx,1:ny,1:Nblks)

    DO v = 1, nVar
        IF( var(v)%name /= 'ref' ) CYCLE
        CALL Fill_Once( var(v) )
        Cref = Blk_C(:,:,1:Nblks)
    END DO
    DO v = 1, nVar
        IF( .NOT. var(v)%active ) CYCLE
        ALLOCATE( var(v)%mlups(trials), var(v)%total(trials) )
        CALL Fill_Once( var(v) )
        var(v)%err = 0
        DO m = 1, nLeaf
            n = LeafByLevel(m)
            var(v)%err = MAX( var(v)%err, MAXVAL( ABS( Blk_C(:,:,n) - Cref(:,:,n) ) ) )
        END DO
    END DO

    CALL SYSTEM_CLOCK( t1, Rate )
    DO t = 1, trials
        DO v = 1, nVar
            IF( .NOT. var(v)%active ) CYCLE
            Blk_C(:,:,1:Nblks) = C0
            tsum = 0
            tcol = 0
            DO s = 1, sweeps
                CALL SYSTEM_CLOCK( t0 )
                CALL Sum_h_All( var(v)%post )
                CALL SYSTEM_CLOCK( t1 )
                tcol = tcol + (t1 - t0)
                CALL Timed_Fill( var(v) )
                CALL SYSTEM_CLOCK( t2 )
                tsum = tsum + (t2 - t1)
            END DO
            var(v)%mlups(t) = REAL(nLeaf,dp)*nx*ny*sweeps / ( REAL(tsum,dp)/Rate ) / 1.e6_dp
            var(v)%total(t) = REAL(nLeaf,dp)*nx*ny*sweeps / ( REAL(tsum+tcol,dp)/Rate ) / 1.e6_dp
        END DO
    END DO

    PRINT '(/A)', '# variant              fill MLUPS(best)  (median)  C=sum(h)+fill(best)   err'
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

    !== C = sum(h) of all leaf blocks, each followed by the variant's post kernel
    SUBROUTINE Sum_h_All( post )
        PROCEDURE(Block_Kernel), POINTER, INTENT(IN) :: post
        INTEGER :: L, m, n
!$OMP PARALLEL PRIVATE( L, m, n )
        DO L = LevelMin, LevelMax
!$OMP DO SCHEDULE(STATIC)
        DO m = LevelStart(L), LevelStart(L+1)-1
            n = LeafByLevel(m)
            CALL Macroscopic_Values_h_AMR( Blk_h(:,:,:,n), Blk_C(:,:,n) )
            IF( ASSOCIATED(post) ) CALL post( n )
        END DO
!$OMP END DO NOWAIT
        END DO
!$OMP END PARALLEL
    END SUBROUTINE Sum_h_All

    SUBROUTINE Timed_Fill( vr )
        TYPE(Variant_t), INTENT(IN) :: vr
        INTEGER :: L, m
!$OMP PARALLEL PRIVATE( L, m )
        IF( ASSOCIATED(vr%pre) ) THEN
            DO L = LevelMin, LevelMax
!$OMP DO SCHEDULE(STATIC)
            DO m = LevelStart(L), LevelStart(L+1)-1
                CALL vr%pre( LeafByLevel(m) )
            END DO
!$OMP END DO NOWAIT
            END DO
!$OMP BARRIER
        END IF
        CALL vr%run( model%ghost_C )
!$OMP END PARALLEL
    END SUBROUTINE Timed_Fill

    !== from the saved state
    SUBROUTINE Fill_Once( vr )
        TYPE(Variant_t), INTENT(IN) :: vr
        Blk_C(:,:,1:Nblks) = C0
        CALL Sum_h_All( vr%post )
        CALL Timed_Fill( vr )
    END SUBROUTINE Fill_Once

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

END PROGRAM Bench_Ghost_C
