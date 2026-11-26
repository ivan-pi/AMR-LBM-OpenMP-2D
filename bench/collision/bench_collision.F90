!================================================
!   Benchmark of the JCP 2016 collision kernel and its variants
!   (Collision_Variants_Mod) on a pool of blocks with a synthetic two-phase
!   state: a lattice of bubbles with a moving fluid, h and g off equilibrium.
!
!     bench_collision [pool_MiB [sweeps [trials]]]
!
!   Kinds of kernels: collide (chemical potential + collision), mom+col (mc_*:
!   with the moments of g, the path of a normal step), moments (mg_*: the
!   moments of g alone, steps with diagnostics).
!
!   pool_MiB : field data of all blocks (default 12: beyond L2, within L3)
!   sweeps   : collisions of all blocks per timed sample (default 20)
!   trials   : timed samples per variant, interleaved (default 15)
!   BENCH_ONLY=name1,name2 restricts the run to these variants.
!
!   Every variant starts each sample from the same saved state. Its result
!   after one collision is compared with the kernel of the model (interior
!   cells of h, g and, for the mc_ variants, P, Ux, Uy) as it was before this
!   study (Model_JCP2016_Ref_Mod; new, mc_new: the model in src/ now): err is the largest
!   difference divided by the largest magnitude of the field.
!================================================
PROGRAM Bench_Collision
    USE OMP_LIB, ONLY: OMP_GET_MAX_THREADS
    USE Kinds_Mod, ONLY: wp, dp, i64
    USE Config_Mod, ONLY: nx, ny, maxLevels
    USE D2Q9_Mod, ONLY: Wa
    USE Physics_Mod, ONLY: W, Sigma, Beta, kappa, taul, tauh, Rhol, Rhoh, Gy, D, M0, wc
    USE Blocks_Mod, ONLY: Blk, Blk_h, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy, Nblks, Init_Block_Pool
    USE Model_Mod, ONLY: Block_Kernel
    USE Model_JCP2016_Mod, ONLY: JCP2016_Model
    USE Model_JCP2016_Ref_Mod, ONLY: JCP2016_Ref_Model
    USE Collision_Variants_Mod
    IMPLICIT NONE

    TYPE :: Variant_t
        CHARACTER(len=20) :: name = ''
        INTEGER :: kind = 1                 ! 1 collision, 2 moments of g + collision, 3 moments of g
        LOGICAL :: al = .FALSE.             ! works on the padded, aligned copy of the fields
        LOGICAL :: active = .TRUE.
        PROCEDURE(Block_Kernel), POINTER, NOPASS :: run => NULL()
        REAL(dp) :: err = 0
        REAL(dp), ALLOCATABLE :: mlups(:)
    END TYPE Variant_t
#include "variant_table.inc"

    REAL(wp), ALLOCATABLE :: h0(:,:,:,:), g0(:,:,:,:), P0(:,:,:), Ux0(:,:,:), Uy0(:,:,:)
    REAL(wp), ALLOCATABLE :: href(:,:,:,:,:), gref(:,:,:,:,:), mref(:,:,:,:,:)
    REAL(dp) :: pool_MiB, bytes
    INTEGER  :: nblk, sweeps, trials, v, t, s, n, length, stat, ir
    INTEGER(i64) :: t1, t2, Rate
    CHARACTER(len=512) :: arg, only
    CHARACTER(len=7), PARAMETER :: kinds(3) = [ 'collide', 'mom+col', 'moments' ]

    pool_MiB = 12
    sweeps   = 20
    trials   = 15
    IF( COMMAND_ARGUMENT_COUNT() >= 1 ) THEN
        CALL GET_COMMAND_ARGUMENT( 1, arg ); READ( arg, * ) pool_MiB
    END IF
    IF( COMMAND_ARGUMENT_COUNT() >= 2 ) THEN
        CALL GET_COMMAND_ARGUMENT( 2, arg ); READ( arg, * ) sweeps
    END IF
    IF( COMMAND_ARGUMENT_COUNT() >= 3 ) THEN
        CALL GET_COMMAND_ARGUMENT( 3, arg ); READ( arg, * ) trials
    END IF
    !== the block size as a run-time value, opaque to the compiler
    arg = ''
    WRITE( arg, * ) nx, ny
    READ( arg, * ) nx_rt, ny_rt

#include "variant_set.inc"
    CALL GET_ENVIRONMENT_VARIABLE( 'BENCH_ONLY', only, length, stat )
    IF( stat == 0 .AND. length > 0 ) THEN
        DO v = 1, nVar
            var(v)%active = INDEX( ','//TRIM(only)//',', ','//TRIM(var(v)%name)//',' ) > 0 &
                            .OR. var(v)%name == 'model' .OR. var(v)%name == 'mc_model' .OR. var(v)%name == 'mg_model'
        END DO
    END IF

    CALL Setup_Physics
    CALL C512_Setup
    IF( collision_avx512_available() == 0 .OR. wp /= dp ) THEN
        DO v = 1, nVar
            IF( INDEX( var(v)%name, 'c512' ) > 0 ) var(v)%active = .FALSE.
        END DO
    END IF
    the_model = JCP2016_Ref_Model()
    new_model = JCP2016_Model()

    bytes = REAL( STORAGE_SIZE(1._wp)/8, dp ) * ( 2*9*(nx+2)*(ny+2) + (nx+4)*(ny+4) + 3*nx*ny )
    nblk  = MAX( NINT( pool_MiB*2._dp**20 / bytes ), 1 )
    CALL Init_Block_Pool
    IF( nblk > SIZE(Blk) ) ERROR STOP 'pool larger than MaxBlks: raise AMR_L0 of the build'
    Nblks = nblk
    CALL Initial_State

    ALLOCATE( h0(0:nx+1,0:ny+1,9,nblk), g0(0:nx+1,0:ny+1,9,nblk), P0(nx,ny,nblk), Ux0(nx,ny,nblk), Uy0(nx,ny,nblk) )
    ALLOCATE( href(nx,ny,9,nblk,3), gref(nx,ny,9,nblk,3), mref(nx,ny,3,nblk,3) )
!$OMP PARALLEL DO SCHEDULE(STATIC)
    DO n = 1, nblk
        h0(:,:,:,n) = Blk_h(:,:,:,n)
        g0(:,:,:,n) = Blk_g(:,:,:,n)
        P0(:,:,n)   = Blk_P(:,:,n)
        Ux0(:,:,n)  = Blk_Ux(:,:,n)
        Uy0(:,:,n)  = Blk_Uy(:,:,n)
    END DO
!$OMP END PARALLEL DO

    CALL Al_Setup( nblk )
    PRINT '(/A,I0,A,I0,A,I0,A,F0.1,A,I0,A,I0,A,I0)', 'bench_collision: block ', nx, 'x', ny, ', ', nblk, ' blocks (', &
        nblk*bytes/2._dp**20, ' MiB), threads ', OMP_GET_MAX_THREADS(), ', sweeps ', sweeps, ', trials ', trials

    !== reference results: the kernels of the model
    DO v = 1, nVar
        IF( var(v)%name /= 'model' .AND. var(v)%name /= 'mc_model' .AND. var(v)%name /= 'mg_model' ) CYCLE
        ir = var(v)%kind
        CALL Restore( .FALSE. )
        CALL Sweep( var(v)%run )
        DO n = 1, nblk
            href(:,:,:,n,ir) = Blk_h(1:nx,1:ny,:,n)
            gref(:,:,:,n,ir) = Blk_g(1:nx,1:ny,:,n)
            mref(:,:,1,n,ir) = Blk_P(:,:,n);  mref(:,:,2,n,ir) = Blk_Ux(:,:,n);  mref(:,:,3,n,ir) = Blk_Uy(:,:,n)
        END DO
    END DO
    DO v = 1, nVar
        IF( .NOT. var(v)%active ) CYCLE
        ALLOCATE( var(v)%mlups(trials) )
        CALL Restore( var(v)%al )
        CALL Sweep( var(v)%run )
        IF( var(v)%al ) THEN
            DO n = 1, nblk
                CALL Al_Get( n )
            END DO
        END IF
        var(v)%err = Difference( var(v)%kind )
    END DO

    CALL SYSTEM_CLOCK( t1, Rate )
    DO t = 1, trials
        DO v = 1, nVar
            IF( .NOT. var(v)%active ) CYCLE
            CALL Restore( var(v)%al )
            CALL SYSTEM_CLOCK( t1 )
            DO s = 1, sweeps
                CALL Sweep( var(v)%run )
            END DO
            CALL SYSTEM_CLOCK( t2 )
            var(v)%mlups(t) = REAL(nblk,dp)*nx*ny*sweeps / ( REAL(t2-t1,dp)/Rate ) / 1.e6_dp
        END DO
    END DO

    PRINT '(/A)', '# variant              kind      MLUPS(best)  MLUPS(median)   err'
    DO v = 1, nVar
        IF( .NOT. var(v)%active ) CYCLE
        CALL Sort( var(v)%mlups )
        PRINT '(A,A22,A8,2F14.2,ES12.2)', 'R ', var(v)%name, kinds(var(v)%kind), &
            var(v)%mlups(trials), var(v)%mlups((trials+1)/2), var(v)%err
    END DO

CONTAINS

    !== parameters of the default rising-bubble case at L0 = 256 (Rising_Bubble_Mod)
    SUBROUTINE Setup_Physics
        REAL(wp), PARAMETER :: Eo = 40, Re = 40, rho_ratio = 1000, mu_ratio = 100
        REAL(wp) :: mu_h
        D     = 256/5._wp
        Rhoh  = 1
        Rhol  = Rhoh/rho_ratio
        Gy    = 5.e-6_wp
        W     = 5
        M0    = 0.03_wp
        wc    = 1._wp/(0.5_wp + M0*3._wp)
        Sigma = Gy*(Rhoh-Rhol)*D**2/Eo
        mu_h  = SQRT( Gy*Rhoh*(Rhoh-Rhol)*D**3 )/Re
        tauh  = 3._wp*(mu_h/Rhoh)
        taul  = 3._wp*(mu_h/mu_ratio/Rhol)
        Beta  = 12._wp*Sigma/W
        kappa = 1.5_wp*Sigma*W
    END SUBROUTINE Setup_Physics

    !== bubbles of radius 20 on a square lattice of period 64 cells (about half
    !   of the cells within the interface), a vortex lattice, h and g of the
    !   model's equilibrium with a deterministic perturbation of 1 %
    SUBROUTINE Initial_State
        REAL(wp), PARAMETER :: Lp = 64, Rb = 20, twopi = 6.283185307179586_wp
        INTEGER  :: n, i, j, k, nbx
        REAL(wp) :: X, Y, Ri, r
        INTEGER(i64) :: seed

        nbx = CEILING( SQRT( REAL(nblk,wp) ) )
        DO n = 1, nblk
            Blk(n)%Level  = maxLevels
            Blk(n)%I      = MOD( n-1, nbx )
            Blk(n)%J      = (n-1)/nbx
            Blk(n)%IsLeaf = .TRUE.
            DO j = -1, ny+2
            DO i = -1, nx+2
                X  = Blk(n)%I*nx + i - 0.5_wp
                Y  = Blk(n)%J*ny + j - 0.5_wp
                Ri = SQRT( (MODULO(X,Lp) - Lp/2)**2 + (MODULO(Y,Lp) - Lp/2)**2 )
                Blk_C(i,j,n) = 0.5_wp + 0.5_wp*TANH( 2._wp*(Ri-Rb)/W )
            END DO
            END DO
            DO j = 1, ny
            DO i = 1, nx
                X = Blk(n)%I*nx + i - 0.5_wp
                Y = Blk(n)%J*ny + j - 0.5_wp
                Blk_Ux(i,j,n) =  0.02_wp*SIN( twopi*X/Lp )*COS( twopi*Y/Lp ) + 0.005_wp
                Blk_Uy(i,j,n) = -0.02_wp*COS( twopi*X/Lp )*SIN( twopi*Y/Lp ) + 0.010_wp
                Blk_P (i,j,n) = -Sigma/Rb*Blk_C(i,j,n) + 1.e-3_wp*COS( twopi*X/Lp )
            END DO
            END DO
            Blk_h(:,:,:,n) = 0
            Blk_g(:,:,:,n) = 0
            CALL the_model%Initialize( n )
            seed = 12345 + n
            DO k = 1, 9
            DO j = 0, ny+1
            DO i = 0, nx+1
                seed = MODULO( seed*1103515245_i64 + 12345_i64, 2147483648_i64 )
                r = REAL( seed, wp )/2147483648._wp - 0.5_wp
                Blk_h(i,j,k,n) = Blk_h(i,j,k,n)*( 1._wp + 0.02_wp*r )
                seed = MODULO( seed*1103515245_i64 + 12345_i64, 2147483648_i64 )
                r = REAL( seed, wp )/2147483648._wp - 0.5_wp
                Blk_g(i,j,k,n) = Blk_g(i,j,k,n)*( 1._wp + 0.02_wp*r ) + 2.e-4_wp*r*Wa(k)
            END DO
            END DO
            END DO
        END DO
    END SUBROUTINE Initial_State

    SUBROUTINE Restore( al )
        LOGICAL, INTENT(IN) :: al
        INTEGER :: n
!$OMP PARALLEL DO SCHEDULE(STATIC)
        DO n = 1, nblk
            IF( al ) CALL Al_Put( n, h0(:,:,:,n), g0(:,:,:,n), P0(:,:,n), Ux0(:,:,n), Uy0(:,:,n) )
            Blk_h(:,:,:,n) = h0(:,:,:,n)
            Blk_g(:,:,:,n) = g0(:,:,:,n)
            Blk_P(:,:,n)   = P0(:,:,n)
            Blk_Ux(:,:,n)  = Ux0(:,:,n)
            Blk_Uy(:,:,n)  = Uy0(:,:,n)
        END DO
!$OMP END PARALLEL DO
    END SUBROUTINE Restore

    SUBROUTINE Sweep( kernel )
        PROCEDURE(Block_Kernel) :: kernel
        INTEGER :: n
!$OMP PARALLEL DO SCHEDULE(STATIC)
        DO n = 1, nblk
            CALL kernel( n )
        END DO
!$OMP END PARALLEL DO
    END SUBROUTINE Sweep

    FUNCTION Difference( ir ) RESULT( err )
        INTEGER, INTENT(IN) :: ir
        REAL(dp) :: err
        INTEGER :: k
        err = MAX( MAXVAL( ABS( Blk_h(1:nx,1:ny,:,1:nblk) - href(:,:,:,:,ir) ) )/MAXVAL( ABS( href(:,:,:,:,ir) ) ), &
                   MAXVAL( ABS( Blk_g(1:nx,1:ny,:,1:nblk) - gref(:,:,:,:,ir) ) )/MAXVAL( ABS( gref(:,:,:,:,ir) ) ) )
        IF( ir >= 2 ) THEN
            err = MAX( err, MAXVAL( ABS( Blk_P (:,:,1:nblk) - mref(:,:,1,:,ir) ) )/MAXVAL( ABS( mref(:,:,1,:,ir) ) ), &
                            MAXVAL( ABS( Blk_Ux(:,:,1:nblk) - mref(:,:,2,:,ir) ) )/MAXVAL( ABS( mref(:,:,2,:,ir) ) ), &
                            MAXVAL( ABS( Blk_Uy(:,:,1:nblk) - mref(:,:,3,:,ir) ) )/MAXVAL( ABS( mref(:,:,3,:,ir) ) ) )
        END IF
        k = 0
    END FUNCTION Difference

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

END PROGRAM Bench_Collision
