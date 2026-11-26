!================================================
!   Template of the ghost-cell fill of C for the benchmark bench_ghost_c:
!   src/Ghost_C_Mod.f90 (as of the study) with options. The including file
!   defines the module name GC_MOD and optionally
!     GC_NG n     the number of ghost layers as a compile-time constant in the strip
!                 routines (the argument ng must equal it)
!     GC_FIXED    (with GC_NG) interpolations with the strip shapes fixed at compile
!                 time, as for h and g, instead of the routines with run-time sizes
!     GC_PUSH     (with GC_NG) no same-level copies in the fill: Push_C( n ) writes the
!                 edges of block n into its neighbours' ghost layers
!     GC_PROBE_NOCOPY, GC_PROBE_NOCF, GC_PROBE_NOWALL   skip the same-level copies,
!                 the coarse/fine fills, the walls of the first pass (wrong results)
!================================================
MODULE GC_MOD
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk, Blk_C, Blk_C2, nScalars, di, dj, facing
    USE BC_Mod, ONLY: Symmetric_C_3, Symmetric_C_4
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Fill_Ghost_Cells_C
#ifdef GC_PUSH
    PUBLIC :: Push_C
#endif
#ifdef GC_NG
#define NGV GC_NG
#else
#define NGV ng
#endif
#ifdef GC_FIXED
#define GC_COARSE_LR( a, b, c, d, e, f ) Coarse_LR( e, f )
#define GC_FINE_LR( a, b, c, d, e, f )   Fine_LR( e, f )
#define GC_COARSE_BT( a, b, c, d, e, f ) Coarse_BT( e, f )
#define GC_FINE_BT( a, b, c, d, e, f )   Fine_BT( e, f )
#define GC_COARSE_CO( a, b, c, d, e, f ) Coarse_CO( e, f )
#define GC_FINE_CO( a, b, c, d, e, f )   Fine_CO( e, f )
#else
#define GC_COARSE_LR Biquadratic_from_Coarse
#define GC_FINE_LR   Biquadratic_from_Fine
#define GC_COARSE_BT Biquadratic_from_Coarse
#define GC_FINE_BT   Biquadratic_from_Fine
#define GC_COARSE_CO Biquadratic_from_Coarse
#define GC_FINE_CO   Biquadratic_from_Fine
#endif

CONTAINS


!== fill ng ghost layers of C on all leaf blocks, in two passes over the
!   levels as for h and g (Fill_Ghost_Cells_hg): copies, restrictions and
!   walls first, then the interpolations from coarser neighbours. The wall
!   ghost cells in the end columns next to an interpolated side copy its
!   ghost cells and are redone then; the collision stencils read them.
SUBROUTINE Fill_Ghost_Cells_C( ng )
    USE Blocks_Mod, ONLY: LevelMin, LevelMax, LeafByLevel, LevelStart
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ng

    INTEGER :: L, k, m, n, X1

    DO L = LevelMin, LevelMax
!$OMP DO schedule(static)
        DO m = LevelStart(L), LevelStart(L+1)-1
            n = LeafByLevel(m)
            CALL Fill( n, L, .FALSE., Blk_C )
#ifndef GC_PROBE_NOWALL
            CALL Walls( n, 1-ng, nx+ng, Blk_C )
#endif
            IF( nScalars == 2 ) THEN
                CALL Fill( n, L, .FALSE., Blk_C2 )
                CALL Walls( n, 1-ng, nx+ng, Blk_C2 )
            END IF
        END DO
!$OMP END DO NOWAIT
    END DO
!$OMP BARRIER

    DO L = LevelMin+1, LevelMax
!$OMP DO schedule(static)
        DO m = LevelStart(L), LevelStart(L+1)-1
            n = LeafByLevel(m)
            IF( .NOT. Blk(n)%Coarser ) CYCLE
            CALL Fill( n, L, .TRUE., Blk_C )
            IF( nScalars == 2 ) CALL Fill( n, L, .TRUE., Blk_C2 )
            DO k = 1, 2
                IF( Blk(n)%nbr(k) <= 0 ) CYCLE
                IF( Blk(Blk(n)%nbr(k))%Level == L ) CYCLE
                X1 = MERGE( 1-ng, nx+1, k == 1 )
                CALL Walls( n, X1, X1+ng-1, Blk_C )
                IF( nScalars == 2 ) CALL Walls( n, X1, X1+ng-1, Blk_C2 )
            END DO
        END DO
!$OMP END DO NOWAIT
    END DO
!$OMP BARRIER

CONTAINS

    !== the sides and corners of block n on level L of the scalar pool C from
    !   same-level and refined neighbours (second = F) or from coarser ones
    !   (second = T)
    SUBROUTINE Fill( n, L, second, C )
        INTEGER, INTENT(IN) :: n, L
        LOGICAL, INTENT(IN) :: second
        REAL(wp), INTENT(INOUT) :: C(-1:nx+2,-1:ny+2,*)
        INTEGER :: k, nb
        DO k = 1, 8
            nb = Blk(n)%nbr(k)
            IF( nb <= 0 ) CYCLE
            IF( (Blk(nb)%Level < L) .NEQV. second ) CYCLE
#if defined(GC_PROBE_NOCOPY) || defined(GC_PUSH)
            IF( Blk(nb)%IsLeaf .AND. Blk(nb)%Level == L ) CYCLE
#endif
#ifdef GC_PROBE_NOCF
            IF( .NOT. ( Blk(nb)%IsLeaf .AND. Blk(nb)%Level == L ) ) CYCLE
#endif
            IF( k <= 2 ) THEN
                CALL Fill_Side_X( n, nb, di(k), ng, C )
            ELSEIF( k <= 4 ) THEN
                CALL Fill_Side_Y( n, nb, dj(k), ng, C )
            ELSE
                CALL Fill_Corner( n, k, ng, C )
            END IF
        END DO
    END SUBROUTINE Fill

    !== symmetry on the bottom and top walls, columns X1:X2, of the scalar pool C
    SUBROUTINE Walls( n, X1, X2, C )
        INTEGER, INTENT(IN) :: n, X1, X2
        REAL(wp), INTENT(INOUT) :: C(-1:nx+2,-1:ny+2,*)
        IF( Blk(n)%nbr(3) < 0 ) CALL Symmetric_C_3( ng, C(:,:,n), X1, X2 )
        IF( Blk(n)%nbr(4) < 0 ) CALL Symmetric_C_4( ng, C(:,:,n), X1, X2 )
    END SUBROUTINE Walls

END SUBROUTINE Fill_Ghost_Cells_C

!******************************************************************************

!== the ng ghost columns of block n towards its left (s = -1) or right (s = +1)
!   neighbour nb. Copies are element loops: as sections of the same pool
!   array they would go through a temporary (no proof that they don't overlap).
SUBROUTINE Fill_Side_X( n, nb, s, ng, Blk_C )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nb, s, ng
    REAL(wp), INTENT(INOUT) :: Blk_C(-1:nx+2,-1:ny+2,*)     ! the scalar pool (C or C2)

    INTEGER  :: ig, i0, c, l, j, Y
    REAL(wp) :: Cf(2+NGV,2*ny), Cc(3,ny/2+2), Cg(NGV,ny)

    ig = MERGE( 0, nx+1, s < 0 )

    IF( Blk(nb)%Level < Blk(n)%Level ) THEN
        i0 = MERGE( nx+1, 0, s < 0 )
        Y  = (ny/2) * MOD(Blk(n)%J,2)
        DO l = 0, 2
            Cc(l+1,:) = Blk_C(i0+s*l,Y:Y+ny/2+1, nb)
        END DO
        CALL GC_COARSE_LR( 3,ny/2+2, NGV,ny, Cc, Cg )
    ELSE
        i0 = MERGE( nx, 1, s < 0 )
        IF( Blk(nb)%IsLeaf ) THEN
            DO j = 1, ny
            DO l = 0, NGV-1
                Blk_C(ig+s*l,j, n) = Blk_C(i0+s*l,j, nb)
            END DO
            END DO
            RETURN
        END IF
        c = Blk(nb)%child + MERGE( 1, 0, s < 0 )      ! the facing children c, c+2
        DO l = 0, 1+NGV
            Cf(l+1,   1:  ny) = Blk_C(i0+s*l,1:ny, c  )
            Cf(l+1,ny+1:2*ny) = Blk_C(i0+s*l,1:ny, c+2)
        END DO
        CALL GC_FINE_LR( 2+NGV,2*ny, NGV,ny, Cf, Cg )
    END IF

    DO j = 1, ny
    DO l = 0, NGV-1
        Blk_C(ig+s*l,j, n) = Cg(l+1,j)
    END DO
    END DO

END SUBROUTINE Fill_Side_X

!******************************************************************************

!== the ng ghost rows of block n towards its bottom (s = -1) or top (s = +1)
!   neighbour nb
SUBROUTINE Fill_Side_Y( n, nb, s, ng, Blk_C )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nb, s, ng
    REAL(wp), INTENT(INOUT) :: Blk_C(-1:nx+2,-1:ny+2,*)     ! the scalar pool (C or C2)

    INTEGER  :: jg, j0, c, l, i, X
    REAL(wp) :: Cf(2*nx,2+NGV), Cc(nx/2+2,3), Cg(nx,NGV)

    jg = MERGE( 0, ny+1, s < 0 )

    IF( Blk(nb)%Level < Blk(n)%Level ) THEN
        j0 = MERGE( ny+1, 0, s < 0 )
        X  = (nx/2) * MOD(Blk(n)%I,2)
        DO l = 0, 2
            Cc(:,l+1) = Blk_C(X:X+nx/2+1,j0+s*l, nb)
        END DO
        CALL GC_COARSE_BT( nx/2+2,3, nx,NGV, Cc, Cg )
    ELSE
        j0 = MERGE( ny, 1, s < 0 )
        IF( Blk(nb)%IsLeaf ) THEN
            DO l = 0, NGV-1
            DO i = 1, nx
                Blk_C(i,jg+s*l, n) = Blk_C(i,j0+s*l, nb)
            END DO
            END DO
            RETURN
        END IF
        c = Blk(nb)%child + MERGE( 2, 0, s < 0 )      ! the facing children c, c+1
        DO l = 0, 1+NGV
            Cf(   1:  nx,l+1) = Blk_C(1:nx,j0+s*l, c  )
            Cf(nx+1:2*nx,l+1) = Blk_C(1:nx,j0+s*l, c+1)
        END DO
        CALL GC_FINE_BT( 2*nx,2+NGV, nx,NGV, Cf, Cg )
    END IF

    DO l = 0, NGV-1
    DO i = 1, nx
        Blk_C(i,jg+s*l, n) = Cg(i,l+1)
    END DO
    END DO

END SUBROUTINE Fill_Side_Y

!******************************************************************************

!== the ng x ng ghost corner of block n towards its neighbour in direction k = 5..8
SUBROUTINE Fill_Corner( n, k, ng, Blk_C )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, k, ng
    REAL(wp), INTENT(INOUT) :: Blk_C(-1:nx+2,-1:ny+2,*)     ! the scalar pool (C or C2)

    INTEGER  :: nb, sx, sy, ig, jg, i0, j0, px, py, c, l, m
    REAL(wp) :: Cf(2+NGV,2+NGV), Cc(3,3), Cg(NGV,NGV)

    nb = Blk(n)%nbr(k)
    sx = di(k)
    sy = dj(k)
    ig = MERGE( 0, nx+1, sx < 0 )
    jg = MERGE( 0, ny+1, sy < 0 )

    IF( Blk(nb)%Level < Blk(n)%Level ) THEN
        px = MOD(Blk(n)%I,2)
        py = MOD(Blk(n)%J,2)
        i0 = MERGE( nx+1-(nx/2)*px, (nx/2)*(1-px), sx < 0 )
        j0 = MERGE( ny+1-(ny/2)*py, (ny/2)*(1-py), sy < 0 )
        DO m = 0, 2
        DO l = 0, 2
            Cc(l+1,m+1) = Blk_C(i0+sx*l,j0+sy*m, nb)
        END DO
        END DO
        CALL GC_COARSE_CO( 3,3, NGV,NGV, Cc, Cg )
    ELSE
        i0 = MERGE( nx, 1, sx < 0 )
        j0 = MERGE( ny, 1, sy < 0 )
        c  = nb
        IF( .NOT. Blk(nb)%IsLeaf ) c = Blk(nb)%child + facing(1,k)     ! the facing child
        DO m = 0, MERGE( NGV-1, 1+NGV, c == nb )
        DO l = 0, MERGE( NGV-1, 1+NGV, c == nb )
            Cf(l+1,m+1) = Blk_C(i0+sx*l,j0+sy*m, c)
        END DO
        END DO
        IF( c == nb ) THEN
            Cg = Cf(1:NGV,1:NGV)
        ELSE
            CALL GC_FINE_CO( 2+NGV,2+NGV, NGV,NGV, Cf, Cg )
        END IF
    END IF

    DO m = 0, NGV-1
    DO l = 0, NGV-1
        Blk_C(ig+sx*l,jg+sy*m, n) = Cg(l+1,m+1)
    END DO
    END DO

END SUBROUTINE Fill_Corner

#ifndef GC_FIXED
#include "biquadratic_runtime.inc"
#endif
#ifdef GC_FIXED
!== interpolations with the strip shapes fixed at compile time (Biquadratic.inc, one component)
#define FINE_NAME    Fine_LR
#define COARSE_NAME  Coarse_LR
#define SIZE_ARGS
#define FINE_SIZES   INTEGER, PARAMETER :: i2 = 2+NGV, j2 = 2*ny, x2 = NGV, y2 = ny, nk = 1
#define COARSE_SIZES INTEGER, PARAMETER :: i2 = 3, j2 = ny/2+2, x2 = NGV, y2 = ny, nk = 1
#include "Biquadratic.inc"
#define FINE_NAME    Fine_BT
#define COARSE_NAME  Coarse_BT
#define SIZE_ARGS
#define FINE_SIZES   INTEGER, PARAMETER :: i2 = 2*nx, j2 = 2+NGV, x2 = nx, y2 = NGV, nk = 1
#define COARSE_SIZES INTEGER, PARAMETER :: i2 = nx/2+2, j2 = 3, x2 = nx, y2 = NGV, nk = 1
#include "Biquadratic.inc"
#define FINE_NAME    Fine_CO
#define COARSE_NAME  Coarse_CO
#define SIZE_ARGS
#define FINE_SIZES   INTEGER, PARAMETER :: i2 = 2+NGV, j2 = 2+NGV, x2 = NGV, y2 = NGV, nk = 1
#define COARSE_SIZES INTEGER, PARAMETER :: i2 = 3, j2 = 3, x2 = NGV, y2 = NGV, nk = 1
#include "Biquadratic.inc"
#endif

#ifdef GC_PUSH
!== the edges of C of leaf block n into the ghost layers of its same-level leaf
!   neighbours (to be called right after C = sum(h) of block n)
SUBROUTINE Push_C( n )
    INTEGER, INTENT(IN) :: n
    INTEGER, PARAMETER :: kopp(8) = [ 2, 1, 4, 3, 8, 7, 6, 5 ]
    INTEGER :: k, nb
    DO k = 1, 8
        nb = Blk(n)%nbr(k)
        IF( nb <= 0 ) CYCLE
        IF( .NOT. Blk(nb)%IsLeaf .OR. Blk(nb)%Level /= Blk(n)%Level ) CYCLE
        CALL Copy_C( kopp(k), Blk_C(:,:,nb), Blk_C(:,:,n) )
    END DO
END SUBROUTINE Push_C

!== ghost layers of Cn in direction k from the edge of its same-level neighbour Cb
SUBROUTINE Copy_C( k, Cn, Cb )
    INTEGER,  INTENT(IN)    :: k
    REAL(wp), INTENT(INOUT) :: Cn(-1:nx+2,-1:ny+2)
    REAL(wp), INTENT(IN)    :: Cb(-1:nx+2,-1:ny+2)
    INTEGER, PARAMETER :: ng = GC_NG
    INTEGER :: i, j, l, m
    SELECT CASE( k )
    CASE( 1 )
        DO j = 1, ny
        DO l = 0, ng-1
            Cn(-l,j) = Cb(nx-l,j)
        END DO
        END DO
    CASE( 2 )
        DO j = 1, ny
        DO l = 0, ng-1
            Cn(nx+1+l,j) = Cb(1+l,j)
        END DO
        END DO
    CASE( 3 )
        DO l = 0, ng-1
        DO i = 1, nx
            Cn(i,-l) = Cb(i,ny-l)
        END DO
        END DO
    CASE( 4 )
        DO l = 0, ng-1
        DO i = 1, nx
            Cn(i,ny+1+l) = Cb(i,1+l)
        END DO
        END DO
    CASE( 5 )
        DO m = 0, ng-1
        DO l = 0, ng-1
            Cn(-l,-m) = Cb(nx-l,ny-m)
        END DO
        END DO
    CASE( 6 )
        DO m = 0, ng-1
        DO l = 0, ng-1
            Cn(-l,ny+1+m) = Cb(nx-l,1+m)
        END DO
        END DO
    CASE( 7 )
        DO m = 0, ng-1
        DO l = 0, ng-1
            Cn(nx+1+l,-m) = Cb(1+l,ny-m)
        END DO
        END DO
    CASE( 8 )
        DO m = 0, ng-1
        DO l = 0, ng-1
            Cn(nx+1+l,ny+1+m) = Cb(1+l,1+m)
        END DO
        END DO
    END SELECT
END SUBROUTINE Copy_C
#endif

END MODULE GC_MOD

#undef GC_MOD
#undef GC_NG
#undef GC_FIXED
#undef GC_PUSH
#undef GC_PROBE_NOCOPY
#undef GC_PROBE_NOCF
#undef GC_PROBE_NOWALL
#undef NGV
#undef GC_COARSE_LR
#undef GC_FINE_LR
#undef GC_COARSE_BT
#undef GC_FINE_BT
#undef GC_COARSE_CO
#undef GC_FINE_CO
