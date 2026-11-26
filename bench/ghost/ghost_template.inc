!================================================
!   Template of the ghost-cell fill of h and g for the benchmark bench_ghost:
!   src/Ghost_hg_Mod.f90 (as of the study) with options. The including file
!   defines the module name GH_MOD and
!     GH_COPY  0  same-level copies as in the solver: loop over the directions,
!                 SELECT CASE( k ) inside, sections of the pool arrays
!              1  SELECT CASE( k ) outside the loops, the blocks as explicit-shape
!                 dummy arguments, the list of directions at run time
!              2  as 1 with one routine per list (all / shift / Lax-Wendroff):
!                 the directions are constants in every branch
!              3  as 2, columns copied cell by cell (directions in the inner loop)
!              4  as 2, pushed: a block writes its edges into the ghost cells of its
!                 neighbours (Push_Block( n ) can run right after the collision of n)
!     GH_CFDIRS   ghost cells interpolated from a coarser neighbour only for the
!                 directions the block's streaming reads (3 or 6 of 8 on a side, 1 or
!                 2 at a corner) unless finer blocks read them (Blk%Finer)
!     GH_PROBE_NOCOPY, GH_PROBE_NOCF, GH_PROBE_NOWALL   skip the same-level copies,
!                 the coarse/fine fills, the walls (wrong results: cost breakdown)
!   All macros are undefined at the end.
!================================================
MODULE GH_MOD
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: nx, ny, maxLevels     ! the fixed block size (strip shapes)
    USE Blocks_Mod, ONLY: Blk, Blk_h, Blk_g, di, dj, facing
    USE BC_Mod, ONLY: Bounce_Back_3, Bounce_Back_4
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Fill_Ghost_Cells_hg
#if GH_COPY == 4
    PUBLIC :: Push_Block
    INTEGER, PARAMETER :: kopp(8) = [ 2, 1, 4, 3, 8, 7, 6, 5 ]
#endif

    !== directions read from the ghost cells in direction k by the streaming:
    !   all moving ones, those of the exact shift, those of Lax-Wendroff
    INTEGER, PARAMETER :: all_dirs(8) = [ 2, 3, 4, 5, 6, 7, 8, 9 ]
    INTEGER, PARAMETER :: nShift(8) = [ 3, 3, 3, 3, 1, 1, 1, 1 ]
    INTEGER, PARAMETER :: shift_dirs(3,8) = RESHAPE( [ 2,6,9, 4,7,8, 3,6,7, 5,8,9, &
                                                       6,0,0, 9,0,0, 7,0,0, 8,0,0 ], [3,8] )
    INTEGER, PARAMETER :: nLW(8) = [ 6, 6, 6, 6, 2, 2, 2, 2 ]
    INTEGER, PARAMETER :: lw_dirs(6,8) = RESHAPE( [ 2,4,6,7,8,9, 2,4,6,7,8,9, 3,5,6,7,8,9, 3,5,6,7,8,9, &
                                                    6,8,0,0,0,0, 7,9,0,0,0,0, 7,9,0,0,0,0, 6,8,0,0,0,0 ], [6,8] )

CONTAINS


!== one ghost layer of h and g on all leaf blocks, in two passes over the
!   levels. The first does the copies from same-level leaves and the
!   restrictions of refined neighbours, which read interior cells only, and
!   the wall conditions. The second interpolates from coarser neighbours,
!   which also reads ghost cells of the coarser block: by the 2:1 balance
!   (diagonal neighbours included) those lie towards blocks of the
!   interpolating block's level or the coarser one, filled by the first
!   pass. So the levels need no barrier between them (NOWAIT), the passes
!   one. The wall conditions run in the first pass only: their values that
!   come from a side ghost cell still to be interpolated (end columns,
!   directions with ex*ey /= 0 pointing out of the block) are read neither
!   by the streaming nor by finer neighbours, which read the wall row on
!   their own side (restricted in the first pass). The static
!   schedules per level are those of the collision and streaming loops:
!   each block stays with one thread.
SUBROUTINE Fill_Ghost_Cells_hg
    USE Blocks_Mod, ONLY: LevelMin, LevelMax, LeafByLevel, LevelStart
    IMPLICIT NONE

    INTEGER :: L, k, m, n, nb

    DO L = LevelMin, LevelMax
!$OMP DO schedule(static)
        DO m = LevelStart(L), LevelStart(L+1)-1
            n = LeafByLevel(m)

            DO k = 1, 8
                nb = Blk(n)%nbr(k)
                IF( nb <= 0 ) CYCLE
                IF( Blk(nb)%Level < L ) CYCLE                 ! second pass
                IF( Blk(nb)%IsLeaf ) THEN
#ifdef GH_PROBE_NOCOPY
                    CYCLE
#endif
#if GH_COPY == 4
                    CALL Push( n, nb, k )
#elif GH_COPY >= 2
                    IF( Blk(n)%Finer ) THEN
                        CALL Copy_All  ( k, Blk_h(:,:,:,n), Blk_h(:,:,:,nb), Blk_g(:,:,:,n), Blk_g(:,:,:,nb) )
                    ELSEIF( L == maxLevels ) THEN
                        CALL Copy_Shift( k, Blk_h(:,:,:,n), Blk_h(:,:,:,nb), Blk_g(:,:,:,n), Blk_g(:,:,:,nb) )
                    ELSE
                        CALL Copy_LW   ( k, Blk_h(:,:,:,n), Blk_h(:,:,:,nb), Blk_g(:,:,:,n), Blk_g(:,:,:,nb) )
                    END IF
#else
                    IF( Blk(n)%Finer ) THEN
                        CALL Copy_Same_Level( n, nb, k, all_dirs )
                    ELSEIF( L == maxLevels ) THEN
                        CALL Copy_Same_Level( n, nb, k, shift_dirs(1:nShift(k),k) )
                    ELSE
                        CALL Copy_Same_Level( n, nb, k, lw_dirs(1:nLW(k),k) )
                    END IF
#endif
#if defined(GH_PROBE_NOCF) || defined(GH_PROBE_NOFINE)
                ELSE
                    CYCLE
#else
                ELSEIF( k <= 2 ) THEN
                    CALL Fill_Side_X( n, nb, di(k) )
                ELSEIF( k <= 4 ) THEN
                    CALL Fill_Side_Y( n, nb, dj(k) )
                ELSE
                    CALL Fill_Corner( n, k )
#endif
                END IF
            END DO

#ifndef GH_PROBE_NOWALL
            CALL Walls( n, 0, nx+1 )
#endif

        END DO
!$OMP END DO NOWAIT
    END DO
!$OMP BARRIER

    DO L = LevelMin+1, LevelMax
!$OMP DO schedule(static)
        DO m = LevelStart(L), LevelStart(L+1)-1
            n = LeafByLevel(m)
            IF( .NOT. Blk(n)%Coarser ) CYCLE
#ifdef GH_PROBE_NOCF
            CYCLE
#endif

            DO k = 1, 8
                nb = Blk(n)%nbr(k)
                IF( nb <= 0 ) CYCLE
                IF( Blk(nb)%Level == L ) CYCLE                ! first pass
#ifdef GH_PROBE_NOCOARSE
                CYCLE
#endif
#ifdef GH_CFDIRS
                !== only the directions the block's streaming reads, unless finer blocks
                !   interpolate from its ghost cells
                IF( .NOT. Blk(n)%Finer ) THEN
                    IF( L == maxLevels ) THEN
                        CALL Fill_From_Coarse( n, nb, k, nShift(k), shift_dirs(1:nShift(k),k) )
                    ELSE
                        CALL Fill_From_Coarse( n, nb, k, nLW(k), lw_dirs(1:nLW(k),k) )
                    END IF
                    CYCLE
                END IF
#endif
                IF( k <= 2 ) THEN
                    CALL Fill_Side_X( n, nb, di(k) )
                ELSEIF( k <= 4 ) THEN
                    CALL Fill_Side_Y( n, nb, dj(k) )
                ELSE
                    CALL Fill_Corner( n, k )
                END IF
            END DO

        END DO
!$OMP END DO NOWAIT
    END DO
!$OMP BARRIER

CONTAINS

    !== bounce-back on the bottom and top walls, columns X1:X2
    SUBROUTINE Walls( n, X1, X2 )
        INTEGER, INTENT(IN) :: n, X1, X2
        IF( Blk(n)%nbr(3) < 0 ) THEN
            CALL Bounce_Back_3( Blk_h(:,:,:,n), X1, X2 )
            CALL Bounce_Back_3( Blk_g(:,:,:,n), X1, X2 )
        END IF
        IF( Blk(n)%nbr(4) < 0 ) THEN
            CALL Bounce_Back_4( Blk_h(:,:,:,n), X1, X2 )
            CALL Bounce_Back_4( Blk_g(:,:,:,n), X1, X2 )
        END IF
    END SUBROUTINE Walls

END SUBROUTINE Fill_Ghost_Cells_hg

!******************************************************************************

#if GH_COPY == 0
!== the directions dirs of the ghost cells of block n from its same-level leaf
!   neighbour nb in direction k: written out with constant sections, the hot
!   path (with computed ones the compiler cannot rule out an overlap of the
!   two sections of the same pool array and copies through a temporary)
SUBROUTINE Copy_Same_Level( n, nb, k, dirs )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nb, k, dirs(:)

    INTEGER :: q, d

    DO q = 1, SIZE(dirs)
        d = dirs(q)
        SELECT CASE( k )
        CASE( 1 )
            Blk_h(0,   1:ny,d, n) = Blk_h(nx,1:ny,d, nb)
            Blk_g(0,   1:ny,d, n) = Blk_g(nx,1:ny,d, nb)
        CASE( 2 )
            Blk_h(nx+1,1:ny,d, n) = Blk_h(1, 1:ny,d, nb)
            Blk_g(nx+1,1:ny,d, n) = Blk_g(1, 1:ny,d, nb)
        CASE( 3 )
            Blk_h(1:nx,0,   d, n) = Blk_h(1:nx,ny,d, nb)
            Blk_g(1:nx,0,   d, n) = Blk_g(1:nx,ny,d, nb)
        CASE( 4 )
            Blk_h(1:nx,ny+1,d, n) = Blk_h(1:nx,1, d, nb)
            Blk_g(1:nx,ny+1,d, n) = Blk_g(1:nx,1, d, nb)
        CASE( 5 )
            Blk_h(0,   0,   d, n) = Blk_h(nx,ny,d, nb)
            Blk_g(0,   0,   d, n) = Blk_g(nx,ny,d, nb)
        CASE( 6 )
            Blk_h(0,   ny+1,d, n) = Blk_h(nx,1, d, nb)
            Blk_g(0,   ny+1,d, n) = Blk_g(nx,1, d, nb)
        CASE( 7 )
            Blk_h(nx+1,0,   d, n) = Blk_h(1, ny,d, nb)
            Blk_g(nx+1,0,   d, n) = Blk_g(1, ny,d, nb)
        CASE( 8 )
            Blk_h(nx+1,ny+1,d, n) = Blk_h(1, 1, d, nb)
            Blk_g(nx+1,ny+1,d, n) = Blk_g(1, 1, d, nb)
        END SELECT
    END DO

END SUBROUTINE Copy_Same_Level
#else
#include "ghost_copy.inc"
#endif
#if GH_COPY == 4
!== push: the edge of block n towards its same-level leaf neighbour nb in direction
!   k goes into the ghost cells of nb in the opposite direction, the directions
!   chosen by what nb reads
SUBROUTINE Push( n, nb, k )
    INTEGER, INTENT(IN) :: n, nb, k
    IF( Blk(nb)%Finer ) THEN
        CALL Copy_All  ( kopp(k), Blk_h(:,:,:,nb), Blk_h(:,:,:,n), Blk_g(:,:,:,nb), Blk_g(:,:,:,n) )
    ELSEIF( Blk(nb)%Level == maxLevels ) THEN
        CALL Copy_Shift( kopp(k), Blk_h(:,:,:,nb), Blk_h(:,:,:,n), Blk_g(:,:,:,nb), Blk_g(:,:,:,n) )
    ELSE
        CALL Copy_LW   ( kopp(k), Blk_h(:,:,:,nb), Blk_h(:,:,:,n), Blk_g(:,:,:,nb), Blk_g(:,:,:,n) )
    END IF
END SUBROUTINE Push

!== all pushes of leaf block n (to be called right after its collision)
SUBROUTINE Push_Block( n )
    INTEGER, INTENT(IN) :: n
    INTEGER :: k, nb
    DO k = 1, 8
        nb = Blk(n)%nbr(k)
        IF( nb <= 0 ) CYCLE
        IF( Blk(nb)%IsLeaf .AND. Blk(nb)%Level == Blk(n)%Level ) CALL Push( n, nb, k )
    END DO
END SUBROUTINE Push_Block
#endif

!******************************************************************************

!== the ghost column of block n from its refined (same level) or coarser
!   left (s = -1) or right (s = +1) neighbour nb
SUBROUTINE Fill_Side_X( n, nb, s )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nb, s

    INTEGER  :: ig, i0, c, l, j, q, Y
    REAL(wp) :: hf(3,2*ny,2:9), gf(3,2*ny,2:9), hc(3,ny/2+2,2:9), gc(3,ny/2+2,2:9)
    REAL(wp) :: hb(1,ny,2:9), gb(1,ny,2:9)

    ig = MERGE( 0, nx+1, s < 0 )

    IF( Blk(nb)%Level < Blk(n)%Level ) THEN
        i0 = MERGE( nx+1, 0, s < 0 )
        Y  = (ny/2) * MOD(Blk(n)%J,2)
        DO q = 2, 9
        DO j = 1, ny/2+2
        DO l = 0, 2
            hc(l+1,j,q) = Blk_h(i0+s*l,Y+j-1,q, nb)
            gc(l+1,j,q) = Blk_g(i0+s*l,Y+j-1,q, nb)
        END DO
        END DO
        END DO
        CALL From_Coarse_LR( hc, hb )
        CALL From_Coarse_LR( gc, gb )
    ELSE
        i0 = MERGE( nx, 1, s < 0 )
        c = Blk(nb)%child + MERGE( 1, 0, s < 0 )      ! the facing children c, c+2
        DO q = 2, 9
        DO j = 1, ny
        DO l = 0, 2
            hf(l+1,   j,q) = Blk_h(i0+s*l,j,q, c  )
            hf(l+1,ny+j,q) = Blk_h(i0+s*l,j,q, c+2)
            gf(l+1,   j,q) = Blk_g(i0+s*l,j,q, c  )
            gf(l+1,ny+j,q) = Blk_g(i0+s*l,j,q, c+2)
        END DO
        END DO
        END DO
        CALL From_Fine_LR( hf, hb )
        CALL From_Fine_LR( gf, gb )
    END IF

    DO q = 2, 9
    DO j = 1, ny
        Blk_h(ig,j,q, n) = hb(1,j,q)
        Blk_g(ig,j,q, n) = gb(1,j,q)
    END DO
    END DO

END SUBROUTINE Fill_Side_X

!******************************************************************************

!== the ghost row of block n from its refined or coarser bottom (s = -1) or
!   top (s = +1) neighbour nb
SUBROUTINE Fill_Side_Y( n, nb, s )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nb, s

    INTEGER  :: jg, j0, c, l, i, q, X
    REAL(wp) :: hf(2*nx,3,2:9), gf(2*nx,3,2:9), hc(nx/2+2,3,2:9), gc(nx/2+2,3,2:9)
    REAL(wp) :: hb(nx,1,2:9), gb(nx,1,2:9)

    jg = MERGE( 0, ny+1, s < 0 )

    IF( Blk(nb)%Level < Blk(n)%Level ) THEN
        j0 = MERGE( ny+1, 0, s < 0 )
        X  = (nx/2) * MOD(Blk(n)%I,2)
        DO l = 0, 2
            hc(:,l+1,:) = Blk_h(X:X+nx/2+1,j0+s*l,2:9, nb)
            gc(:,l+1,:) = Blk_g(X:X+nx/2+1,j0+s*l,2:9, nb)
        END DO
        CALL From_Coarse_BT( hc, hb )
        CALL From_Coarse_BT( gc, gb )
    ELSE
        j0 = MERGE( ny, 1, s < 0 )
        c = Blk(nb)%child + MERGE( 2, 0, s < 0 )      ! the facing children c, c+1
        DO l = 0, 2
            hf(   1:  nx,l+1,:) = Blk_h(1:nx,j0+s*l,2:9, c  )
            hf(nx+1:2*nx,l+1,:) = Blk_h(1:nx,j0+s*l,2:9, c+1)
            gf(   1:  nx,l+1,:) = Blk_g(1:nx,j0+s*l,2:9, c  )
            gf(nx+1:2*nx,l+1,:) = Blk_g(1:nx,j0+s*l,2:9, c+1)
        END DO
        CALL From_Fine_BT( hf, hb )
        CALL From_Fine_BT( gf, gb )
    END IF

    DO q = 2, 9
    DO i = 1, nx
        Blk_h(i,jg,q, n) = hb(i,1,q)
        Blk_g(i,jg,q, n) = gb(i,1,q)
    END DO
    END DO

END SUBROUTINE Fill_Side_Y

!******************************************************************************

!== the ghost corner of block n from its refined or coarser neighbour in
!   direction k = 5..8
SUBROUTINE Fill_Corner( n, k )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, k

    INTEGER  :: nb, sx, sy, ig, jg, i0, j0, px, py, c, l, m, q
    REAL(wp) :: hf(3,3,2:9), gf(3,3,2:9), hb(1,1,2:9), gb(1,1,2:9)

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
        DO q = 2, 9
        DO m = 0, 2
        DO l = 0, 2
            hf(l+1,m+1,q) = Blk_h(i0+sx*l,j0+sy*m,q, nb)
            gf(l+1,m+1,q) = Blk_g(i0+sx*l,j0+sy*m,q, nb)
        END DO
        END DO
        END DO
        CALL From_Coarse_Corner( hf, hb )
        CALL From_Coarse_Corner( gf, gb )
    ELSE
        i0 = MERGE( nx, 1, sx < 0 )
        j0 = MERGE( ny, 1, sy < 0 )
        c  = Blk(nb)%child + facing(1,k)             ! the facing child
        DO q = 2, 9
        DO m = 0, 2
        DO l = 0, 2
            hf(l+1,m+1,q) = Blk_h(i0+sx*l,j0+sy*m,q, c)
            gf(l+1,m+1,q) = Blk_g(i0+sx*l,j0+sy*m,q, c)
        END DO
        END DO
        END DO
        CALL From_Fine_Corner( hf, hb )
        CALL From_Fine_Corner( gf, gb )
    END IF

    Blk_h(ig,jg,2:9, n) = hb(1,1,:)
    Blk_g(ig,jg,2:9, n) = gb(1,1,:)

END SUBROUTINE Fill_Corner

!******************************************************************************

!== the strip interpolations for the three strip shapes of the 8 moving
!   directions, with the shape fixed at compile time (from Biquadratic.inc):
!   left/right edges (LR), bottom/top edges (BT) and corners. With fixed
!   shapes the compiler unrolls the short loops: ~27% less time in the
!   ghost cells of h and g than with the sizes as arguments.
#define FINE_NAME    From_Fine_LR
#define COARSE_NAME  From_Coarse_LR
#define SIZE_ARGS
#define FINE_SIZES   INTEGER, PARAMETER :: i2 = 3, j2 = 2*ny, x2 = 1, y2 = ny, nk = 8
#define COARSE_SIZES INTEGER, PARAMETER :: i2 = 3, j2 = ny/2+2, x2 = 1, y2 = ny, nk = 8
#include "Biquadratic.inc"

!******************************************************************************

#define FINE_NAME    From_Fine_BT
#define COARSE_NAME  From_Coarse_BT
#define SIZE_ARGS
#define FINE_SIZES   INTEGER, PARAMETER :: i2 = 2*nx, j2 = 3, x2 = nx, y2 = 1, nk = 8
#define COARSE_SIZES INTEGER, PARAMETER :: i2 = nx/2+2, j2 = 3, x2 = nx, y2 = 1, nk = 8
#include "Biquadratic.inc"

!******************************************************************************

#define FINE_NAME    From_Fine_Corner
#define COARSE_NAME  From_Coarse_Corner
#define SIZE_ARGS
#define FINE_SIZES   INTEGER, PARAMETER :: i2 = 3, j2 = 3, x2 = 1, y2 = 1, nk = 8
#define COARSE_SIZES INTEGER, PARAMETER :: i2 = 3, j2 = 3, x2 = 1, y2 = 1, nk = 8
#include "Biquadratic.inc"

#ifdef GH_CFDIRS
#include "ghost_from_coarse.inc"
#endif

END MODULE GH_MOD

#undef GH_MOD
#undef GH_COPY
#undef GH_PROBE_NOCOPY
#undef GH_PROBE_NOCF
#undef GH_PROBE_NOWALL
#undef GH_PROBE_NOCOARSE
#undef GH_PROBE_NOFINE
#undef GH_CFDIRS
