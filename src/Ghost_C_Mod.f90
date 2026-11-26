!================================================
!   Ghost cells of the phase field C (ng = 1 or 2 layers), and of the
!   second scalar C2 of two-scalar models (Blocks_Mod nScalars = 2)
!
!   As for h and g (Ghost_hg_Mod, where the strip layout is described): a
!   copy from a same-level leaf, the restriction of the facing children of a
!   refined neighbour, or the interpolation of a coarser one, with one
!   routine per strip shape. The strips across an edge are 2+ng fine or 3
!   coarse cells deep. The strip routines are the template
!   Ghost_C_Strips.inc, instantiated for one and for two layers: with the
!   number of layers and so every strip shape fixed at compile time the fill
!   is 2.3-3.4 times faster than with one set of routines for any ng
!   (docs/benchmarks/2026-11-26-kernel-performance.md).
!   Pushing the copies from the streaming loop, as for h and
!   g, gains nothing more on 8x8 blocks and is not done.
!================================================
MODULE Ghost_C_Mod
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk, Blk_C, Blk_C2, nScalars, di, dj, facing
    USE BC_Mod, ONLY: Symmetric_C_3, Symmetric_C_4
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Fill_Ghost_Cells_C

CONTAINS


!== fill ng (1 or 2) ghost layers of C on all leaf blocks, in two passes over
!   the levels as for h and g (Fill_Ghost_Cells_hg): copies, restrictions and
!   walls first, then the interpolations from coarser neighbours. The wall
!   ghost cells in the end columns next to an interpolated side copy its
!   ghost cells and are redone then; the collision stencils read them.
SUBROUTINE Fill_Ghost_Cells_C( ng )
    USE Blocks_Mod, ONLY: LevelMin, LevelMax, LeafByLevel, LevelStart
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ng

    INTEGER :: L, k, m, n, X1

    IF( ng /= 1 .AND. ng /= 2 ) ERROR STOP 'Fill_Ghost_Cells_C: one or two ghost layers'

    DO L = LevelMin, LevelMax
!$OMP DO schedule(static)
        DO m = LevelStart(L), LevelStart(L+1)-1
            n = LeafByLevel(m)
            CALL Fill( n, L, .FALSE., Blk_C )
            CALL Walls( n, 1-ng, nx+ng, Blk_C )
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
            IF( Blk(nb)%IsLeaf .AND. .NOT. second ) THEN        ! same-level leaf
                IF( ng == 1 ) THEN
                    CALL Copy_1( k, C(:,:,n), C(:,:,nb) )
                ELSE
                    CALL Copy_2( k, C(:,:,n), C(:,:,nb) )
                END IF
            ELSEIF( ng == 1 ) THEN
                IF( k <= 2 ) THEN
                    CALL Side_X_1( n, nb, di(k), C )
                ELSEIF( k <= 4 ) THEN
                    CALL Side_Y_1( n, nb, dj(k), C )
                ELSE
                    CALL Corner_1( n, k, C )
                END IF
            ELSE
                IF( k <= 2 ) THEN
                    CALL Side_X_2( n, nb, di(k), C )
                ELSEIF( k <= 4 ) THEN
                    CALL Side_Y_2( n, nb, dj(k), C )
                ELSE
                    CALL Corner_2( n, k, C )
                END IF
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

!== the copies and strip fills for one and for two ghost layers
#define NGC          1
#define COPY_NAME    Copy_1
#define SIDE_X_NAME  Side_X_1
#define SIDE_Y_NAME  Side_Y_1
#define CORNER_NAME  Corner_1
#include "Ghost_C_Strips.inc"

!******************************************************************************

#define NGC          2
#define COPY_NAME    Copy_2
#define SIDE_X_NAME  Side_X_2
#define SIDE_Y_NAME  Side_Y_2
#define CORNER_NAME  Corner_2
#include "Ghost_C_Strips.inc"

END MODULE Ghost_C_Mod
