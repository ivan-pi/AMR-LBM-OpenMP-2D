!================================================
!   Neighbour links of newly created child blocks
!================================================
MODULE Neighbors_Mod
    USE Blocks_Mod, ONLY: Neighbor_of, Link_Neighbor, Copy_Neighbor_Link
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Find_Neighbors

CONTAINS


SUBROUTINE Find_Neighbors( par, n, I, J )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: par, n, I, J

    INTEGER :: child

    child = 1 + I + 2*J

    CALL Neighbor( par, n, 1+I, child+1-2*I )   ! Left   & Right neighbors
    CALL Neighbor( par, n, 3+J, child+2-4*J )   ! Bottom & Top   neighbors

    CALL Diagonal_Neighbor( par, n, child )

CONTAINS

    SUBROUTINE Diagonal_Neighbor( p, n, c )
        INTEGER, INTENT(IN) :: p, n, c

        !== boundary codes of the corner neighbours
        INTEGER, PARAMETER :: SW = -5
        INTEGER, PARAMETER :: NW = -6
        INTEGER, PARAMETER :: SE = -7
        INTEGER, PARAMETER :: NE = -8

        IF( c == 1 ) THEN
            IF( Neighbor_of(p, 5) < 0 ) THEN
                CALL Copy_Neighbor_Link( n, p, 5 )
            ELSE
                CALL Link_Neighbor( n, 5, Neighbor_of(p, 5), 4 )
            END IF
            IF( Neighbor_of(p, 1) < 0 ) THEN
                CALL Link_Neighbor( n, 6, NW, 0 )
            ELSE
                CALL Link_Neighbor( n, 6, Neighbor_of(p, 1), 4 )
            END IF
            IF( Neighbor_of(p, 3) < 0 ) THEN
                CALL Link_Neighbor( n, 7, SE, 0 )
            ELSE
                CALL Link_Neighbor( n, 7, Neighbor_of(p, 3), 4 )
            END IF
            CALL Link_Neighbor( n, 8, p, 4 )
        ELSEIF( c == 2 ) THEN
            IF( Neighbor_of(p, 3) < 0 ) THEN
                CALL Link_Neighbor( n, 5, SW, 0 )
            ELSE
                CALL Link_Neighbor( n, 5, Neighbor_of(p, 3), 3 )
            END IF
            CALL Link_Neighbor( n, 6, p, 3 )
            IF( Neighbor_of(p, 7) < 0 ) THEN
                CALL Copy_Neighbor_Link( n, p, 7 )
            ELSE
                CALL Link_Neighbor( n, 7, Neighbor_of(p, 7), 3 )
            END IF
            IF( Neighbor_of(p, 2) < 0 ) THEN
                CALL Link_Neighbor( n, 8, NE, 0 )
            ELSE
                CALL Link_Neighbor( n, 8, Neighbor_of(p, 2), 3 )
            END IF
        ELSEIF( c == 3 ) THEN
            IF( Neighbor_of(p, 1) < 0 ) THEN
                CALL Link_Neighbor( n, 5, SW, 0 )
            ELSE
                CALL Link_Neighbor( n, 5, Neighbor_of(p, 1), 2 )
            END IF
            IF( Neighbor_of(p, 6) < 0 ) THEN
                CALL Copy_Neighbor_Link( n, p, 6 )
            ELSE
                CALL Link_Neighbor( n, 6, Neighbor_of(p, 6), 2 )
            END IF
            CALL Link_Neighbor( n, 7, p, 2 )
            IF( Neighbor_of(p, 4) < 0 ) THEN
                CALL Link_Neighbor( n, 8, NE, 0 )
            ELSE
                CALL Link_Neighbor( n, 8, Neighbor_of(p, 4), 2 )
            END IF
        ELSEIF( c == 4 ) THEN
            CALL Link_Neighbor( n, 5, p, 1 )
            IF( Neighbor_of(p, 4) < 0 ) THEN
                CALL Link_Neighbor( n, 6, NW, 0 )
            ELSE
                CALL Link_Neighbor( n, 6, Neighbor_of(p, 4), 1 )
            END IF
            IF( Neighbor_of(p, 2) < 0 ) THEN
                CALL Link_Neighbor( n, 7, SE, 0 )
            ELSE
                CALL Link_Neighbor( n, 7, Neighbor_of(p, 2), 1 )
            END IF
            IF( Neighbor_of(p, 8) < 0 ) THEN
                CALL Copy_Neighbor_Link( n, p, 8 )
            ELSE
                CALL Link_Neighbor( n, 8, Neighbor_of(p, 8), 1 )
            END IF
        END IF

    END SUBROUTINE Diagonal_Neighbor

END SUBROUTINE Find_Neighbors

!******************************************************************************

SUBROUTINE Neighbor( p, n, n2, c )
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: p, n, n2, c  ! n2: neighbor inside the adjacent block

    INTEGER :: n1   ! neighbor within the current block

    IF( Neighbor_of(p, n2) < 0 ) THEN
        CALL Copy_Neighbor_Link( n, p, n2 )
    ELSE
        CALL Link_Neighbor( n, n2, Neighbor_of(p, n2), c )
    END IF

    n1 = n2+1 - 2*MOD(n2-1,2)   ! opposite neighbor

    CALL Link_Neighbor( n, n1, p, c )

END


END MODULE Neighbors_Mod
