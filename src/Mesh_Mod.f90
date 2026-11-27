!================================================
!   Block (quadtree) mesh: the domain of base blocks, children, leaf lists,
!   neighbour table
!================================================
MODULE Mesh_Mod
    USE Kinds_Mod, ONLY: i64
    USE Blocks_Mod, ONLY: Block_t, Blk, Nblks
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Create_Domain, Add_Child_Blocks, Array_of_Leaf_Blocks
    PUBLIC :: t_nbr_table, n_nbr_rebuilds

    !== time spent in Rebuild_Neighbor_Table (SYSTEM_CLOCK ticks) and number of calls
    INTEGER(i64) :: t_nbr_table = 0
    INTEGER    :: n_nbr_rebuilds = 0

CONTAINS

!== the domain of nBaseX x nBaseY base blocks (slots 1, 2, ... row by row from
!   the bottom), refined uniformly to the coarsest level minLevels
SUBROUTINE Create_Domain
    USE Config_Mod, ONLY: minLevels, nBaseX, nBaseY
    USE Blocks_Mod, ONLY: Init_Block_Pool
    IMPLICIT NONE

    INTEGER :: L, n, I, J

    CALL Init_Block_Pool

    DO J = 0, nBaseY-1
    DO I = 0, nBaseX-1
        Nblks = Nblks + 1
        Blk(Nblks) = Block_t( Level=0, I=I, J=J, IsLeaf=.TRUE. )
    END DO
    END DO

    DO L = 1, minLevels
        DO n = 1, Nblks     ! Nblks grows inside Add_Child_Blocks
            IF( Blk(n)%Level == L-1 ) CALL Add_Child_Blocks( n, Nblks+1 )
        END DO
    END DO

    CALL Array_of_Leaf_Blocks

END SUBROUTINE Create_Domain

!******************************************************************************

!== children of leaf par in the slots c1:c1+3 (from the free queue, or c1 = Nblks+1)
SUBROUTINE Add_Child_Blocks( par, c1 )
    USE Blocks_Mod, ONLY: Reserve_Blocks, Clear_Slot
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: par, c1

    INTEGER :: c, I, J

    CALL Reserve_Blocks( c1+3 )

    Blk(par)%child  = c1
    Blk(par)%IsLeaf = .FALSE.

    c = c1
    DO J = 0, 1
    DO I = 0, 1
        Blk(c) = Block_t( Level=Blk(par)%Level+1, I=2*Blk(par)%I+I, J=2*Blk(par)%J+J, &
                          Parent=par, IsLeaf=.TRUE. )
        CALL Clear_Slot( c )
        c = c + 1
    END DO
    END DO

    Nblks = MAX( Nblks, c1+3 )

END SUBROUTINE Add_Child_Blocks

!******************************************************************************

!== the leaf lists and the neighbour table of the current mesh
SUBROUTINE Array_of_Leaf_Blocks
    USE Config_Mod, ONLY: maxLevels, nBaseX, nBaseY
    USE Blocks_Mod, ONLY: nLeaf, LeafByLevel, LevelStart, LevelMin, LevelMax
    IMPLICIT NONE

    INTEGER :: n, i, L, nz
    INTEGER :: zorder(Nblks), cnt(0:maxLevels+1)

    !== the leaves in Z-order (Morton order): a depth-first walk of the
    !   quadtree, as the children 1..4 of a block are (0,0), (1,0), (0,1), (1,1)
    nz = 0
    DO n = 1, nBaseX*nBaseY         ! the base blocks, row by row from the bottom
        CALL Visit( n )
    END DO
    nLeaf = nz

    !== counting sort of the leaves by level, stable: Z-order within a level, so
    !   that a static share of a level is a compact patch whose neighbours are
    !   mostly blocks of the same thread
    cnt = 0
    DO i = 1, nLeaf
        L = Blk(zorder(i))%Level
        cnt(L) = cnt(L) + 1
    END DO
    LevelStart(0) = 1
    DO L = 0, maxLevels
        LevelStart(L+1) = LevelStart(L) + cnt(L)
    END DO
    LevelMin = FINDLOC( cnt(0:maxLevels) > 0, .TRUE., DIM=1 ) - 1
    LevelMax = FINDLOC( cnt(0:maxLevels) > 0, .TRUE., DIM=1, BACK=.TRUE. ) - 1

    IF( ALLOCATED(LeafByLevel) ) DEALLOCATE( LeafByLevel )
    ALLOCATE( LeafByLevel(nLeaf) )
    cnt(0:maxLevels) = LevelStart(0:maxLevels)
    DO i = 1, nLeaf
        n = zorder(i)
        L = Blk(n)%Level
        LeafByLevel(cnt(L)) = n
        cnt(L) = cnt(L) + 1
    END DO

    CALL Rebuild_Neighbor_Table

CONTAINS

    RECURSIVE SUBROUTINE Visit( b )
        INTEGER, INTENT(IN) :: b
        INTEGER :: c
        IF( Blk(b)%IsLeaf ) THEN
            nz = nz + 1
            zorder(nz) = b
        ELSE
            DO c = 0, 3
                CALL Visit( Blk(b)%child + c )
            END DO
        END IF
    END SUBROUTINE Visit

END

!******************************************************************************

!== Blk(:)%nbr of every live block from the block positions alone: the block
!   at the neighbouring position on the same level, else the block one level
!   coarser that covers it (the 2:1 balance guarantees one of them exists).
!   Periodic in x (and in y with periodicY); otherwise above the top or below
!   the bottom the neighbour is the boundary code -k. Called after every change
!   of the mesh; the block positions are the only neighbour information that
!   is stored.
SUBROUTINE Rebuild_Neighbor_Table
    USE Config_Mod, ONLY: MaxBlks, nBaseX, nBaseY, periodicY
    USE Blocks_Mod, ONLY: di, dj
    IMPLICIT NONE

    !== block at position (I,J) of level L: at( Slot(L,I,J) ), 0 if none
    !   (kept all 0 between calls, so only the entries of live blocks are touched)
    INTEGER, SAVE :: at(MaxBlks) = 0

    INTEGER :: n, k, L, I, J, b
    INTEGER(i64) :: t0, t1

    CALL SYSTEM_CLOCK( t0 )

    DO n = 1, Nblks
        IF( Blk(n)%Level >= 0 ) at( Slot( Blk(n)%Level, Blk(n)%I, Blk(n)%J ) ) = n
    END DO

    DO n = 1, Nblks
        L = Blk(n)%Level
        IF( L < 0 ) CYCLE
        DO k = 1, 8
            J = Blk(n)%J + dj(k)
            IF( periodicY ) J = MODULO( J, nBaseY*2**L )
            IF( J < 0 .OR. J >= nBaseY*2**L ) THEN
                Blk(n)%nbr(k) = -k
                CYCLE
            END IF
            I = MODULO( Blk(n)%I + di(k), nBaseX*2**L )
            b = at( Slot(L, I, J) )
            IF( b == 0 .AND. L > 0 ) b = at( Slot(L-1, I/2, J/2) )
            IF( b == 0 ) THEN
                PRINT '(A,6I7)', ' n, k, L, I, J:', n, k, L, Blk(n)%I, Blk(n)%J
                ERROR STOP 'Rebuild_Neighbor_Table: no neighbour (2:1 balance violated?)'
            END IF
            Blk(n)%nbr(k) = b
        END DO
    END DO

    DO n = 1, Nblks
        IF( Blk(n)%Level >= 0 ) at( Slot( Blk(n)%Level, Blk(n)%I, Blk(n)%J ) ) = 0
    END DO

    !== blocks next to a refined same-level block (the finer blocks beyond it
    !   interpolate from their ghost cells) and next to a coarser block
    !   (their ghost cells there are interpolated)
    DO n = 1, Nblks
        IF( Blk(n)%Level < 0 ) CYCLE
        Blk(n)%Finer   = .FALSE.
        Blk(n)%Coarser = .FALSE.
        DO k = 1, 8
            b = Blk(n)%nbr(k)
            IF( b <= 0 ) CYCLE
            IF( Blk(b)%Level == Blk(n)%Level .AND. .NOT. Blk(b)%IsLeaf ) Blk(n)%Finer = .TRUE.
            IF( Blk(b)%Level <  Blk(n)%Level ) Blk(n)%Coarser = .TRUE.
        END DO
    END DO

    CALL SYSTEM_CLOCK( t1 )
    t_nbr_table    = t_nbr_table + (t1 - t0)
    n_nbr_rebuilds = n_nbr_rebuilds + 1

CONTAINS

    !== slot of position (I,J) on level L (levels stored one after the other,
    !   nBaseX*nBaseY*4**L positions on level L)
    PURE INTEGER FUNCTION Slot( L, I, J )
        INTEGER, INTENT(IN) :: L, I, J
        Slot = 1 + nBaseX*nBaseY*(4**L - 1)/3 + I + J*nBaseX*2**L
    END FUNCTION Slot

END SUBROUTINE Rebuild_Neighbor_Table


END MODULE Mesh_Mod
