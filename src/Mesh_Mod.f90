!================================================
!   Block (quadtree) mesh: the domain of three base blocks, children,
!   leaf lists
!================================================
MODULE Mesh_Mod
    USE Kinds_Mod, ONLY: wp, dp
    USE Neighbors_Mod, ONLY: Find_Neighbors
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Create_Domain, Add_Child_Blocks, Array_of_Leaf_Blocks, Determine_Min_Max_Levels

CONTAINS

!== the domain L0 x 3*L0 as three base blocks stacked vertically, refined
!   uniformly to the coarsest level minLevels
SUBROUTINE Create_Domain
    USE Config_Mod, ONLY: minLevels
    USE Blocks_Mod, ONLY: LevelMin, LevelMax, Init_Block_Pool
    IMPLICIT NONE

    CALL Init_Block_Pool

    LevelMin = minLevels
    LevelMax = minLevels

    CALL Initialize_Base_Block( 1 ) ! lower  part of the domain
    CALL Initialize_Base_Block( 2 ) ! middle part of the domain
    CALL Initialize_Base_Block( 3 ) ! upper  part of the domain
    CALL Fix_Initial_Neighbors_3Blocks

    CALL Initial_Mesh( LevelMin )

END SUBROUTINE Create_Domain

!******************************************************************************

SUBROUTINE Initialize_Base_Block( N )
    USE Config_Mod, ONLY: L0, nx
    USE Blocks_Mod, ONLY: Blk, Nblks, Link_Neighbor
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N

    INTEGER :: I

    Nblks = N

    Blk(N)%Level = 0
    Blk(N)%dx    = L0 / nx
    Blk(N)%I     = 0
    Blk(N)%J     = N-1  !N=1: the lower block

    Blk(N)%child(:)= -1         ! Children will be determined soon
    Blk(N)%Parent  = -1         ! This is the base block and has no parent

    !============ This is the numbering for the neighboring blocks
    !   6   4   8
    !   1       2
    !   5   3   7
    !============ This is the numbering for the neighboring blocks

    DO I = 1, 8
        CALL Link_Neighbor( N, I, -I, 0 )   !External Boundaries
    END DO
    CALL Link_Neighbor( N, 1, N, 0 )        !Periodic
    CALL Link_Neighbor( N, 2, N, 0 )        !Periodic

    Blk(N)%IsLeaf  = .True.
    Blk(N)%hasLeaf = .False.
    Blk(N)%Refin   = .False.
    Blk(N)%Deref   = .False.

END SUBROUTINE Initialize_Base_Block

!******************************************************************************

SUBROUTINE Fix_Initial_Neighbors_3Blocks
    USE Blocks_Mod, ONLY: Link_Neighbor
    IMPLICIT NONE

    CALL Link_Neighbor( 1, 4, 2, 0 )
    CALL Link_Neighbor( 1, 6, 2, 0 )
    CALL Link_Neighbor( 1, 8, 2, 0 )

    CALL Link_Neighbor( 2, 3, 1, 0 )
    CALL Link_Neighbor( 2, 4, 3, 0 )
    CALL Link_Neighbor( 2, 5, 1, 0 )
    CALL Link_Neighbor( 2, 6, 3, 0 )
    CALL Link_Neighbor( 2, 7, 1, 0 )
    CALL Link_Neighbor( 2, 8, 3, 0 )

    CALL Link_Neighbor( 3, 3, 2, 0 )
    CALL Link_Neighbor( 3, 5, 2, 0 )
    CALL Link_Neighbor( 3, 7, 2, 0 )

END

!******************************************************************************

SUBROUTINE Initial_Mesh( nLevels )
    USE Blocks_Mod, ONLY: Blk, Nblks
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nLevels

    INTEGER :: L, n

    DO L = 1, nLevels
        DO n = 1, Nblks     ! Nblks will change inside Sub:Add_Child_Blocks
            IF( Blk(n)%Level == L-1 ) THEN
                CALL Add_Child_Blocks( n, Nblks+1 )
            END IF
        END DO
    END DO

    CALL Array_of_Leaf_Blocks

END SUBROUTINE Initial_Mesh

!******************************************************************************

SUBROUTINE Add_Child_Blocks( par, FreeBlk )
    USE Blocks_Mod, ONLY: Blk, Blk_h, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy, Nblks, Reserve_Blocks
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: par, FreeBlk

    INTEGER :: N, grandpar, I, J

    N = FreeBlk - 1

    CALL Reserve_Blocks( N+4 )

    DO i = 1, 4
        Blk(par)%child(i) = N + i
    END DO
    Blk(par)%IsLeaf  = .False.
    Blk(par)%hasLeaf = .True.

    grandpar = Blk(par)%Parent
    IF( grandpar > 0 ) Blk(grandpar)%hasLeaf = .False.

    DO J = 0, 1
    DO I = 0, 1

        N = N + 1

        Blk(N)%Level = Blk(par)%Level + 1
        Blk(N)%dx    = Blk(par)%dx / 2

        Blk(N)%I = 2 * Blk(par)%I + I
        Blk(N)%J = 2 * Blk(par)%J + J

        Blk(N)%Parent = par

        Blk(N)%child(:) = N

        CALL Find_Neighbors( par, N, I, J )

        Blk(N)%IsLeaf  = .True.
        Blk(N)%hasLeaf = .False.
        Blk(N)%Refin   = .False.
        Blk(N)%Deref   = .False.

        Blk_h(:,:,:,N)  = 0
        Blk_g(:,:,:,N)  = 0
        Blk_C(:,:,N)  = 0
        Blk_P(:,:,N)  = 0
        Blk_Ux(:,:,N) = 0
        Blk_Uy(:,:,N) = 0

    END DO
    END DO

    Nblks = MAX(Nblks,N)    ! the slots may have been reused from the free queue

END SUBROUTINE Add_Child_Blocks

!******************************************************************************

SUBROUTINE Array_of_Leaf_Blocks
    USE Config_Mod, ONLY: maxLevels
    USE Blocks_Mod, ONLY: Blk, Nblks, LeafBlks, nLeaf, LeafByLevel, LevelStart
    IMPLICIT NONE

    INTEGER :: n, i, L
    INTEGER :: tmpArray(Nblks), cnt(0:maxLevels+1)

    i = 0
    DO n = 1, Nblks
        IF( Blk(n)%IsLeaf ) THEN
            i = i + 1
            tmpArray(i) = n
        END IF
    END DO

    nLeaf = i
    IF( ALLOCATED(LeafBlks) ) DEALLOCATE( LeafBlks )
    ALLOCATE( LeafBlks(i) )

    LeafBlks(:) = tmpArray(1:i)

    !== counting sort of the leaves by level (stable: keeps block order per level)
    cnt = 0
    DO i = 1, nLeaf
        L = Blk(LeafBlks(i))%Level
        cnt(L) = cnt(L) + 1
    END DO
    LevelStart(0) = 1
    DO L = 0, maxLevels
        LevelStart(L+1) = LevelStart(L) + cnt(L)
    END DO

    IF( ALLOCATED(LeafByLevel) ) DEALLOCATE( LeafByLevel )
    ALLOCATE( LeafByLevel(nLeaf) )
    cnt(0:maxLevels) = LevelStart(0:maxLevels)
    DO i = 1, nLeaf
        n = LeafBlks(i)
        L = Blk(n)%Level
        LeafByLevel(cnt(L)) = n
        cnt(L) = cnt(L) + 1
    END DO

END

!******************************************************************************

SUBROUTINE Determine_Min_Max_Levels
    USE Config_Mod, ONLY: maxLevels
    USE Blocks_Mod, ONLY: LevelMin, LevelMax, LevelStart
    IMPLICIT NONE

    INTEGER :: L

    LevelMin =  100
    LevelMax = -100
    DO L = 0, maxLevels
        IF( LevelStart(L+1) > LevelStart(L) ) THEN
            LevelMin = MIN( LevelMin, L )
            LevelMax = MAX( LevelMax, L )
        END IF
    END DO

END SUBROUTINE Determine_Min_Max_Levels


END MODULE Mesh_Mod
