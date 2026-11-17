!================================================
!   Block storage of the quadtree mesh
!
!   Every block, leaf or not, lives in Blk(1:Nblks); slots of derefined
!   blocks are kept in a free list (Root ... emptyBlk, nUnused slots) and
!   reused by the next refinement. The leaves are listed in LeafBlks and,
!   sorted by level, in LeafByLevel.
!================================================
MODULE Blocks_Mod
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: maxLevels, MaxBlks
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Blk, Nblks, nLeaf, LeafBlks, LeafByLevel, LevelStart, LevelMin, LevelMax, &
              nUnused, Root, emptyBlk, Linked_Free_Blocks

    TYPE Ptr
        INTEGER, POINTER :: n
    END TYPE Ptr
    TYPE BLOCK
        INTEGER :: Level, dx
        INTEGER :: I, J
        INTEGER :: Parent
        INTEGER, POINTER :: child(:)
        TYPE(Ptr)       :: nghbr(1:8)
        LOGICAL :: IsLeaf
        LOGICAL :: hasLeaf
        LOGICAL :: Refin
        LOGICAL :: Deref
        REAL(wp), DIMENSION(:, :, :), ALLOCATABLE :: h, g
        REAL(wp), DIMENSION(:, :),   ALLOCATABLE :: C, P, Ux, Uy
    END TYPE BLOCK
    TYPE(BLOCK) :: Blk(MaxBlks)

    INTEGER :: Nblks, nLeaf

    !== coarsest and finest level that currently has leaves
    INTEGER :: LevelMin, LevelMax

    INTEGER, ALLOCATABLE :: LeafBlks(:)

    !== leaf blocks sorted by refinement level (rebuilt in Array_of_Leaf_Blocks):
    !   the leaves on level L are LeafByLevel( LevelStart(L) : LevelStart(L+1)-1 )
    INTEGER, ALLOCATABLE :: LeafByLevel(:)
    INTEGER :: LevelStart(0:maxLevels+1)

    !== free list of the slots of derefined blocks (groups of 4 siblings)
    TYPE Linked_Free_Blocks
        INTEGER :: N
        TYPE(Linked_Free_Blocks), POINTER :: next
    END TYPE Linked_Free_Blocks
    TYPE(Linked_Free_Blocks), POINTER :: Root, emptyBlk
    INTEGER :: nUnused = 0

END MODULE Blocks_Mod
