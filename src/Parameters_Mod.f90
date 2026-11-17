
MODULE Parameters_Mod
    USE Kinds_Mod, ONLY: wp

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Nblks, nLeaf, T, KK, LevelMin, LevelMax, nx, ny, L0, maxLevels, minLevels, &
              MaxBlks, LeafBlks, LeafByLevel, LevelStart, Blk, Root, emptyBlk, &
              Linked_Free_Blocks

    INTEGER :: T, Nblks, nLeaf

    INTEGER :: KK = 0

    INTEGER :: LevelMin, LevelMax

    !== case size, set by CMake (-DAMR_BLOCK_SIZE, -DAMR_L0, -DAMR_MIN_LEVEL):
    !   square blocks of nx x ny cells, a domain of L0 x 3*L0 lattice units at
    !   the finest level, and a uniform background refinement to minLevels
#if !defined(AMR_BLOCK_SIZE) || !defined(AMR_L0) || !defined(AMR_MIN_LEVEL)
#error "AMR_BLOCK_SIZE, AMR_L0 and AMR_MIN_LEVEL must be defined (see CMakeLists.txt)"
#endif
    INTEGER, PARAMETER :: nx = AMR_BLOCK_SIZE
    INTEGER, PARAMETER :: ny = AMR_BLOCK_SIZE
    INTEGER, PARAMETER :: L0 = AMR_L0
    INTEGER, PARAMETER :: maxLevels = BIT_SIZE(L0) - 1 - LEADZ(L0/nx)   ! log2(L0/nx)
    INTEGER, PARAMETER :: minLevels = AMR_MIN_LEVEL

    INTEGER, PARAMETER :: MaxBlks = 6 * (L0/nx)**2

    INTEGER, ALLOCATABLE :: LeafBlks(:)

    !== leaf blocks sorted by refinement level (rebuilt in Array_of_Leaf_Blocks):
    !   the leaves on level L are LeafByLevel( LevelStart(L) : LevelStart(L+1)-1 )
    INTEGER, ALLOCATABLE :: LeafByLevel(:)
    INTEGER :: LevelStart(0:maxLevels+1)

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

    TYPE Linked_Free_Blocks
        INTEGER :: N
        TYPE(Linked_Free_Blocks), POINTER :: next
    END TYPE Linked_Free_Blocks
    TYPE(Linked_Free_Blocks), POINTER :: Root, emptyBlk

END MODULE Parameters_Mod
