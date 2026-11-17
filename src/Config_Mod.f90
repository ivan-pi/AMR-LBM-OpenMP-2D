!================================================
!   Compile-time configuration of the solver
!
!   Set by CMake (-DAMR_BLOCK_SIZE, -DAMR_L0, -DAMR_MIN_LEVEL): square blocks
!   of nx x ny cells, a domain of L0 x 3*L0 lattice units at the finest
!   level, and a uniform background refinement to minLevels. As constants
!   they let the compiler specialise every block kernel for the block size.
!================================================
MODULE Config_Mod
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: nx, ny, L0, maxLevels, minLevels, MaxBlks

#if !defined(AMR_BLOCK_SIZE) || !defined(AMR_L0) || !defined(AMR_MIN_LEVEL)
#error "AMR_BLOCK_SIZE, AMR_L0 and AMR_MIN_LEVEL must be defined (see CMakeLists.txt)"
#endif
    INTEGER, PARAMETER :: nx = AMR_BLOCK_SIZE
    INTEGER, PARAMETER :: ny = AMR_BLOCK_SIZE
    INTEGER, PARAMETER :: L0 = AMR_L0
    INTEGER, PARAMETER :: maxLevels = BIT_SIZE(L0) - 1 - LEADZ(L0/nx)   ! log2(L0/nx)
    INTEGER, PARAMETER :: minLevels = AMR_MIN_LEVEL

    !== upper bound of the number of blocks (all levels of the quadtree)
    INTEGER, PARAMETER :: MaxBlks = 6 * (L0/nx)**2

END MODULE Config_Mod
