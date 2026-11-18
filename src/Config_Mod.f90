!================================================
!   Compile-time configuration of the solver
!
!   Set by CMake (-DAMR_BLOCK_SIZE, -DAMR_L0, -DAMR_MIN_LEVEL): square blocks
!   of nx x ny cells, a domain of L0 x 3*L0 lattice units at the finest
!   level, and a uniform background refinement to minLevels. As constants
!   they let the compiler specialise every block kernel for the block size.
!================================================
MODULE Config_Mod
    USE Kinds_Mod, ONLY: wp
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: nx, ny, L0, maxLevels, minLevels, MaxBlks, inv_dx_level

#if !defined(AMR_BLOCK_SIZE) || !defined(AMR_L0) || !defined(AMR_MIN_LEVEL)
#error "AMR_BLOCK_SIZE, AMR_L0 and AMR_MIN_LEVEL must be defined (see CMakeLists.txt)"
#endif
    INTEGER, PARAMETER :: nx = AMR_BLOCK_SIZE
    INTEGER, PARAMETER :: ny = AMR_BLOCK_SIZE
    INTEGER, PARAMETER :: L0 = AMR_L0
    INTEGER, PARAMETER :: maxLevels = BIT_SIZE(L0) - 1 - LEADZ(L0/nx)   ! log2(L0/nx)
    INTEGER, PARAMETER :: minLevels = AMR_MIN_LEVEL

    !== number of blocks of the fully refined quadtree (three base blocks,
    !   4**L blocks below each on level L): an upper bound of the live blocks
    INTEGER, PARAMETER :: MaxBlks = 4**(maxLevels+1) - 1

    !== 1/dx of the cells on level L, in finest-level lattice units. dx is
    !   (L0/nx)/2**L, a power of two, so the reciprocals are exact. Kernels take
    !   these instead of the integer cell size (no conversions, no divisions).
    INTEGER :: lev      ! implied-do variable of the table
    REAL(wp), PARAMETER :: inv_dx_level(0:maxLevels) = &
        [( REAL(nx,wp) / REAL(L0,wp) * 2._wp**lev, lev = 0, maxLevels )]

END MODULE Config_Mod
