!================================================
!   Compile-time configuration of the solver
!
!   Set by CMake (-DAMR_BLOCK_SIZE, -DAMR_L0, -DAMR_MIN_LEVEL,
!   -DAMR_DOMAIN_HEIGHT): square blocks of nx x ny cells, base blocks of
!   L0 x L0 lattice units at the finest level, and a uniform background
!   refinement to minLevels. As constants they let the compiler specialise
!   every block kernel for the block size. The domain is nBaseX x nBaseY
!   base blocks (L0 x nBaseY*L0, by default L0 x 3*L0), periodic in x,
!   with walls at the bottom and top, or periodic in y as well
!   (-DAMR_PERIODIC_Y=ON).
!================================================
MODULE Config_Mod
    USE Kinds_Mod, ONLY: wp
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: nx, ny, L0, nBaseX, nBaseY, periodicY, maxLevels, minLevels, MaxBlks, dx_level, inv_dx_level

#if !defined(AMR_BLOCK_SIZE) || !defined(AMR_L0) || !defined(AMR_MIN_LEVEL) || !defined(AMR_DOMAIN_HEIGHT)
#error "AMR_BLOCK_SIZE, AMR_L0, AMR_MIN_LEVEL and AMR_DOMAIN_HEIGHT must be defined (see CMakeLists.txt)"
#endif
    INTEGER, PARAMETER :: nx = AMR_BLOCK_SIZE
    INTEGER, PARAMETER :: ny = AMR_BLOCK_SIZE
    INTEGER, PARAMETER :: L0 = AMR_L0
    INTEGER, PARAMETER :: maxLevels = BIT_SIZE(L0) - 1 - LEADZ(L0/nx)   ! log2(L0/nx)
    INTEGER, PARAMETER :: minLevels = AMR_MIN_LEVEL
    INTEGER, PARAMETER :: nBaseX = 1, nBaseY = AMR_DOMAIN_HEIGHT
#ifdef AMR_PERIODIC_Y
    LOGICAL, PARAMETER :: periodicY = .TRUE.
#else
    LOGICAL, PARAMETER :: periodicY = .FALSE.
#endif

    !== number of blocks of the fully refined quadtrees (4**L blocks below
    !   each base block on level L): an upper bound of the live blocks
    INTEGER, PARAMETER :: MaxBlks = nBaseX*nBaseY * (4**(maxLevels+1) - 1) / 3

    !== cell size dx = (L0/nx)/2**L of level L in finest-level lattice units,
    !   and 1/dx. dx is a power of two, so the reciprocals are exact. Kernels
    !   take 1/dx instead of the integer cell size (no conversions, no divisions).
    INTEGER :: lev      ! implied-do variable of the tables
    INTEGER,  PARAMETER :: dx_level(0:maxLevels) = [( (L0/nx) / 2**lev, lev = 0, maxLevels )]
    REAL(wp), PARAMETER :: inv_dx_level(0:maxLevels) = &
        [( REAL(nx,wp) / REAL(L0,wp) * 2._wp**lev, lev = 0, maxLevels )]

END MODULE Config_Mod
