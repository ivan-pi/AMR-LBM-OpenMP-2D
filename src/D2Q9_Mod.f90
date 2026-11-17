!================================================
!   D2Q9 lattice: velocities and weights
!================================================
MODULE D2Q9_Mod
    USE Kinds_Mod, ONLY: wp
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: ex, ey, Wa

    INTEGER,PARAMETER :: ex(9) = [0, 1, 0,-1, 0, 1,-1,-1, 1]
    INTEGER,PARAMETER :: ey(9) = [0, 0, 1, 0,-1, 1, 1,-1,-1]
    REAL(wp),PARAMETER :: Wa(9) = [16,4, 4, 4, 4, 1, 1, 1, 1] / 36._wp

END MODULE D2Q9_Mod
