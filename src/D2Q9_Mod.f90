!================================================
!   D2Q9 lattice: velocities and weights
!================================================
MODULE D2Q9_Mod
    USE Kinds_Mod, ONLY: wp
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: ex, ey, opp, Wa, exr, eyr, w0, ws, wd

    INTEGER,PARAMETER :: ex(9) = [0, 1, 0,-1, 0, 1,-1,-1, 1]
    INTEGER,PARAMETER :: ey(9) = [0, 0, 1, 0,-1, 1, 1,-1,-1]
    INTEGER,PARAMETER :: opp(9) = [1, 4, 5, 2, 3, 8, 9, 6, 7]    ! the opposite direction
    REAL(wp),PARAMETER :: Wa(9) = [16,4, 4, 4, 4, 1, 1, 1, 1] / 36._wp

    !== the directions as reals (for arithmetic: no integer conversions in loops)
    REAL(wp),PARAMETER :: exr(9) = REAL(ex,wp)
    REAL(wp),PARAMETER :: eyr(9) = REAL(ey,wp)
    !== weights of the rest, straight and diagonal directions
    REAL(wp),PARAMETER :: w0 = 4._wp/9._wp, ws = 1._wp/9._wp, wd = 1._wp/36._wp

END MODULE D2Q9_Mod
