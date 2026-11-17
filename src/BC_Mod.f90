!================================================
!   Boundary conditions (bounce-back walls, symmetric phase field)
!================================================
MODULE BC_Mod
    USE Kinds_Mod, ONLY: wp
    USE Parameters_Mod, ONLY: nx, ny     ! the fixed block size
    USE D2Q9_Mod, ONLY: ex, ey
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Bounce_Back_3, Bounce_Back_4, Symmetric_C_3, Symmetric_C_4

CONTAINS


SUBROUTINE Bounce_Back_3( f )
    IMPLICIT NONE
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)

    INTEGER :: X, Y, i, J, Xn

    Y = 0
    DO X = 0, nx+1
        f(X,Y,1) = f(X,Y+1,1)
        DO i = 2, 9
            Xn = MIN( MAX(0,X+ex(i)*ey(i)), nx+1 )
            J = (i-2) + 4*MOD(i/2,2)
            f(X,Y,i) = f(Xn,Y+1,J)
        END DO
    END DO

END

!******************************************************************************

SUBROUTINE Bounce_Back_4( f )
    IMPLICIT NONE
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)

    INTEGER :: X, Y, i, J, Xn

    Y = ny+1
    DO X = 0, nx+1
        f(X,Y,1) = f(X,Y-1,1)
        DO i = 2, 9
            Xn = MIN( MAX(0,X-ex(i)*ey(i)), nx+1 )
            J = (i-2) + 4*MOD(i/2,2)
            f(X,Y,i) = f(Xn,Y-1,J)
        END DO
    END DO

END

!******************************************************************************

SUBROUTINE Symmetric_C_3( ng, C )
    !== symmetry (zero-gradient) condition for C on the bottom wall, ng layers
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: ng
    REAL(wp), INTENT(INOUT) :: C(-1:nx+2,-1:ny+2)

    INTEGER :: l

    DO l = 1, ng
        C(1-ng:nx+ng,1-l) = C(1-ng:nx+ng,l)
    END DO

END

!******************************************************************************

SUBROUTINE Symmetric_C_4( ng, C )
    !== symmetry (zero-gradient) condition for C on the top wall, ng layers
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: ng
    REAL(wp), INTENT(INOUT) :: C(-1:nx+2,-1:ny+2)

    INTEGER :: l

    DO l = 1, ng
        C(1-ng:nx+ng,ny+l) = C(1-ng:nx+ng,ny+1-l)
    END DO

END

END MODULE BC_Mod
