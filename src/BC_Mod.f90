!================================================
!   Boundary conditions (bounce-back walls, symmetric phase field)
!
!   Each fills the ghost column(s) of a block at the left (1) or right (2)
!   wall in the rows i1:i2, or the ghost row(s) at the bottom (3) or top
!   (4) wall in the columns i1:i2, from the first interior column(s) or
!   row(s) including their ghost cells. The ghost fills call them for the
!   whole side once the copies and restrictions are done; for C again at
!   the ends next to a side interpolated from a coarser block
!   (Ghost_C_Mod). The routines of the four sides are mirror images: the
!   left and right ones are those of the bottom and top with x and y
!   swapped (ex*ey is symmetric in x and y).
!   Boundaries_Mod assigns them to the sides of the domain.
!================================================
MODULE BC_Mod
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: nx, ny     ! the fixed block size
    USE D2Q9_Mod, ONLY: ex, ey, opp
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Bounce_Back_1, Bounce_Back_2, Bounce_Back_3, Bounce_Back_4
    PUBLIC :: Symmetric_C_1, Symmetric_C_2, Symmetric_C_3, Symmetric_C_4

CONTAINS


SUBROUTINE Bounce_Back_1( f, Y1, Y2 )
    IMPLICIT NONE
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)
    INTEGER,  INTENT(IN)    :: Y1, Y2

    INTEGER :: X, Y, i, Yn

    X = 0
    DO Y = Y1, Y2
        f(X,Y,1) = f(X+1,Y,1)
        DO i = 2, 9
            Yn = MIN( MAX(0,Y+ex(i)*ey(i)), ny+1 )
            f(X,Y,i) = f(X+1,Yn,opp(i))
        END DO
    END DO

END

!******************************************************************************

SUBROUTINE Bounce_Back_2( f, Y1, Y2 )
    IMPLICIT NONE
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)
    INTEGER,  INTENT(IN)    :: Y1, Y2

    INTEGER :: X, Y, i, Yn

    X = nx+1
    DO Y = Y1, Y2
        f(X,Y,1) = f(X-1,Y,1)
        DO i = 2, 9
            Yn = MIN( MAX(0,Y-ex(i)*ey(i)), ny+1 )
            f(X,Y,i) = f(X-1,Yn,opp(i))
        END DO
    END DO

END

!******************************************************************************

SUBROUTINE Bounce_Back_3( f, X1, X2 )
    IMPLICIT NONE
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)
    INTEGER,  INTENT(IN)    :: X1, X2

    INTEGER :: X, Y, i, Xn

    Y = 0
    DO X = X1, X2
        f(X,Y,1) = f(X,Y+1,1)
        DO i = 2, 9
            Xn = MIN( MAX(0,X+ex(i)*ey(i)), nx+1 )
            f(X,Y,i) = f(Xn,Y+1,opp(i))
        END DO
    END DO

END

!******************************************************************************

SUBROUTINE Bounce_Back_4( f, X1, X2 )
    IMPLICIT NONE
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)
    INTEGER,  INTENT(IN)    :: X1, X2

    INTEGER :: X, Y, i, Xn

    Y = ny+1
    DO X = X1, X2
        f(X,Y,1) = f(X,Y-1,1)
        DO i = 2, 9
            Xn = MIN( MAX(0,X-ex(i)*ey(i)), nx+1 )
            f(X,Y,i) = f(Xn,Y-1,opp(i))
        END DO
    END DO

END

!******************************************************************************

SUBROUTINE Symmetric_C_1( ng, C, Y1, Y2 )
    !== symmetry (zero-gradient) condition for C on the left wall, ng layers
    IMPLICIT NONE
    INTEGER,  INTENT(IN)    :: ng, Y1, Y2
    REAL(wp), INTENT(INOUT) :: C(-1:nx+2,-1:ny+2)

    INTEGER :: l

    DO l = 1, ng
        C(1-l,Y1:Y2) = C(l,Y1:Y2)
    END DO

END

!******************************************************************************

SUBROUTINE Symmetric_C_2( ng, C, Y1, Y2 )
    !== symmetry (zero-gradient) condition for C on the right wall, ng layers
    IMPLICIT NONE
    INTEGER,  INTENT(IN)    :: ng, Y1, Y2
    REAL(wp), INTENT(INOUT) :: C(-1:nx+2,-1:ny+2)

    INTEGER :: l

    DO l = 1, ng
        C(nx+l,Y1:Y2) = C(nx+1-l,Y1:Y2)
    END DO

END

!******************************************************************************

SUBROUTINE Symmetric_C_3( ng, C, X1, X2 )
    !== symmetry (zero-gradient) condition for C on the bottom wall, ng layers
    IMPLICIT NONE
    INTEGER,  INTENT(IN)    :: ng, X1, X2
    REAL(wp), INTENT(INOUT) :: C(-1:nx+2,-1:ny+2)

    INTEGER :: l

    DO l = 1, ng
        C(X1:X2,1-l) = C(X1:X2,l)
    END DO

END

!******************************************************************************

SUBROUTINE Symmetric_C_4( ng, C, X1, X2 )
    !== symmetry (zero-gradient) condition for C on the top wall, ng layers
    IMPLICIT NONE
    INTEGER,  INTENT(IN)    :: ng, X1, X2
    REAL(wp), INTENT(INOUT) :: C(-1:nx+2,-1:ny+2)

    INTEGER :: l

    DO l = 1, ng
        C(X1:X2,ny+l) = C(X1:X2,ny+1-l)
    END DO

END

END MODULE BC_Mod
