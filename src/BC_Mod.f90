!================================================
!   Boundary conditions (bounce-back walls, symmetric phase field)
!
!   Each fills the ghost row(s) of a block at the bottom (3) or top (4)
!   wall in the columns X1:X2, from the first interior row(s) including
!   their side ghost cells. The ghost fills call them for the whole row once
!   the copies and restrictions are done; for C again in the end columns
!   next to a side interpolated from a coarser block (Ghost_C_Mod).
!================================================
MODULE BC_Mod
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: nx, ny     ! the fixed block size
    USE D2Q9_Mod, ONLY: ex, ey, opp
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Bounce_Back_3, Bounce_Back_4, Symmetric_C_3, Symmetric_C_4

CONTAINS


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
