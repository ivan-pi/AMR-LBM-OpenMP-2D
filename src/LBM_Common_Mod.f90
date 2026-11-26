!================================================
!   Lattice Boltzmann kernels shared by all models: streaming of the
!   distribution functions and the phase-field moment C = sum(h)
!================================================
MODULE LBM_Common_Mod
    USE Kinds_Mod, ONLY: wp
    USE D2Q9_Mod, ONLY: ex, ey
    USE Config_Mod, ONLY: nx, ny     ! the fixed block size
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Perfect_Shift_AMR, Streaming_LaxWendroff_AMR, Macroscopic_Values_h_AMR

CONTAINS

!== Both streaming kernels compute a direction into a temporary of one
!   nx x ny plane and copy it back before the next direction. (With one
!   temporary for all directions and a single array assignment at the end,
!   gfortran 15 vectorizes the copy across the directions with shuffles,
!   2.4 times slower; the plane also stays in L1.
!   20-180 % faster with gfortran, ifx and flang, same results.)
SUBROUTINE Perfect_Shift_AMR( f )
    IMPLICIT NONE
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)

    INTEGER :: i, X, Y
    REAL(wp) :: fnew(nx,ny)

    DO i = 2, 9
        DO Y = 1, ny
        DO X = 1, nx
            fnew(X,Y) = f(X-ex(i),Y-ey(i),i)
        END DO
        END DO
        DO Y = 1, ny
        DO X = 1, nx
            f(X,Y,i) = fnew(X,Y)
        END DO
        END DO
    END DO

END SUBROUTINE Perfect_Shift_AMR

!**********************************************************************

SUBROUTINE Streaming_LaxWendroff_AMR( cfl, f )
    IMPLICIT NONE
    REAL(wp), INTENT(IN)    :: cfl
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)

    INTEGER :: i, X, Y
    REAL(wp) :: fnew(nx,ny), cU, c0, cD

    !== Lax-Wendroff weights of the upwind, central and downwind values
    cU = 0.5_wp*cfl*(1._wp+cfl)
    c0 = 1._wp - cfl*cfl
    cD = 0.5_wp*cfl*(1._wp-cfl)

    DO i = 2, 9
        DO Y = 1, ny
        DO X = 1, nx
            fnew(X,Y) = cU * f(X-ex(i),Y-ey(i),i) &
                      + c0 * f(X,Y,i) &
                      - cD * f(X+ex(i),Y+ey(i),i)
        END DO
        END DO
        DO Y = 1, ny
        DO X = 1, nx
            f(X,Y,i) = fnew(X,Y)
        END DO
        END DO
    END DO

END SUBROUTINE Streaming_LaxWendroff_AMR

!**********************************************************************

SUBROUTINE Macroscopic_Values_h_AMR( h, C )
    IMPLICIT NONE
    REAL(wp), INTENT(IN)    :: h(0:nx+1,0:ny+1,9)
    REAL(wp), INTENT(INOUT) :: C(  -1:nx+2,-1:ny+2)

    INTEGER :: k

    C(1:nx,1:ny) = h(1:nx,1:ny,1)
    DO k = 2, 9
        C(1:nx,1:ny) = C(1:nx,1:ny) + h(1:nx,1:ny,k)
    END DO

END SUBROUTINE Macroscopic_Values_h_AMR

END MODULE LBM_Common_Mod
