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

SUBROUTINE Perfect_Shift_AMR( f )
    IMPLICIT NONE
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)

    INTEGER :: i, X, Y
    REAL(wp) :: fnew(nx,ny,2:9)

    DO i = 2, 9
    DO Y = 1, ny
    DO X = 1, nx

        fnew(X,Y,i) = f(X-ex(i),Y-ey(i),i)

    END DO
    END DO
    END DO

    f(1:nx,1:ny,2:9) = fnew

END SUBROUTINE Perfect_Shift_AMR

!**********************************************************************

SUBROUTINE Streaming_LaxWendroff_AMR( cfl, f )
    IMPLICIT NONE
    REAL(wp), INTENT(IN)    :: cfl
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)

    INTEGER :: i, X, Y
    REAL(wp) :: fnew(nx,ny,2:9), cU, c0, cD

    !== Lax-Wendroff weights of the upwind, central and downwind values
    cU = cfl/2*(1+cfl)
    c0 = 1 - cfl**2
    cD = cfl/2*(1-cfl)

    DO i = 2, 9
    DO Y = 1, ny
    DO X = 1, nx

        fnew(X,Y,i) = cU * f(X-ex(i),Y-ey(i),i) &
                    + c0 * f(X,Y,i) &
                    - cD * f(X+ex(i),Y+ey(i),i)

    END DO
    END DO
    END DO

    f(1:nx,1:ny,2:9) = fnew

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
