!================================================
!   Lattice Boltzmann kernels shared by all models: streaming of the
!   distribution functions, the phase-field moment C = sum(h), and the
!   isotropic finite differences (gradient, chemical potential) of C
!================================================
MODULE LBM_Common_Mod
    USE Kinds_Mod, ONLY: wp
    USE D2Q9_Mod, ONLY: ex, ey, Beta, kappa
    USE Parameters_Mod, ONLY: nx, ny     ! the fixed block size
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Perfect_Shift_AMR, Streaming_LaxWendroff_AMR, Macroscopic_Values_h_AMR, &
              Chemical_Potential_AMR, Gradient_AMR

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

!**********************************************************************

SUBROUTINE Chemical_Potential_AMR( dx, C, mu )
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: dx
    REAL(wp), INTENT(IN)  :: C(-1:nx+2,-1:ny+2)
    REAL(wp), INTENT(OUT) :: mu(nx,ny)

    INTEGER :: i, j
    REAL(wp) :: D2C(nx,ny)

    CALL Laplacian_AMR( C, D2C )

    DO j = 1, ny
    DO i = 1, nx
        mu(i,j) = 4._wp*Beta*C(i,j)*(C(i,j)-1._wp)*(C(i,j)-0.5_wp) - kappa*D2C(i,j)/dx**2
    END DO
    END DO

CONTAINS

    SUBROUTINE Laplacian_AMR( C, D2C )
        REAL(wp), INTENT(IN)  :: C(-1:nx+2,-1:ny+2)
        REAL(wp), INTENT(OUT) :: D2C(nx,ny)

        INTEGER :: i, j

        DO j = 1, ny
        DO i = 1, nx

            D2C(i,j) = ( C(i-1,j-1)+C(i+1,j-1)+C(i-1,j+1)+C(i+1,j+1) &
                     +4*(C(i  ,j-1)+C(i-1,j  )+C(i+1,j  )+C(i  ,j+1)) - 20*C(i,j) )/6

        END DO
        END DO

    END SUBROUTINE Laplacian_AMR

END SUBROUTINE Chemical_Potential_AMR

!**********************************************************************

SUBROUTINE Gradient_AMR( dx, C, DcDx, DcDy )
    !== isotropic gradient of C, in units of the finest lattice
    INTEGER, INTENT(IN)  :: dx
    REAL(wp), INTENT(IN)  :: C(-1:nx+2,-1:ny+2)
    REAL(wp), INTENT(OUT) :: DcDx(nx,ny), DcDy(nx,ny)

    INTEGER :: i, j

    DO j = 1, ny
    DO i = 1, nx

        DcDx(i,j) = ( (C(i+1,j) - C(i-1,j))/3._wp  + &
                    (C(i+1,j+1) + C(i+1,j-1) - C(i-1,j+1) - C(i-1,j-1))/12._wp )/REAL(dx,wp)
        DcDy(i,j) = ( (C(i,j+1) - C(i,j-1))/3._wp + &
                    (C(i+1,j+1) + C(i-1,j+1) - C(i+1,j-1) - C(i-1,j-1))/12._wp )/REAL(dx,wp)

    END DO
    END DO

END SUBROUTINE Gradient_AMR

END MODULE LBM_Common_Mod
