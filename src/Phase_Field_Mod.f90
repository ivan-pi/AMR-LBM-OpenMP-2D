!================================================
!   Finite differences of the phase field C on one block (isotropic
!   stencils, C with ghost layers): chemical potential, gradient, and the
!   largest gradient used by the refinement criteria
!================================================
MODULE Phase_Field_Mod
    USE Kinds_Mod, ONLY: wp
    USE Physics_Mod, ONLY: Beta, kappa
    USE Config_Mod, ONLY: nx, ny     ! the fixed block size
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Chemical_Potential_AMR, Gradient_AMR, Max_Gradient

CONTAINS

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

!**********************************************************************

!== largest |grad C| of a block (central differences), in units of the
!   finest lattice: the interface indicator of the refinement criteria
FUNCTION Max_Gradient( dx, C ) RESULT( err )
    INTEGER, INTENT(IN ) :: dx
    REAL(wp), INTENT(IN ) :: C(-1:nx+2,-1:ny+2)
    REAL(wp) :: err

    INTEGER :: I, J
    REAL(wp) :: DC

    err = 0

    DO J = 1, ny
    DO I = 1, nx

        DC = (C(I+1,J) - C(I-1,J))**2 + (C(I,J+1) - C(I,J-1))**2

        err = MAX( err, DC )

    END DO
    END DO

    err = SQRT(err)/(2._wp*dx)

END FUNCTION Max_Gradient

END MODULE Phase_Field_Mod
