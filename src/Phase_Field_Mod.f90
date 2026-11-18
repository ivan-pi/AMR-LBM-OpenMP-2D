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
    PUBLIC :: Chemical_Potential_AMR, Gradient_AMR, Max_Gradient, grad2_min

    !== |grad C|**2 at or below which a cell is treated as away from any interface:
    !   the interface normal is zero there (|grad C| <= 1e-16, i.e. round-off)
    REAL(wp), PARAMETER :: grad2_min = 1.e-32_wp

CONTAINS

!**********************************************************************

SUBROUTINE Chemical_Potential_AMR( inv_dx, C, mu )
    !== mu = 4 beta C (C-1) (C-1/2) - kappa Laplacian(C), isotropic 9-point Laplacian
    IMPLICIT NONE
    REAL(wp), INTENT(IN)  :: inv_dx         ! 1/dx in finest-level lattice units
    REAL(wp), INTENT(IN)  :: C(-1:nx+2,-1:ny+2)
    REAL(wp), INTENT(OUT) :: mu(nx,ny)

    INTEGER :: i, j
    REAL(wp) :: c4B, cK, S

    c4B = 4._wp*Beta
    cK  = kappa*inv_dx*inv_dx/6._wp         ! kappa/dx**2 and the 1/6 of the stencil

    DO j = 1, ny
    DO i = 1, nx
        S = ( C(i-1,j-1) + C(i+1,j-1) + C(i-1,j+1) + C(i+1,j+1) ) &
          + 4._wp*( C(i,j-1) + C(i-1,j) + C(i+1,j) + C(i,j+1) ) - 20._wp*C(i,j)
        mu(i,j) = c4B*C(i,j)*(C(i,j)-1._wp)*(C(i,j)-0.5_wp) - cK*S
    END DO
    END DO

END SUBROUTINE Chemical_Potential_AMR

!**********************************************************************

SUBROUTINE Gradient_AMR( inv_dx, C, DcDx, DcDy )
    !== isotropic gradient of C, in units of the finest lattice
    REAL(wp), INTENT(IN)  :: inv_dx         ! 1/dx in finest-level lattice units
    REAL(wp), INTENT(IN)  :: C(-1:nx+2,-1:ny+2)
    REAL(wp), INTENT(OUT) :: DcDx(nx,ny), DcDy(nx,ny)

    INTEGER :: i, j
    REAL(wp) :: c1, c2

    c1 = inv_dx/3._wp       ! weights of the straight and the diagonal neighbours
    c2 = inv_dx/12._wp

    DO j = 1, ny
    DO i = 1, nx
        DcDx(i,j) = c1*( C(i+1,j) - C(i-1,j) ) &
                  + c2*( C(i+1,j+1) + C(i+1,j-1) - C(i-1,j+1) - C(i-1,j-1) )
        DcDy(i,j) = c1*( C(i,j+1) - C(i,j-1) ) &
                  + c2*( C(i+1,j+1) + C(i-1,j+1) - C(i+1,j-1) - C(i-1,j-1) )
    END DO
    END DO

END SUBROUTINE Gradient_AMR

!**********************************************************************

!== largest |grad C| of a block (central differences), in units of the
!   finest lattice: the interface indicator of the refinement criteria.
!   The maximum is kept per cell of a row (one SIMD lane each) and reduced
!   at the end: a single running maximum chains every cell, and MAXVAL
!   carries the IEEE NaN handling, both slow
FUNCTION Max_Gradient( inv_dx, C ) RESULT( err )
    REAL(wp), INTENT(IN ) :: inv_dx         ! 1/dx in finest-level lattice units
    REAL(wp), INTENT(IN ) :: C(-1:nx+2,-1:ny+2)
    REAL(wp) :: err

    INTEGER  :: I, J
    REAL(wp) :: rowmax(nx)

    rowmax = 0
    DO J = 1, ny
    DO I = 1, nx
        rowmax(I) = MAX( rowmax(I), (C(I+1,J) - C(I-1,J))**2 + (C(I,J+1) - C(I,J-1))**2 )
    END DO
    END DO

    err = 0
    DO I = 1, nx
        err = MAX( err, rowmax(I) )
    END DO
    err = SQRT(err)*(0.5_wp*inv_dx)

END FUNCTION Max_Gradient

END MODULE Phase_Field_Mod
