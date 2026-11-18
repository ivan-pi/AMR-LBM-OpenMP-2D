!================================================
!   Biquadratic/bilinear interpolation between coarse and fine levels:
!   strips for the ghost cells (Biquadratic_*), whole blocks for the
!   refinement (Bilinear_from_Coarse*)
!================================================
MODULE Interpolations_Mod
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: nx, ny     ! block size, for the whole-block prolongation
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Biquadratic_from_Fine, Biquadratic_from_Coarse, &
              Bilinear_from_Coarse, Bilinear_from_Coarse_f

CONTAINS


SUBROUTINE Biquadratic_from_Fine( i2, j2, X2, Y2, Cin_F, Cout_C )
    !== single field: the multi-component version with one component
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: i2, j2, X2, Y2
    REAL(wp), INTENT(IN ) :: Cin_F(i2,j2)
    REAL(wp), INTENT(OUT) :: Cout_C(X2,Y2)

    CALL Biquadratic_from_Fine_f( i2, j2, X2, Y2, 1, Cin_F, Cout_C )

END

!******************************************************************************

!== the strip interpolations for any strip shape (sizes as arguments)
#define FINE_NAME    Biquadratic_from_Fine_f
#define COARSE_NAME  Biquadratic_from_Coarse_f
#define SIZE_ARGS    i2, j2, x2, y2, nk,
#define FINE_SIZES   INTEGER, INTENT(IN) :: i2, j2, x2, y2, nk
#define COARSE_SIZES INTEGER, INTENT(IN) :: i2, j2, x2, y2, nk
#include "Biquadratic.inc"

!******************************************************************************

SUBROUTINE Biquadratic_from_Coarse( I2, J2, x2, y2, Cin_C, Cout_F )
    !== single field: the multi-component version with one component
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: I2, J2, x2, y2
    REAL(wp), INTENT(IN ) :: Cin_C(I2,J2)
    REAL(wp), INTENT(OUT) :: Cout_F(x2,y2)

    CALL Biquadratic_from_Coarse_f( I2, J2, x2, y2, 1, Cin_C, Cout_F )

END

!******************************************************************************

SUBROUTINE Bilinear_from_Coarse_f( fin_C, fout_F )
    !== one block (9 directions) to the 2x2 children, bilinear, extrapolated
    !   at the block edges (no ghost cells)
    REAL(wp), INTENT(IN ) :: fin_C (nx,ny,9)
    REAL(wp), INTENT(OUT) :: fout_F(2*nx,2*ny,9)

    INTEGER :: i, X, a, b
    REAL(wp) :: tmp(2*nx,ny,9)

    DO i = 1, 2*nx
        X = MAX( 1, MIN(nx-1,i/2) )
        a = 2*i - 4*X + 1
        b = 4 - a

        tmp(i,:,:) = b * fin_C(X,:,:) + a * fin_C(X+1,:,:)

    END DO
    DO i = 1, 2*ny
        X = MAX( 1, MIN(ny-1,i/2) )
        a = 2*i - 4*X + 1
        b = 4 - a

        fout_F(:,i,:) = b * tmp(:,X,:) + a * tmp(:,X+1,:)

    END DO
    fout_F = fout_F / 16._wp

END SUBROUTINE Bilinear_from_Coarse_f

!******************************************************************************

SUBROUTINE Bilinear_from_Coarse( Cin_C, Cout_F )
    !== one block of C with one ghost layer to the 2x2 children, bilinear
    REAL(wp), INTENT(IN ) :: Cin_C (0:  nx+1,0:  ny+1)
    REAL(wp), INTENT(OUT) :: Cout_F(0:2*nx+1,0:2*ny+1)

    INTEGER :: i, X, a, b
    REAL(wp) :: tmp(0:2*nx+1,0:ny+1)

    DO i = 0, 2*nx+1
        X = i/2
        a = 2*i - 4*X + 1
        b = 4 - a

        tmp(i,:) = b * Cin_C(X,:) + a * Cin_C(X+1,:)

    END DO
    DO i = 0, 2*ny+1
        X = i/2
        a = 2*i - 4*X + 1
        b = 4 - a

        Cout_F(:,i) = b * tmp(:,X) + a * tmp(:,X+1)

    END DO
    Cout_F = Cout_F / 16._wp

END SUBROUTINE Bilinear_from_Coarse

END MODULE Interpolations_Mod
