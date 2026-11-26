!== the strip interpolations of C with the sizes at run time, as Interpolations_Mod had
!   them before the study (used by the variants without GC_FIXED)
SUBROUTINE Biquadratic_from_Fine( i2, j2, X2, Y2, Cin_F, Cout_C )
    !== single field: the multi-component version with one component
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: i2, j2, X2, Y2
    REAL(wp), INTENT(IN ) :: Cin_F(i2,j2)
    REAL(wp), INTENT(OUT) :: Cout_C(X2,Y2)

    CALL Biquadratic_from_Fine_f( i2, j2, X2, Y2, 1, Cin_F, Cout_C )

END SUBROUTINE

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

END SUBROUTINE

!******************************************************************************

