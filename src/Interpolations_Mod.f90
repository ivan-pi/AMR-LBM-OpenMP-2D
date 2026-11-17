!================================================
!   Biquadratic/bilinear interpolation between coarse and fine levels:
!   strips for the ghost cells (Biquadratic_*), whole blocks for the
!   refinement and derefinement (Bilinear_*)
!================================================
MODULE Interpolations_Mod
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: nx, ny     ! block size, for the whole-block prolongation
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Biquadratic_from_Fine, Biquadratic_from_Fine_f, Biquadratic_from_Coarse, &
              Biquadratic_from_Coarse_f, Bilinear_from_Fine, Bilinear_from_Fine_f, &
              Bilinear_from_Coarse, Bilinear_from_Coarse_f

CONTAINS


SUBROUTINE Biquadratic_from_Fine( i2, j2, X2, Y2, Cin_F, Cout_C )
    !== single field: the multi-component version with one component
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: i2, j2, X2, Y2
    REAL(wp), INTENT(IN ) :: Cin_F(i2,j2)
    REAL(wp), INTENT(OUT) :: Cout_C(X2,Y2)

    CALL Biquadratic_from_Fine_f( i2, j2, X2, Y2, Cin_F, Cout_C, 1 )

END

!******************************************************************************

SUBROUTINE Biquadratic_from_Fine_f( i2, j2, X2, Y2, fin_F, fout_C, nk )
    !== nk components (directions) stored last

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: i2, j2, X2, Y2, nk
    REAL(wp), INTENT(IN ) :: fin_F (i2,j2,nk)
    REAL(wp), INTENT(OUT) :: fout_C(X2,Y2,nk)

    INTEGER :: X, J, i, a, b, c
    REAL(wp) :: tmp(i2,Y2,nk)

    b = 6
    DO X = 1, Y2
        J = MOD(X+1,2)
        i = 2*X-1 - J
        a = 3 - 4*J
        c = 2 - a

        tmp(:,X,:) = a * fin_F(:,i,:) + b * fin_F(:,i+1,:) + c * fin_F(:,i+2,:)

    END DO
    DO X = 1, X2
        J = MOD(X+1,2)
        i = 2*X-1 - J
        a = 3 - 4*J
        c = 2 - a

        fout_C(X,:,:) = a * tmp(i,:,:) + b * tmp(i+1,:,:) + c * tmp(i+2,:,:)

    END DO
    fout_C = fout_C / 64._wp

END

!******************************************************************************

SUBROUTINE Biquadratic_from_Coarse( I2, J2, x2, y2, Cin_C, Cout_F )
    !== single field: the multi-component version with one component
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: I2, J2, x2, y2
    REAL(wp), INTENT(IN ) :: Cin_C(I2,J2)
    REAL(wp), INTENT(OUT) :: Cout_F(x2,y2)

    CALL Biquadratic_from_Coarse_f( I2, J2, x2, y2, Cin_C, Cout_F, 1 )

END

!******************************************************************************

SUBROUTINE Biquadratic_from_Coarse_f( I2, J2, x2, y2, fin_C, fout_F, nk )
    !== nk components (directions) stored last

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: I2, J2, x2, y2, nk
    REAL(wp), INTENT(IN ) :: fin_C (I2,J2,nk)
    REAL(wp), INTENT(OUT) :: fout_F(x2,y2,nk)

    INTEGER :: x, I, a, b, c
    REAL(wp) :: tmp(I2,y2,nk)

    b = 30
    DO x = 1, y2
        I = (x+1)/2
        a = 5 - 8*MOD(x+1,2)
        c = 2-a

        tmp(:,x,:) = a * fin_C(:,I,:) + b * fin_C(:,I+1,:) + c * fin_C(:,I+2,:)

    END DO
    DO x = 1, x2
        I = (x+1)/2
        a = 5 - 8*MOD(x+1,2)
        c = 2-a

        fout_F(x,:,:) = a * tmp(I,:,:) + b * tmp(I+1,:,:) + c * tmp(I+2,:,:)

    END DO
    fout_F = fout_F / 1024._wp

END

!******************************************************************************

SUBROUTINE Bilinear_from_Fine( i2, j2, X2, Y2, Cin_F, Cout_C )
    !== single field: the multi-component version with one component
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: i2, j2, X2, Y2
    REAL(wp), INTENT(IN ) :: Cin_F(i2,j2)
    REAL(wp), INTENT(OUT) :: Cout_C(X2,Y2)

    CALL Bilinear_from_Fine_f( i2, j2, X2, Y2, Cin_F, Cout_C, 1 )

END

!******************************************************************************

SUBROUTINE Bilinear_from_Fine_f( i2, j2, X2, Y2, hin_F, hout_C, nk )
    !== nk components (directions) stored last

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: i2, j2, X2, Y2, nk
    REAL(wp), INTENT(IN ) :: hin_F (i2,j2,nk)
    REAL(wp), INTENT(OUT) :: hout_C(X2,Y2,nk)

    INTEGER :: X, i
    REAL(wp) :: tmp(i2,Y2,nk)

    DO X = 1, Y2
        i = 2*X-1

        tmp(:,X,:) = 0.5_wp * ( hin_F(:,i,:) + hin_F(:,i+1,:) )

    END DO
    DO X = 1, X2
        i = 2*X-1

        hout_C(X,:,:) = 0.5_wp * ( tmp(i,:,:) + tmp(i+1,:,:) )

    END DO

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
