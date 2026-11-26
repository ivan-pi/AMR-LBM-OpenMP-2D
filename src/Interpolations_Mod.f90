!================================================
!   Bilinear interpolation of whole blocks from a coarse to a fine level,
!   for the refinement (Bilinear_from_Coarse*). (The biquadratic strip
!   interpolations of the ghost cells are instantiated from Biquadratic.inc
!   in Ghost_hg_Mod and Ghost_C_Mod.)
!================================================
MODULE Interpolations_Mod
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: nx, ny     ! block size, for the whole-block prolongation
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Bilinear_from_Coarse, Bilinear_from_Coarse_f

CONTAINS


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
