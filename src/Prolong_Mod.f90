!================================================
!   Prolongation of block data to new children (refinement)
!================================================
MODULE Prolong_Mod
    USE Kinds_Mod, ONLY: wp
    USE Model_Mod, ONLY: Block_Kernel
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Prolong_Bilinear

CONTAINS

!******************************************************************************

!== interpolate block n to its new children c..c+3; Moments_g (the model's
!   kernel) derives their pressure and velocity from the interpolated g
SUBROUTINE Prolong_Bilinear( n, c, Moments_g )
    USE Parameters_Mod, ONLY: Blk, nx, ny
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, c
    PROCEDURE(Block_Kernel) :: Moments_g

    REAL(wp) :: hnew_F(2*nx,2*ny,9), Hnew_C(nx,ny,9)
    REAL(wp) :: gnew_F(2*nx,2*ny,9), Gnew_C(nx,ny,9)
    REAL(wp) :: C_F(0:2*nx+1,0:2*ny+1), C_C(0:nx+1,0:ny+1)

    Hnew_C = Blk(n)%h(1:nx,1:ny,:)
    Gnew_C = Blk(n)%g(1:nx,1:ny,:)
    C_C    = Blk(n)%C(0:nx+1,0:ny+1)

    DEALLOCATE( Blk(n)%h, Blk(n)%g, Blk(n)%C, Blk(n)%P, Blk(n)%Ux, Blk(n)%Uy )

    CALL Bilinear_Interpolate_Extrapolate( Hnew_C, hnew_F )
    CALL Bilinear_Interpolate_Extrapolate( Gnew_C, gnew_F )
    CALL Bilinear_Interpolation          ( C_C,    C_F    )

!== now fill children data
    CALL Fill_Children_Data( c  , 0 , 0 , gnew_F, hnew_F, C_F )
    CALL Fill_Children_Data( c+1, nx, 0 , gnew_F, hnew_F, C_F )
    CALL Fill_Children_Data( c+2, 0 , ny, gnew_F, hnew_F, C_F )
    CALL Fill_Children_Data( c+3, nx, ny, gnew_F, hnew_F, C_F )

CONTAINS

    SUBROUTINE Bilinear_Interpolate_Extrapolate( fin_C, fout_F )
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

    END SUBROUTINE Bilinear_Interpolate_Extrapolate

    SUBROUTINE Bilinear_Interpolation( Cin_C, Cout_F )
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

    END SUBROUTINE Bilinear_Interpolation

    SUBROUTINE Fill_Children_Data( c, I_idx, J_idx, gnew_F, hnew_F, C_F )
        USE Parameters_Mod, ONLY: Blk, nx, ny
        INTEGER, INTENT(IN) :: c, I_idx, J_idx
        REAL(wp), INTENT(IN) :: gnew_F(2*nx,2*ny,9)
        REAL(wp), INTENT(IN) :: hnew_F(2*nx,2*ny,9)
        REAL(wp), INTENT(IN) :: C_F(0:2*nx+1,0:2*ny+1)

        INTEGER :: I, J

        DO J = 1, ny
        DO I = 1, nx
            Blk(c)%g(I,J,:) = gnew_F(I+I_idx,J+J_idx,:)
            Blk(c)%h(I,J,:) = hnew_F(I+I_idx,J+J_idx,:)
        END DO
        END DO
        DO J = 0, ny+1
        DO I = 0, nx+1
            Blk(c)%C(I,J) = C_F(I+I_idx,J+J_idx)
        END DO
        END DO

        CALL Moments_g( c )

    END SUBROUTINE Fill_Children_Data

END SUBROUTINE Prolong_Bilinear


END MODULE Prolong_Mod
