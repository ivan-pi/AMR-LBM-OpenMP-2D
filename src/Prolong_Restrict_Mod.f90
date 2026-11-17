!================================================
!   Transfer of block data between levels: prolongation to new children
!   (refinement) and restriction of children to their parent (derefinement)
!================================================
MODULE Prolong_Restrict_Mod
    USE Kinds_Mod, ONLY: wp
    USE Model_Mod, ONLY: Block_Kernel
    USE Interpolations_Mod, ONLY: Bilinear_from_Coarse, Bilinear_from_Coarse_f, &
                                  Bilinear_from_Fine, Bilinear_from_Fine_f
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Prolong_Bilinear, Restrict_Bilinear

CONTAINS

!******************************************************************************

!== interpolate block n to its new children c..c+3; Moments_g (the model's
!   kernel) derives their pressure and velocity from the interpolated g
SUBROUTINE Prolong_Bilinear( n, c, Moments_g )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk
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

    CALL Bilinear_from_Coarse_f( Hnew_C, hnew_F )
    CALL Bilinear_from_Coarse_f( Gnew_C, gnew_F )
    CALL Bilinear_from_Coarse  ( C_C,    C_F    )

!== now fill children data
    CALL Fill_Children_Data( c  , 0 , 0 , gnew_F, hnew_F, C_F )
    CALL Fill_Children_Data( c+1, nx, 0 , gnew_F, hnew_F, C_F )
    CALL Fill_Children_Data( c+2, 0 , ny, gnew_F, hnew_F, C_F )
    CALL Fill_Children_Data( c+3, nx, ny, gnew_F, hnew_F, C_F )

CONTAINS

    SUBROUTINE Fill_Children_Data( c, I_idx, J_idx, gnew_F, hnew_F, C_F )
        USE Config_Mod, ONLY: nx, ny
        USE Blocks_Mod, ONLY: Blk
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

!******************************************************************************


!******************************************************************************

SUBROUTINE Restrict_Bilinear( n )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    INTEGER :: child, i, j, i1, i2, j1, j2

    REAL(wp), DIMENSION(2*nx,2*ny,9) :: h, g
    REAL(wp), DIMENSION(   2*nx,2*ny) :: C, P, Ux, Uy
    REAL(wp) :: hc(nx,ny,9), gc(nx,ny,9), Cc(nx,ny)     ! contiguous results for h, g, C

!== newly added (02/26/2016)
    ALLOCATE( Blk(n)%h(0:nx+1,0:ny+1,9), Blk(n)%g(0:nx+1,0:ny+1,9) )
    ALLOCATE( Blk(n)%C(  -1:nx+2,-1:ny+2) )
    ALLOCATE( Blk(n)%P(nx,ny), Blk(n)%Ux(nx,ny), Blk(n)%Uy(nx,ny) )
!-- New Correction! (Added 2016-08-30):
    Blk(n)%h  = 0
    Blk(n)%g  = 0
    Blk(n)%C  = 0
    Blk(n)%P  = 0
    Blk(n)%Ux = 0
    Blk(n)%Uy = 0
!== newly added (02/26/2016) - Correction (2016-08-30)

    child = Blk(n)%child(1)

    DO j = 0, 1
        j1 = 1  + ny*j
        j2 = ny + ny*j
        DO i = 0, 1
            i1 = 1  + nx*i
            i2 = nx + nx*i

            h (i1:i2,j1:j2,:) = Blk(child)%h (1:nx,1:ny,:)
            g (i1:i2,j1:j2,:) = Blk(child)%g (1:nx,1:ny,:)
            C (   i1:i2,j1:j2) = Blk(child)%C (   1:nx,1:ny)
            P (   i1:i2,j1:j2) = Blk(child)%P (   1:nx,1:ny)
            Ux(   i1:i2,j1:j2) = Blk(child)%Ux(   1:nx,1:ny)
            Uy(   i1:i2,j1:j2) = Blk(child)%Uy(   1:nx,1:ny)

            child = child + 1

        END DO
    END DO

    CALL Bilinear_from_Fine_f( 2*nx,2*ny, nx,ny, h , hc, 9 )
    CALL Bilinear_from_Fine_f( 2*nx,2*ny, nx,ny, g , gc, 9 )
    CALL Bilinear_from_Fine  ( 2*nx,2*ny, nx,ny, C , Cc )
    CALL Bilinear_from_Fine  ( 2*nx,2*ny, nx,ny, P , Blk(n)%P  )
    CALL Bilinear_from_Fine  ( 2*nx,2*ny, nx,ny, Ux, Blk(n)%Ux )
    CALL Bilinear_from_Fine  ( 2*nx,2*ny, nx,ny, Uy, Blk(n)%Uy )
    Blk(n)%h(1:nx,1:ny,:) = hc
    Blk(n)%g(1:nx,1:ny,:) = gc
    Blk(n)%C(1:nx,1:ny)   = Cc

END

END MODULE Prolong_Restrict_Mod
