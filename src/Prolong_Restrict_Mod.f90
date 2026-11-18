!================================================
!   Transfer of block data between levels: prolongation to new children
!   (refinement) and restriction of children to their parent (derefinement)
!================================================
MODULE Prolong_Restrict_Mod
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk, Blk_h, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy, Clear_Slot
    USE Interpolations_Mod, ONLY: Bilinear_from_Coarse, Bilinear_from_Coarse_f
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Prolong_Bilinear, Restrict_Bilinear

CONTAINS

!******************************************************************************

!== interpolate h, g and C (with one ghost layer) of block n to its new
!   children c..c+3; their pressure and velocity are left to the model
SUBROUTINE Prolong_Bilinear( n, c )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, c

    INTEGER  :: q, ox, oy
    REAL(wp) :: hF(2*nx,2*ny,9), gF(2*nx,2*ny,9), CF(0:2*nx+1,0:2*ny+1)

    CALL Bilinear_from_Coarse_f( Blk_h(1:nx,1:ny,:,n), hF )
    CALL Bilinear_from_Coarse_f( Blk_g(1:nx,1:ny,:,n), gF )
    CALL Bilinear_from_Coarse  ( Blk_C(0:nx+1,0:ny+1,n), CF )

    !== child c+q covers the quadrant (MOD(q,2), q/2) of the parent
    DO q = 0, 3
        ox = nx * MOD(q,2)
        oy = ny * (q/2)
        Blk_h(1:nx,1:ny,:,c+q)       = hF(ox+1:ox+nx, oy+1:oy+ny, :)
        Blk_g(1:nx,1:ny,:,c+q)       = gF(ox+1:ox+nx, oy+1:oy+ny, :)
        Blk_C(0:nx+1,0:ny+1,c+q)     = CF(ox:ox+nx+1, oy:oy+ny+1)
    END DO

END SUBROUTINE Prolong_Bilinear

!******************************************************************************

!== the parent n of four leaves gets their 2x2 averages, child c+q into its
!   quadrant (MOD(q,2), q/2); its ghost cells are cleared
SUBROUTINE Restrict_Bilinear( n )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    INTEGER :: c, q, i1, i2, j1, j2

    c = Blk(n)%child
    CALL Clear_Slot( n )

    DO q = 0, 3
        i1 = 1 + (nx/2) * MOD(q,2)
        j1 = 1 + (ny/2) * (q/2)
        i2 = i1 + nx/2 - 1
        j2 = j1 + ny/2 - 1
        Blk_h(i1:i2,j1:j2,:,n)     = Average_2x2( Blk_h(1:nx,1:ny,:,c+q) )
        Blk_g(i1:i2,j1:j2,:,n)     = Average_2x2( Blk_g(1:nx,1:ny,:,c+q) )
        Blk_C(i1:i2,j1:j2,n:n)     = Average_2x2( Blk_C(1:nx,1:ny,c+q:c+q) )
        Blk_P(i1:i2,j1:j2,n:n)     = Average_2x2( Blk_P(:,:,c+q:c+q) )
        Blk_Ux(i1:i2,j1:j2,n:n)    = Average_2x2( Blk_Ux(:,:,c+q:c+q) )
        Blk_Uy(i1:i2,j1:j2,n:n)    = Average_2x2( Blk_Uy(:,:,c+q:c+q) )
    END DO

CONTAINS

    !== mean of each 2x2 group of cells (pairs in y first, then in x)
    PURE FUNCTION Average_2x2( f ) RESULT( fc )
        REAL(wp), INTENT(IN) :: f(:,:,:)
        REAL(wp) :: fc(SIZE(f,1)/2,SIZE(f,2)/2,SIZE(f,3))

        fc = 0.5_wp * ( 0.5_wp * ( f(1::2,1::2,:) + f(1::2,2::2,:) ) &
                      + 0.5_wp * ( f(2::2,1::2,:) + f(2::2,2::2,:) ) )

    END FUNCTION Average_2x2

END SUBROUTINE Restrict_Bilinear

END MODULE Prolong_Restrict_Mod
