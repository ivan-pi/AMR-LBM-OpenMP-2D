!================================================
!   Restriction of children data to the parent (derefinement)
!================================================
MODULE Restrict_Mod
    USE Kinds_Mod, ONLY: wp
    USE Interpolations_Mod, ONLY: Bilinear_from_Fine, Bilinear_from_Fine_f
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Restrict_Bilinear

CONTAINS

!******************************************************************************

SUBROUTINE Restrict_Bilinear( n )
    USE Parameters_Mod, ONLY: Blk, nx, ny
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

END MODULE Restrict_Mod
