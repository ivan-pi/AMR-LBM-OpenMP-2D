!== ghost cells of block n in direction k from its coarser neighbour nb, for the
!   nd directions dirs only: the strip is gathered per direction and interpolated
!   with the routine for that number of directions (shape fixed at compile time)
SUBROUTINE Fill_From_Coarse( n, nb, k, nd, dirs )
    INTEGER, INTENT(IN) :: n, nb, k, nd, dirs(nd)

    INTEGER  :: sx, sy, ig, jg, i0, j0, px, py, X, Y, l, m, i, j, q, d
    REAL(wp) :: hx(3,ny/2+2,6), gx(3,ny/2+2,6), hbx(1,ny,6), gbx(1,ny,6)
    REAL(wp) :: hy(nx/2+2,3,6), gy(nx/2+2,3,6), hby(nx,1,6), gby(nx,1,6)
    REAL(wp) :: hc(3,3,2), gc(3,3,2), hbc(1,1,2), gbc(1,1,2)

    sx = di(k)
    sy = dj(k)
    px = MOD(Blk(n)%I,2)
    py = MOD(Blk(n)%J,2)
    IF( k <= 2 ) THEN
        ig = MERGE( 0, nx+1, sx < 0 )
        i0 = MERGE( nx+1, 0, sx < 0 )
        Y  = (ny/2)*py
        DO q = 1, nd
            d = dirs(q)
            DO j = 1, ny/2+2
            DO l = 0, 2
                hx(l+1,j,q) = Blk_h(i0+sx*l,Y+j-1,d, nb)
                gx(l+1,j,q) = Blk_g(i0+sx*l,Y+j-1,d, nb)
            END DO
            END DO
        END DO
        IF( nd == 3 ) THEN
            CALL From_Coarse_LR3( hx, hbx );  CALL From_Coarse_LR3( gx, gbx )
        ELSE
            CALL From_Coarse_LR6( hx, hbx );  CALL From_Coarse_LR6( gx, gbx )
        END IF
        DO q = 1, nd
            d = dirs(q)
            DO j = 1, ny
                Blk_h(ig,j,d, n) = hbx(1,j,q)
                Blk_g(ig,j,d, n) = gbx(1,j,q)
            END DO
        END DO
    ELSEIF( k <= 4 ) THEN
        jg = MERGE( 0, ny+1, sy < 0 )
        j0 = MERGE( ny+1, 0, sy < 0 )
        X  = (nx/2)*px
        DO q = 1, nd
            d = dirs(q)
            DO l = 0, 2
            DO i = 1, nx/2+2
                hy(i,l+1,q) = Blk_h(X+i-1,j0+sy*l,d, nb)
                gy(i,l+1,q) = Blk_g(X+i-1,j0+sy*l,d, nb)
            END DO
            END DO
        END DO
        IF( nd == 3 ) THEN
            CALL From_Coarse_BT3( hy, hby );  CALL From_Coarse_BT3( gy, gby )
        ELSE
            CALL From_Coarse_BT6( hy, hby );  CALL From_Coarse_BT6( gy, gby )
        END IF
        DO q = 1, nd
            d = dirs(q)
            DO i = 1, nx
                Blk_h(i,jg,d, n) = hby(i,1,q)
                Blk_g(i,jg,d, n) = gby(i,1,q)
            END DO
        END DO
    ELSE
        ig = MERGE( 0, nx+1, sx < 0 )
        jg = MERGE( 0, ny+1, sy < 0 )
        i0 = MERGE( nx+1-(nx/2)*px, (nx/2)*(1-px), sx < 0 )
        j0 = MERGE( ny+1-(ny/2)*py, (ny/2)*(1-py), sy < 0 )
        DO q = 1, nd
            d = dirs(q)
            DO m = 0, 2
            DO l = 0, 2
                hc(l+1,m+1,q) = Blk_h(i0+sx*l,j0+sy*m,d, nb)
                gc(l+1,m+1,q) = Blk_g(i0+sx*l,j0+sy*m,d, nb)
            END DO
            END DO
        END DO
        IF( nd == 1 ) THEN
            CALL From_Coarse_Corner1( hc, hbc );  CALL From_Coarse_Corner1( gc, gbc )
        ELSE
            CALL From_Coarse_Corner2( hc, hbc );  CALL From_Coarse_Corner2( gc, gbc )
        END IF
        DO q = 1, nd
            Blk_h(ig,jg,dirs(q), n) = hbc(1,1,q)
            Blk_g(ig,jg,dirs(q), n) = gbc(1,1,q)
        END DO
    END IF

END SUBROUTINE Fill_From_Coarse

#define FINE_NAME    Unused_Fine_LR3
#define COARSE_NAME  From_Coarse_LR3
#define SIZE_ARGS
#define FINE_SIZES   INTEGER, PARAMETER :: i2 = 3, j2 = 2*ny, x2 = 1, y2 = ny, nk = 3
#define COARSE_SIZES INTEGER, PARAMETER :: i2 = 3, j2 = ny/2+2, x2 = 1, y2 = ny, nk = 3
#include "Biquadratic.inc"
#define FINE_NAME    Unused_Fine_LR6
#define COARSE_NAME  From_Coarse_LR6
#define SIZE_ARGS
#define FINE_SIZES   INTEGER, PARAMETER :: i2 = 3, j2 = 2*ny, x2 = 1, y2 = ny, nk = 6
#define COARSE_SIZES INTEGER, PARAMETER :: i2 = 3, j2 = ny/2+2, x2 = 1, y2 = ny, nk = 6
#include "Biquadratic.inc"
#define FINE_NAME    Unused_Fine_BT3
#define COARSE_NAME  From_Coarse_BT3
#define SIZE_ARGS
#define FINE_SIZES   INTEGER, PARAMETER :: i2 = 2*nx, j2 = 3, x2 = nx, y2 = 1, nk = 3
#define COARSE_SIZES INTEGER, PARAMETER :: i2 = nx/2+2, j2 = 3, x2 = nx, y2 = 1, nk = 3
#include "Biquadratic.inc"
#define FINE_NAME    Unused_Fine_BT6
#define COARSE_NAME  From_Coarse_BT6
#define SIZE_ARGS
#define FINE_SIZES   INTEGER, PARAMETER :: i2 = 2*nx, j2 = 3, x2 = nx, y2 = 1, nk = 6
#define COARSE_SIZES INTEGER, PARAMETER :: i2 = nx/2+2, j2 = 3, x2 = nx, y2 = 1, nk = 6
#include "Biquadratic.inc"
#define FINE_NAME    Unused_Fine_Corner1
#define COARSE_NAME  From_Coarse_Corner1
#define SIZE_ARGS
#define FINE_SIZES   INTEGER, PARAMETER :: i2 = 3, j2 = 3, x2 = 1, y2 = 1, nk = 1
#define COARSE_SIZES INTEGER, PARAMETER :: i2 = 3, j2 = 3, x2 = 1, y2 = 1, nk = 1
#include "Biquadratic.inc"
#define FINE_NAME    Unused_Fine_Corner2
#define COARSE_NAME  From_Coarse_Corner2
#define SIZE_ARGS
#define FINE_SIZES   INTEGER, PARAMETER :: i2 = 3, j2 = 3, x2 = 1, y2 = 1, nk = 2
#define COARSE_SIZES INTEGER, PARAMETER :: i2 = 3, j2 = 3, x2 = 1, y2 = 1, nk = 2
#include "Biquadratic.inc"
