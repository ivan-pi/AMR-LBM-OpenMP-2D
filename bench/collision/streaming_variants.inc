!**********************************************************************
!== Streaming of h and g of one block (kind 'stream', variants st_*), with
!   the CFL number st_cfl of a coarse level for the interpolating schemes:
!     st_shift    Perfect_Shift_AMR of the solver (finest level)
!     st_lw       Streaming_LaxWendroff_AMR of the solver: 3 points along e_k
!     st_9pt_lw   general 9-point stencil, 9 run-time weights per direction; the
!                 weights are those of st_lw (six zeros), so the results must agree
!     st_9pt_tp   the same kernel with the tensor product of the 1-D Lax-Wendroff
!                 weights in x and y (9 non-zero weights for the diagonals)
!     st_9pt_ax3  9 points for the diagonal directions, 3 for the axis directions
!                 (where the tensor product has only 3 non-zero weights)
!     st_sep_tp   the tensor-product stencil as two 1-D passes (x, then y)
SUBROUTINE St_Weights( tensor )
    LOGICAL, INTENT(IN) :: tensor
    INTEGER :: k, a, b
    REAL(wp) :: w1(-1:1), wx(-1:1), wy(-1:1)

    w1 = [ 0.5_wp*st_cfl*(1._wp+st_cfl), 1._wp - st_cfl*st_cfl, -0.5_wp*st_cfl*(1._wp-st_cfl) ]   ! upwind, central, downwind
    st_w9 = 0
    DO k = 2, 9
        wx = [ 0._wp, 1._wp, 0._wp ]
        wy = [ 0._wp, 1._wp, 0._wp ]
        IF( ex(k) /= 0 ) wx = [( w1(-a*ex(k)), a = -1, 1 )]     ! the value at x - e is the upwind one
        IF( ey(k) /= 0 ) wy = [( w1(-b*ey(k)), b = -1, 1 )]
        st_wx(:,k) = wx
        st_wy(:,k) = wy
        IF( tensor ) THEN
            DO b = -1, 1
                st_w9(:,b,k) = wx*wy(b)
            END DO
        ELSE
            st_w9(-ex(k),-ey(k),k) = w1(-1)
            st_w9(0,0,k)           = w1(0)
            st_w9( ex(k), ey(k),k) = w1(1)
        END IF
    END DO
END SUBROUTINE St_Weights

SUBROUTINE V_st_shift( n )
    INTEGER, INTENT(IN) :: n
    CALL Perfect_Shift_AMR( Blk_h(:,:,:,n) )
    CALL Perfect_Shift_AMR( Blk_g(:,:,:,n) )
END SUBROUTINE V_st_shift

SUBROUTINE V_st_lw( n )
    INTEGER, INTENT(IN) :: n
    CALL Streaming_LaxWendroff_AMR( st_cfl, Blk_h(:,:,:,n) )
    CALL Streaming_LaxWendroff_AMR( st_cfl, Blk_g(:,:,:,n) )
END SUBROUTINE V_st_lw

SUBROUTINE V_st_9pt_lw( n )
    INTEGER, INTENT(IN) :: n
    IF( st_tensor /= 1 ) THEN
!$OMP CRITICAL
        IF( st_tensor /= 1 ) THEN
            CALL St_Weights( .FALSE. );  st_tensor = 1
        END IF
!$OMP END CRITICAL
    END IF
    CALL K_st_9pt( st_w9, Blk_h(:,:,:,n) )
    CALL K_st_9pt( st_w9, Blk_g(:,:,:,n) )
END SUBROUTINE V_st_9pt_lw

SUBROUTINE St_Use_Tensor
    IF( st_tensor /= 2 ) THEN
!$OMP CRITICAL
        IF( st_tensor /= 2 ) THEN
            CALL St_Weights( .TRUE. );  st_tensor = 2
        END IF
!$OMP END CRITICAL
    END IF
END SUBROUTINE St_Use_Tensor

SUBROUTINE V_st_9pt_tp( n )
    INTEGER, INTENT(IN) :: n
    CALL St_Use_Tensor
    CALL K_st_9pt( st_w9, Blk_h(:,:,:,n) )
    CALL K_st_9pt( st_w9, Blk_g(:,:,:,n) )
END SUBROUTINE V_st_9pt_tp

SUBROUTINE V_st_9pt_ax3( n )
    INTEGER, INTENT(IN) :: n
    CALL St_Use_Tensor
    CALL K_st_9pt_ax3( st_w9, Blk_h(:,:,:,n) )
    CALL K_st_9pt_ax3( st_w9, Blk_g(:,:,:,n) )
END SUBROUTINE V_st_9pt_ax3

SUBROUTINE V_st_sep_tp( n )
    INTEGER, INTENT(IN) :: n
    CALL St_Use_Tensor
    CALL K_st_sep( st_wx, st_wy, Blk_h(:,:,:,n) )
    CALL K_st_sep( st_wx, st_wy, Blk_g(:,:,:,n) )
END SUBROUTINE V_st_sep_tp

!== fnew(x,k) = sum over the 3x3 neighbours (a,b) of w(a,b,k) f(x+(a,b),k)
SUBROUTINE K_st_9pt( w, f )
    REAL(wp), INTENT(IN)    :: w(-1:1,-1:1,2:9)
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)
    INTEGER :: k, X, Y
    REAL(wp) :: fnew(nx,ny,2:9)

    DO k = 2, 9
    DO Y = 1, ny
    DO X = 1, nx
        fnew(X,Y,k) = ( w(-1,-1,k)*f(X-1,Y-1,k) + w(0,-1,k)*f(X,Y-1,k) + w(1,-1,k)*f(X+1,Y-1,k) ) &
                    + ( w(-1, 0,k)*f(X-1,Y  ,k) + w(0, 0,k)*f(X,Y  ,k) + w(1, 0,k)*f(X+1,Y  ,k) ) &
                    + ( w(-1, 1,k)*f(X-1,Y+1,k) + w(0, 1,k)*f(X,Y+1,k) + w(1, 1,k)*f(X+1,Y+1,k) )
    END DO
    END DO
    END DO
    f(1:nx,1:ny,2:9) = fnew
END SUBROUTINE K_st_9pt

!== 3 points along the axis for directions 2-5, 9 points for the diagonals 6-9
SUBROUTINE K_st_9pt_ax3( w, f )
    REAL(wp), INTENT(IN)    :: w(-1:1,-1:1,2:9)
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)
    INTEGER :: k, X, Y
    REAL(wp) :: fnew(nx,ny,2:9)

    DO k = 2, 5
    DO Y = 1, ny
    DO X = 1, nx
        fnew(X,Y,k) = w(-ex(k),-ey(k),k)*f(X-ex(k),Y-ey(k),k) + w(0,0,k)*f(X,Y,k) + w(ex(k),ey(k),k)*f(X+ex(k),Y+ey(k),k)
    END DO
    END DO
    END DO
    DO k = 6, 9
    DO Y = 1, ny
    DO X = 1, nx
        fnew(X,Y,k) = ( w(-1,-1,k)*f(X-1,Y-1,k) + w(0,-1,k)*f(X,Y-1,k) + w(1,-1,k)*f(X+1,Y-1,k) ) &
                    + ( w(-1, 0,k)*f(X-1,Y  ,k) + w(0, 0,k)*f(X,Y  ,k) + w(1, 0,k)*f(X+1,Y  ,k) ) &
                    + ( w(-1, 1,k)*f(X-1,Y+1,k) + w(0, 1,k)*f(X,Y+1,k) + w(1, 1,k)*f(X+1,Y+1,k) )
    END DO
    END DO
    END DO
    f(1:nx,1:ny,2:9) = fnew
END SUBROUTINE K_st_9pt_ax3

!== tensor-product stencil as a pass in x (rows 0:ny+1) and a pass in y
SUBROUTINE K_st_sep( wx, wy, f )
    REAL(wp), INTENT(IN)    :: wx(-1:1,2:9), wy(-1:1,2:9)
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)
    INTEGER :: k, X, Y
    REAL(wp) :: tmp(nx,0:ny+1)

    DO k = 2, 9
        DO Y = 0, ny+1
        DO X = 1, nx
            tmp(X,Y) = wx(-1,k)*f(X-1,Y,k) + wx(0,k)*f(X,Y,k) + wx(1,k)*f(X+1,Y,k)
        END DO
        END DO
        DO Y = 1, ny
        DO X = 1, nx
            f(X,Y,k) = wy(-1,k)*tmp(X,Y-1) + wy(0,k)*tmp(X,Y) + wy(1,k)*tmp(X,Y+1)
        END DO
        END DO
    END DO
END SUBROUTINE K_st_sep

!**********************************************************************
!== the loop over the directions of the Lax-Wendroff kernel: as written,
!   unrolled by directive, unrolling forbidden by directive
#define ST_V V_st_lw_plain
#define ST_K K_st_lw_plain
#define ST_UNROLL 0
#include "streaming_lw.inc"
#define ST_V V_st_lw_unroll
#define ST_K K_st_lw_unroll
#define ST_UNROLL 1
#include "streaming_lw.inc"
#define ST_V V_st_lw_nounroll
#define ST_K K_st_lw_nounroll
#define ST_UNROLL 2
#include "streaming_lw.inc"

!**********************************************************************
!== one temporary per direction instead of fnew(nx,ny,2:9) and one copy at the end
#define ST_V V_st_shift_2d
#define ST_K K_st_shift_2d
#include "streaming_2d.inc"
#define ST_V V_st_lw_2d
#define ST_K K_st_lw_2d
#define ST_LW
#include "streaming_2d.inc"

!**********************************************************************
!== rank-1 planes: linearized index with nested loops (lin), one flat loop
!   over the plane including the ghost columns in between (flat)
#define ST_V V_st_shift_lin
#define ST_K K_st_shift_lin
#include "streaming_1d.inc"
#define ST_V V_st_lw_lin
#define ST_K K_st_lw_lin
#define ST_LW
#include "streaming_1d.inc"
#define ST_V V_st_shift_flat
#define ST_K K_st_shift_flat
#define ST_FLAT
#include "streaming_1d.inc"
#define ST_V V_st_lw_flat
#define ST_K K_st_lw_flat
#define ST_LW
#define ST_FLAT
#include "streaming_1d.inc"

!== perfect shift in place, no temporary: each plane is moved by its offset,
!   from the end for a positive offset (flat, ghost columns overwritten)
SUBROUTINE V_st_shift_move( n )
    INTEGER, INTENT(IN) :: n
    CALL K_st_shift_move( Blk_h(:,:,:,n) )
    CALL K_st_shift_move( Blk_g(:,:,:,n) )
END SUBROUTINE V_st_shift_move

SUBROUTINE K_st_shift_move( f )
    REAL(wp), INTENT(INOUT) :: f((nx+2)*(ny+2),9)
    INTEGER, PARAMETER :: ld = nx+2, m1 = ld + 2, m2 = ld*ny + nx + 1
    INTEGER :: i, m, off

    DO i = 2, 9
        off = ex(i) + ld*ey(i)
        IF( off > 0 ) THEN
            DO m = m2, m1, -1
                f(m,i) = f(m-off,i)
            END DO
        ELSE
            DO m = m1, m2
                f(m,i) = f(m-off,i)
            END DO
        END IF
    END DO
END SUBROUTINE K_st_shift_move
