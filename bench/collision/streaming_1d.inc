!== streaming with the planes of f as rank-1 arrays (sequence association),
!   one temporary per direction. The shift of a direction is the offset
!   off = ex + (nx+2) ey in the plane.
!     ST_FLAT not defined: nested loops over the interior, index linearized by hand
!     ST_FLAT defined:     one loop from cell (1,1) to cell (nx,ny) of the plane,
!                          ghost columns in between included (they are overwritten)
!     ST_LW defined: Lax-Wendroff, else perfect shift
SUBROUTINE ST_V( n )
    INTEGER, INTENT(IN) :: n
    CALL ST_K( st_cfl, Blk_h(:,:,:,n) )
    CALL ST_K( st_cfl, Blk_g(:,:,:,n) )
END SUBROUTINE ST_V

SUBROUTINE ST_K( cfl, f )
    REAL(wp), INTENT(IN)    :: cfl
    REAL(wp), INTENT(INOUT) :: f((nx+2)*(ny+2),9)
    INTEGER, PARAMETER :: ld = nx+2, m1 = ld + 2, m2 = ld*ny + nx + 1     ! cells (1,1) and (nx,ny)
    INTEGER :: i, X, Y, m, off
    REAL(wp) :: fnew(m1:m2), cU, c0, cD

    cU = 0.5_wp*cfl*(1._wp+cfl)
    c0 = 1._wp - cfl*cfl
    cD = 0.5_wp*cfl*(1._wp-cfl)

    DO i = 2, 9
        off = ex(i) + ld*ey(i)
#ifdef ST_FLAT
        DO m = m1, m2
#else
        DO Y = 1, ny
        DO X = 1, nx
            m = X + 1 + ld*Y
#endif
#ifdef ST_LW
            fnew(m) = cU * f(m-off,i) + c0 * f(m,i) - cD * f(m+off,i)
#else
            fnew(m) = f(m-off,i)
#endif
#ifdef ST_FLAT
        END DO
        DO m = m1, m2
            f(m,i) = fnew(m)
        END DO
#else
        END DO
        END DO
        DO Y = 1, ny
        DO X = 1, nx
            m = X + 1 + ld*Y
            f(m,i) = fnew(m)
        END DO
        END DO
#endif
    END DO
END SUBROUTINE ST_K
#undef ST_V
#undef ST_K
#undef ST_LW
#undef ST_FLAT
