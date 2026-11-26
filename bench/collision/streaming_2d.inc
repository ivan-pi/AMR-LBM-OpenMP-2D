!== streaming with one nx x ny temporary per direction, copied back inside
!   the loop over the directions. ST_LW defined: Lax-Wendroff, else perfect shift
SUBROUTINE ST_V( n )
    INTEGER, INTENT(IN) :: n
    CALL ST_K( st_cfl, Blk_h(:,:,:,n) )
    CALL ST_K( st_cfl, Blk_g(:,:,:,n) )
END SUBROUTINE ST_V

SUBROUTINE ST_K( cfl, f )
    REAL(wp), INTENT(IN)    :: cfl
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)
    INTEGER :: i, X, Y
    REAL(wp) :: fnew(nx,ny), cU, c0, cD

    cU = 0.5_wp*cfl*(1._wp+cfl)
    c0 = 1._wp - cfl*cfl
    cD = 0.5_wp*cfl*(1._wp-cfl)

    DO i = 2, 9
        DO Y = 1, ny
        DO X = 1, nx
#ifdef ST_LW
            fnew(X,Y) = cU * f(X-ex(i),Y-ey(i),i) + c0 * f(X,Y,i) - cD * f(X+ex(i),Y+ey(i),i)
#else
            fnew(X,Y) = f(X-ex(i),Y-ey(i),i)
#endif
        END DO
        END DO
        DO Y = 1, ny
        DO X = 1, nx
            f(X,Y,i) = fnew(X,Y)
        END DO
        END DO
    END DO
END SUBROUTINE ST_K
#undef ST_V
#undef ST_K
#undef ST_LW
