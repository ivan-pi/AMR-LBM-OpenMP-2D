!== Streaming_LaxWendroff_AMR of the solver, with a directive on the loop over
!   the directions: ST_UNROLL 0 none, 1 unroll all 8, 2 do not unroll
SUBROUTINE ST_V( n )
    INTEGER, INTENT(IN) :: n
    CALL ST_K( st_cfl, Blk_h(:,:,:,n) )
    CALL ST_K( st_cfl, Blk_g(:,:,:,n) )
END SUBROUTINE ST_V

SUBROUTINE ST_K( cfl, f )
    REAL(wp), INTENT(IN)    :: cfl
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)
    INTEGER :: i, X, Y
    REAL(wp) :: fnew(nx,ny,2:9), cU, c0, cD

    cU = 0.5_wp*cfl*(1._wp+cfl)
    c0 = 1._wp - cfl*cfl
    cD = 0.5_wp*cfl*(1._wp-cfl)

#if ST_UNROLL == 1
!DIR$ UNROLL(8)
!GCC$ unroll 8
#elif ST_UNROLL == 2
!DIR$ NOUNROLL
!GCC$ unroll 0
#endif
    DO i = 2, 9
    DO Y = 1, ny
    DO X = 1, nx
        fnew(X,Y,i) = cU * f(X-ex(i),Y-ey(i),i) + c0 * f(X,Y,i) - cD * f(X+ex(i),Y+ey(i),i)
    END DO
    END DO
    END DO
    f(1:nx,1:ny,2:9) = fnew
END SUBROUTINE ST_K
#undef ST_V
#undef ST_K
#undef ST_UNROLL
