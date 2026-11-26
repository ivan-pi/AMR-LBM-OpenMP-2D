    SELECT CASE( k )
    CASE( 1 )
#undef GH_K
#define GH_K 1
#if GH_COPY == 3
        DO j = 1, ny
        DO q = 1, GH_ND
            d = GH_DIR
            hn(0,j,d) = hb(nx,j,d)
            gn(0,j,d) = gb(nx,j,d)
        END DO
        END DO
#else
        DO q = 1, GH_ND
            d = GH_DIR
            DO j = 1, ny
                hn(0,j,d) = hb(nx,j,d)
            END DO
            DO j = 1, ny
                gn(0,j,d) = gb(nx,j,d)
            END DO
        END DO
#endif
    CASE( 2 )
#undef GH_K
#define GH_K 2
#if GH_COPY == 3
        DO j = 1, ny
        DO q = 1, GH_ND
            d = GH_DIR
            hn(nx+1,j,d) = hb(1,j,d)
            gn(nx+1,j,d) = gb(1,j,d)
        END DO
        END DO
#else
        DO q = 1, GH_ND
            d = GH_DIR
            DO j = 1, ny
                hn(nx+1,j,d) = hb(1,j,d)
            END DO
            DO j = 1, ny
                gn(nx+1,j,d) = gb(1,j,d)
            END DO
        END DO
#endif
    CASE( 3 )
#undef GH_K
#define GH_K 3
        DO q = 1, GH_ND
            d = GH_DIR
            DO i = 1, nx
                hn(i,0,d) = hb(i,ny,d)
            END DO
            DO i = 1, nx
                gn(i,0,d) = gb(i,ny,d)
            END DO
        END DO
    CASE( 4 )
#undef GH_K
#define GH_K 4
        DO q = 1, GH_ND
            d = GH_DIR
            DO i = 1, nx
                hn(i,ny+1,d) = hb(i,1,d)
            END DO
            DO i = 1, nx
                gn(i,ny+1,d) = gb(i,1,d)
            END DO
        END DO
    CASE( 5 )
#undef GH_K
#define GH_K 5
        DO q = 1, GH_ND
            d = GH_DIR
            hn(0,0,d) = hb(nx,ny,d)
            gn(0,0,d) = gb(nx,ny,d)
        END DO
    CASE( 6 )
#undef GH_K
#define GH_K 6
        DO q = 1, GH_ND
            d = GH_DIR
            hn(0,ny+1,d) = hb(nx,1,d)
            gn(0,ny+1,d) = gb(nx,1,d)
        END DO
    CASE( 7 )
#undef GH_K
#define GH_K 7
        DO q = 1, GH_ND
            d = GH_DIR
            hn(nx+1,0,d) = hb(1,ny,d)
            gn(nx+1,0,d) = gb(1,ny,d)
        END DO
    CASE( 8 )
#undef GH_K
#define GH_K 8
        DO q = 1, GH_ND
            d = GH_DIR
            hn(nx+1,ny+1,d) = hb(1,1,d)
            gn(nx+1,ny+1,d) = gb(1,1,d)
        END DO
    END SELECT
#undef GH_K
#undef GH_ND
#undef GH_DIR
