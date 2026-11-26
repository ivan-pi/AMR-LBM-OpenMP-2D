!================================================
!   Template: copy of the edge of a same-level block fb, gb into the ghost
!   cells hn, gn of its neighbour, for one list of directions. Ghost_hg_Mod
!   instantiates it for the three lists (all moving directions, those read by
!   the exact shift, those read by Lax-Wendroff):
!
!     COPY_NAME   name of the routine
!     COPY_ND     number of directions of the list for the ghost cells in
!     COPY_DIR    direction COPY_K (1..8), and its q-th direction
!
!   k is selected outside the loops and the blocks are explicit-shape dummy
!   arguments, so in every branch the directions and all sections are
!   constants (bench/ghost: 5-30 % faster than one routine with the list as
!   an argument and SELECT CASE inside the loop). Macros undefined at the end.
!================================================
SUBROUTINE COPY_NAME( k, hn, hb, gn, gb )
    IMPLICIT NONE
    INTEGER,  INTENT(IN)    :: k                    ! direction of the ghost cells, seen from hn
    REAL(wp), INTENT(INOUT) :: hn(0:nx+1,0:ny+1,9), gn(0:nx+1,0:ny+1,9)
    REAL(wp), INTENT(IN)    :: hb(0:nx+1,0:ny+1,9), gb(0:nx+1,0:ny+1,9)

    INTEGER :: q, d, i, j

    SELECT CASE( k )
    CASE( 1 )
#undef COPY_K
#define COPY_K 1
        DO q = 1, COPY_ND
            d = COPY_DIR
            DO j = 1, ny
                hn(0,j,d) = hb(nx,j,d)
            END DO
            DO j = 1, ny
                gn(0,j,d) = gb(nx,j,d)
            END DO
        END DO
    CASE( 2 )
#undef COPY_K
#define COPY_K 2
        DO q = 1, COPY_ND
            d = COPY_DIR
            DO j = 1, ny
                hn(nx+1,j,d) = hb(1,j,d)
            END DO
            DO j = 1, ny
                gn(nx+1,j,d) = gb(1,j,d)
            END DO
        END DO
    CASE( 3 )
#undef COPY_K
#define COPY_K 3
        DO q = 1, COPY_ND
            d = COPY_DIR
            DO i = 1, nx
                hn(i,0,d) = hb(i,ny,d)
            END DO
            DO i = 1, nx
                gn(i,0,d) = gb(i,ny,d)
            END DO
        END DO
    CASE( 4 )
#undef COPY_K
#define COPY_K 4
        DO q = 1, COPY_ND
            d = COPY_DIR
            DO i = 1, nx
                hn(i,ny+1,d) = hb(i,1,d)
            END DO
            DO i = 1, nx
                gn(i,ny+1,d) = gb(i,1,d)
            END DO
        END DO
    CASE( 5 )
#undef COPY_K
#define COPY_K 5
        DO q = 1, COPY_ND
            d = COPY_DIR
            hn(0,0,d) = hb(nx,ny,d)
            gn(0,0,d) = gb(nx,ny,d)
        END DO
    CASE( 6 )
#undef COPY_K
#define COPY_K 6
        DO q = 1, COPY_ND
            d = COPY_DIR
            hn(0,ny+1,d) = hb(nx,1,d)
            gn(0,ny+1,d) = gb(nx,1,d)
        END DO
    CASE( 7 )
#undef COPY_K
#define COPY_K 7
        DO q = 1, COPY_ND
            d = COPY_DIR
            hn(nx+1,0,d) = hb(1,ny,d)
            gn(nx+1,0,d) = gb(1,ny,d)
        END DO
    CASE( 8 )
#undef COPY_K
#define COPY_K 8
        DO q = 1, COPY_ND
            d = COPY_DIR
            hn(nx+1,ny+1,d) = hb(1,1,d)
            gn(nx+1,ny+1,d) = gb(1,1,d)
        END DO
    END SELECT

END SUBROUTINE COPY_NAME

#undef COPY_NAME
#undef COPY_ND
#undef COPY_DIR
#undef COPY_K
