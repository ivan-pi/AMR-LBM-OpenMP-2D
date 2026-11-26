#if GH_COPY == 1
!== SELECT CASE( k ) outside; hn, gn: block n, hb, gb: its neighbour
SUBROUTINE Copy_Same_Level( n, nb, k, dirs )
    INTEGER, INTENT(IN) :: n, nb, k, dirs(:)
    CALL Copy_K( k, SIZE(dirs), dirs, Blk_h(:,:,:,n), Blk_h(:,:,:,nb), Blk_g(:,:,:,n), Blk_g(:,:,:,nb) )
END SUBROUTINE Copy_Same_Level

SUBROUTINE Copy_K( k, nd, dirs, hn, hb, gn, gb )
    INTEGER,  INTENT(IN)    :: k, nd, dirs(nd)
    REAL(wp), INTENT(INOUT) :: hn(0:nx+1,0:ny+1,9), gn(0:nx+1,0:ny+1,9)
    REAL(wp), INTENT(IN)    :: hb(0:nx+1,0:ny+1,9), gb(0:nx+1,0:ny+1,9)
    INTEGER :: q, d, i, j
#define GH_ND    nd
#define GH_DIR   dirs(q)
#include "ghost_copy_cases.inc"
END SUBROUTINE Copy_K
#else
!== one routine per list of directions: constants in every branch of SELECT CASE( k )
SUBROUTINE Copy_All( k, hn, hb, gn, gb )
    INTEGER,  INTENT(IN)    :: k
    REAL(wp), INTENT(INOUT) :: hn(0:nx+1,0:ny+1,9), gn(0:nx+1,0:ny+1,9)
    REAL(wp), INTENT(IN)    :: hb(0:nx+1,0:ny+1,9), gb(0:nx+1,0:ny+1,9)
    INTEGER :: q, d, i, j
#define GH_ND    8
#define GH_DIR   all_dirs(q)
#include "ghost_copy_cases.inc"
END SUBROUTINE Copy_All

SUBROUTINE Copy_Shift( k, hn, hb, gn, gb )
    INTEGER,  INTENT(IN)    :: k
    REAL(wp), INTENT(INOUT) :: hn(0:nx+1,0:ny+1,9), gn(0:nx+1,0:ny+1,9)
    REAL(wp), INTENT(IN)    :: hb(0:nx+1,0:ny+1,9), gb(0:nx+1,0:ny+1,9)
    INTEGER :: q, d, i, j
#define GH_ND    nShift(GH_K)
#define GH_DIR   shift_dirs(q,GH_K)
#include "ghost_copy_cases.inc"
END SUBROUTINE Copy_Shift

SUBROUTINE Copy_LW( k, hn, hb, gn, gb )
    INTEGER,  INTENT(IN)    :: k
    REAL(wp), INTENT(INOUT) :: hn(0:nx+1,0:ny+1,9), gn(0:nx+1,0:ny+1,9)
    REAL(wp), INTENT(IN)    :: hb(0:nx+1,0:ny+1,9), gb(0:nx+1,0:ny+1,9)
    INTEGER :: q, d, i, j
#define GH_ND    nLW(GH_K)
#define GH_DIR   lw_dirs(q,GH_K)
#include "ghost_copy_cases.inc"
END SUBROUTINE Copy_LW
#endif
