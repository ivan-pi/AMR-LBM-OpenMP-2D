!================================================
!   Kind parameters
!
!   wp : working precision of all floating-point fields and
!        computations of the solver: real64 by default, real32 when
!        configured with -DAMR_PRECISION=4 (defines WP_SINGLE)
!   dp : double precision, used for diagnostics that accumulate
!        many terms (e.g. the total mass) independent of wp.
!   i64: 64-bit integer, for clock counts and hashes
!================================================
MODULE Kinds_Mod
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: REAL32, REAL64, INT64
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: wp, dp, i64

#ifdef WP_SINGLE
    INTEGER, PARAMETER :: wp = REAL32
#else
    INTEGER, PARAMETER :: wp = REAL64
#endif
    INTEGER, PARAMETER :: dp = REAL64
    INTEGER, PARAMETER :: i64 = INT64

END MODULE Kinds_Mod
