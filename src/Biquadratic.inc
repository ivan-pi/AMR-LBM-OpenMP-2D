!================================================
!   Template: biquadratic interpolation of a strip of nk components (stored
!   last) from a fine to a coarse level and from a coarse to a fine level.
!   The including module sets the two routine names and the sizes, either
!   as dummy arguments (any strip) or as PARAMETERs (one fixed strip shape,
!   which lets the compiler unroll the short loops); see Interpolations_Mod
!   and Ghost_hg_Mod. Input i2 x j2 values per component, output x2 x y2.
!   With COARSE_ONLY defined only the coarse-to-fine routine is generated.
!   The macros are undefined at the end, ready for the next instantiation.
!================================================

#ifndef COARSE_ONLY
SUBROUTINE FINE_NAME( SIZE_ARGS fin_F, fout_C )
    IMPLICIT NONE
    FINE_SIZES
    REAL(wp), INTENT(IN ) :: fin_F (i2,j2,nk)
    REAL(wp), INTENT(OUT) :: fout_C(x2,y2,nk)

    INTEGER :: X, J, i, a, b, c
    REAL(wp) :: tmp(i2,y2,nk), wa, wb, wc

    !== weights (a, b, c)/8 per direction; the second pass carries both
    !   factors 1/8 (1/64 is a power of two, so this equals dividing at the end)
    b = 6
    DO X = 1, y2
        J = MOD(X+1,2)
        i = 2*X-1 - J
        a = 3 - 4*J
        c = 2 - a

        tmp(:,X,:) = a * fin_F(:,i,:) + b * fin_F(:,i+1,:) + c * fin_F(:,i+2,:)

    END DO
    wb = b / 64._wp
    DO X = 1, x2
        J = MOD(X+1,2)
        i = 2*X-1 - J
        a = 3 - 4*J
        c = 2 - a
        wa = a / 64._wp
        wc = c / 64._wp

        fout_C(X,:,:) = wa * tmp(i,:,:) + wb * tmp(i+1,:,:) + wc * tmp(i+2,:,:)

    END DO

END SUBROUTINE FINE_NAME

!******************************************************************************
#endif

SUBROUTINE COARSE_NAME( SIZE_ARGS fin_C, fout_F )
    IMPLICIT NONE
    COARSE_SIZES
    REAL(wp), INTENT(IN ) :: fin_C (i2,j2,nk)
    REAL(wp), INTENT(OUT) :: fout_F(x2,y2,nk)

    INTEGER :: x, I, a, b, c
    REAL(wp) :: tmp(i2,y2,nk), wa, wb, wc

    !== weights (a, b, c)/32 per direction; the second pass carries both
    !   factors 1/32 (1/1024 is a power of two, so this equals dividing at the end)
    b = 30
    DO x = 1, y2
        I = (x+1)/2
        a = 5 - 8*MOD(x+1,2)
        c = 2-a

        tmp(:,x,:) = a * fin_C(:,I,:) + b * fin_C(:,I+1,:) + c * fin_C(:,I+2,:)

    END DO
    wb = b / 1024._wp
    DO x = 1, x2
        I = (x+1)/2
        a = 5 - 8*MOD(x+1,2)
        c = 2-a
        wa = a / 1024._wp
        wc = c / 1024._wp

        fout_F(x,:,:) = wa * tmp(I,:,:) + wb * tmp(I+1,:,:) + wc * tmp(I+2,:,:)

    END DO

END SUBROUTINE COARSE_NAME

#undef FINE_NAME
#undef COARSE_NAME
#undef SIZE_ARGS
#undef FINE_SIZES
#undef COARSE_SIZES
#undef COARSE_ONLY
