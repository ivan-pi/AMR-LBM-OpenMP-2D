SUBROUTINE MG_V( n )
    INTEGER, INTENT(IN) :: n
    REAL(wp) :: mu(nx,ny), inv_dx
    inv_dx = inv_dx_level( Blk(n)%Level )
#ifdef MG_MU_PASS
    CALL Chemical_Potential_AMR( inv_dx, Blk_C(:,:,n), mu )
#endif
    CALL MG_K( inv_dx, Blk_g(:,:,:,n), Blk_C(:,:,n), mu, Blk_P(:,:,n), Blk_Ux(:,:,n), Blk_Uy(:,:,n) )
END SUBROUTINE MG_V

SUBROUTINE MG_K( inv_dx, g, C, mu, P, Ux, Uy )
    REAL(wp), INTENT(IN)  :: inv_dx
    REAL(wp), INTENT(IN)  :: g(0:nx+1,0:ny+1,9), C(-1:nx+2,-1:ny+2), mu(nx,ny)
    REAL(wp), INTENT(OUT) :: P (nx,ny), Ux(nx,ny), Uy(nx,ny)

    INTEGER :: i, j
    REAL(wp) :: dRho, dRho6, c1, c2, c4B, cK, hGy
    REAL(wp) :: C0, S0, Sx, Sy, Rho, rRho, xmu, DcDx, DcDy, vx, vy, dpp, dpm, s24, s35, s68, s79, d68, d79

    dRho  = Rhoh - Rhol
    dRho6 = dRho*c16
    c1    = inv_dx/3._wp
    c2    = inv_dx/12._wp
    c4B   = 4._wp*Beta
    cK    = kappa*inv_dx*inv_dx/6._wp
    hGy   = 0.5_wp*Gy

    DO j = 1, ny
#ifndef MG_NOSIMD
!$OMP SIMD PRIVATE( C0, S0, Sx, Sy, Rho, rRho, xmu, DcDx, DcDy, vx, vy, dpp, dpm, s24, s35, s68, s79, d68, d79 )
#endif
    DO i = 1, nx
        C0 = C(i,j)
        !== isotropic gradient of C (as Gradient_AMR)
        dpp  = C(i+1,j+1) - C(i-1,j-1)
        dpm  = C(i+1,j-1) - C(i-1,j+1)
        DcDx = c1*( C(i+1,j) - C(i-1,j) ) + c2*( dpp + dpm )
        DcDy = c1*( C(i,j+1) - C(i,j-1) ) + c2*( dpp - dpm )
#ifdef MG_MU_PASS
        xmu = mu(i,j)
#else
        !== chemical potential (as Chemical_Potential_AMR)
        xmu = c4B*C0*(C0-1._wp)*(C0-0.5_wp) - cK*( ( C(i-1,j-1) + C(i+1,j-1) + C(i-1,j+1) + C(i+1,j+1) ) &
            + 4._wp*( C(i,j-1) + C(i-1,j) + C(i+1,j) + C(i,j+1) ) - 20._wp*C0 )
#endif
#ifdef MG_SYM
        !== sums and differences of the opposite pairs, shared by the three moments
        s24 = g(i,j,2) + g(i,j,4);  s35 = g(i,j,3) + g(i,j,5)
        s68 = g(i,j,6) + g(i,j,8);  s79 = g(i,j,7) + g(i,j,9)
        d68 = g(i,j,6) - g(i,j,8);  d79 = g(i,j,7) - g(i,j,9)
        S0 = g(i,j,1) + ( (s24 + s35) + (s68 + s79) )
        Sx = ( g(i,j,2) - g(i,j,4) ) + ( d68 - d79 )
        Sy = ( g(i,j,3) - g(i,j,5) ) + ( d68 + d79 )
#else
        S0 = g(i,j,1) + ( (g(i,j,2) + g(i,j,3)) + (g(i,j,4) + g(i,j,5)) ) &
                      + ( (g(i,j,6) + g(i,j,7)) + (g(i,j,8) + g(i,j,9)) )
        Sx = ( g(i,j,2) - g(i,j,4) ) + ( (g(i,j,6) - g(i,j,7)) + (g(i,j,9) - g(i,j,8)) )
        Sy = ( g(i,j,3) - g(i,j,5) ) + ( (g(i,j,6) + g(i,j,7)) - (g(i,j,8) + g(i,j,9)) )
#endif
        Rho  = Rhol + C0*dRho
#ifdef MG_NODIV
        rRho = 1._wp + 1.e-3_wp*Rho
#else
        rRho = 1._wp/Rho
#endif
        vx = ( 3._wp*Sx + 0.5_wp*( xmu*DcDx ) )*rRho
        vy = ( 3._wp*Sy + 0.5_wp*( xmu*DcDy ) + (Rhoh - Rho)*hGy )*rRho
        Ux(i,j) = vx
        Uy(i,j) = vy
        P (i,j) = S0 + dRho6*( vx*DcDx + vy*DcDy )
    END DO
    END DO
END SUBROUTINE MG_K
#undef MG_V
#undef MG_K
#undef MG_MU_PASS
#undef MG_NOSIMD
#undef MG_SYM
#undef MG_NODIV
