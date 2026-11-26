!================================================
!   Template: collision of h and g of the JCP 2016 model for one block, cell
!   by cell in a single pass. Model_JCP2016_Mod instantiates it twice:
!
!     KERNEL_NAME           name of the routine
!     WITH_MOMENTS defined  P, Ux, Uy are computed first from g (zeroth and
!                           first moments, with the forces) and stored;
!           not defined     P, Ux, Uy are input
!
!   The macros are undefined at the end, ready for the next instantiation.
!
!   Everything a cell needs is computed in the one loop body, the chemical
!   potential and (WITH_MOMENTS) the moments included: the differences of C
!   are shared, no work array is passed between loops, and h, g are loaded
!   once. (bench/collision: +30-65 % over three separate passes.)
!
!   The directions come in opposite pairs (2,4), (3,5), (6,8), (7,9),
!   e_k' = -e_k, with weight w and e.u = +-v. Per pair
!
!       Gamma+- = Gs +- Ga,   Gs = w (base0 + 4.5 v^2),   Ga = 3 w v
!       A+-     = Gamma+- amu - w dRho/3        (amu = mu + dRho/3)
!       f+-     = (+-ey - vy) Gamma+- Fb        (Fb: buoyancy)
!       geq+-   = (Pc - rho/3) w + Gamma+- (rho/3 - f-coefficient/2)
!                 - A+- (+-eC - u.grad C)/2
!
!   h: BGK relaxation to C (Gamma + sharpening term) with rate wc.
!
!   g: MRT g <- g - M^-1 S M (g - geq) with S = diag(1,...,1,s9,s9), plus the
!   forcing of the mixed differences eM. All moments but the two stresses
!   relax with rate 1, so
!
!       g+- = geq+- + A+- (eM+- - u.grad_m C) + f+-  +- m
!           = Pw + Gamma+- r+- + A+- (B+- + D0)      +- m
!
!   with  B+- = eM+- -+ eC/2 = (4 C(x+-e) - 3 C(x) - C(x+-2e))/(4 dx)
!         D0  = u.grad C/2 - u.grad_m C,   r+- = rho/3 + (+-ey - vy) Fb/2,
!         Pw  = (Pc - rho/3) w
!   and m = +-m8 for the straight pairs, +-m9 for the diagonal ones: (1-s9)/4
!   times the stress moments of g - geq,
!
!       m8 ~ (dg2 + dg4) - (dg3 + dg5),     m9 ~ (dg6 + dg8) - (dg7 + dg9).
!
!   Only the sums of geq over a pair enter them. In these the odd parts
!   cancel, and so do Pw and the constant part of Gs between the two pairs,
!   which leaves closed expressions in u, grad C and the buoyancy (below).
!   The individual geq_k and g_k - geq_k are never formed: about a quarter
!   fewer operations than relaxing every direction against its geq, and each
!   g is stored as soon as it is computed (no values held to the end of the
!   body). The results differ from that form by round-off only (1e-14).
!================================================
SUBROUTINE KERNEL_NAME( inv_dx, C, P, Ux, Uy, h, g )
    IMPLICIT NONE
    REAL(wp), INTENT(IN)    :: inv_dx      ! 1/dx in finest-level lattice units
    REAL(wp), INTENT(IN)    :: C(-1:nx+2,-1:ny+2)
#ifdef WITH_MOMENTS
    REAL(wp), INTENT(OUT)   :: P (nx,ny), Ux(nx,ny), Uy(nx,ny)
#else
    REAL(wp), INTENT(IN)    :: P (nx,ny), Ux(nx,ny), Uy(nx,ny)
#endif
    REAL(wp), INTENT(INOUT) :: h(0:nx+1,0:ny+1,9), g(0:nx+1,0:ny+1,9)

    INTEGER :: i, j
    !== per block
    REAL(wp) :: cC, cM, dRho, dRho3, dRho6, invTaul, dInvTau, tws, twd, omwc, w0d, wsd, wdd
    REAL(wp) :: c4B, cK, hGy, ws3, wd3, ws45, wd18, wd12
    !== per cell
    REAL(wp) :: C0, C3, eC2, eC3, eC6, eC7, B2, B3, B4, B5, B6, B7, B8, B9, t68, t79
    REAL(wp) :: DcDx, DcDy, DmDx, DmDy, g2, rn, ni, nj, xmu, vx, vy, vxx, vyy, base0, uDc, D0
    REAL(wp) :: Rho, rho3, hF, Pc, Pr, amu, cH, tCH, t, km, m8, m9, r0, rp, rm, Q
    REAL(wp) :: v, Gs, Ga, Gp, Gm, sq, Pw
    REAL(wp) :: S0, Sx, Sy, rRho

    cC      = 0.5_wp *inv_dx               ! central difference (C(x+e)-C(x-e))/(2 dx)
    cM      = 0.25_wp*inv_dx               ! biased difference over 4 dx
    dRho    = Rhoh - Rhol
    dRho3   = dRho*c13
    dRho6   = dRho*c16
    invTaul = 1._wp/taul                    ! 1/tau is linear in C between
    dInvTau = 1._wp/tauh - invTaul          ! the light (C=0) and heavy (C=1) fluid
    tws     = 12._wp*M0/W * ws              ! sharpening: 12 M0 (1-C)/W, times the weight
    twd     = 12._wp*M0/W * wd
    omwc    = 1._wp - wc
    w0d     = w0*dRho3
    wsd     = ws*dRho3
    wdd     = wd*dRho3
    c4B     = 4._wp*Beta                    ! chemical potential: bulk part and
    cK      = kappa*inv_dx*inv_dx/6._wp     ! kappa/dx**2 with the 1/6 of the Laplacian
    hGy     = 0.5_wp*Gy
    ws3     = 3._wp*ws
    wd3     = 3._wp*wd
    ws45    = 4.5_wp*ws
    wd18    = 18._wp*wd
    wd12    = 12._wp*wd

    DO j = 1, ny
!$OMP SIMD PRIVATE( C0, C3, eC2, eC3, eC6, eC7, B2, B3, B4, B5, B6, B7, B8, B9, t68, t79, &
!$OMP   DcDx, DcDy, DmDx, DmDy, g2, rn, ni, nj, xmu, vx, vy, vxx, vyy, base0, uDc, D0, &
!$OMP   Rho, rho3, hF, Pc, Pr, amu, cH, tCH, t, km, m8, m9, r0, rp, rm, Q, &
!$OMP   v, Gs, Ga, Gp, Gm, sq, Pw, S0, Sx, Sy, rRho )
    DO i = 1, nx
        C0 = C(i,j)

        !== central differences along e2 = (1,0), e3 = (0,1), e6 = (1,1),
        !   e7 = (-1,1); the opposite directions have the negative values
        eC2 = ( C(i+1,j  ) - C(i-1,j  ) )*cC
        eC3 = ( C(i  ,j+1) - C(i  ,j-1) )*cC
        eC6 = ( C(i+1,j+1) - C(i-1,j-1) )*cC
        eC7 = ( C(i-1,j+1) - C(i+1,j-1) )*cC

        !== biased differences B_k = ( 4 C(x+e_k) - 3 C(x) - C(x+2e_k) )/(4 dx): the
        !   mixed difference eM_k = ( 5 C(x+e) - C(x-e) - 3 C(x) - C(x+2e) )/(4 dx)
        !   minus half the central one
        C3 = 3._wp*C0
        B2 = ( 4._wp*C(i+1,j  ) - C3 - C(i+2,j  ) )*cM
        B4 = ( 4._wp*C(i-1,j  ) - C3 - C(i-2,j  ) )*cM
        B3 = ( 4._wp*C(i  ,j+1) - C3 - C(i  ,j+2) )*cM
        B5 = ( 4._wp*C(i  ,j-1) - C3 - C(i  ,j-2) )*cM
        B6 = ( 4._wp*C(i+1,j+1) - C3 - C(i+2,j+2) )*cM
        B8 = ( 4._wp*C(i-1,j-1) - C3 - C(i-2,j-2) )*cM
        B7 = ( 4._wp*C(i-1,j+1) - C3 - C(i-2,j+2) )*cM
        B9 = ( 4._wp*C(i+1,j-1) - C3 - C(i+2,j-2) )*cM

        !== isotropic gradients: sum_k w_k e_k (differences along e_k)/cs^2;
        !   for the mixed differences eM_k - eM_k' = (B_k - B_k') + eC_k
        DcDx = c23*eC2 + c16*( eC6 - eC7 )
        DcDy = c23*eC3 + c16*( eC6 + eC7 )
        t68  = ( B6 - B8 ) + eC6
        t79  = ( B7 - B9 ) + eC7
        DmDx = c13*( ( B2 - B4 ) + eC2 ) + c112*( t68 - t79 )
        DmDy = c13*( ( B3 - B5 ) + eC3 ) + c112*( t68 + t79 )

        !== unit normal; zero away from the interface (|grad C|**2 <= grad2_min),
        !   branch-free (a mask over a clamped root) so that the loop vectorizes
        g2 = DcDx*DcDx + DcDy*DcDy
        rn = MERGE( 1._wp, 0._wp, g2 > grad2_min ) / SQRT( MAX( g2, grad2_min ) )
        ni = DcDx*rn
        nj = DcDy*rn

        !== chemical potential mu = 4 beta C (C-1) (C-1/2) - kappa Laplacian(C),
        !   isotropic 9-point Laplacian (as Chemical_Potential_AMR)
        xmu = c4B*C0*(C0-1._wp)*(C0-0.5_wp) - cK*( ( C(i-1,j-1) + C(i+1,j-1) + C(i-1,j+1) + C(i+1,j+1) ) &
            + 4._wp*( C(i,j-1) + C(i-1,j) + C(i+1,j) + C(i,j+1) ) - 20._wp*C0 )

        Rho = Rhol + C0*dRho
#ifdef WITH_MOMENTS
        !== pressure and velocity from g and the forces (as Macroscopic_Values_g_AMR):
        !   u = ( sum(g e)/cs^2 + F/2 )/rho
        S0 = g(i,j,1) + ( (g(i,j,2) + g(i,j,3)) + (g(i,j,4) + g(i,j,5)) ) &
                      + ( (g(i,j,6) + g(i,j,7)) + (g(i,j,8) + g(i,j,9)) )
        Sx = ( g(i,j,2) - g(i,j,4) ) + ( (g(i,j,6) - g(i,j,7)) + (g(i,j,9) - g(i,j,8)) )
        Sy = ( g(i,j,3) - g(i,j,5) ) + ( (g(i,j,6) + g(i,j,7)) - (g(i,j,8) + g(i,j,9)) )
        rRho = 1._wp/Rho
        vx = ( 3._wp*Sx + 0.5_wp*( xmu*DcDx ) )*rRho
        vy = ( 3._wp*Sy + 0.5_wp*( xmu*DcDy + (Rhoh - Rho)*Gy ) )*rRho
        Pc = S0 + dRho6*( vx*DcDx + vy*DcDy )
        Ux(i,j) = vx
        Uy(i,j) = vy
        P (i,j) = Pc
#else
        vx = Ux(i,j)
        vy = Uy(i,j)
        Pc = P(i,j)
#endif

        !== per-cell quantities
        vxx   = vx*vx
        vyy   = vy*vy
        base0 = 1._wp - 1.5_wp*( vxx + vyy )
        uDc   = vx*DcDx + vy*DcDy
        D0    = 0.5_wp*uDc - ( vx*DmDx + vy*DmDy )
        rho3  = Rho*c13
        hF    = (Rhoh - Rho)*hGy                ! half the buoyancy
        amu   = dRho3 + xmu                     ! A_k = Gamma_k amu - w_k dRho3
        Pr    = Pc - rho3
        r0    = rho3 - hF*vy                    ! r for ey = 0, +1, -1
        rp    = r0 + hF
        rm    = r0 - hF
        cH    = wc*C0                           ! h relaxes towards C (Gamma + sharpening)
        tCH   = 1._wp - C0

        !== (1 - s9)/4 for the relaxation rate s9 = 1/tau with tau = 1/t + 1/2, where
        !   t = 1/tau_fluid is interpolated in C (C clamped, as it may overshoot
        !   [0,1] at the interface): 1 - s9 = (2 - t)/(2 + t)
        t  = MIN( MAX( C0, 0._wp ), 1._wp ) * dInvTau + invTaul
        km = 0.25_wp*( 2._wp - t )/( 2._wp + t )

        !== stress moments of g - geq, relaxed. Sums of geq over a pair:
        !     geq+ + geq- = 2 Pw + Gamma+ q+ + Gamma- q- - Aa eC + As u.grad C,
        !   q+- = rho/3 - (+-ey - vy) Fb/2, As = Gs amu - w dRho3, Aa = Ga amu
        Q  = 2._wp*( rho3 + hF*vy ) + amu*uDc
        m8 = km*( ( (g(i,j,2) + g(i,j,4)) - (g(i,j,3) + g(i,j,5)) ) &
                - ( ws45*(vxx - vyy)*Q + ws3*( 2._wp*vy*hF - amu*( vx*eC2 - vy*eC3 ) ) ) )
        m9 = km*( ( (g(i,j,6) + g(i,j,8)) - (g(i,j,7) + g(i,j,9)) ) &
                - ( wd18*vx*vy*Q - wd12*hF*vx - wd3*amu*( (vx + vy)*eC6 - (vy - vx)*eC7 ) ) )

        !== rest direction (e = 0): no stress moments
        Gp = w0*base0
        h(i,j,1) = omwc*h(i,j,1) + cH*Gp
        g(i,j,1) = Pr*w0 + Gp*r0 + ( Gp*amu - w0d )*D0

        !== pair (2,4): e = (+-1, 0)
        Gs = ws*( base0 + 4.5_wp*vxx )
        Ga = ws3*vx
        Gp = Gs + Ga
        Gm = Gs - Ga
        sq = tCH*tws*ni
        h(i,j,2) = omwc*h(i,j,2) + cH*( Gp + sq )
        h(i,j,4) = omwc*h(i,j,4) + cH*( Gm - sq )
        Pw = Pr*ws + m8
        g(i,j,2) = Pw + Gp*r0 + ( Gp*amu - wsd )*( B2 + D0 )
        g(i,j,4) = Pw + Gm*r0 + ( Gm*amu - wsd )*( B4 + D0 )

        !== pair (3,5): e = (0, +-1)
        Gs = ws*( base0 + 4.5_wp*vyy )
        Ga = ws3*vy
        Gp = Gs + Ga
        Gm = Gs - Ga
        sq = tCH*tws*nj
        h(i,j,3) = omwc*h(i,j,3) + cH*( Gp + sq )
        h(i,j,5) = omwc*h(i,j,5) + cH*( Gm - sq )
        Pw = Pr*ws - m8
        g(i,j,3) = Pw + Gp*rp + ( Gp*amu - wsd )*( B3 + D0 )
        g(i,j,5) = Pw + Gm*rm + ( Gm*amu - wsd )*( B5 + D0 )

        !== pair (6,8): e = +-(1, 1)
        v  = vx + vy
        Gs = wd*( base0 + 4.5_wp*v*v )
        Ga = wd3*v
        Gp = Gs + Ga
        Gm = Gs - Ga
        sq = tCH*twd*( ni + nj )
        h(i,j,6) = omwc*h(i,j,6) + cH*( Gp + sq )
        h(i,j,8) = omwc*h(i,j,8) + cH*( Gm - sq )
        Pw = Pr*wd + m9
        g(i,j,6) = Pw + Gp*rp + ( Gp*amu - wdd )*( B6 + D0 )
        g(i,j,8) = Pw + Gm*rm + ( Gm*amu - wdd )*( B8 + D0 )

        !== pair (7,9): e = +-(-1, 1)
        v  = vy - vx
        Gs = wd*( base0 + 4.5_wp*v*v )
        Ga = wd3*v
        Gp = Gs + Ga
        Gm = Gs - Ga
        sq = tCH*twd*( nj - ni )
        h(i,j,7) = omwc*h(i,j,7) + cH*( Gp + sq )
        h(i,j,9) = omwc*h(i,j,9) + cH*( Gm - sq )
        Pw = Pr*wd - m9
        g(i,j,7) = Pw + Gp*rp + ( Gp*amu - wdd )*( B7 + D0 )
        g(i,j,9) = Pw + Gm*rm + ( Gm*amu - wdd )*( B9 + D0 )
    END DO
    END DO

END SUBROUTINE KERNEL_NAME

#undef KERNEL_NAME
#undef WITH_MOMENTS
