!================================================
!   Improved-locality conservative phase-field LBM (PRE 2017)
!
!   A. Fakhari, T. Mitchell, C. Leonardi, D. Bolster, "Improved locality
!   of the phase-field lattice-Boltzmann model for immiscible fluids at
!   high density ratios", Phys. Rev. E 96, 053301 (2017);
!   uniform-grid reference code (GPL-3.0):
!   https://github.com/AbbasHub/Improved_Phase-Field_LBM_2D
!
!   Differences to Model_JCP2016_Mod:
!     - interface (h): the sharpening term 4C(1-C)/W e.n is a source term
!       of the h equation instead of part of its equilibrium
!     - hydrodynamics (g): g is normalised by rho*cs^2 (sum g = 3p/rho);
!       surface tension mu*grad C, the pressure force -3p/rho*(Rhoh-Rhol)/3
!       *grad C, the viscous force from the non-equilibrium stress of g and
!       buoyancy enter as a forcing term. Only the isotropic gradient and
!       Laplacian of C are needed: one ghost layer of C instead of two.
!     - collision of g: MRT with S = diag(1,...,1,s,s), s = 1/(tau+1/2), as
!       in the paper (Eqs. 27-29; the reference code uses BGK); the viscous
!       force uses the stress relaxed with the same rates (Eq. 31)
!     - relaxation time from the linearly interpolated dynamic viscosity,
!       tau = mu/(rho cs^2) (Eqs. 24-25)
!   Blk(n)%P holds the physical pressure p, as in the other model.
!   In the AMR setting derivatives are scaled with the block spacing dx
!   (finest-level lattice units).
!================================================
MODULE Model_PRE2017_Mod
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: nx, ny, inv_dx_level     ! the fixed block size, 1/dx per level
    USE D2Q9_Mod, ONLY: Wa, exr, eyr, w0, ws, wd
    USE Physics_Mod, ONLY: W, wc, Gy, taul, tauh, Rhol, Rhoh
    USE Model_Mod, ONLY: LBM_Model
    USE Phase_Field_Mod, ONLY: Chemical_Potential_AMR, Gradient_AMR, Max_Gradient, grad2_min
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: PRE2017_Model

    REAL(wp), PARAMETER :: c13 = 1._wp/3._wp

CONTAINS

FUNCTION PRE2017_Model() RESULT( model )
    TYPE(LBM_Model) :: model

    model%name       =  'improved-locality phase-field LBM (Fakhari et al., PRE 2017)'
    model%ghost_C    =  1       ! isotropic differences only
    model%Initialize => Initialize_Block
    model%Collide    => Collide_Block
    model%Moments_g  => Moments_g_Block
    model%Moments_and_Collide => Moments_and_Collide_Block
    model%Refinement_Error    => Interface_Light_Phase_Error

END FUNCTION PRE2017_Model

!**********************************************************************

!== refinement criterion: the interface (largest |grad C|) and in addition the
!   light phase (C < 1) of block b, max( |grad C|, max(1-C) ). Keeps the
!   interior of bubbles on the finest level, which this model needs for
!   stability (g is normalised by rho/3).
FUNCTION Interface_Light_Phase_Error( b ) RESULT( err )
    USE Blocks_Mod, ONLY: Blk, Blk_C
    INTEGER, INTENT(IN) :: b
    REAL(wp) :: err

    err = MAX( Max_Gradient( inv_dx_level(Blk(b)%Level), Blk_C(:,:,b) ), MAXVAL( 1._wp - Blk_C(1:nx,1:ny, b) ) )

END FUNCTION Interface_Light_Phase_Error

!**********************************************************************

SUBROUTINE Collide_Block( n )
    USE Blocks_Mod, ONLY: Blk, Blk_h, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    REAL(wp), DIMENSION(nx,ny) :: mu, DcDx, DcDy

    CALL Chemical_Potential_AMR( inv_dx_level(Blk(n)%Level), Blk_C(:,:,n), mu )
    CALL Gradient_AMR( inv_dx_level(Blk(n)%Level), Blk_C(:,:,n), DcDx, DcDy )
    CALL Collision_h_g( Blk_C(:,:,n), mu, DcDx, DcDy, Blk_P(:,:,n), &
                        Blk_Ux(:,:,n), Blk_Uy(:,:,n), Blk_h(:,:,:,n), Blk_g(:,:,:,n) )

END SUBROUTINE Collide_Block

!**********************************************************************

SUBROUTINE Moments_g_Block( n )
    USE Blocks_Mod, ONLY: Blk, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    REAL(wp), DIMENSION(nx,ny) :: mu, DcDx, DcDy

    CALL Chemical_Potential_AMR( inv_dx_level(Blk(n)%Level), Blk_C(:,:,n), mu )
    CALL Gradient_AMR( inv_dx_level(Blk(n)%Level), Blk_C(:,:,n), DcDx, DcDy )
    CALL Moments_g( Blk_g(:,:,:,n), Blk_C(:,:,n), mu, DcDx, DcDy, Blk_P(:,:,n), Blk_Ux(:,:,n), Blk_Uy(:,:,n) )

END SUBROUTINE Moments_g_Block

!**********************************************************************

!== moments of g and collision of one block with one chemical potential
!   and gradient of C
SUBROUTINE Moments_and_Collide_Block( n )
    USE Blocks_Mod, ONLY: Blk, Blk_h, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    REAL(wp), DIMENSION(nx,ny) :: mu, DcDx, DcDy

    CALL Chemical_Potential_AMR( inv_dx_level(Blk(n)%Level), Blk_C(:,:,n), mu )
    CALL Gradient_AMR( inv_dx_level(Blk(n)%Level), Blk_C(:,:,n), DcDx, DcDy )
    CALL Moments_g( Blk_g(:,:,:,n), Blk_C(:,:,n), mu, DcDx, DcDy, Blk_P(:,:,n), Blk_Ux(:,:,n), Blk_Uy(:,:,n) )
    CALL Collision_h_g( Blk_C(:,:,n), mu, DcDx, DcDy, Blk_P(:,:,n), &
                        Blk_Ux(:,:,n), Blk_Uy(:,:,n), Blk_h(:,:,:,n), Blk_g(:,:,:,n) )

END SUBROUTINE Moments_and_Collide_Block

!**********************************************************************

!== equilibrium h and g of block n from C, P, Ux, Uy
SUBROUTINE Initialize_Block( n )
    USE Blocks_Mod, ONLY: Blk, Blk_h, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    INTEGER :: i, j, k
    REAL(wp) :: eU, eN, Gam, tmp, g2, c4W
    REAL(wp), DIMENSION(nx,ny) :: DcDx, DcDy, ni, nj, Pn, U2, aC

    c4W = 4._wp/W

    !== C: the interior cells only (an associated section starts at index 1)
    ASSOCIATE( C => Blk_C(1:nx,1:ny,n), Ux => Blk_Ux(:,:,n), Uy => Blk_Uy(:,:,n) )

        CALL Gradient_AMR( inv_dx_level(Blk(n)%Level), Blk_C(:,:,n), DcDx, DcDy )

        DO j = 1, ny
        DO i = 1, nx
            !== unit normal; zero away from the interface (|grad C|**2 <= grad2_min),
            !   branch-free (a mask over a clamped root) so that the loop vectorizes
            g2 = DcDx(i,j)**2 + DcDy(i,j)**2
            tmp = MERGE( 1._wp, 0._wp, g2 > grad2_min ) / SQRT( MAX( g2, grad2_min ) )
            ni(i,j) = DcDx(i,j) * tmp
            nj(i,j) = DcDy(i,j) * tmp
            Pn(i,j) = 3._wp*Blk_P(i,j, n) / ( Rhol + C(i,j) * (Rhoh - Rhol) )
            U2(i,j) = Ux(i,j)*Ux(i,j) + Uy(i,j)*Uy(i,j)
            aC(i,j) = c4W*C(i,j)*(1._wp - C(i,j))
        END DO
        END DO

        DO k = 1, 9
            DO j = 1, ny
            DO i = 1, nx
                eU  = exr(k) * Ux(i,j) + eyr(k) * Uy(i,j)
                eN  = exr(k) * ni(i,j) + eyr(k) * nj(i,j)
                Gam = Wa(k) * ( 1._wp + eU*(3._wp + 4.5_wp*eU) - 1.5_wp*U2(i,j) )

                Blk_h(i,j,k, n) = C(i,j)*Gam - 0.5_wp * Wa(k)*aC(i,j)*eN
                Blk_g(i,j,k, n) = Pn(i,j)*Wa(k) + (Gam - Wa(k))
            END DO
            END DO
        END DO

    END ASSOCIATE

END SUBROUTINE Initialize_Block

!**********************************************************************

!== collision of h and g for one block, cell by cell in a single pass. The
!   directions come in opposite pairs (2,4), (3,5), (6,8), (7,9), e_k' = -e_k,
!   sharing e.u, e.n and e.F (the odd parts change sign). Per cell:
!   1) equilibria, h relaxation with the sharpening source, g - geq
!   2) the stress of g - geq and the total force (as in Moments_g)
!   3) reduced MRT of g with the forcing term gS = 3 w e.F/rho: the stress
!      moments m8 = Pxx-Pyy, m9 = Pxy relax with s = 1/tau, all others with 1.
!      gS is odd within a pair and m8, m9 are even, so the shift of geq by
!      -gS/2 cancels in m8, m9 and they come from g - geq directly.
SUBROUTINE Collision_h_g( C, mu, DcDx, DcDy, P, Ux, Uy, h, g )
    REAL(wp), INTENT(IN)    :: C(-1:nx+2,-1:ny+2)
    REAL(wp), INTENT(IN), DIMENSION(nx,ny) :: mu, DcDx, DcDy    ! chemical potential, grad C
    REAL(wp), INTENT(IN)    :: P (nx,ny), Ux(nx,ny), Uy(nx,ny)
    REAL(wp), INTENT(INOUT) :: h(0:nx+1,0:ny+1,9), g(0:nx+1,0:ny+1,9)

    INTEGER :: i, j
    !== per block
    REAL(wp) :: c4W, dRho, dRho3, RT0, dRT, omwc, cHs
    !== per cell
    REAL(wp) :: C0, mu0, dcx, dcy, g2, rn, ni, nj, vx, vy, base0, Rho, r3, pw, wcC, aCc
    REAL(wp) :: Cc, tau, s, sm, sxx, syy, sxy, txx, tyy, txy, vis, prs, Fx, Fy, m8, m9
    REAL(wp) :: v, q, base, Gp, Gm, hS, gS2, gS3, gS6, gS7
    REAL(wp) :: ge2, ge3, ge4, ge5, ge6, ge7, ge8, ge9, dg2, dg3, dg4, dg5, dg6, dg7, dg8, dg9

    c4W   = 4._wp/W                     ! sharpening: 4 C (1-C)/W
    dRho  = Rhoh - Rhol
    dRho3 = dRho*c13
    RT0   = Rhol*taul                   ! tau: mu = rho tau/3 linear in C
    dRT   = Rhoh*tauh - Rhol*taul
    omwc  = 1._wp - wc
    cHs   = 1._wp - 0.5_wp*wc           ! the source enters h with (1 - wc/2)

    DO j = 1, ny
!$OMP SIMD PRIVATE( C0, mu0, dcx, dcy, g2, rn, ni, nj, vx, vy, base0, Rho, r3, pw, wcC, aCc, &
!$OMP   Cc, tau, s, sm, sxx, syy, sxy, txx, tyy, txy, vis, prs, Fx, Fy, m8, m9, &
!$OMP   v, q, base, Gp, Gm, hS, gS2, gS3, gS6, gS7, &
!$OMP   ge2, ge3, ge4, ge5, ge6, ge7, ge8, ge9, dg2, dg3, dg4, dg5, dg6, dg7, dg8, dg9 )
    DO i = 1, nx
        C0  = C(i,j)
        mu0 = mu(i,j)
        dcx = DcDx(i,j)
        dcy = DcDy(i,j)

        !== unit normal; zero away from the interface (|grad C|**2 <= grad2_min),
        !   branch-free (a mask over a clamped root) so that the loop vectorizes
        g2 = dcx*dcx + dcy*dcy
        rn = MERGE( 1._wp, 0._wp, g2 > grad2_min ) / SQRT( MAX( g2, grad2_min ) )
        ni = dcx*rn
        nj = dcy*rn

        !== per-cell quantities; relaxation time (with the 1/2) from the dynamic
        !   viscosity interpolated linearly in C (see Moments_g)
        vx    = Ux(i,j)
        vy    = Uy(i,j)
        base0 = 1._wp - 1.5_wp*( vx*vx + vy*vy )
        Rho   = Rhol + C0*dRho
        r3    = 3._wp/Rho
        pw    = P(i,j)*r3 - 1._wp           ! geq_k = w_k (Pn - 1) + Gamma_k
        wcC   = wc*C0
        aCc   = c4W*C0*(1._wp - C0)
        Cc    = MIN( MAX( C0, 0._wp ), 1._wp )
        tau   = ( RT0 + Cc*dRT ) / ( Rhol + Cc*dRho ) + 0.5_wp
        s     = 1._wp/tau

        !== rest direction: no source, no stress, no forcing
        Gp = w0*base0
        h(i,j,1) = omwc*h(i,j,1) + wcC*Gp
        g(i,j,1) = w0*pw + Gp

        !== pair (2,4): e = (+-1, 0)
        v = vx;  q = ni
        base = base0 + 4.5_wp*v*v
        Gp = ws*( base + 3._wp*v )
        Gm = ws*( base - 3._wp*v )
        hS = ws*aCc*q*cHs
        h(i,j,2) = omwc*h(i,j,2) + wcC*Gp + hS
        h(i,j,4) = omwc*h(i,j,4) + wcC*Gm - hS
        ge2 = ws*pw + Gp
        ge4 = ws*pw + Gm
        dg2 = g(i,j,2) - ge2
        dg4 = g(i,j,4) - ge4

        !== pair (3,5): e = (0, +-1)
        v = vy;  q = nj
        base = base0 + 4.5_wp*v*v
        Gp = ws*( base + 3._wp*v )
        Gm = ws*( base - 3._wp*v )
        hS = ws*aCc*q*cHs
        h(i,j,3) = omwc*h(i,j,3) + wcC*Gp + hS
        h(i,j,5) = omwc*h(i,j,5) + wcC*Gm - hS
        ge3 = ws*pw + Gp
        ge5 = ws*pw + Gm
        dg3 = g(i,j,3) - ge3
        dg5 = g(i,j,5) - ge5

        !== pair (6,8): e = +-(1, 1)
        v = vx + vy;  q = ni + nj
        base = base0 + 4.5_wp*v*v
        Gp = wd*( base + 3._wp*v )
        Gm = wd*( base - 3._wp*v )
        hS = wd*aCc*q*cHs
        h(i,j,6) = omwc*h(i,j,6) + wcC*Gp + hS
        h(i,j,8) = omwc*h(i,j,8) + wcC*Gm - hS
        ge6 = wd*pw + Gp
        ge8 = wd*pw + Gm
        dg6 = g(i,j,6) - ge6
        dg8 = g(i,j,8) - ge8

        !== pair (7,9): e = +-(-1, 1)
        v = vy - vx;  q = nj - ni
        base = base0 + 4.5_wp*v*v
        Gp = wd*( base + 3._wp*v )
        Gm = wd*( base - 3._wp*v )
        hS = wd*aCc*q*cHs
        h(i,j,7) = omwc*h(i,j,7) + wcC*Gp + hS
        h(i,j,9) = omwc*h(i,j,9) + wcC*Gm - hS
        ge7 = wd*pw + Gp
        ge9 = wd*pw + Gm
        dg7 = g(i,j,7) - ge7
        dg9 = g(i,j,9) - ge9

        !== stress of g - geq; viscous force with the MRT-relaxed stress
        !   (deviatoric part with s, trace with 1), pressure, surface tension, buoyancy
        sxx = ( dg2 + dg4 ) + ( (dg6 + dg8) + (dg7 + dg9) )
        syy = ( dg3 + dg5 ) + ( (dg6 + dg8) + (dg7 + dg9) )
        sxy = ( dg6 + dg8 ) - ( dg7 + dg9 )
        txx = 0.5_wp*( sxx + syy + s*(sxx - syy) )
        tyy = 0.5_wp*( sxx + syy - s*(sxx - syy) )
        txy = s*sxy
        vis = (tau - 0.5_wp)*dRho
        prs = mu0 - (pw + 1._wp)*dRho3
        Fx  = prs*dcx - vis*( txx*dcx + txy*dcy )
        Fy  = prs*dcy - vis*( txy*dcx + tyy*dcy ) + (Rhoh - Rho)*Gy

        !== forcing term gS = 3 w (e.F)/rho, entering as geq - gS/2 + gS
        gS2 = 0.5_wp*ws*r3*Fx
        gS3 = 0.5_wp*ws*r3*Fy
        gS6 = 0.5_wp*wd*r3*( Fx + Fy )
        gS7 = 0.5_wp*wd*r3*( Fy - Fx )

        !== reduced MRT: the stress moments relax with s, all others with 1
        sm = 0.25_wp*( 1._wp - s )
        m8 = ( (dg2 + dg4) - (dg3 + dg5) )*sm
        m9 = ( (dg6 + dg8) - (dg7 + dg9) )*sm
        g(i,j,2) = ge2 + gS2 + m8
        g(i,j,4) = ge4 - gS2 + m8
        g(i,j,3) = ge3 + gS3 - m8
        g(i,j,5) = ge5 - gS3 - m8
        g(i,j,6) = ge6 + gS6 + m9
        g(i,j,8) = ge8 - gS6 + m9
        g(i,j,7) = ge7 + gS7 - m9
        g(i,j,9) = ge9 - gS7 - m9
    END DO
    END DO

END SUBROUTINE Collision_h_g

!**********************************************************************

!== normalised pressure and velocity from g, cell by cell in one pass. The
!   force (as in Collision_h_g) needs the stress of g - geq, with geq from the
!   new pressure and the velocity of the previous step (in Ux, Uy). Only the
!   sums g_k + g_k' of opposite directions enter the stress, and in them the
!   odd parts of the equilibrium cancel: per pair g_k + g_k' - 2 w (Pn - 1
!   + 1 - 3/2 u.u + 9/2 (e.u)**2).
!   Relaxation time (with the 1/2), Eqs. (24)-(25): the dynamic viscosity is
!   interpolated linearly in C (clamped to [0,1]) between the light and heavy
!   fluid and tau - 1/2 = mu/(rho cs^2), with taul, tauh = 3 mu/rho of the
!   pure phases. Force: surface tension mu grad C, pressure -Pn (Rhoh-Rhol)/3
!   grad C, the viscous force -(tau-1/2) (S P) . grad C (Rhoh-Rhol) with the
!   non-equilibrium stress P relaxed with the MRT rates (deviatoric part with
!   s = 1/tau, trace with 1), and buoyancy.
SUBROUTINE Moments_g( g, C, mu, DcDx, DcDy, P, Ux, Uy )
    REAL(wp), INTENT(IN)    :: g(0:nx+1,0:ny+1,9), C(-1:nx+2,-1:ny+2)
    REAL(wp), INTENT(IN), DIMENSION(nx,ny) :: mu, DcDx, DcDy    ! chemical potential, grad C
    REAL(wp), INTENT(OUT)   :: P(nx,ny)
    REAL(wp), INTENT(INOUT) :: Ux(nx,ny), Uy(nx,ny)

    INTEGER :: i, j
    !== per block
    REAL(wp) :: dRho, dRho3, RT0, dRT
    !== per cell
    REAL(wp) :: S0, Sx, Sy, vx, vy, base0, pb, dx2, dy2, dd6, dd7, sxx, syy, sxy
    REAL(wp) :: C0, Rho, Cc, tau, s, txx, tyy, txy, vis, prs, Fx, Fy, hR

    dRho  = Rhoh - Rhol
    dRho3 = dRho*c13
    RT0   = Rhol*taul
    dRT   = Rhoh*tauh - Rhol*taul

    DO j = 1, ny
!$OMP SIMD PRIVATE( S0, Sx, Sy, vx, vy, base0, pb, dx2, dy2, dd6, dd7, sxx, syy, sxy, &
!$OMP   C0, Rho, Cc, tau, s, txx, tyy, txy, vis, prs, Fx, Fy, hR )
    DO i = 1, nx
        !== zeroth and first moments: Pn = sum g, sum g e
        S0 = g(i,j,1) + ( (g(i,j,2) + g(i,j,3)) + (g(i,j,4) + g(i,j,5)) ) &
                      + ( (g(i,j,6) + g(i,j,7)) + (g(i,j,8) + g(i,j,9)) )
        Sx = ( g(i,j,2) - g(i,j,4) ) + ( (g(i,j,6) - g(i,j,7)) + (g(i,j,9) - g(i,j,8)) )
        Sy = ( g(i,j,3) - g(i,j,5) ) + ( (g(i,j,6) + g(i,j,7)) - (g(i,j,8) + g(i,j,9)) )

        !== pair sums of g - geq (previous velocity), and the stress
        vx    = Ux(i,j)
        vy    = Uy(i,j)
        base0 = 1._wp - 1.5_wp*( vx*vx + vy*vy )
        pb    = S0 - 1._wp + base0                      ! Pn - 1 + base0
        dx2 = ( g(i,j,2) + g(i,j,4) ) - 2._wp*ws*( pb + 4.5_wp*vx*vx )
        dy2 = ( g(i,j,3) + g(i,j,5) ) - 2._wp*ws*( pb + 4.5_wp*vy*vy )
        dd6 = ( g(i,j,6) + g(i,j,8) ) - 2._wp*wd*( pb + 4.5_wp*(vx + vy)**2 )
        dd7 = ( g(i,j,7) + g(i,j,9) ) - 2._wp*wd*( pb + 4.5_wp*(vy - vx)**2 )
        sxx = dx2 + ( dd6 + dd7 )
        syy = dy2 + ( dd6 + dd7 )
        sxy = dd6 - dd7

        !== density, relaxation time, force
        C0  = C(i,j)
        Rho = Rhol + C0*dRho
        Cc  = MIN( MAX( C0, 0._wp ), 1._wp )
        tau = ( RT0 + Cc*dRT ) / ( Rhol + Cc*dRho ) + 0.5_wp
        s   = 1._wp/tau
        txx = 0.5_wp*( sxx + syy + s*(sxx - syy) )
        tyy = 0.5_wp*( sxx + syy - s*(sxx - syy) )
        txy = s*sxy
        vis = (tau - 0.5_wp)*dRho
        prs = mu(i,j) - S0*dRho3
        Fx  = prs*DcDx(i,j) - vis*( txx*DcDx(i,j) + txy*DcDy(i,j) )
        Fy  = prs*DcDy(i,j) - vis*( txy*DcDx(i,j) + tyy*DcDy(i,j) ) + (Rhoh - Rho)*Gy

        !== velocity u = sum g e + F/(2 rho), pressure p = Pn rho/3
        hR = 0.5_wp/Rho
        Ux(i,j) = Sx + Fx*hR
        Uy(i,j) = Sy + Fy*hR
        P (i,j) = S0*Rho*c13
    END DO
    END DO

END SUBROUTINE Moments_g

END MODULE Model_PRE2017_Mod
