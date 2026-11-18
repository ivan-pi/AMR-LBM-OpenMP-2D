!================================================
!   Mass-conserving phase-field LBM (JCP 2016)
!
!   A. Fakhari, M. Geier, T. Lee, "A mass-conserving lattice Boltzmann
!   method with dynamic grid refinement for immiscible two-phase flows",
!   J. Comput. Phys. 315 (2016) 434-457
!
!   Conservative Allen-Cahn equation for C (h) with the interface
!   sharpening built into the equilibrium; pressure-evolution equation (g)
!   whose forcing uses central and mixed directional differences of C
!   (two ghost layers of C); reduced MRT collision of g.
!================================================
MODULE Model_JCP2016_Mod
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: nx, ny, inv_dx_level     ! the fixed block size, 1/dx per level
    USE D2Q9_Mod, ONLY: Wa, ex, ey, exr, eyr, w0, ws, wd
    USE Physics_Mod, ONLY: W, M0, wc, Gy, taul, tauh, Rhol, Rhoh
    USE Model_Mod, ONLY: LBM_Model
    USE Phase_Field_Mod, ONLY: Chemical_Potential_AMR, Gradient_AMR, Interface_Error, grad2_min
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: JCP2016_Model

    REAL(wp), PARAMETER :: c13 = 1._wp/3._wp, c112 = 1._wp/12._wp, c23 = 2._wp/3._wp, c16 = 1._wp/6._wp

CONTAINS

FUNCTION JCP2016_Model() RESULT( model )
    TYPE(LBM_Model) :: model

    model%name       =  'mass-conserving phase-field LBM (Fakhari, Geier & Lee, JCP 2016)'
    model%ghost_C    =  2       ! mixed differences reach two cells
    model%Initialize => Initialize_Block
    model%Collide    => Collide_Block
    model%Moments_g  => Moments_g_Block
    model%Moments_and_Collide => Moments_and_Collide_Block
    model%Refinement_Error    => Interface_Error

END FUNCTION JCP2016_Model

!**********************************************************************

SUBROUTINE Collide_Block( n )
    USE Blocks_Mod, ONLY: Blk, Blk_h, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    REAL(wp) :: mu(nx,ny), inv_dx

    inv_dx = inv_dx_level( Blk(n)%Level )
    CALL Chemical_Potential_AMR( inv_dx, Blk_C(:,:,n), mu )
    CALL Collision_h_g_AMR( inv_dx, Blk_C(:,:,n), mu, Blk_P(:,:,n), &
                            Blk_Ux(:,:,n), Blk_Uy(:,:,n), Blk_h(:,:,:,n), Blk_g(:,:,:,n) )

END SUBROUTINE Collide_Block

!**********************************************************************

SUBROUTINE Moments_g_Block( n )
    USE Blocks_Mod, ONLY: Blk, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    REAL(wp) :: mu(nx,ny), inv_dx

    inv_dx = inv_dx_level( Blk(n)%Level )
    CALL Chemical_Potential_AMR( inv_dx, Blk_C(:,:,n), mu )
    CALL Macroscopic_Values_g_AMR( inv_dx, Blk_g(:,:,:,n), Blk_C(:,:,n), mu, &
                                   Blk_P(:,:,n), Blk_Ux(:,:,n), Blk_Uy(:,:,n) )

END SUBROUTINE Moments_g_Block

!**********************************************************************

!== moments of g and collision of one block with one chemical potential
SUBROUTINE Moments_and_Collide_Block( n )
    USE Blocks_Mod, ONLY: Blk, Blk_h, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    REAL(wp) :: mu(nx,ny), inv_dx

    inv_dx = inv_dx_level( Blk(n)%Level )
    CALL Chemical_Potential_AMR( inv_dx, Blk_C(:,:,n), mu )
    CALL Macroscopic_Values_g_AMR( inv_dx, Blk_g(:,:,:,n), Blk_C(:,:,n), mu, &
                                   Blk_P(:,:,n), Blk_Ux(:,:,n), Blk_Uy(:,:,n) )
    CALL Collision_h_g_AMR( inv_dx, Blk_C(:,:,n), mu, Blk_P(:,:,n), &
                            Blk_Ux(:,:,n), Blk_Uy(:,:,n), Blk_h(:,:,:,n), Blk_g(:,:,:,n) )

END SUBROUTINE Moments_and_Collide_Block

!**********************************************************************

!== equilibrium h and g of block n from C, P, Ux, Uy (no body force);
!   runs once per block at the start, so it is written for clarity
SUBROUTINE Initialize_Block( n )
    USE Blocks_Mod, ONLY: Blk, Blk_h, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    INTEGER :: i, j, k
    REAL(wp) :: inv_dx, Rho, dRho3, cT, ni, nj, DcDx, DcDy, g2, rn, U, V
    REAL(wp) :: eDc(9), hlp(9), Gamma(9), eU(9), eN(9), mu(nx,ny)

    inv_dx = inv_dx_level( Blk(n)%Level )
    dRho3  = (Rhoh - Rhol)*c13
    cT     = 12._wp*M0/W

    CALL Chemical_Potential_AMR( inv_dx, Blk_C(:,:,n), mu )

    DO j = 1, ny
    DO i = 1, nx
        !== central differences along e_k and the isotropic gradient
        DO k = 1, 9
            eDc(k) = ( Blk_C(i+ex(k),j+ey(k), n) - Blk_C(i-ex(k),j-ey(k), n) )*(0.5_wp*inv_dx)
        END DO
        DcDx = c13*( eDc(2)-eDc(4) ) + c112*( eDc(6)-eDc(7)-eDc(8)+eDc(9) )
        DcDy = c13*( eDc(3)-eDc(5) ) + c112*( eDc(6)+eDc(7)-eDc(8)-eDc(9) )

        !== unit normal (zero away from the interface)
        g2 = DcDx*DcDx + DcDy*DcDy
        rn = MERGE( 1._wp, 0._wp, g2 > grad2_min ) / SQRT( MAX( g2, grad2_min ) )
        ni = DcDx*rn
        nj = DcDy*rn

        U   = Blk_Ux(i,j, n)
        V   = Blk_Uy(i,j, n)
        Rho = Rhol + Blk_C(i,j, n) * (Rhoh - Rhol)

        eU(:) = exr(:)*U  + eyr(:)*V
        eN(:) = exr(:)*ni + eyr(:)*nj
        Gamma(:) = Wa(:) * ( 1._wp + eU(:)*(3._wp + 4.5_wp*eU(:)) - 1.5_wp*(U*U + V*V) )

        Blk_h(i,j,:, n) = Blk_C(i,j, n) * ( Gamma(:) + cT*(1._wp - Blk_C(i,j, n)) * Wa(:) * eN(:) )

        hlp(:) = ( (Gamma(:)-Wa(:))*dRho3 + Gamma(:)*mu(i,j) ) * ( eDc(:) - (U*DcDx + V*DcDy) )
        Blk_g(i,j,:, n) = Blk_P(i,j, n)*Wa(:) + (Gamma(:)-Wa(:))*(Rho*c13) - 0.5_wp*hlp(:)
    END DO
    END DO

END SUBROUTINE Initialize_Block

!**********************************************************************

SUBROUTINE Collision_h_g_AMR( inv_dx, C, mu, P, Ux, Uy, h, g )
    !== Collision of h and g for one block, cell by cell in a single pass.
    !   The directions come in opposite pairs (2,4), (3,5), (6,8), (7,9),
    !   e_k' = -e_k: a pair shares e.u, e.n and the central difference of C
    !   (which changes sign), so its equilibria are a common part plus/minus
    !   an odd part.
    !   h: BGK relaxation to C (Gamma + sharpening term) with rate wc.
    !   g: MRT g <- g - M^-1 S M (g - geq) with S = diag(1,...,1,s9,s9) reduces
    !   to g = geq plus (1-s9) times the stress moments m8, m9 of g - geq,
    !   mapped back as +-m8/4 (directions 2-5) and +-m9/4 (6-9); plus the
    !   forcing term hM of the mixed differences.
    IMPLICIT NONE
    REAL(wp), INTENT(IN)    :: inv_dx      ! 1/dx in finest-level lattice units
    REAL(wp), INTENT(IN)    :: C(-1:nx+2,-1:ny+2), mu(nx,ny)    ! mu: chemical potential
    REAL(wp), INTENT(IN)    :: P (nx,ny), Ux(nx,ny), Uy(nx,ny)
    REAL(wp), INTENT(INOUT) :: h(0:nx+1,0:ny+1,9), g(0:nx+1,0:ny+1,9)

    INTEGER :: i, j
    !== per block
    REAL(wp) :: cC, cM, dRho, dRho3, invTaul, dInvTau, cT, omwc, w0d, wsd, wdd
    !== per cell
    REAL(wp) :: C0, C3, eC2, eC3, eC6, eC7, eM2, eM3, eM4, eM5, eM6, eM7, eM8, eM9
    REAL(wp) :: DcDx, DcDy, DmDx, DmDy, g2, rn, ni, nj, vx, vy, base0, uDc, uDm
    REAL(wp) :: Rho, rho3, Fb, Pc, amu, cH, tCH, t, s9, m8, m9
    REAL(wp) :: v, base, Gp, Gm, sq, Apl, Ami, fp, fm, geqp, geqm
    REAL(wp) :: gb2, gb3, gb4, gb5, gb6, gb7, gb8, gb9, dg2, dg3, dg4, dg5, dg6, dg7, dg8, dg9

    cC      = 0.5_wp *inv_dx               ! central difference (C(x+e)-C(x-e))/(2 dx)
    cM      = 0.25_wp*inv_dx               ! mixed difference over 4 dx
    dRho    = Rhoh - Rhol
    dRho3   = dRho*c13
    invTaul = 1._wp/taul                    ! 1/tau is linear in C between
    dInvTau = 1._wp/tauh - invTaul          ! the light (C=0) and heavy (C=1) fluid
    cT      = 12._wp*M0/W                   ! sharpening: tC = 12 M0 (1-C)/W
    omwc    = 1._wp - wc
    w0d     = w0*dRho3
    wsd     = ws*dRho3
    wdd     = wd*dRho3

    DO j = 1, ny
!$OMP SIMD PRIVATE( C0, C3, eC2, eC3, eC6, eC7, eM2, eM3, eM4, eM5, eM6, eM7, eM8, eM9, &
!$OMP   DcDx, DcDy, DmDx, DmDy, g2, rn, ni, nj, vx, vy, base0, uDc, uDm, &
!$OMP   Rho, rho3, Fb, Pc, amu, cH, tCH, t, s9, m8, m9, &
!$OMP   v, base, Gp, Gm, sq, Apl, Ami, fp, fm, geqp, geqm, &
!$OMP   gb2, gb3, gb4, gb5, gb6, gb7, gb8, gb9, dg2, dg3, dg4, dg5, dg6, dg7, dg8, dg9 )
    DO i = 1, nx
        C0 = C(i,j)

        !== central differences along e2 = (1,0), e3 = (0,1), e6 = (1,1),
        !   e7 = (-1,1); the opposite directions have the negative values
        eC2 = ( C(i+1,j  ) - C(i-1,j  ) )*cC
        eC3 = ( C(i  ,j+1) - C(i  ,j-1) )*cC
        eC6 = ( C(i+1,j+1) - C(i-1,j-1) )*cC
        eC7 = ( C(i-1,j+1) - C(i+1,j-1) )*cC

        !== mixed differences ( -C(x-e) - 3 C(x) + 5 C(x+e) - C(x+2e) )/(4 dx)
        C3  = 3._wp*C0
        eM2 = ( 5._wp*C(i+1,j  ) - C(i-1,j  ) - C3 - C(i+2,j  ) )*cM
        eM4 = ( 5._wp*C(i-1,j  ) - C(i+1,j  ) - C3 - C(i-2,j  ) )*cM
        eM3 = ( 5._wp*C(i  ,j+1) - C(i  ,j-1) - C3 - C(i  ,j+2) )*cM
        eM5 = ( 5._wp*C(i  ,j-1) - C(i  ,j+1) - C3 - C(i  ,j-2) )*cM
        eM6 = ( 5._wp*C(i+1,j+1) - C(i-1,j-1) - C3 - C(i+2,j+2) )*cM
        eM8 = ( 5._wp*C(i-1,j-1) - C(i+1,j+1) - C3 - C(i-2,j-2) )*cM
        eM7 = ( 5._wp*C(i-1,j+1) - C(i+1,j-1) - C3 - C(i-2,j+2) )*cM
        eM9 = ( 5._wp*C(i+1,j-1) - C(i-1,j+1) - C3 - C(i+2,j-2) )*cM

        !== isotropic gradients: sum_k w_k e_k (differences along e_k)/cs^2
        DcDx = c23*eC2 + c16*( eC6 - eC7 )
        DcDy = c23*eC3 + c16*( eC6 + eC7 )
        DmDx = c13*( eM2 - eM4 ) + c112*( (eM6 - eM8) - (eM7 - eM9) )
        DmDy = c13*( eM3 - eM5 ) + c112*( (eM6 - eM8) + (eM7 - eM9) )

        !== unit normal; zero away from the interface (|grad C|**2 <= grad2_min),
        !   branch-free (a mask over a clamped root) so that the loop vectorizes
        g2 = DcDx*DcDx + DcDy*DcDy
        rn = MERGE( 1._wp, 0._wp, g2 > grad2_min ) / SQRT( MAX( g2, grad2_min ) )
        ni = DcDx*rn
        nj = DcDy*rn

        !== per-cell quantities
        vx    = Ux(i,j)
        vy    = Uy(i,j)
        base0 = 1._wp - 1.5_wp*( vx*vx + vy*vy )
        uDc   = vx*DcDx + vy*DcDy
        uDm   = vx*DmDx + vy*DmDy
        Rho   = Rhol + C0*dRho
        rho3  = Rho*c13
        Fb    = (Rhoh - Rho)*Gy
        Pc    = P(i,j)
        amu    = dRho3 + mu(i,j)                 ! A_k = Gamma_k amu - w_k dRho3
        cH    = wc*C0                           ! h relaxes towards C (Gamma + tC w e.n)
        tCH   = cT*(1._wp - C0)

        !== relaxation rate s9 = 1/tau with tau = 1/x + 1/2, where x = 1/tau_fluid is
        !   interpolated in C (C clamped, as it may overshoot [0,1] at the interface)
        t  = MIN( MAX( C0, 0._wp ), 1._wp ) * dInvTau + invTaul
        s9 = t / ( 1._wp + 0.5_wp*t )

        !== rest direction (e = 0): no stress moments
        Gp   = w0*base0
        h(i,j,1) = omwc*h(i,j,1) + cH*Gp
        Apl   = Gp*amu - w0d
        fp   = -vy*Gp*Fb
        geqp = Pc*w0 + (Gp - w0)*rho3 + 0.5_wp*( Apl*uDc - fp )
        g(i,j,1) = geqp - Apl*uDm + fp

        !== pair (2,4): e = (+-1, 0)
        v    = vx
        base = base0 + 4.5_wp*v*v
        Gp   = ws*( base + 3._wp*v )
        Gm   = ws*( base - 3._wp*v )
        sq   = tCH*ws*ni
        h(i,j,2) = omwc*h(i,j,2) + cH*( Gp + sq )
        h(i,j,4) = omwc*h(i,j,4) + cH*( Gm - sq )
        Apl   = Gp*amu - wsd
        Ami   = Gm*amu - wsd
        fp   = -vy*Gp*Fb
        fm   = -vy*Gm*Fb
        geqp = Pc*ws + (Gp - ws)*rho3 - 0.5_wp*( Apl*( eC2 - uDc) + fp )
        geqm = Pc*ws + (Gm - ws)*rho3 - 0.5_wp*( Ami*(-eC2 - uDc) + fm )
        dg2  = g(i,j,2) - geqp
        dg4  = g(i,j,4) - geqm
        gb2  = geqp + Apl*( eM2 - uDm ) + fp
        gb4  = geqm + Ami*( eM4 - uDm ) + fm

        !== pair (3,5): e = (0, +-1)
        v    = vy
        base = base0 + 4.5_wp*v*v
        Gp   = ws*( base + 3._wp*v )
        Gm   = ws*( base - 3._wp*v )
        sq   = tCH*ws*nj
        h(i,j,3) = omwc*h(i,j,3) + cH*( Gp + sq )
        h(i,j,5) = omwc*h(i,j,5) + cH*( Gm - sq )
        Apl   = Gp*amu - wsd
        Ami   = Gm*amu - wsd
        fp   = ( 1._wp - vy)*Gp*Fb
        fm   = (-1._wp - vy)*Gm*Fb
        geqp = Pc*ws + (Gp - ws)*rho3 - 0.5_wp*( Apl*( eC3 - uDc) + fp )
        geqm = Pc*ws + (Gm - ws)*rho3 - 0.5_wp*( Ami*(-eC3 - uDc) + fm )
        dg3  = g(i,j,3) - geqp
        dg5  = g(i,j,5) - geqm
        gb3  = geqp + Apl*( eM3 - uDm ) + fp
        gb5  = geqm + Ami*( eM5 - uDm ) + fm

        !== pair (6,8): e = +-(1, 1)
        v    = vx + vy
        base = base0 + 4.5_wp*v*v
        Gp   = wd*( base + 3._wp*v )
        Gm   = wd*( base - 3._wp*v )
        sq   = tCH*wd*( ni + nj )
        h(i,j,6) = omwc*h(i,j,6) + cH*( Gp + sq )
        h(i,j,8) = omwc*h(i,j,8) + cH*( Gm - sq )
        Apl   = Gp*amu - wdd
        Ami   = Gm*amu - wdd
        fp   = ( 1._wp - vy)*Gp*Fb
        fm   = (-1._wp - vy)*Gm*Fb
        geqp = Pc*wd + (Gp - wd)*rho3 - 0.5_wp*( Apl*( eC6 - uDc) + fp )
        geqm = Pc*wd + (Gm - wd)*rho3 - 0.5_wp*( Ami*(-eC6 - uDc) + fm )
        dg6  = g(i,j,6) - geqp
        dg8  = g(i,j,8) - geqm
        gb6  = geqp + Apl*( eM6 - uDm ) + fp
        gb8  = geqm + Ami*( eM8 - uDm ) + fm

        !== pair (7,9): e = +-(-1, 1)
        v    = vy - vx
        base = base0 + 4.5_wp*v*v
        Gp   = wd*( base + 3._wp*v )
        Gm   = wd*( base - 3._wp*v )
        sq   = tCH*wd*( nj - ni )
        h(i,j,7) = omwc*h(i,j,7) + cH*( Gp + sq )
        h(i,j,9) = omwc*h(i,j,9) + cH*( Gm - sq )
        Apl   = Gp*amu - wdd
        Ami   = Gm*amu - wdd
        fp   = ( 1._wp - vy)*Gp*Fb
        fm   = (-1._wp - vy)*Gm*Fb
        geqp = Pc*wd + (Gp - wd)*rho3 - 0.5_wp*( Apl*( eC7 - uDc) + fp )
        geqm = Pc*wd + (Gm - wd)*rho3 - 0.5_wp*( Ami*(-eC7 - uDc) + fm )
        dg7  = g(i,j,7) - geqp
        dg9  = g(i,j,9) - geqm
        gb7  = geqp + Apl*( eM7 - uDm ) + fp
        gb9  = geqm + Ami*( eM9 - uDm ) + fm

        !== relax the stress moments of g - geq with s9 (all others with 1)
        m8 = ( (dg2 + dg4) - (dg3 + dg5) ) * (1._wp - s9) * 0.25_wp
        m9 = ( (dg6 + dg8) - (dg7 + dg9) ) * (1._wp - s9) * 0.25_wp
        g(i,j,2) = gb2 + m8
        g(i,j,4) = gb4 + m8
        g(i,j,3) = gb3 - m8
        g(i,j,5) = gb5 - m8
        g(i,j,6) = gb6 + m9
        g(i,j,8) = gb8 + m9
        g(i,j,7) = gb7 - m9
        g(i,j,9) = gb9 - m9
    END DO
    END DO

END SUBROUTINE Collision_h_g_AMR

!**********************************************************************

SUBROUTINE Macroscopic_Values_g_AMR( inv_dx, g, C, mu, P, Ux, Uy )
    !== pressure and velocity from g (zeroth and first moments) and the forces
    IMPLICIT NONE
    REAL(wp), INTENT(IN)  :: inv_dx         ! 1/dx in finest-level lattice units
    REAL(wp), INTENT(IN)  :: g(0:nx+1,0:ny+1,9), C(-1:nx+2,-1:ny+2), mu(nx,ny)
    REAL(wp), INTENT(OUT) :: P (nx,ny), Ux(nx,ny), Uy(nx,ny)

    INTEGER :: i, j
    REAL(wp) :: dRho, dRho6, S0, Sx, Sy, Rho, rRho, Fx, Fy
    REAL(wp) :: DcDx(nx,ny), DcDy(nx,ny)

    dRho  = Rhoh - Rhol
    dRho6 = dRho*c16

    CALL Gradient_AMR( inv_dx, C, DcDx, DcDy )

    DO j = 1, ny
!$OMP SIMD PRIVATE( S0, Sx, Sy, Rho, rRho, Fx, Fy )
    DO i = 1, nx
        S0 = g(i,j,1) + ( (g(i,j,2) + g(i,j,3)) + (g(i,j,4) + g(i,j,5)) ) &
                      + ( (g(i,j,6) + g(i,j,7)) + (g(i,j,8) + g(i,j,9)) )
        Sx = ( g(i,j,2) - g(i,j,4) ) + ( (g(i,j,6) - g(i,j,7)) + (g(i,j,9) - g(i,j,8)) )
        Sy = ( g(i,j,3) - g(i,j,5) ) + ( (g(i,j,6) + g(i,j,7)) - (g(i,j,8) + g(i,j,9)) )

        Rho  = Rhol + C(i,j)*dRho
        rRho = 1._wp/Rho

        Fx = mu(i,j)*DcDx(i,j)
        Fy = mu(i,j)*DcDy(i,j) + (Rhoh - Rho)*Gy

        !== u = ( sum(g e)/cs^2 + F/2 )/rho, from rho u = sum(f e) + F/2 and
        !   sum(g e) = cs^2 sum(f e). (Eq. (35) of the paper prints the force term
        !   without the factor cs^2; this form is consistent with its Eq. (24).)
        Ux(i,j) = ( 3._wp*Sx + 0.5_wp*Fx )*rRho
        Uy(i,j) = ( 3._wp*Sy + 0.5_wp*Fy )*rRho

        P (i,j) = S0 + dRho6*( Ux(i,j)*DcDx(i,j) + Uy(i,j)*DcDy(i,j) )
    END DO
    END DO

END SUBROUTINE Macroscopic_Values_g_AMR

END MODULE Model_JCP2016_Mod
