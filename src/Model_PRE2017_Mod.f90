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
    USE D2Q9_Mod, ONLY: Wa, exr, eyr
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

!== collision of h and g for one block:
!   1) block buffers: gradient, chemical potential, normal, rho, tau
!   2) direction loop: h relaxation with the sharpening source, g equilibria,
!      and the stress of the non-equilibrium part g - geq (viscous force)
!   3) forces, then direction loops for the reduced-MRT relaxation of g
!      with the forcing term
SUBROUTINE Collision_h_g( C, mu, DcDx, DcDy, P, Ux, Uy, h, g )
    REAL(wp), INTENT(IN)    :: C(-1:nx+2,-1:ny+2)
    REAL(wp), INTENT(IN), DIMENSION(nx,ny) :: mu, DcDx, DcDy    ! chemical potential, grad C
    REAL(wp), INTENT(IN)    :: P (nx,ny), Ux(nx,ny), Uy(nx,ny)
    REAL(wp), INTENT(INOUT) :: h(0:nx+1,0:ny+1,9), g(0:nx+1,0:ny+1,9)

    INTEGER :: i, j, k
    REAL(wp) :: tmp, eU, eN, Gam, hS, dg, g2, c4W
    !== signs of g_k - geq_k in the stress moments m8 = Pxx-Pyy and m9 = Pxy
    REAL(wp), PARAMETER :: w8(9) = [0, 1,-1, 1,-1, 0, 0, 0, 0]
    REAL(wp), PARAMETER :: w9(9) = [0, 0, 0, 0, 0, 1,-1, 1,-1]
    REAL(wp), DIMENSION(nx,ny,9) :: geq, gS
    REAL(wp), DIMENSION(nx,ny) :: ni, nj, Rho, r3Rho, Pn, tau, U2, aC, &
                                  sxx, sxy, syy, Fx, Fy, m8, m9, sm

    c4W = 4._wp/W

    DO j = 1, ny
    DO i = 1, nx
        !== unit normal; zero away from the interface (|grad C|**2 <= grad2_min),
        !   branch-free (a mask over a clamped root) so that the loop vectorizes
        g2 = DcDx(i,j)**2 + DcDy(i,j)**2
        tmp = MERGE( 1._wp, 0._wp, g2 > grad2_min ) / SQRT( MAX( g2, grad2_min ) )
        ni (i,j) = DcDx(i,j) * tmp
        nj (i,j) = DcDy(i,j) * tmp
        Rho(i,j) = Rhol + C(i,j) * (Rhoh - Rhol)
        r3Rho(i,j) = 3._wp/Rho(i,j)
        Pn (i,j) = P(i,j) * r3Rho(i,j)
        tau(i,j) = Relaxation_Time( C(i,j) )
        U2 (i,j) = Ux(i,j)*Ux(i,j) + Uy(i,j)*Uy(i,j)
        aC (i,j) = c4W*C(i,j)*(1._wp - C(i,j))
    END DO
    END DO

    sxx = 0._wp
    sxy = 0._wp
    syy = 0._wp
    DO k = 1, 9
        DO j = 1, ny
        DO i = 1, nx
            eU  = exr(k) * Ux(i,j) + eyr(k) * Uy(i,j)
            eN  = exr(k) * ni(i,j) + eyr(k) * nj(i,j)
            Gam = Wa(k) * ( 1._wp + eU*(3._wp + 4.5_wp*eU) - 1.5_wp*U2(i,j) )

            !== h: equilibrium C*Gamma, sharpening source hS (half in the equilibrium)
            hS = Wa(k) * aC(i,j) * eN
            h(i,j,k) = h(i,j,k)*(1._wp-wc) + ( C(i,j)*Gam - 0.5_wp*hS )*wc + hS

            !== g: equilibrium and the stress of the non-equilibrium part
            geq(i,j,k) = Pn(i,j)*Wa(k) + (Gam - Wa(k))
            dg = g(i,j,k) - geq(i,j,k)
            sxx(i,j) = sxx(i,j) + dg*exr(k)*exr(k)
            sxy(i,j) = sxy(i,j) + dg*exr(k)*eyr(k)
            syy(i,j) = syy(i,j) + dg*eyr(k)*eyr(k)
        END DO
        END DO
    END DO

    CALL Total_Force( DcDx, DcDy, mu, Rho, Pn, tau, sxx, sxy, syy, Fx, Fy )

    !== forcing term gS, the equilibrium shifted by -gS/2, and the stress
    !   moments of g - geq (g holds the shifted equilibrium afterwards)
    m8 = 0._wp
    m9 = 0._wp
    DO k = 1, 9
        DO j = 1, ny
        DO i = 1, nx
            gS (i,j,k) = Wa(k) * ( exr(k)*Fx(i,j) + eyr(k)*Fy(i,j) ) * r3Rho(i,j)
            geq(i,j,k) = geq(i,j,k) - 0.5_wp*gS(i,j,k)
            m8(i,j) = m8(i,j) + w8(k)*( g(i,j,k) - geq(i,j,k) )
            m9(i,j) = m9(i,j) + w9(k)*( g(i,j,k) - geq(i,j,k) )
        END DO
        END DO
    END DO

    !== reduced MRT: the stress moments relax with s = 1/tau, all others with 1
    sm = 0.25_wp*( 1._wp - 1._wp/tau )        ! one division per cell for both moments
    m8 = m8 * sm
    m9 = m9 * sm
    DO k = 1, 9
        g(1:nx,1:ny,k) = geq(:,:,k) + w8(k)*m8 + w9(k)*m9 + gS(:,:,k)
    END DO

END SUBROUTINE Collision_h_g

!**********************************************************************

!== normalised pressure and velocity from g; the viscous force uses the
!   stress of g - geq with the velocity of the previous step (in Ux, Uy)
SUBROUTINE Moments_g( g, C, mu, DcDx, DcDy, P, Ux, Uy )
    REAL(wp), INTENT(IN)    :: g(0:nx+1,0:ny+1,9), C(-1:nx+2,-1:ny+2)
    REAL(wp), INTENT(IN), DIMENSION(nx,ny) :: mu, DcDx, DcDy    ! chemical potential, grad C
    REAL(wp), INTENT(OUT)   :: P(nx,ny)
    REAL(wp), INTENT(INOUT) :: Ux(nx,ny), Uy(nx,ny)

    INTEGER :: i, j, k
    REAL(wp) :: eU, Gam, dg
    REAL(wp), DIMENSION(nx,ny) :: Rho, hRho, Pn, tau, U2, &
                                  S0, Sx, Sy, sxx, sxy, syy, Fx, Fy

    !== zeroth and first moments, one fused loop over the directions
    !   (OMP SIMD on the cell loops of both direction sums: with the block size
    !   known at compile time ifx otherwise unrolls k and i completely and
    !   vectorizes j with gathers, ~30% slower collision; gfortran unaffected)
    S0 = 0._wp
    Sx = 0._wp
    Sy = 0._wp
    DO k = 1, 9
        DO j = 1, ny
!$OMP SIMD
        DO i = 1, nx
            S0(i,j) = S0(i,j) + g(i,j,k)
            Sx(i,j) = Sx(i,j) + g(i,j,k)*exr(k)
            Sy(i,j) = Sy(i,j) + g(i,j,k)*eyr(k)
        END DO
        END DO
    END DO

    DO j = 1, ny
    DO i = 1, nx
        Rho(i,j) = Rhol + C(i,j) * (Rhoh - Rhol)
        Pn (i,j) = S0(i,j)
        tau(i,j) = Relaxation_Time( C(i,j) )
        U2 (i,j) = Ux(i,j)*Ux(i,j) + Uy(i,j)*Uy(i,j)
    END DO
    END DO

    !== stress of g - geq (geq with the new pressure and the previous velocity)
    sxx = 0._wp
    sxy = 0._wp
    syy = 0._wp
    DO k = 1, 9
        DO j = 1, ny
!$OMP SIMD PRIVATE( eU, Gam, dg )
        DO i = 1, nx
            eU  = exr(k) * Ux(i,j) + eyr(k) * Uy(i,j)
            Gam = Wa(k) * ( 1._wp + eU*(3._wp + 4.5_wp*eU) - 1.5_wp*U2(i,j) )
            dg  = g(i,j,k) - ( Pn(i,j)*Wa(k) + (Gam - Wa(k)) )
            sxx(i,j) = sxx(i,j) + dg*exr(k)*exr(k)
            sxy(i,j) = sxy(i,j) + dg*exr(k)*eyr(k)
            syy(i,j) = syy(i,j) + dg*eyr(k)*eyr(k)
        END DO
        END DO
    END DO

    CALL Total_Force( DcDx, DcDy, mu, Rho, Pn, tau, sxx, sxy, syy, Fx, Fy )

    hRho = 0.5_wp/Rho                           ! one division per cell for both components
    Ux = Sx + Fx*hRho
    Uy = Sy + Fy*hRho
    P  = Pn*Rho*c13

END SUBROUTINE Moments_g

!**********************************************************************

!== surface tension, pressure, viscous and buoyancy force; the viscous force
!   -(tau-1/2) (S P) . grad C (Rhoh-Rhol) uses the non-equilibrium stress P
!   relaxed with the MRT rates: deviatoric part with s = 1/tau, trace with 1
SUBROUTINE Total_Force( DcDx, DcDy, mu, Rho, Pn, tau, sxx, sxy, syy, Fx, Fy )
    REAL(wp), INTENT(IN), DIMENSION(nx,ny)  :: DcDx, DcDy, mu, Rho, Pn, tau, sxx, sxy, syy
    REAL(wp), INTENT(OUT), DIMENSION(nx,ny) :: Fx, Fy

    INTEGER :: i, j
    REAL(wp) :: dRho3, s, txx, tyy, txy

    dRho3 = (Rhoh - Rhol)*c13

    DO j = 1, ny
    DO i = 1, nx
        s   = 1._wp/tau(i,j)
        txx = 0.5_wp*( sxx(i,j) + syy(i,j) + s*(sxx(i,j) - syy(i,j)) )
        tyy = 0.5_wp*( sxx(i,j) + syy(i,j) - s*(sxx(i,j) - syy(i,j)) )
        txy = s*sxy(i,j)
        Fx(i,j) = ( mu(i,j) - Pn(i,j)*dRho3 )*DcDx(i,j) &
                - (tau(i,j)-0.5_wp) * ( txx*DcDx(i,j) + txy*DcDy(i,j) ) * (Rhoh-Rhol)
        Fy(i,j) = ( mu(i,j) - Pn(i,j)*dRho3 )*DcDy(i,j) &
                - (tau(i,j)-0.5_wp) * ( txy*DcDx(i,j) + tyy*DcDy(i,j) ) * (Rhoh-Rhol) &
                + (Rhoh - Rho(i,j))*Gy
    END DO
    END DO

END SUBROUTINE Total_Force

!**********************************************************************

!== relaxation time (including the 1/2) of the hydrodynamics, Eqs. (24)-(25):
!   the dynamic viscosity mu is interpolated linearly in C between the light
!   (C=0) and heavy (C=1) fluid and tau - 1/2 = mu/(rho cs^2); with
!   taul, tauh = 3 mu/rho of the pure phases. C is clamped to [0,1].
ELEMENTAL FUNCTION Relaxation_Time( C ) RESULT( tau )
    REAL(wp), INTENT(IN) :: C
    REAL(wp) :: tau

    REAL(wp) :: Cc

    Cc  = MIN( MAX( C, 0._wp ), 1._wp )
    tau = ( Rhol*taul + Cc*(Rhoh*tauh - Rhol*taul) ) / ( Rhol + Cc*(Rhoh - Rhol) ) + 0.5_wp

END FUNCTION Relaxation_Time

END MODULE Model_PRE2017_Mod
