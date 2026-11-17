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
    USE Config_Mod, ONLY: nx, ny     ! the fixed block size
    USE D2Q9_Mod, ONLY: Wa, ex, ey
    USE Physics_Mod, ONLY: W, M0, wc, Gy, taul, tauh, Rhol, Rhoh
    USE Model_Mod, ONLY: LBM_Model
    USE Phase_Field_Mod, ONLY: Chemical_Potential_AMR, Gradient_AMR, Max_Gradient
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: JCP2016_Model

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

!== refinement criterion: the interface, i.e. the largest |grad C| of block b
FUNCTION Interface_Error( b ) RESULT( err )
    USE Blocks_Mod, ONLY: Blk, Blk_C
    INTEGER, INTENT(IN) :: b
    REAL(wp) :: err

    err = Max_Gradient( Blk(b)%dx, Blk_C(:,:,b) )

END FUNCTION Interface_Error

!**********************************************************************

SUBROUTINE Collide_Block( n )
    USE Blocks_Mod, ONLY: Blk, Blk_h, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    REAL(wp) :: mu(nx,ny)

    CALL Chemical_Potential_AMR( Blk(n)%dx, Blk_C(:,:,n), mu )
    CALL Collision_h_g_AMR( Blk(n)%dx, Blk_C(:,:,n), mu, Blk_P(:,:,n), &
                            Blk_Ux(:,:,n), Blk_Uy(:,:,n), Blk_h(:,:,:,n), Blk_g(:,:,:,n) )

END SUBROUTINE Collide_Block

!**********************************************************************

SUBROUTINE Moments_g_Block( n )
    USE Blocks_Mod, ONLY: Blk, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    REAL(wp) :: mu(nx,ny)

    CALL Chemical_Potential_AMR( Blk(n)%dx, Blk_C(:,:,n), mu )
    CALL Macroscopic_Values_g_AMR( Blk(n)%dx, Blk_g(:,:,:,n), Blk_C(:,:,n), mu, &
                                   Blk_P(:,:,n), Blk_Ux(:,:,n), Blk_Uy(:,:,n) )

END SUBROUTINE Moments_g_Block

!**********************************************************************

!== moments of g and collision of one block with one chemical potential
SUBROUTINE Moments_and_Collide_Block( n )
    USE Blocks_Mod, ONLY: Blk, Blk_h, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    REAL(wp) :: mu(nx,ny)

    CALL Chemical_Potential_AMR( Blk(n)%dx, Blk_C(:,:,n), mu )
    CALL Macroscopic_Values_g_AMR( Blk(n)%dx, Blk_g(:,:,:,n), Blk_C(:,:,n), mu, &
                                   Blk_P(:,:,n), Blk_Ux(:,:,n), Blk_Uy(:,:,n) )
    CALL Collision_h_g_AMR( Blk(n)%dx, Blk_C(:,:,n), mu, Blk_P(:,:,n), &
                            Blk_Ux(:,:,n), Blk_Uy(:,:,n), Blk_h(:,:,:,n), Blk_g(:,:,:,n) )

END SUBROUTINE Moments_and_Collide_Block

!**********************************************************************

!== equilibrium h and g of block n from C, P, Ux, Uy (no body force)
SUBROUTINE Initialize_Block( n )
    USE Blocks_Mod, ONLY: Blk, Blk_h, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    INTEGER :: dx, i, j
    REAL(wp) :: Fb, Rho, dRho3, ni, nj, DcDx, DcDy
    REAL(wp) :: eDc(9), hlp(9), Gamma(9), mu(nx,ny)

    dx    = Blk(n)%dx
    dRho3 = (Rhoh - Rhol)/3._wp

    CALL Chemical_Potential_AMR( dx, Blk_C(:,:,n), mu )

    DO j = 1, ny
    DO i = 1, nx

        CALL ED_C_AMR( dx, Blk_C(i-1:i+1,j-1:j+1, n), eDc(:) )

        CALL Gradient_ED_AMR( eDc(:), DcDx, DcDy )

        CALL normal_AMR( DcDx, DcDy, ni, nj )

        Rho = Rhol + Blk_C(i,j, n) * (Rhoh - Rhol)

        Fb = 0  !(Rhoh - Rho) * Gy

        CALL Equilibrium_h_g_AMR( Blk_C(i,j, n), Blk_Ux(i,j, n), Blk_Uy(i,j, n), &
                                  ni, nj, Gamma(:), Blk_h(i,j,:, n) )

        hlp(:) = ( (Gamma(:)-Wa(:))*dRho3 + Gamma(:)*mu(i,j) ) * ( eDc(:) &
                - (Blk_Ux(i,j, n)*DcDx + Blk_Uy(i,j, n)*DcDy) ) + (ey(:)-Blk_Uy(i,j, n))*Gamma(:) * Fb

        Blk_g(i,j,:, n) = Blk_P(i,j, n)*Wa(:) &
                         + (Gamma(:)-Wa(:))*Rho/3._wp - 0.5_wp*hlp(:)

    END DO
    END DO

CONTAINS

    SUBROUTINE Equilibrium_h_g_AMR( C, U, V, ni, nj, Gamma, heq )
        REAL(wp), INTENT(IN)  :: C, U, V, ni, nj
        REAL(wp), INTENT(OUT) :: Gamma(9), heq(:)     ! assumed shape: h(i,j,:) is strided

        REAL(wp) :: U2, eU(9), eN(9), tmp

        tmp = ( 12*(1-C) )/W * M0

        U2 = U*U + V*V

        eU(:) = ex(:) * U  + ey(:) * V
        eN(:) = ex(:) * ni + ey(:) * nj

        Gamma(:) = Wa(:) * ( 1._wp + eU(:)*(3._wp + 4.5_wp*eU(:)) - 1.5_wp*U2 )

        heq(:) = C * ( Gamma(:) + tmp * Wa(:) * eN(:) )

    END SUBROUTINE Equilibrium_h_g_AMR

    SUBROUTINE normal_AMR( DcDx, DcDy, ni, nj )
        REAL(wp), INTENT(IN)  :: DcDx, DcDy
        REAL(wp), INTENT(OUT) :: ni, nj

        REAL(wp) :: tmp

        tmp = SQRT( DcDx**2 + DcDy**2 + 1.e-32_wp )

        ni = DcDx / tmp
        nj = DcDy / tmp

    END SUBROUTINE normal_AMR

    SUBROUTINE ED_C_AMR( dx, C, EdotDelC )
        INTEGER, INTENT(IN)  :: dx
        REAL(wp), INTENT(IN)  :: C(-1:,-1:)
        REAL(wp), INTENT(OUT) :: EdotDelC(9)

        INTEGER :: i

        DO i = 1, 9
            EdotDelC(i) = C(ex(i),ey(i)) - C(-ex(i),-ey(i))
        END DO
        EdotDelC(:) = EdotDelC(:) * 0.5_wp/REAL(dx,wp)

    END SUBROUTINE ED_C_AMR

    SUBROUTINE Gradient_ED_AMR( ED, Dx, Dy )
        REAL(wp),INTENT(IN) :: ED(9)
        REAL(wp),INTENT(OUT):: Dx, Dy

        Dx = ( ED(2)-ED(4) )/3._wp + ( ED(6)-ED(7)-ED(8)+ED(9) )/12._wp
        Dy = ( ED(3)-ED(5) )/3._wp + ( ED(6)+ED(7)-ED(8)-ED(9) )/12._wp

    END SUBROUTINE Gradient_ED_AMR

END SUBROUTINE Initialize_Block

!**********************************************************************

SUBROUTINE Collision_h_g_AMR( dx, C, mu, P, Ux, Uy, h, g )
    !== Collision of h and g for one block:
    !   1) block buffers: directional stencils, gradients, normals, rho, s9
    !   2) a direction loop with cells innermost does equilibria and relaxation.
    !   MRT: g <- g - M^-1 S M (g - geq) with S = diag(1,...,1,s9,s9). As all
    !   other moments relax with rate 1, this reduces to g = geq plus (1-s9)
    !   times the two stress moments m8, m9 of g - geq, mapped back as +-m8/4
    !   (directions 2-5) and +-m9/4 (directions 6-9).
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: dx
    REAL(wp), INTENT(IN)    :: C(-1:nx+2,-1:ny+2), mu(nx,ny)    ! mu: chemical potential
    REAL(wp), INTENT(IN)    :: P (nx,ny), Ux(nx,ny), Uy(nx,ny)
    REAL(wp), INTENT(INOUT) :: h(0:nx+1,0:ny+1,9), g(0:nx+1,0:ny+1,9)

    INTEGER :: i, j, k, ix, iy
    REAL(wp) :: dRho3, invTaul, dInvTau, tmp, eU, eN, Gam, A, hlpC, geq
    !== signs of g_k - geq_k in the stress moments m8 (k=2..5) and m9 (k=6..9)
    REAL(wp), PARAMETER :: w8(9) = [0, 1,-1, 1,-1, 0, 0, 0, 0]
    REAL(wp), PARAMETER :: w9(9) = [0, 0, 0, 0, 0, 1,-1, 1,-1]
    REAL(wp), DIMENSION(nx,ny,9) :: eC, eM, hM
    REAL(wp), DIMENSION(nx,ny) :: DcDx, DcDy, DmDx, DmDy, ni, nj, Rho, Fb, s9, &
                                  U2, uDc, uDm, tC, m8, m9

    dRho3   = (Rhoh - Rhol)/3._wp
    invTaul = 1._wp/taul                    ! 1/tau is linear in C between
    dInvTau = 1._wp/tauh - 1._wp/taul       ! the light (C=0) and heavy (C=1) fluid

    !== directional differences along e_k: central (eC) and mixed (eM);
    !   Cm, Cp, Cpp are C shifted by -e_k, +e_k, +2e_k, C0 the block interior
    eC(:,:,1) = 0._wp
    eM(:,:,1) = 0._wp
    DO k = 2, 9
        ix = ex(k); iy = ey(k)
        ASSOCIATE( C0  => C(1     :nx     ,1     :ny     ), &
                   Cm  => C(1-ix  :nx-ix  ,1-iy  :ny-iy  ), &
                   Cp  => C(1+ix  :nx+ix  ,1+iy  :ny+iy  ), &
                   Cpp => C(1+2*ix:nx+2*ix,1+2*iy:ny+2*iy) )
            eC(:,:,k) = ( Cp - Cm ) * 0.5_wp/REAL(dx,wp)
            eM(:,:,k) = ( - Cm - 3._wp*C0 + 5._wp*Cp - Cpp ) * 0.25_wp/REAL(dx,wp)
        END ASSOCIATE
    END DO

    !== per-cell quantities: gradients, normal, density, forcing, relaxation rate
    !   (one fused loop; as separate array statements this is ~3% slower)
    DO j = 1, ny
!$OMP SIMD PRIVATE( tmp )
    DO i = 1, nx
        DcDx(i,j) = ( eC(i,j,2)-eC(i,j,4) )/3._wp + ( eC(i,j,6)-eC(i,j,7)-eC(i,j,8)+eC(i,j,9) )/12._wp
        DcDy(i,j) = ( eC(i,j,3)-eC(i,j,5) )/3._wp + ( eC(i,j,6)+eC(i,j,7)-eC(i,j,8)-eC(i,j,9) )/12._wp
        DmDx(i,j) = ( eM(i,j,2)-eM(i,j,4) )/3._wp + ( eM(i,j,6)-eM(i,j,7)-eM(i,j,8)+eM(i,j,9) )/12._wp
        DmDy(i,j) = ( eM(i,j,3)-eM(i,j,5) )/3._wp + ( eM(i,j,6)+eM(i,j,7)-eM(i,j,8)-eM(i,j,9) )/12._wp

        tmp = SQRT( DcDx(i,j)**2 + DcDy(i,j)**2 + 1.e-32_wp )
        ni(i,j) = DcDx(i,j) / tmp
        nj(i,j) = DcDy(i,j) / tmp

        tC  (i,j) = ( 12*(1-C(i,j)) )/W * M0
        Rho (i,j) = Rhol + C(i,j) * (Rhoh - Rhol)
        Fb  (i,j) = (Rhoh - Rho(i,j)) * Gy
        U2  (i,j) = Ux(i,j)*Ux(i,j) + Uy(i,j)*Uy(i,j)
        uDc (i,j) = Ux(i,j)*DcDx(i,j) + Uy(i,j)*DcDy(i,j)
        uDm (i,j) = Ux(i,j)*DmDx(i,j) + Uy(i,j)*DmDy(i,j)

        !== relaxation rate s9 = 1/tau with tau = 1/x + 1/2, where x = 1/tau_fluid is
        !   interpolated in C (C clamped, as it may overshoot [0,1] at the interface)
        tmp = MIN( MAX( C(i,j), 0._wp ), 1._wp ) * dInvTau + invTaul
        s9(i,j) = tmp / ( 1._wp + 0.5_wp*tmp )
    END DO
    END DO

    !== equilibria, h relaxation, g equilibrium and forcing, direction by direction;
    !   g temporarily holds geq, hM the forcing term of the mixed difference
    m8 = 0._wp
    m9 = 0._wp
    DO k = 1, 9
        DO j = 1, ny
!$OMP SIMD PRIVATE( A, Gam, eN, eU, geq, hlpC )
        DO i = 1, nx
            eU  = ex(k) * Ux(i,j) + ey(k) * Uy(i,j)
            eN  = ex(k) * ni(i,j) + ey(k) * nj(i,j)
            Gam = Wa(k) * ( 1._wp + eU*(3._wp + 4.5_wp*eU) - 1.5_wp*U2(i,j) )

            h(i,j,k) = h(i,j,k) * (1._wp-wc) + ( C(i,j) * ( Gam + tC(i,j) * Wa(k) * eN ) )*wc

            A    = (Gam-Wa(k))*dRho3 + Gam*mu(i,j)
            hlpC = A * ( eC(i,j,k) - uDc(i,j) ) + (ey(k)-Uy(i,j)) * Gam * Fb(i,j)
            hM(i,j,k) = A * ( eM(i,j,k) - uDm(i,j) ) + (ey(k)-Uy(i,j)) * Gam * Fb(i,j)
            geq  = P(i,j)*Wa(k) + (Gam-Wa(k))*Rho(i,j)/3._wp - 0.5_wp * hlpC

            m8(i,j) = m8(i,j) + w8(k)*( g(i,j,k) - geq )
            m9(i,j) = m9(i,j) + w9(k)*( g(i,j,k) - geq )
            g(i,j,k) = geq
        END DO
        END DO
    END DO

    !== relax the stress moments and add the forcing term
    m8 = m8 * (1._wp-s9) * 0.25_wp
    m9 = m9 * (1._wp-s9) * 0.25_wp
    DO k = 1, 9
        g(1:nx,1:ny,k) = ( g(1:nx,1:ny,k) + w8(k)*m8 + w9(k)*m9 ) + hM(:,:,k)
    END DO

END SUBROUTINE Collision_h_g_AMR

!**********************************************************************

SUBROUTINE Macroscopic_Values_g_AMR( dx, g, C, mu, P, Ux, Uy )
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: dx
    REAL(wp), INTENT(IN)  :: g(0:nx+1,0:ny+1,9), C(-1:nx+2,-1:ny+2), mu(nx,ny)
    REAL(wp), INTENT(OUT) :: P (nx,ny), Ux(nx,ny), Uy(nx,ny)

    INTEGER :: i, j, k
    REAL(wp) :: Rho, dRho6, Fx, Fy, Fb
    REAL(wp) :: DcDx(nx,ny), DcDy(nx,ny)
    REAL(wp) :: Sx(nx,ny), Sy(nx,ny), S0(nx,ny)

    dRho6 = (Rhoh-Rhol)/6._wp

    CALL Gradient_AMR( dx, C, DcDx, DcDy )

    !== zeroth and first moments, accumulated direction by direction in one fused
    !   loop (as three array statements this is ~8% slower)
    DO j = 1, ny
    DO i = 1, nx
        Sx(i,j) = g(i,j,2)*ex(2)
        Sy(i,j) = g(i,j,2)*ey(2)
        S0(i,j) = g(i,j,1) + g(i,j,2)
    END DO
    END DO
    DO k = 3, 9
        DO j = 1, ny
        DO i = 1, nx
            Sx(i,j) = Sx(i,j) + g(i,j,k)*ex(k)
            Sy(i,j) = Sy(i,j) + g(i,j,k)*ey(k)
            S0(i,j) = S0(i,j) + g(i,j,k)
        END DO
        END DO
    END DO

    DO j = 1, ny
    DO i = 1, nx

        Rho = Rhol + C(i,j) * (Rhoh - Rhol)

        Fb = (Rhoh - Rho) * Gy

        Fx = mu(i,j) * DcDx(i,j)
        Fy = mu(i,j) * DcDy(i,j) + Fb

        !== u = ( sum(g e)/cs^2 + F/2 )/rho, from rho u = sum(f e) + F/2 and
        !   sum(g e) = cs^2 sum(f e). (Eq. (35) of the paper prints the force term
        !   without the factor cs^2; this form is consistent with its Eq. (24).)
        Ux(i,j) = ( 3._wp*Sx(i,j) + 0.5_wp*Fx )/Rho
        Uy(i,j) = ( 3._wp*Sy(i,j) + 0.5_wp*Fy )/Rho

        P (i,j) = S0(i,j) + dRho6*( Ux(i,j)*DcDx(i,j) + Uy(i,j)*DcDy(i,j) )

    END DO
    END DO

END SUBROUTINE Macroscopic_Values_g_AMR

END MODULE Model_JCP2016_Mod
