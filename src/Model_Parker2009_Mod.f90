!================================================
!   Multicomponent LBM with chemical-potential driven diffusion
!
!   J. M. Parker, PhD thesis, Oregon State University (2008);
!   J. M. Parker, G. N. Jovanovic, Chem. Biochem. Eng. Q. 23 (2009) 399-409.
!   Analysis and the surface-tension relation: docs/parker_model/README.md
!
!   Two species on D2Q9: species 1 in h (C = sum h), species 2 in g
!   (C2 = sum g). Each species relaxes towards the equilibrium of a fluid
!   moving with u + u'_i, u the mixture velocity and
!       u'_i = -D grad ln gamma_i  -  sum_j x_j (-D grad ln gamma_j)
!   (mass-average reference), ln gamma_i from Margules' equation. The u'u'
!   terms of the second moments are scaled by kappa: the mixture carries the
!   capillary stress T = kappa sum_i c_i u'_i u'_i and the surface tension is
!   sigma = kappa int sum_i c_i u'_(i,n)^2 dn.
!
!   Collision in central Hermite moments about u (Galilean invariant):
!   a0 = c_i, a1 = c_i u'_i, a2 = kappa c_i u'_i u'_i, a3 = a4 = 0
!   (Parker's equilibrium at rest); rates: a1, a3 omega_D, deviatoric a2
!   omega_nu, trace of a2 and a4 omega_E. D = cs^2 (1/omega_D - 1/2),
!   nu = cs^2 (1/omega_nu - 1/2).
!
!   Blk_P holds the pressure c cs^2, Blk_Ux/Uy the mixture velocity.
!   Derivatives are scaled with the block spacing (finest-level units); all
!   levels use the same time step, so D, nu and kappa are level-independent.
!================================================
MODULE Model_Parker2009_Mod
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: nx, ny, inv_dx_level
    USE D2Q9_Mod, ONLY: Wa, exr, eyr
    USE Model_Mod, ONLY: LBM_Model
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Parker2009_Model, Diff, kappa_s, nu, omega_E, A12, A21, Ln_Gamma, x_eps

    !== parameters (lattice units), set by the case before the initialization
    REAL(wp) :: Diff    = 0.2_wp        ! diffusivity D
    REAL(wp) :: kappa_s = 1._wp         ! scaling of the u'u' terms (surface tension)
    REAL(wp) :: nu      = 0.02_wp       ! kinematic viscosity
    REAL(wp) :: omega_E = 1._wp         ! rate of the trace of the stress and of a4
    REAL(wp) :: A12 = 2.7627_wp, A21 = 2.4474_wp    ! Margules (Parker's periodic.f90)

    !== mole fractions are clamped to [x_eps, 1-x_eps] in ln gamma
    REAL(wp), PARAMETER :: x_eps = 1.e-9_wp
    REAL(wp), PARAMETER :: cs2 = 1._wp/3._wp, cs4 = 1._wp/9._wp
    !== weights of the isotropic gradient: straight and diagonal neighbours
    REAL(wp), PARAMETER :: g1 = 1._wp/3._wp, g2 = 1._wp/12._wp

CONTAINS

FUNCTION Parker2009_Model() RESULT( model )
    TYPE(LBM_Model) :: model

    model%name       =  'two-species LBM, chemical-potential diffusion (Parker & Jovanovic 2009)'
    model%ghost_C    =  1
    model%n_scalars  =  2
    model%Initialize => Initialize_Block
    model%Collide    => Collide_Block
    model%Moments_g  => Moments_Block
    model%Moments_and_Collide => Moments_and_Collide_Block
    model%Refinement_Error    => Activity_Error
    model%refine_tol =  0.02_wp
    model%deref_tol  =  0.01_wp
    model%max_island_blocks = 16    ! droplet interiors: no coarse/fine edges inside a phase

END FUNCTION Parker2009_Model

!**********************************************************************

!== ln gamma_1, ln gamma_2 (Margules) of the mole fraction x1
ELEMENTAL SUBROUTINE Ln_Gamma( x1, lg1, lg2 )
    REAL(wp), INTENT(IN)  :: x1
    REAL(wp), INTENT(OUT) :: lg1, lg2

    REAL(wp) :: xa, xb

    xa = MIN( MAX( x1, x_eps ), 1._wp - x_eps )
    xb = 1._wp - xa
    lg1 = ( A12 + 2._wp*(A21 - A12)*xa ) * xb*xb
    lg2 = ( A21 + 2._wp*(A12 - A21)*xb ) * xa*xa

END SUBROUTINE Ln_Gamma

!**********************************************************************

!== ln gamma of both species on the block with one ghost layer
SUBROUTINE Block_Ln_Gamma( C1, C2, lg1, lg2 )
    REAL(wp), INTENT(IN)  :: C1(-1:nx+2,-1:ny+2), C2(-1:nx+2,-1:ny+2)
    REAL(wp), INTENT(OUT) :: lg1(0:nx+1,0:ny+1), lg2(0:nx+1,0:ny+1)

    CALL Ln_Gamma( C1(0:nx+1,0:ny+1) / ( C1(0:nx+1,0:ny+1) + C2(0:nx+1,0:ny+1) ), lg1, lg2 )

END SUBROUTINE Block_Ln_Gamma

!**********************************************************************

!== isotropic gradient of a field with one ghost layer, finest-level units
PURE SUBROUTINE Gradient( inv_dx, q, i, j, qx, qy )
    REAL(wp), INTENT(IN)  :: inv_dx, q(0:nx+1,0:ny+1)
    INTEGER,  INTENT(IN)  :: i, j
    REAL(wp), INTENT(OUT) :: qx, qy

    qx = inv_dx*( g1*( q(i+1,j) - q(i-1,j) ) &
                + g2*( q(i+1,j+1) + q(i+1,j-1) - q(i-1,j+1) - q(i-1,j-1) ) )
    qy = inv_dx*( g1*( q(i,j+1) - q(i,j-1) ) &
                + g2*( q(i+1,j+1) + q(i-1,j+1) - q(i+1,j-1) - q(i-1,j-1) ) )

END SUBROUTINE Gradient

!**********************************************************************

!== refinement criterion: the largest |grad ln gamma_i| of block b (finest-level
!   units). It is zero in the bulk phases and largest in the interfaces (at
!   equilibrium it equals |grad ln x_i|); grad mu_i = grad ln a_i is not used,
!   as it vanishes across an equilibrium interface.
FUNCTION Activity_Error( b ) RESULT( err )
    USE Blocks_Mod, ONLY: Blk, Blk_C, Blk_C2
    INTEGER, INTENT(IN) :: b
    REAL(wp) :: err

    INTEGER  :: i, j
    REAL(wp) :: inv_dx, ax, ay, bx, by
    REAL(wp), DIMENSION(0:nx+1,0:ny+1) :: lg1, lg2

    inv_dx = inv_dx_level( Blk(b)%Level )
    CALL Block_Ln_Gamma( Blk_C(:,:,b), Blk_C2(:,:,b), lg1, lg2 )
    err = 0
    DO j = 1, ny
    DO i = 1, nx
        CALL Gradient( inv_dx, lg1, i, j, ax, ay )
        CALL Gradient( inv_dx, lg2, i, j, bx, by )
        err = MAX( err, ax*ax + ay*ay, bx*bx + by*by )
    END DO
    END DO
    err = SQRT( err )

END FUNCTION Activity_Error

!**********************************************************************

!== mixture velocity and pressure c cs^2 of block n
SUBROUTINE Moments_Block( n )
    USE Blocks_Mod, ONLY: Blk_h, Blk_g, Blk_C, Blk_C2, Blk_P, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    CALL Mixture_Moments( Blk_h(:,:,:,n), Blk_g(:,:,:,n), Blk_C(:,:,n), Blk_C2(:,:,n), &
                          Blk_P(:,:,n), Blk_Ux(:,:,n), Blk_Uy(:,:,n) )

END SUBROUTINE Moments_Block

SUBROUTINE Mixture_Moments( h, g, C1, C2, P, Ux, Uy )
    REAL(wp), INTENT(IN)  :: h(0:nx+1,0:ny+1,9), g(0:nx+1,0:ny+1,9)
    REAL(wp), INTENT(IN)  :: C1(-1:nx+2,-1:ny+2), C2(-1:nx+2,-1:ny+2)
    REAL(wp), INTENT(OUT) :: P(nx,ny), Ux(nx,ny), Uy(nx,ny)

    INTEGER  :: i, j, k
    REAL(wp) :: c, jx, jy

    DO j = 1, ny
    DO i = 1, nx
        c  = C1(i,j) + C2(i,j)
        jx = 0
        jy = 0
        DO k = 2, 9
            jx = jx + exr(k)*( h(i,j,k) + g(i,j,k) )
            jy = jy + eyr(k)*( h(i,j,k) + g(i,j,k) )
        END DO
        Ux(i,j) = jx / c
        Uy(i,j) = jy / c
        P (i,j) = c * cs2
    END DO
    END DO

END SUBROUTINE Mixture_Moments

!**********************************************************************

SUBROUTINE Collide_Block( n )
    USE Blocks_Mod, ONLY: Blk, Blk_h, Blk_g, Blk_C, Blk_C2, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    CALL Collision( inv_dx_level(Blk(n)%Level), Blk_C(:,:,n), Blk_C2(:,:,n), &
                    Blk_Ux(:,:,n), Blk_Uy(:,:,n), Blk_h(:,:,:,n), Blk_g(:,:,:,n) )

END SUBROUTINE Collide_Block

SUBROUTINE Moments_and_Collide_Block( n )
    INTEGER, INTENT(IN) :: n

    CALL Moments_Block( n )
    CALL Collide_Block( n )

END SUBROUTINE Moments_and_Collide_Block

!**********************************************************************

!== collision of both species of one block
SUBROUTINE Collision( inv_dx, C1, C2, Ux, Uy, h, g )
    REAL(wp), INTENT(IN)    :: inv_dx
    REAL(wp), INTENT(IN)    :: C1(-1:nx+2,-1:ny+2), C2(-1:nx+2,-1:ny+2)
    REAL(wp), INTENT(IN)    :: Ux(nx,ny), Uy(nx,ny)
    REAL(wp), INTENT(INOUT) :: h(0:nx+1,0:ny+1,9), g(0:nx+1,0:ny+1,9)

    INTEGER  :: i, j
    REAL(wp) :: wD, wV, c1c, c2c, x1, x2, p1x, p1y, p2x, p2y, qx, qy
    REAL(wp), DIMENSION(0:nx+1,0:ny+1) :: lg1, lg2

    wD = 1._wp / ( Diff/cs2 + 0.5_wp )
    wV = 1._wp / ( nu/cs2 + 0.5_wp )

    CALL Block_Ln_Gamma( C1, C2, lg1, lg2 )

    DO j = 1, ny
    DO i = 1, nx
        c1c = C1(i,j)
        c2c = C2(i,j)
        x1  = c1c / ( c1c + c2c )
        x2  = 1._wp - x1

        !== diffusion velocities, mass-average reference
        CALL Gradient( inv_dx, lg1, i, j, p1x, p1y )
        CALL Gradient( inv_dx, lg2, i, j, p2x, p2y )
        p1x = -Diff*p1x;  p1y = -Diff*p1y
        p2x = -Diff*p2x;  p2y = -Diff*p2y
        qx = x1*p1x + x2*p2x
        qy = x1*p1y + x2*p2y

        CALL Collide_Species( h(i,j,:), c1c, Ux(i,j), Uy(i,j), p1x - qx, p1y - qy, wD, wV )
        CALL Collide_Species( g(i,j,:), c2c, Ux(i,j), Uy(i,j), p2x - qx, p2y - qy, wD, wV )
    END DO
    END DO

END SUBROUTINE Collision

!**********************************************************************

!== central Hermite moment collision of one species at one cell
PURE SUBROUTINE Collide_Species( f, c, ux, uy, px, py, wD, wV )
    REAL(wp), INTENT(INOUT) :: f(9)
    REAL(wp), INTENT(IN)    :: c, ux, uy, px, py, wD, wV

    INTEGER  :: k
    REAL(wp) :: dx, dy, dx2, dy2, fk
    REAL(wp) :: k00, k10, k01, k20, k02, k11, k21, k12, k22
    REAL(wp) :: h10, h01, h20, h02, h11, h21, h12, h22, tr, dv, trq, dvq
    REAL(wp) :: m10, m01, m20, m02, m11, m21, m12, m22

    !== central moments about u
    k00 = 0; k10 = 0; k01 = 0; k20 = 0; k02 = 0; k11 = 0; k21 = 0; k12 = 0; k22 = 0
    DO k = 1, 9
        fk  = f(k)
        dx  = exr(k) - ux
        dy  = eyr(k) - uy
        dx2 = dx*dx
        dy2 = dy*dy
        k00 = k00 + fk
        k10 = k10 + fk*dx
        k01 = k01 + fk*dy
        k20 = k20 + fk*dx2
        k02 = k02 + fk*dy2
        k11 = k11 + fk*dx*dy
        k21 = k21 + fk*dx2*dy
        k12 = k12 + fk*dx*dy2
        k22 = k22 + fk*dx2*dy2
    END DO

    !== central Hermite moments, relaxed (c = k00 is conserved)
    h10 = k10 + wD*( c*px - k10 )
    h01 = k01 + wD*( c*py - k01 )
    h21 = (1._wp - wD)*( k21 - cs2*k01 )
    h12 = (1._wp - wD)*( k12 - cs2*k10 )
    h11 = k11 + wV*( kappa_s*c*px*py - k11 )
    tr  = k20 + k02 - 2._wp*cs2*c
    dv  = k20 - k02
    trq = kappa_s*c*( px*px + py*py )
    dvq = kappa_s*c*( px*px - py*py )
    tr  = tr + omega_E*( trq - tr )
    dv  = dv + wV*( dvq - dv )
    h20 = 0.5_wp*( tr + dv )
    h02 = 0.5_wp*( tr - dv )
    h22 = (1._wp - omega_E)*( k22 - cs2*(k20 + k02) + cs4*c )

    !== back to central moments
    k10 = h10
    k01 = h01
    k20 = h20 + cs2*c
    k02 = h02 + cs2*c
    k11 = h11
    k21 = h21 + cs2*k01
    k12 = h12 + cs2*k10
    k22 = h22 + cs2*( k20 + k02 ) - cs4*c

    !== raw moments (binomial shift by u)
    m10 = k10 + ux*c
    m01 = k01 + uy*c
    m20 = k20 + 2._wp*ux*k10 + ux*ux*c
    m02 = k02 + 2._wp*uy*k01 + uy*uy*c
    m11 = k11 + ux*k01 + uy*k10 + ux*uy*c
    m21 = k21 + 2._wp*ux*k11 + uy*k20 + ux*ux*k01 + 2._wp*ux*uy*k10 + ux*ux*uy*c
    m12 = k12 + 2._wp*uy*k11 + ux*k02 + uy*uy*k10 + 2._wp*ux*uy*k01 + ux*uy*uy*c
    m22 = k22 + 2._wp*ux*k12 + 2._wp*uy*k21 + ux*ux*k02 + uy*uy*k20 + 4._wp*ux*uy*k11 &
        + 2._wp*ux*uy*uy*k10 + 2._wp*ux*ux*uy*k01 + ux*ux*uy*uy*c

    !== populations
    f(1) = c - m20 - m02 + m22
    f(2) = 0.5_wp *(  m10 + m20 - m12 - m22 )
    f(3) = 0.5_wp *(  m01 + m02 - m21 - m22 )
    f(4) = 0.5_wp *( -m10 + m20 + m12 - m22 )
    f(5) = 0.5_wp *( -m01 + m02 + m21 - m22 )
    f(6) = 0.25_wp*(  m11 + m21 + m12 + m22 )
    f(7) = 0.25_wp*( -m11 + m21 - m12 + m22 )
    f(8) = 0.25_wp*(  m11 - m21 - m12 + m22 )
    f(9) = 0.25_wp*( -m11 - m21 + m12 + m22 )

END SUBROUTINE Collide_Species

!**********************************************************************

!== h and g of block n: equilibria of both species for C, C2, Ux, Uy (u' = 0)
SUBROUTINE Initialize_Block( n )
    USE Blocks_Mod, ONLY: Blk_h, Blk_g, Blk_C, Blk_C2, Blk_P, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    INTEGER  :: i, j, k
    REAL(wp) :: eu, uu, G

    DO k = 1, 9
        DO j = 1, ny
        DO i = 1, nx
            eu = exr(k)*Blk_Ux(i,j,n) + eyr(k)*Blk_Uy(i,j,n)
            uu = Blk_Ux(i,j,n)**2 + Blk_Uy(i,j,n)**2
            G  = Wa(k)*( 1._wp + 3._wp*eu + 4.5_wp*eu*eu - 1.5_wp*uu )
            Blk_h(i,j,k,n) = Blk_C (i,j,n)*G
            Blk_g(i,j,k,n) = Blk_C2(i,j,n)*G
        END DO
        END DO
    END DO
    Blk_P(:,:,n) = cs2*( Blk_C(1:nx,1:ny,n) + Blk_C2(1:nx,1:ny,n) )

END SUBROUTINE Initialize_Block

END MODULE Model_Parker2009_Mod
