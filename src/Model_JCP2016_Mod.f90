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
!
!   The collision kernel is the template Collision_JCP2016.inc; its form is
!   the outcome of the measurements in
!   docs/benchmarks/2026-11-26-kernel-performance.md.
!================================================
MODULE Model_JCP2016_Mod
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: nx, ny, inv_dx_level     ! the fixed block size, 1/dx per level
    USE D2Q9_Mod, ONLY: Wa, ex, ey, exr, eyr, w0, ws, wd
    USE Physics_Mod, ONLY: W, M0, wc, Gy, taul, tauh, Rhol, Rhoh, Beta, kappa
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

    CALL Collision_h_g_AMR( inv_dx_level( Blk(n)%Level ), Blk_C(:,:,n), Blk_P(:,:,n), &
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

!== moments of g and collision of one block in one pass over its cells
SUBROUTINE Moments_and_Collide_Block( n )
    USE Blocks_Mod, ONLY: Blk, Blk_h, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    CALL Moments_Collision_h_g_AMR( inv_dx_level( Blk(n)%Level ), Blk_C(:,:,n), Blk_P(:,:,n), &
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

!== the collision kernel (src/Collision_JCP2016.inc), with P, Ux, Uy as input
!   and with the moments of g computed in the same pass
#define KERNEL_NAME Collision_h_g_AMR
#include "Collision_JCP2016.inc"

!**********************************************************************

#define KERNEL_NAME Moments_Collision_h_g_AMR
#define WITH_MOMENTS
#include "Collision_JCP2016.inc"

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
