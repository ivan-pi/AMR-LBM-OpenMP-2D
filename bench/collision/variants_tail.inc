!**********************************************************************
!== the kernels of the reference model (the frozen copy Model_JCP2016_Ref_Mod)
!   and of the current model in src/, through their procedure pointers
SUBROUTINE V_new( n )
    INTEGER, INTENT(IN) :: n
    CALL new_model%Collide( n )
END SUBROUTINE V_new

SUBROUTINE V_mc_new( n )
    INTEGER, INTENT(IN) :: n
    CALL new_model%Moments_and_Collide( n )
END SUBROUTINE V_mc_new

SUBROUTINE V_model( n )
    INTEGER, INTENT(IN) :: n
    CALL the_model%Collide( n )
END SUBROUTINE V_model

SUBROUTINE V_mc_model( n )
    INTEGER, INTENT(IN) :: n
    CALL the_model%Moments_and_Collide( n )
END SUBROUTINE V_mc_model

!**********************************************************************
!== block size passed at run time (the original algebra)
SUBROUTINE V_rt_size( n )
    INTEGER, INTENT(IN) :: n
    REAL(wp) :: mu(nx,ny), inv_dx
    inv_dx = inv_dx_level( Blk(n)%Level )
    CALL Chemical_Potential_AMR( inv_dx, Blk_C(:,:,n), mu )
    CALL K_rt_size( nx_rt, ny_rt, inv_dx, Blk_C(:,:,n), mu, Blk_P(:,:,n), &
                    Blk_Ux(:,:,n), Blk_Uy(:,:,n), Blk_h(:,:,:,n), Blk_g(:,:,:,n) )
END SUBROUTINE V_rt_size

SUBROUTINE K_rt_size( nxr, nyr, inv_dx, C, mu, P, Ux, Uy, h, g )
    INTEGER,  INTENT(IN)    :: nxr, nyr
    REAL(wp), INTENT(IN)    :: inv_dx
    REAL(wp), INTENT(IN)    :: C(-1:nxr+2,-1:nyr+2), mu(nxr,nyr)
    REAL(wp), INTENT(IN)    :: P (nxr,nyr), Ux(nxr,nyr), Uy(nxr,nyr)
    REAL(wp), INTENT(INOUT) :: h(0:nxr+1,0:nyr+1,9), g(0:nxr+1,0:nyr+1,9)
#include "ref_decl.inc"
#include "ref_setup.inc"
#include "acc_std.inc"
    DO j = 1, nyr
#include "ref_simd.inc"
    DO i = 1, nxr
#include "ref_cell_a.inc"
#include "state_load.inc"
#include "ref_cell_b.inc"
    END DO
    END DO
#include "acc_undef.inc"
END SUBROUTINE K_rt_size

!**********************************************************************
!== copy of Macroscopic_Values_g_AMR of the model (private there)
SUBROUTINE Moments_g_Copy( inv_dx, g, C, mu, P, Ux, Uy )
    REAL(wp), INTENT(IN)  :: inv_dx
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
        Ux(i,j) = ( 3._wp*Sx + 0.5_wp*Fx )*rRho
        Uy(i,j) = ( 3._wp*Sy + 0.5_wp*Fy )*rRho
        P (i,j) = S0 + dRho6*( Ux(i,j)*DcDx(i,j) + Uy(i,j)*DcDy(i,j) )
    END DO
    END DO

END SUBROUTINE Moments_g_Copy

!**********************************************************************
!== What-if data layout (variants al_*): rows padded to ldp = nx+8 values and
!   placed so that the interior cells of every row of every field start at a
!   multiple of 64 bytes (a full vector never straddles two cache lines).
!   The fields live in a pool of their own, filled from / copied back to the
!   block pool by the benchmark.
SUBROUTINE Al_Setup( nblk )
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LOC, C_INTPTR_T
    INTEGER, INTENT(IN) :: nblk

    ALLOCATE( A_h(sz_h*nblk+8), A_g(sz_h*nblk+8), A_C(sz_C*nblk+8), A_P(sz_P*nblk+8), A_Ux(sz_P*nblk+8), A_Uy(sz_P*nblk+8) )
    o_h  = First( A_h, 1 )      ! the dummy h starts at i = 0: one value before the aligned one
    o_g  = First( A_g, 1 )
    o_C  = First( A_C, 2 )      ! C starts at i = -1
    o_P  = First( A_P, 0 )
    o_Ux = First( A_Ux, 0 )
    o_Uy = First( A_Uy, 0 )

CONTAINS
    !== first index k with a(k+shift) on a 64-byte boundary
    INTEGER FUNCTION First( a, shift )
        REAL(wp), INTENT(IN), TARGET :: a(:)
        INTEGER,  INTENT(IN) :: shift
        DO First = 1, 64
            IF( MOD( TRANSFER( C_LOC( a(First+shift) ), 0_C_INTPTR_T ), 64_C_INTPTR_T ) == 0 ) RETURN
        END DO
    END FUNCTION First
END SUBROUTINE Al_Setup

SUBROUTINE Al_Put( n, hs, gs, Ps, Uxs, Uys )
    INTEGER,  INTENT(IN) :: n
    REAL(wp), INTENT(IN) :: hs(0:nx+1,0:ny+1,9), gs(0:nx+1,0:ny+1,9), Ps(nx,ny), Uxs(nx,ny), Uys(nx,ny)
    CALL Put( A_h(o_h+(n-1)*sz_h), A_g(o_g+(n-1)*sz_h), A_C(o_C+(n-1)*sz_C), Blk_C(:,:,n) )
    A_P (o_P +(n-1)*sz_P:o_P +n*sz_P-1) = RESHAPE( Ps,  [sz_P] )
    A_Ux(o_Ux+(n-1)*sz_P:o_Ux+n*sz_P-1) = RESHAPE( Uxs, [sz_P] )
    A_Uy(o_Uy+(n-1)*sz_P:o_Uy+n*sz_P-1) = RESHAPE( Uys, [sz_P] )
CONTAINS
    SUBROUTINE Put( h, g, C, Cs )
        REAL(wp), INTENT(OUT) :: h(0:ldp-1,0:ny+1,9), g(0:ldp-1,0:ny+1,9), C(-1:ldp-2,-1:ny+2)
        REAL(wp), INTENT(IN)  :: Cs(-1:nx+2,-1:ny+2)
        h = 0;  g = 0;  C = 0
        h(0:nx+1,:,:) = hs
        g(0:nx+1,:,:) = gs
        C(-1:nx+2,:)  = Cs
    END SUBROUTINE Put
END SUBROUTINE Al_Put

!== copy block n back to the block pool
SUBROUTINE Al_Get( n )
    INTEGER, INTENT(IN) :: n
    CALL Get( A_h(o_h+(n-1)*sz_h), A_g(o_g+(n-1)*sz_h) )
    Blk_P (:,:,n) = RESHAPE( A_P (o_P +(n-1)*sz_P:o_P +n*sz_P-1), [nx,ny] )
    Blk_Ux(:,:,n) = RESHAPE( A_Ux(o_Ux+(n-1)*sz_P:o_Ux+n*sz_P-1), [nx,ny] )
    Blk_Uy(:,:,n) = RESHAPE( A_Uy(o_Uy+(n-1)*sz_P:o_Uy+n*sz_P-1), [nx,ny] )
CONTAINS
    SUBROUTINE Get( h, g )
        REAL(wp), INTENT(IN) :: h(0:ldp-1,0:ny+1,9), g(0:ldp-1,0:ny+1,9)
        Blk_h(:,:,:,n) = h(0:nx+1,:,:)
        Blk_g(:,:,:,n) = g(0:nx+1,:,:)
    END SUBROUTINE Get
END SUBROUTINE Al_Get

!**********************************************************************
!== the kernel written with AVX-512 intrinsics (collision_avx512.c): exact
!   divisions (c512) or reciprocal estimates with Newton steps (c512_rcp)
#define C512_WRAP(VN, KN) \
SUBROUTINE VN( n ); \
    INTEGER, INTENT(IN) :: n; \
    CALL KN( inv_dx_level( Blk(n)%Level ), c512_par, Blk_C(:,:,n), Blk_P(:,:,n), \
             Blk_Ux(:,:,n), Blk_Uy(:,:,n), Blk_h(:,:,:,n), Blk_g(:,:,:,n) ); \
END SUBROUTINE VN
C512_WRAP( V_c512,        collision_jcp2016_avx512 )
C512_WRAP( V_mc_c512,     moments_collision_jcp2016_avx512 )
C512_WRAP( V_c512_rcp,    collision_jcp2016_avx512_rcp )
C512_WRAP( V_mc_c512_rcp, moments_collision_jcp2016_avx512_rcp )
C512_WRAP( V_x_c512_nodiv,    collision_jcp2016_avx512_nodiv )
C512_WRAP( V_mc_x_c512_nodiv, moments_collision_jcp2016_avx512_nodiv )
#undef C512_WRAP

!== the physical parameters as the C kernel takes them (call once, after the setup)
SUBROUTINE C512_Setup
    c512_par = [ Rhol, Rhoh, Gy, taul, tauh, M0, W, wc, Beta, kappa, grad2_min ]
END SUBROUTINE C512_Setup
#include "moments_variants.inc"
