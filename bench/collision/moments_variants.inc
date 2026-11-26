!**********************************************************************
!== Moments of g alone (pressure, velocity): Moments_g_Block of the model,
!   called on steps with diagnostics and for new blocks. Variants mg_*:
!     mg_model      the reference: chemical potential pass, gradient pass, moments loop
!     mg_new        the model in src/ now
!     mg_fuse_grad  gradient in the moments loop (mu from its own pass)
!     mg_fused      chemical potential and gradient in the moments loop: one pass
!     mg_fused_nosimd   without the SIMD directive
!     mg_fused_sym  pair sums and differences of g shared between the moments
!     mg_x_nodiv    1/rho replaced by a multiply-add (wrong results: cost of the division)
SUBROUTINE V_mg_model( n )
    INTEGER, INTENT(IN) :: n
    CALL the_model%Moments_g( n )
END SUBROUTINE V_mg_model

SUBROUTINE V_mg_new( n )
    INTEGER, INTENT(IN) :: n
    CALL new_model%Moments_g( n )
END SUBROUTINE V_mg_new

#define MG_V V_mg_fuse_grad
#define MG_K K_mg_fuse_grad
#define MG_MU_PASS
#include "moments_kernel.inc"

#define MG_V V_mg_fused
#define MG_K K_mg_fused
#include "moments_kernel.inc"

#define MG_V V_mg_fused_nosimd
#define MG_K K_mg_fused_nosimd
#define MG_NOSIMD
#include "moments_kernel.inc"

#define MG_V V_mg_fused_sym
#define MG_K K_mg_fused_sym
#define MG_SYM
#include "moments_kernel.inc"

#define MG_V V_mg_x_nodiv
#define MG_K K_mg_x_nodiv
#define MG_SYM
#define MG_NODIV
#include "moments_kernel.inc"
