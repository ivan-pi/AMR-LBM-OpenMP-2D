!================================================
!   An LBM model as a set of block kernels
!
!   The drivers select a model; the time loop, the mesh adaptation and
!   the initialization call the model only through these kernels, which
!   act on the leaf block with the given index. The model also chooses
!   its refinement criterion, as the blocks it needs on the finest level
!   depend on the model.
!================================================
MODULE Model_Mod
    USE Kinds_Mod, ONLY: wp
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Block_Kernel, Block_Error, LBM_Model

    ABSTRACT INTERFACE
        SUBROUTINE Block_Kernel( n )
            INTEGER, INTENT(IN) :: n        ! block index
        END SUBROUTINE Block_Kernel

        !== refinement criterion: an error indicator of block b, compared
        !   with the refinement / derefinement tolerances
        FUNCTION Block_Error( b ) RESULT( err )
            IMPORT :: wp
            INTEGER, INTENT(IN) :: b        ! block index
            REAL(wp) :: err
        END FUNCTION Block_Error
    END INTERFACE

    TYPE :: LBM_Model
        CHARACTER(len=120) :: name = ''
        !== ghost layers of C read by Collide (filled before the collision)
        INTEGER :: ghost_C = 1
        !== scalars with ghost layers: 1 (C = sum h) or 2 (also C2 = sum g)
        INTEGER :: n_scalars = 1
        !== h and g from C (with ghost layers), P, Ux, Uy (initial state)
        PROCEDURE(Block_Kernel), POINTER, NOPASS :: Initialize => NULL()
        !== collision of h and g
        PROCEDURE(Block_Kernel), POINTER, NOPASS :: Collide    => NULL()
        !== pressure and velocity from g (needs C with one ghost layer)
        PROCEDURE(Block_Kernel), POINTER, NOPASS :: Moments_g  => NULL()
        !== Moments_g followed by Collide, sharing the derivatives of C
        PROCEDURE(Block_Kernel), POINTER, NOPASS :: Moments_and_Collide => NULL()
        !== refinement criterion: blocks with Refinement_Error > refine_tol are
        !   refined, groups of siblings all below deref_tol are derefined
        PROCEDURE(Block_Error),  POINTER, NOPASS :: Refinement_Error => NULL()
        REAL(wp) :: refine_tol = 0.002_wp
        REAL(wp) :: deref_tol  = 0.001_wp
        !== also refine coarse islands of at most this many leaf blocks enclosed by
        !   finer blocks (e.g. droplet interiors: no coarse/fine edges inside a
        !   phase); 0: off
        INTEGER  :: max_island_blocks = 0
    END TYPE LBM_Model

END MODULE Model_Mod
