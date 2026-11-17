!================================================
!   An LBM model as a set of block kernels
!
!   The drivers select a model; the time loop, the mesh adaptation and
!   the initialization call the model only through these kernels, which
!   act on the leaf block with the given index.
!================================================
MODULE Model_Mod
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Block_Kernel, LBM_Model

    ABSTRACT INTERFACE
        SUBROUTINE Block_Kernel( n )
            INTEGER, INTENT(IN) :: n        ! block index
        END SUBROUTINE Block_Kernel
    END INTERFACE

    TYPE :: LBM_Model
        CHARACTER(len=80) :: name = ''
        !== ghost layers of C read by Collide (filled before the collision)
        INTEGER :: ghost_C = 1
        !== h and g from C (with ghost layers), P, Ux, Uy (initial state)
        PROCEDURE(Block_Kernel), POINTER, NOPASS :: Initialize => NULL()
        !== collision of h and g
        PROCEDURE(Block_Kernel), POINTER, NOPASS :: Collide    => NULL()
        !== pressure and velocity from g (needs C with one ghost layer)
        PROCEDURE(Block_Kernel), POINTER, NOPASS :: Moments_g  => NULL()
        !== Moments_g followed by Collide, sharing the derivatives of C
        PROCEDURE(Block_Kernel), POINTER, NOPASS :: Moments_and_Collide => NULL()
    END TYPE LBM_Model

END MODULE Model_Mod
