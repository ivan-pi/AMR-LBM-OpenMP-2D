!================================================
!   A test case for the AMR-LBM solver
!
!   The drivers select a case and a model (Model_Mod); the solver calls the
!   case only through these procedures. Setup sets the physical parameters
!   (Physics_Mod) and returns the end time t* (in units of sqrt(D/Gy));
!   Initial_State sets C (with its ghost layers), P, Ux and Uy of a leaf
!   block, from which the model computes h and g; Diagnostics records the
!   case's own quantities at time t* (a negative time closes its files).
!================================================
MODULE Case_Mod
    USE Kinds_Mod, ONLY: wp
    USE Model_Mod, ONLY: Block_Kernel
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: AMR_Case

    ABSTRACT INTERFACE
        SUBROUTINE Case_Setup( t_end )
            IMPORT :: wp
            REAL(wp), INTENT(OUT) :: t_end
        END SUBROUTINE Case_Setup

        SUBROUTINE Case_Diagnostics( Time )
            IMPORT :: wp
            REAL(wp), INTENT(IN) :: Time
        END SUBROUTINE Case_Diagnostics
    END INTERFACE

    TYPE :: AMR_Case
        CHARACTER(len=80) :: name = ''
        PROCEDURE(Case_Setup),       POINTER, NOPASS :: Setup         => NULL()
        PROCEDURE(Block_Kernel),     POINTER, NOPASS :: Initial_State => NULL()
        PROCEDURE(Case_Diagnostics), POINTER, NOPASS :: Diagnostics   => NULL()
    END TYPE AMR_Case

END MODULE Case_Mod
