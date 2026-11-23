!================================================
!   A test case for the AMR-LBM solver
!
!   The drivers select a case and a model (Model_Mod); the solver calls the
!   case only through these procedures. Setup sets the physical parameters
!   (Physics_Mod) and returns the end time t_end in the case's own time unit
!   and that unit's length per lattice time step (all times passed to the
!   output and the diagnostics are in that unit); Initial_State sets C (with
!   its ghost layers), P, Ux and Uy of a leaf block, from which the model
!   computes h and g; Diagnostics records the case's own quantities at a
!   time; Finish (optional) closes the case's files after the last step;
!   Adjust_Model (optional), called after Setup, may change the model's mesh
!   policy (refine_tol, deref_tol, max_island_blocks), not its kernels.
!================================================
MODULE Case_Mod
    USE Kinds_Mod, ONLY: wp
    USE Model_Mod, ONLY: Block_Kernel, LBM_Model
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: AMR_Case

    ABSTRACT INTERFACE
        SUBROUTINE Case_Setup( t_end, time_per_step )
            IMPORT :: wp
            REAL(wp), INTENT(OUT) :: t_end, time_per_step
        END SUBROUTINE Case_Setup

        SUBROUTINE Case_Diagnostics( Time )
            IMPORT :: wp
            REAL(wp), INTENT(IN) :: Time
        END SUBROUTINE Case_Diagnostics

        SUBROUTINE Case_Finish
        END SUBROUTINE Case_Finish

        SUBROUTINE Case_Adjust_Model( model )
            IMPORT :: LBM_Model
            TYPE(LBM_Model), INTENT(INOUT) :: model
        END SUBROUTINE Case_Adjust_Model
    END INTERFACE

    TYPE :: AMR_Case
        CHARACTER(len=80) :: name = ''
        PROCEDURE(Case_Setup),       POINTER, NOPASS :: Setup         => NULL()
        PROCEDURE(Block_Kernel),     POINTER, NOPASS :: Initial_State => NULL()
        PROCEDURE(Case_Diagnostics), POINTER, NOPASS :: Diagnostics   => NULL()
        PROCEDURE(Case_Finish),      POINTER, NOPASS :: Finish        => NULL()
        PROCEDURE(Case_Adjust_Model), POINTER, NOPASS :: Adjust_Model => NULL()
    END TYPE AMR_Case

END MODULE Case_Mod
