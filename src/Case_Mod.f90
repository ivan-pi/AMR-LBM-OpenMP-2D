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
!   time; Finish (optional) closes the case's files after the last step.
!   side(1:4) are the left, right, bottom and top sides of the domain
!   (Boundaries_Mod: Periodic() or Wall(k), the default periodic); the
!   namelist &boundaries of the case file can override them.
!================================================
MODULE Case_Mod
    USE Kinds_Mod, ONLY: wp
    USE Model_Mod, ONLY: Block_Kernel
    USE Boundaries_Mod, ONLY: Side_BC
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
    END INTERFACE

    TYPE :: AMR_Case
        CHARACTER(len=80) :: name = ''
        PROCEDURE(Case_Setup),       POINTER, NOPASS :: Setup         => NULL()
        PROCEDURE(Block_Kernel),     POINTER, NOPASS :: Initial_State => NULL()
        PROCEDURE(Case_Diagnostics), POINTER, NOPASS :: Diagnostics   => NULL()
        PROCEDURE(Case_Finish),      POINTER, NOPASS :: Finish        => NULL()
        TYPE(Side_BC) :: side(4)
    END TYPE AMR_Case

END MODULE Case_Mod
