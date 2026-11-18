!================================================
!   The sides of the domain: periodic, or a boundary with a ghost-cell
!   condition per field
!
!   Side k = 1, 2, 3, 4 is the left, right, bottom and top side, the
!   neighbour directions 1..4 of Blocks_Mod. A periodic side wraps the
!   neighbour table around (Mesh_Mod); the opposite side must be periodic
!   too. Any other side gives the blocks along it the boundary code -k as
!   neighbour, and the ghost fills (Ghost_hg_Mod, Ghost_C_Mod) call the
!   side's procedures for h, g and C on those blocks:
!     h, g( f, i1, i2 )      ghost column (k = 1, 2) or row (k = 3, 4) of a
!                            distribution function at the positions i1:i2
!                            along the side
!     C( ng, C, i1, i2 )     ng ghost layers of the phase field
!   They are called once per block and step, so a procedure pointer costs
!   nothing measurable; the kernels themselves stay plain loops (BC_Mod).
!   The conditions are separate per field so that a variant can replace one
!   of them, e.g. a wetting condition for C on a wall:
!       side = Wall(3);  side%C => My_Wetting_C_3
!
!   The case chooses the sides (AMR_Case%side); the namelist &boundaries of
!   the case file can override each one by name ('periodic' or 'wall').
!   Set_Boundaries checks them and makes them the sides of the run.
!================================================
MODULE Boundaries_Mod
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: nx, ny
    USE BC_Mod
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Side_BC, Periodic, Wall, Set_Boundaries, Side, Is_Periodic, side_name

    ABSTRACT INTERFACE
        SUBROUTINE Ghost_f( f, i1, i2 )
            IMPORT :: wp, nx, ny
            REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)
            INTEGER,  INTENT(IN)    :: i1, i2
        END SUBROUTINE Ghost_f

        SUBROUTINE Ghost_C( ng, C, i1, i2 )
            IMPORT :: wp, nx, ny
            INTEGER,  INTENT(IN)    :: ng, i1, i2
            REAL(wp), INTENT(INOUT) :: C(-1:nx+2,-1:ny+2)
        END SUBROUTINE Ghost_C
    END INTERFACE

    !== one side of the domain (the default: periodic)
    TYPE :: Side_BC
        CHARACTER(len=32) :: name = 'periodic'
        LOGICAL :: periodic = .TRUE.
        PROCEDURE(Ghost_f), POINTER, NOPASS :: h => NULL()
        PROCEDURE(Ghost_f), POINTER, NOPASS :: g => NULL()
        PROCEDURE(Ghost_C), POINTER, NOPASS :: C => NULL()
    END TYPE Side_BC

    CHARACTER(len=6), PARAMETER :: side_name(4) = [ 'left  ', 'right ', 'bottom', 'top   ' ]

    !== the sides of the run (set by Set_Boundaries)
    TYPE(Side_BC), PROTECTED :: Side(4)

CONTAINS

!== a periodic side
FUNCTION Periodic() RESULT( s )
    TYPE(Side_BC) :: s
    s = Side_BC()
END FUNCTION Periodic

!******************************************************************************

!== a no-slip wall on side k: halfway bounce-back of h (no flux of the phase
!   field) and g (no-slip), zero normal gradient of C
FUNCTION Wall( k ) RESULT( s )
    INTEGER, INTENT(IN) :: k
    TYPE(Side_BC) :: s

    s%name     = 'wall'
    s%periodic = .FALSE.
    SELECT CASE( k )
    CASE( 1 )
        s%h => Bounce_Back_1;  s%g => Bounce_Back_1;  s%C => Symmetric_C_1
    CASE( 2 )
        s%h => Bounce_Back_2;  s%g => Bounce_Back_2;  s%C => Symmetric_C_2
    CASE( 3 )
        s%h => Bounce_Back_3;  s%g => Bounce_Back_3;  s%C => Symmetric_C_3
    CASE( 4 )
        s%h => Bounce_Back_4;  s%g => Bounce_Back_4;  s%C => Symmetric_C_4
    CASE DEFAULT
        ERROR STOP 'Wall: the side must be 1, 2, 3 or 4'
    END SELECT

END FUNCTION Wall

!******************************************************************************

!== the side of kind `name` for side k
FUNCTION Side_by_Name( name, k ) RESULT( s )
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER,          INTENT(IN) :: k
    TYPE(Side_BC) :: s

    SELECT CASE( name )
    CASE( 'periodic' )
        s = Periodic()
    CASE( 'wall' )
        s = Wall( k )
    CASE DEFAULT
        PRINT '(5A)', ' &boundaries: unknown kind ''', name, ''' for the ', TRIM(side_name(k)), &
                      ' side (known: periodic, wall)'
        ERROR STOP '&boundaries: unknown kind of side'
    END SELECT

END FUNCTION Side_by_Name

!******************************************************************************

!== the sides of the run: the case's choice, overridden by &boundaries of the
!   case file; periodic sides must come in opposite pairs, the others need
!   all three ghost-cell procedures
SUBROUTINE Set_Boundaries( case_sides )
    USE Case_File_Mod, ONLY: Open_Case_File, Close_Case_File, To_Lower
    IMPLICIT NONE
    TYPE(Side_BC), INTENT(IN) :: case_sides(4)

    CHARACTER(len=32)  :: left, right, bottom, top, names(4)
    CHARACTER(len=256) :: msg
    INTEGER :: k, u, stat
    LOGICAL :: found

    NAMELIST /boundaries/ left, right, bottom, top

    Side = case_sides

    left = '';  right = '';  bottom = '';  top = ''
    CALL Open_Case_File( u, found )
    IF( found ) THEN
        msg = ''
        READ( u, NML=boundaries, IOSTAT=stat, IOMSG=msg )
        CALL Close_Case_File( u, stat, msg, 'boundaries', required=.FALSE. )
    END IF
    names = [ left, right, bottom, top ]
    DO k = 1, 4
        IF( LEN_TRIM(names(k)) > 0 ) Side(k) = Side_by_Name( To_Lower(TRIM(ADJUSTL(names(k)))), k )
    END DO

    IF( (Side(1)%periodic .NEQV. Side(2)%periodic) .OR. (Side(3)%periodic .NEQV. Side(4)%periodic) ) &
        ERROR STOP 'Set_Boundaries: opposite sides must both be periodic or both not'
    DO k = 1, 4
        IF( Side(k)%periodic ) CYCLE
        IF( .NOT. ( ASSOCIATED(Side(k)%h) .AND. ASSOCIATED(Side(k)%g) .AND. ASSOCIATED(Side(k)%C) ) ) &
            ERROR STOP 'Set_Boundaries: a non-periodic side needs ghost-cell procedures for h, g and C'
    END DO

    PRINT '(/A,4(2X,A,": ",A))', 'Boundaries:', ( TRIM(side_name(k)), TRIM(Side(k)%name), k = 1, 4 )


END SUBROUTINE Set_Boundaries

!******************************************************************************

!== side k (1..4) is periodic
PURE LOGICAL FUNCTION Is_Periodic( k )
    INTEGER, INTENT(IN) :: k
    Is_Periodic = Side(k)%periodic
END FUNCTION Is_Periodic

END MODULE Boundaries_Mod
