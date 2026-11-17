!================================================
!   Output writer, selected once at initialization
!
!   The environment variable AMR_OUTPUT chooses one writer:
!     vtk      VTK XML series in vtk/ with .pvd time-series files (default)
!     vtkhdf   two VTKHDF time-series files (needs -DAMR_WITH_FVTKHDF=ON)
!     tecplot  TecPlot ASCII files Rising_Bubble_2D_NN.dat
!     none     no field output
!   Every writer has the interface Write_Step( Time ); writers that keep
!   files open between steps also provide a closing routine.
!================================================
MODULE Output_Mod
    USE Kinds_Mod, ONLY: wp
    USE VTK_Output_Mod, ONLY: Write_VTU_Step
#ifdef USE_FVTKHDF
    USE VTK_Output_Mod, ONLY: Write_VTKHDF_Step, Close_VTKHDF
#endif
    USE TecPlot_Mod, ONLY: TecPlot_CellCenter
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Select_Writer, Write_Output, Close_Output

    ABSTRACT INTERFACE
        SUBROUTINE Write_Step( Time )
            IMPORT :: wp
            REAL(wp), INTENT(IN) :: Time    ! physical time of the output
        END SUBROUTINE Write_Step
        SUBROUTINE Close_Step
        END SUBROUTINE Close_Step
    END INTERFACE

    !== the selected writer; both stay NULL for "none", Closer also for
    !   writers without open files
    PROCEDURE(Write_Step), POINTER :: Writer => NULL()
    PROCEDURE(Close_Step), POINTER :: Closer => NULL()

CONTAINS

!******************************************************************************

    SUBROUTINE Select_Writer
        CHARACTER(len=32) :: val
        INTEGER :: length, stat

        CALL GET_ENVIRONMENT_VARIABLE( 'AMR_OUTPUT', val, length, stat )
        IF( stat /= 0 .OR. length == 0 ) val = 'vtk'
        CALL To_Lower( val )

        SELECT CASE( TRIM(ADJUSTL(val)) )
        CASE( 'vtk' )
            Writer => Write_VTU_Step
        CASE( 'vtkhdf' )
#ifdef USE_FVTKHDF
            Writer => Write_VTKHDF_Step
            Closer => Close_VTKHDF
#else
            ERROR STOP 'AMR_OUTPUT=vtkhdf: built without fVTKHDF (configure with -DAMR_WITH_FVTKHDF=ON)'
#endif
        CASE( 'tecplot' )
            Writer => TecPlot_CellCenter
        CASE( 'none' )
        CASE DEFAULT
            PRINT '(3A)', ' AMR_OUTPUT: unknown writer "', TRIM(val), '" (vtk, vtkhdf, tecplot, none)'
            ERROR STOP
        END SELECT

        PRINT '(/2A)', 'Output writer: ', TRIM(ADJUSTL(val))

    CONTAINS

        SUBROUTINE To_Lower( s )
            CHARACTER(len=*), INTENT(INOUT) :: s
            INTEGER :: k

            DO k = 1, LEN(s)
                IF( s(k:k) >= 'A' .AND. s(k:k) <= 'Z' ) s(k:k) = ACHAR( IACHAR(s(k:k)) + 32 )
            END DO

        END SUBROUTINE To_Lower

    END SUBROUTINE Select_Writer

!******************************************************************************

    SUBROUTINE Write_Output( Time )
        REAL(wp), INTENT(IN) :: Time

        IF( ASSOCIATED(Writer) ) CALL Writer( Time )

    END SUBROUTINE Write_Output

!******************************************************************************

    SUBROUTINE Close_Output

        IF( ASSOCIATED(Closer) ) CALL Closer

    END SUBROUTINE Close_Output

END MODULE Output_Mod
