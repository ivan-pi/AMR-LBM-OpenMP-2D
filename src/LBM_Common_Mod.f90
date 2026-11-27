!================================================
!   Lattice Boltzmann kernels shared by all models: streaming of the
!   distribution functions and the phase-field moment C = sum(h)
!================================================
MODULE LBM_Common_Mod
    USE Kinds_Mod, ONLY: wp
    USE D2Q9_Mod, ONLY: ex, ey
    USE Config_Mod, ONLY: nx, ny     ! the fixed block size
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Perfect_Shift_AMR, Streaming_LaxWendroff_AMR, Streaming_LaxWendroff9Pt_AMR, &
              Macroscopic_Values_h_AMR, Streaming_Kernel, Streaming_Coarse_AMR, coarseStencil, &
              Select_Coarse_Streaming

    ABSTRACT INTERFACE
        SUBROUTINE Streaming_Kernel( cfl, f )
            IMPORT :: wp, nx, ny
            REAL(wp), INTENT(IN)    :: cfl
            REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)
        END SUBROUTINE Streaming_Kernel
    END INTERFACE

    !== the streaming of the levels below the finest and the number of points
    !   of its stencil (3 or 9), set by Select_Coarse_Streaming; the ghost cells
    !   of h and g follow the stencil (Ghost_hg_Mod)
    PROCEDURE(Streaming_Kernel), POINTER, PROTECTED :: Streaming_Coarse_AMR => Streaming_LaxWendroff_AMR
    INTEGER, PROTECTED :: coarseStencil = 3

CONTAINS

!== Both streaming kernels compute a direction into a temporary of one
!   nx x ny plane and copy it back before the next direction. (With one
!   temporary for all directions and a single array assignment at the end,
!   gfortran 15 vectorizes the copy across the directions with shuffles,
!   2.4 times slower; the plane also stays in L1.
!   20-180 % faster with gfortran, ifx and flang, same results.)
SUBROUTINE Perfect_Shift_AMR( f )
    IMPLICIT NONE
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)

    INTEGER :: i, X, Y
    REAL(wp) :: fnew(nx,ny)

    DO i = 2, 9
        DO Y = 1, ny
        DO X = 1, nx
            fnew(X,Y) = f(X-ex(i),Y-ey(i),i)
        END DO
        END DO
        DO Y = 1, ny
        DO X = 1, nx
            f(X,Y,i) = fnew(X,Y)
        END DO
        END DO
    END DO

END SUBROUTINE Perfect_Shift_AMR

!**********************************************************************

!== Lax-Wendroff streaming with CFL number cfl < 1 (the levels below the
!   finest), three points along the direction e_i (upwind, central,
!   downwind): the 1-D scheme on the characteristic
SUBROUTINE Streaming_LaxWendroff_AMR( cfl, f )
    IMPLICIT NONE
    REAL(wp), INTENT(IN)    :: cfl
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)

    INTEGER :: i, X, Y
    REAL(wp) :: fnew(nx,ny), cU, c0, cD

    !== Lax-Wendroff weights of the upwind, central and downwind values
    cU = 0.5_wp*cfl*(1._wp+cfl)
    c0 = 1._wp - cfl*cfl
    cD = 0.5_wp*cfl*(1._wp-cfl)

    DO i = 2, 9
        DO Y = 1, ny
        DO X = 1, nx
            fnew(X,Y) = cU * f(X-ex(i),Y-ey(i),i) &
                      + c0 * f(X,Y,i) &
                      - cD * f(X+ex(i),Y+ey(i),i)
        END DO
        END DO
        DO Y = 1, ny
        DO X = 1, nx
            f(X,Y,i) = fnew(X,Y)
        END DO
        END DO
    END DO

END SUBROUTINE Streaming_LaxWendroff_AMR

!**********************************************************************

!== Lax-Wendroff streaming with nine points: for the diagonal directions the
!   tensor product of the 1-D schemes in x and in y, done as a pass in x
!   (into a plane with the two ghost rows) and a pass in y. It is the 2-D
!   Lax-Wendroff scheme (4-point cross derivative) plus terms of third and
!   fourth order. The axis directions (2:5) are those of the 3-point scheme.
!   It reads all four diagonal directions in the ghost corners (Ghost_hg_Mod)
!   and costs about the same as the 3-point scheme.
!   tools/check_lax_wendroff.py checks the order and the stability of both
!   schemes (Taylor and von Neumann analysis).
SUBROUTINE Streaming_LaxWendroff9Pt_AMR( cfl, f )
    IMPLICIT NONE
    REAL(wp), INTENT(IN)    :: cfl
    REAL(wp), INTENT(INOUT) :: f(0:nx+1,0:ny+1,9)

    INTEGER :: i, X, Y
    REAL(wp) :: fnew(nx,0:ny+1), cU, c0, cD

    cU = 0.5_wp*cfl*(1._wp+cfl)
    c0 = 1._wp - cfl*cfl
    cD = 0.5_wp*cfl*(1._wp-cfl)

    DO i = 2, 5
        DO Y = 1, ny
        DO X = 1, nx
            fnew(X,Y) = cU * f(X-ex(i),Y-ey(i),i) &
                      + c0 * f(X,Y,i) &
                      - cD * f(X+ex(i),Y+ey(i),i)
        END DO
        END DO
        DO Y = 1, ny
        DO X = 1, nx
            f(X,Y,i) = fnew(X,Y)
        END DO
        END DO
    END DO

    DO i = 6, 9
        DO Y = 0, ny+1
        DO X = 1, nx
            fnew(X,Y) = cU * f(X-ex(i),Y,i) &
                      + c0 * f(X,Y,i) &
                      - cD * f(X+ex(i),Y,i)
        END DO
        END DO
        DO Y = 1, ny
        DO X = 1, nx
            f(X,Y,i) = cU * fnew(X,Y-ey(i)) &
                     + c0 * fnew(X,Y) &
                     - cD * fnew(X,Y+ey(i))
        END DO
        END DO
    END DO

END SUBROUTINE Streaming_LaxWendroff9Pt_AMR

!**********************************************************************

!== the streaming of the levels below the finest, from the environment
!   variable AMR_COARSE_STENCIL: 3 (default) or 9 points
SUBROUTINE Select_Coarse_Streaming
    IMPLICIT NONE
    CHARACTER(len=32) :: val
    INTEGER :: length, stat

    CALL GET_ENVIRONMENT_VARIABLE( 'AMR_COARSE_STENCIL', val, length, stat )
    IF( stat /= 0 .OR. length == 0 ) val = '3'

    SELECT CASE( TRIM(ADJUSTL(val)) )
    CASE( '3' )
        coarseStencil = 3
        Streaming_Coarse_AMR => Streaming_LaxWendroff_AMR
    CASE( '9' )
        coarseStencil = 9
        Streaming_Coarse_AMR => Streaming_LaxWendroff9Pt_AMR
    CASE DEFAULT
        PRINT '(3A)', ' AMR_COARSE_STENCIL: unknown stencil "', TRIM(val), '" (3, 9)'
        ERROR STOP
    END SELECT

    PRINT '(/A,I0,A)', 'Streaming on the coarse levels: Lax-Wendroff, ', coarseStencil, ' points'

END SUBROUTINE Select_Coarse_Streaming

!**********************************************************************

SUBROUTINE Macroscopic_Values_h_AMR( h, C )
    IMPLICIT NONE
    REAL(wp), INTENT(IN)    :: h(0:nx+1,0:ny+1,9)
    REAL(wp), INTENT(INOUT) :: C(  -1:nx+2,-1:ny+2)

    INTEGER :: k

    C(1:nx,1:ny) = h(1:nx,1:ny,1)
    DO k = 2, 9
        C(1:nx,1:ny) = C(1:nx,1:ny) + h(1:nx,1:ny,k)
    END DO

END SUBROUTINE Macroscopic_Values_h_AMR

END MODULE LBM_Common_Mod
