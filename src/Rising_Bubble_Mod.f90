!================================================
!
!       2D Two-Phase AMR-LBM Solver
!   Conservative Phase-Field + Pressure Evolution
!               Rising Bubble
!
!       based on the following paper:
!   A. Fakhari et al., "A mass-conserving lattice Boltzmann method
!   with dynamic grid refinement for immiscible two-phase flows"
!   Journal of Computational Physics 315 (2016) 434-457
!
!   Written by Abbas Fakhari    06/05/2014
!   updated                     08/30/2016
!   updated by Ivan Pribec      11/15/2026
!     (modules, OpenMP, VTK/VTKHDF output, CMake build)
!================================================

!================================================
!   Rising-bubble test case (Case_Mod): parameters, initial state and the
!   rise of the bubble; the solver (Solver_Mod) runs it with an LBM model.
!
!   Bubble of diameter D = L0/5 at (L0/2, L0/2) in an L0 x 3*L0 channel,
!   periodic sides, walls at the bottom and top (Fakhari, Geier & Lee,
!   J. Comput. Phys. 315 (2016), Sec. 5.1). The dimensionless numbers are
!   read from the namelist &rising_bubble in the case file (Case_File_Mod,
!   see cases/); without one the defaults below are used. The sides of the
!   channel are the case's default boundaries (&boundaries can change
!   them). The rise of the bubble is written to bubble_rise.dat.
!================================================
MODULE Rising_Bubble_Mod
    USE Kinds_Mod, ONLY: wp, dp
    USE Case_Mod, ONLY: AMR_Case
    USE Boundaries_Mod, ONLY: Periodic, Wall
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Rising_Bubble_Case

    !== case parameters (namelist &rising_bubble), defaults: the original case
    REAL(wp) :: Eo        = 40          ! Eotvos number Gy (rho_h-rho_l) D^2 / sigma
    REAL(wp) :: Re        = 40          ! gravity Reynolds number sqrt(Gy rho_h (rho_h-rho_l) D^3) / mu_h
    REAL(wp) :: rho_ratio = 1000        ! rho_h / rho_l
    REAL(wp) :: mu_ratio  = 100         ! mu_h / mu_l
    REAL(wp) :: Gy_lu     = 5.e-6_wp    ! gravity (lattice units)
    REAL(wp) :: Mobility  = 0.03_wp     ! used if Pe <= 0
    REAL(wp) :: Pe        = 0           ! Peclet number sqrt(Gy D) W / M, if > 0
    REAL(wp) :: Width     = 5           ! interface width W, used if Cn <= 0
    REAL(wp) :: Cn        = 0           ! Cahn number W / L0, if > 0
    REAL(wp) :: t_end     = 8           ! end time t* = t sqrt(Gy/D)

    REAL(wp) :: mu_h                        ! dynamic viscosity of the heavy fluid
    REAL(wp) :: R                           ! bubble radius
    INTEGER  :: X0, Y0                      ! bubble centre

    !== bubble_rise.dat
    INTEGER :: u
    LOGICAL :: IsOpen = .FALSE.

CONTAINS

FUNCTION Rising_Bubble_Case() RESULT( problem )
    TYPE(AMR_Case) :: problem

    problem%name          =  'rising bubble (Fakhari, Geier & Lee, JCP 2016, Sec. 5.1)'
    problem%Setup         => Setup
    problem%Initial_State => Initial_State
    problem%Diagnostics   => Bubble_Diagnostics
    problem%Finish        => Close_Diagnostics
    problem%side          =  [ Periodic(), Periodic(), Wall(3), Wall(4) ]

END FUNCTION Rising_Bubble_Case

!******************************************************************************

SUBROUTINE Setup( t_end_out, time_per_step )
    USE Config_Mod, ONLY: L0
    USE Physics_Mod, ONLY: W, Sigma, Beta, kappa, taul, tauh, Rhol, Rhoh, Gy, D, M0, wc
    USE Case_File_Mod, ONLY: Open_Case_File, Close_Case_File
    IMPLICIT NONE
    REAL(wp), INTENT(OUT) :: t_end_out      ! end time t*
    REAL(wp), INTENT(OUT) :: time_per_step  ! t* of one lattice time step: sqrt(Gy/D)

    NAMELIST /rising_bubble/ Eo, Re, rho_ratio, mu_ratio, Gy_lu, Mobility, Pe, Width, Cn, t_end

    CHARACTER(len=256) :: msg
    INTEGER  :: stat, u
    LOGICAL  :: found
    REAL(wp) :: mul

    !== bubble of diameter D = L0/5 at (L0/2, L0/2)
    R  = L0/10._wp
    D  = 2._wp*R
    X0 = L0/2
    Y0 = L0/2

    CALL Open_Case_File( u, found )
    IF( found ) THEN
        msg = ''
        READ( u, NML=rising_bubble, IOSTAT=stat, IOMSG=msg )
        CALL Close_Case_File( u, stat, msg, 'rising_bubble', required=.TRUE. )
    END IF

    Rhoh = 1
    Rhol = Rhoh / rho_ratio

    Gy = Gy_lu

    IF( Cn > 0 ) THEN
        W = Cn * L0
    ELSE
        W = Width
    END IF
    IF( Pe > 0 ) THEN
        M0 = SQRT(Gy*D) * W / Pe
    ELSE
        M0 = Mobility
    END IF
    wc = 1._wp/(0.5_wp + M0*3._wp)

    Sigma = Gy * (Rhoh-Rhol) * D**2 / Eo
    mu_h  = SQRT(Gy*Rhoh*(Rhoh-Rhol)*D**3) / Re

    mul   = mu_h/mu_ratio

    tauh = 3._wp*(mu_h/Rhoh)
    taul = 3._wp*(mul/Rhol)

    Beta  = 12.0_wp * Sigma/W
    kappa =  1.5_wp * Sigma*W

    PRINT '(/A,5F12.3)', 'Eo, Re_h, muh, mul, Pe =', Eo, Re, mu_h, mul, SQRT(Gy*D)*W/M0

    t_end_out     = t_end
    time_per_step = SQRT(Gy/D)

END SUBROUTINE Setup

!******************************************************************************

!== bubble at rest in leaf block n: C (with ghost layers), P, Ux, Uy
SUBROUTINE Initial_State( n )
    USE Config_Mod, ONLY: nx, ny, dx_level
    USE Blocks_Mod, ONLY: Blk, Blk_C, Blk_P, Blk_Ux, Blk_Uy
    USE Physics_Mod, ONLY: W, Sigma
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    INTEGER  :: dx, i, j
    REAL(wp) :: Ri, Xc, Yc

    dx = dx_level( Blk(n)%Level )

    !== phase field, including the ghost layers
    DO j = -1, ny+2
    DO i = -1, nx+2

        Xc = Blk(n)%I * dx * nx + dx * (i-0.5_wp)
        Yc = Blk(n)%J * dx * ny + dx * (j-0.5_wp)

        Ri = SQRT((Xc-X0)**2 + (Yc-Y0)**2)

        Blk_C(i,j, n) = 0.5_wp + 0.5_wp*TANH(2._wp*(Ri-R)/W)

    END DO
    END DO
    !== at rest, with the Laplace pressure jump
    DO j = 1, ny
    DO i = 1, nx

        Blk_Ux(i,j, n) = 0
        Blk_Uy(i,j, n) = 0

        Blk_P(i,j, n) = -Sigma/R * Blk_C(i,j, n)

    END DO
    END DO

END SUBROUTINE Initial_State

!******************************************************************************

!== rise of the bubble: centroid y_c and rise velocity u_r of the light phase
!   (weight 1-C), and the rise Reynolds number Re_r = rho_h u_r D / mu_h,
!   appended to bubble_rise.dat.
!   Cells with C >= 0.9 are left out: faint tails of C slightly below 1 in
!   coarse cells far from the bubble (wake, coarse-fine edges) carry a large
!   area and change whenever blocks are refined or derefined, which makes
!   y_c jump.
SUBROUTINE Bubble_Diagnostics( Time )
    USE Config_Mod, ONLY: nx, ny, dx_level
    USE Blocks_Mod, ONLY: Blk, nLeaf, LeafByLevel, Blk_C, Blk_Uy
    USE Physics_Mod, ONLY: Rhoh, D
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: Time

    INTEGER  :: m, n, i, j, dx
    REAL(dp) :: A, Ay, Av, wgt

    IF( .NOT. IsOpen ) THEN
        OPEN( NEWUNIT=u, FILE='bubble_rise.dat', STATUS='REPLACE', ACTION='WRITE' )
        IsOpen = .True.
        WRITE( u, '(A)' ) '#      t*           y_c          u_r          Re_r'
    END IF

    A  = 0
    Ay = 0
    Av = 0
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED( Blk, Blk_C, Blk_Uy, LeafByLevel, nLeaf ) REDUCTION(+:A,Ay,Av)
    DO m = 1, nLeaf
        n  = LeafByLevel(m)
        dx = dx_level( Blk(n)%Level )
        DO j = 1, ny
        DO i = 1, nx
            IF( Blk_C(i,j, n) >= 0.9_wp ) CYCLE
            wgt = REAL( 1._wp - MAX( Blk_C(i,j, n), 0._wp ), dp ) * dx**2
            A   = A  + wgt
            Ay  = Ay + wgt * ( Blk(n)%J*dx*ny + dx*(j-0.5_dp) )
            Av  = Av + wgt * Blk_Uy(i,j, n)
        END DO
        END DO
    END DO
!$OMP END PARALLEL DO

    WRITE( u, '(4ES13.5)' ) Time, Ay/A, Av/A, Rhoh*(Av/A)*D/mu_h
    FLUSH( u )

END SUBROUTINE Bubble_Diagnostics

!******************************************************************************

SUBROUTINE Close_Diagnostics
    IF( IsOpen ) CLOSE( u )
    IsOpen = .False.
END SUBROUTINE Close_Diagnostics

END MODULE Rising_Bubble_Mod
