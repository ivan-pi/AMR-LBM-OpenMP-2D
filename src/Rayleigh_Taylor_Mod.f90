!================================================
!   Rayleigh-Taylor instability test case (Case_Mod)
!
!   A heavy fluid (C = 1) above a light fluid (C = 0) in an L0 x nBaseY*L0
!   channel (periodic sides, walls at the bottom and top). The interface
!   lies at mid-height, perturbed by a single cosine mode,
!       y_I(x) = H/2 + A L0 cos(2 pi x / L0),
!   as in Fakhari, Mitchell, Leonardi & Bolster, Phys. Rev. E 96, 053301
!   (2017), Sec. IV B (domain L x 4L: build with -DAMR_DOMAIN_HEIGHT=4).
!   With L = L0, U0 = sqrt(g L) and t0 = sqrt(L / (g At)) the dimensionless
!   numbers are
!       At = (rho_h - rho_l) / (rho_h + rho_l),  Re = rho_h U0 L / mu_h,
!       Ca = mu_h U0 / sigma,                    Pe = U0 L / M,
!   and the time is t* = t / t0. They are read from the namelist
!   &rayleigh_taylor in the case file (Case_File_Mod, see cases/); without
!   one the defaults below (the paper's At = 0.5 case at L0 = 256) are used.
!   The sides of the channel are the case's default boundaries (&boundaries
!   can change them).
!
!   The models apply gravity as the buoyancy (rho_h - rho) g relative to the
!   heavy fluid, so the heavy fluid stays at rest under P = 0 and the light
!   fluid is initialised with the hydrostatic pressure of that force.
!
!   The bubble front (light fluid rising at x = 0) and the spike tip (heavy
!   fluid falling at x = L0/2), measured from the unperturbed interface in
!   units of L0, are written to rt_fronts.dat.
!================================================
MODULE Rayleigh_Taylor_Mod
    USE Kinds_Mod, ONLY: wp
    USE Case_Mod, ONLY: AMR_Case
    USE Boundaries_Mod, ONLY: Periodic, Wall
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Rayleigh_Taylor_Case

    !== case parameters (namelist &rayleigh_taylor), defaults: PRE 2017, Fig. 5-6
    REAL(wp) :: At        = 0.5_wp      ! Atwood number
    REAL(wp) :: Re        = 3000        ! Reynolds number rho_h sqrt(g L0) L0 / mu_h
    REAL(wp) :: mu_ratio  = 1           ! mu_h / mu_l
    REAL(wp) :: Ca        = 0.26_wp     ! capillary number mu_h sqrt(g L0) / sigma
    REAL(wp) :: Gy_lu     = 2.e-6_wp    ! gravity (lattice units); t0 = 16000 steps at L0 = 256
    REAL(wp) :: Amplitude = 0.1_wp      ! perturbation amplitude / L0
    REAL(wp) :: Mobility  = 0.03_wp     ! used if Pe <= 0
    REAL(wp) :: Pe        = 1000        ! Peclet number sqrt(g L0) L0 / M, if > 0
    REAL(wp) :: Width     = 5           ! interface width W, used if Cn <= 0
    REAL(wp) :: Cn        = 0           ! Cahn number W / L0, if > 0
    REAL(wp) :: t_end     = 2.5_wp      ! end time t* = t / t0 (3.2 in the paper, needs height 4)

    NAMELIST /rayleigh_taylor/ At, Re, mu_ratio, Ca, Gy_lu, Amplitude, Mobility, Pe, Width, Cn, t_end

    REAL(wp) :: Y0                          ! unperturbed interface height (H/2)

    !== phase field along a vertical line at finest-level resolution
    REAL(wp), ALLOCATABLE :: column(:)

    !== rt_fronts.dat
    INTEGER :: u
    LOGICAL :: IsOpen = .FALSE.

    REAL(wp), PARAMETER :: pi = 4*ATAN(1._wp)

CONTAINS

FUNCTION Rayleigh_Taylor_Case() RESULT( problem )
    TYPE(AMR_Case) :: problem

    problem%name          =  'Rayleigh-Taylor instability (Fakhari et al., PRE 2017, Sec. IV B)'
    problem%Setup         => Setup
    problem%Initial_State => Initial_State
    problem%Diagnostics   => Front_Diagnostics
    problem%Finish        => Close_Diagnostics
    problem%side          =  [ Periodic(), Periodic(), Wall(3), Wall(4) ]

END FUNCTION Rayleigh_Taylor_Case

!******************************************************************************

!== the namelist &rayleigh_taylor from the case file
SUBROUTINE Read_Parameters
    USE Case_File_Mod, ONLY: Open_Case_File, Close_Case_File
    CHARACTER(len=256) :: msg
    INTEGER :: stat, unit
    LOGICAL :: found

    CALL Open_Case_File( unit, found )
    IF( .NOT. found ) RETURN
    msg = ''
    READ( unit, NML=rayleigh_taylor, IOSTAT=stat, IOMSG=msg )
    CALL Close_Case_File( unit, stat, msg, 'rayleigh_taylor', required=.TRUE. )

END SUBROUTINE Read_Parameters

!******************************************************************************

SUBROUTINE Setup( t_end_out, time_per_step )
    USE Config_Mod, ONLY: L0, nBaseY
    USE Physics_Mod, ONLY: W, Sigma, Beta, kappa, taul, tauh, Rhol, Rhoh, Gy, D, M0, wc
    IMPLICIT NONE
    REAL(wp), INTENT(OUT) :: t_end_out      ! end time t*
    REAL(wp), INTENT(OUT) :: time_per_step  ! t* of one lattice time step: 1/t0 = sqrt(g At / L0)

    REAL(wp) :: U0, mu_h, mu_l

    CALL Read_Parameters
    IF( At <= 0 .OR. At >= 1 ) ERROR STOP 'Rayleigh-Taylor: At must lie in (0,1)'

    D  = L0                         ! reference length
    Y0 = 0.5_wp * nBaseY * L0

    Rhoh = 1
    Rhol = Rhoh * (1 - At) / (1 + At)

    Gy = Gy_lu
    U0 = SQRT(Gy*D)

    IF( Cn > 0 ) THEN
        W = Cn * L0
    ELSE
        W = Width
    END IF
    IF( Pe > 0 ) THEN
        M0 = U0 * D / Pe
    ELSE
        M0 = Mobility
    END IF
    wc = 1._wp/(0.5_wp + M0*3._wp)

    mu_h  = Rhoh * U0 * D / Re
    mu_l  = mu_h / mu_ratio
    Sigma = mu_h * U0 / Ca

    tauh = 3._wp*(mu_h/Rhoh)
    taul = 3._wp*(mu_l/Rhol)

    Beta  = 12.0_wp * Sigma/W
    kappa =  1.5_wp * Sigma*W

    PRINT '(/A,F7.4,F9.1,F7.3,F9.1,A,2F9.5,A,F9.1)', 'At, Re, Ca, Pe =', At, Re, Ca, U0*D/M0, &
        '   tau_h, tau_l =', tauh + 0.5_wp, taul + 0.5_wp, '   t0 =', SQRT(D/(Gy*At))
    IF( nBaseY < 4 .AND. t_end > 2.5_wp ) PRINT '(A)', &
        ' Note: the spike reaches the bottom wall at about t* = 3 unless AMR_DOMAIN_HEIGHT >= 4'

    ALLOCATE( column(0:nBaseY*L0-1) )

    t_end_out     = t_end
    time_per_step = SQRT(Gy*At/D)

END SUBROUTINE Setup

!******************************************************************************

!== fluids at rest in leaf block n: C (with ghost layers), P, Ux, Uy.
!   P is hydrostatic for the buoyancy force (rho_h - rho) g of the models:
!   dP/dy = (rho_h - rho_l) g (1 - C), P = 0 in the heavy fluid. With
!   1 - C = 1/(1 + exp(4 (y - y_I)/W)) this integrates to
!   P = -(rho_h - rho_l) g (W/4) ln(1 + exp(-4 (y - y_I)/W)).
SUBROUTINE Initial_State( n )
    USE Config_Mod, ONLY: nx, ny, L0, dx_level
    USE Blocks_Mod, ONLY: Blk, Blk_C, Blk_P, Blk_Ux, Blk_Uy
    USE Physics_Mod, ONLY: W, Gy, Rhol, Rhoh
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    INTEGER  :: dx, i, j
    REAL(wp) :: Xc, Yc, Yi, s

    dx = dx_level( Blk(n)%Level )

    !== phase field, including the ghost layers
    DO j = -1, ny+2
    DO i = -1, nx+2
        Xc = Blk(n)%I * dx * nx + dx * (i-0.5_wp)
        Yc = Blk(n)%J * dx * ny + dx * (j-0.5_wp)
        Yi = Y0 + Amplitude * L0 * COS( 2*pi*Xc/L0 )
        Blk_C(i,j, n) = 0.5_wp + 0.5_wp*TANH( 2._wp*(Yc-Yi)/W )
    END DO
    END DO

    !== at rest, hydrostatic pressure (ln(1+exp(s)) written without overflow)
    DO j = 1, ny
    DO i = 1, nx
        Xc = Blk(n)%I * dx * nx + dx * (i-0.5_wp)
        Yc = Blk(n)%J * dx * ny + dx * (j-0.5_wp)
        Yi = Y0 + Amplitude * L0 * COS( 2*pi*Xc/L0 )
        s  = -4._wp*(Yc-Yi)/W

        Blk_Ux(i,j, n) = 0
        Blk_Uy(i,j, n) = 0
        Blk_P (i,j, n) = -(Rhoh-Rhol) * Gy * 0.25_wp*W * ( MAX(s,0._wp) + LOG(1 + EXP(-ABS(s))) )
    END DO
    END DO

END SUBROUTINE Initial_State

!******************************************************************************

!== bubble front (x = 0) and spike tip (x = L0/2): the lowest height at which
!   C crosses 1/2 on that vertical line, relative to the unperturbed
!   interface and in units of L0, appended to rt_fronts.dat
SUBROUTINE Front_Diagnostics( Time )
    USE Config_Mod, ONLY: L0
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: Time

    IF( .NOT. IsOpen ) THEN
        OPEN( NEWUNIT=u, FILE='rt_fronts.dat', STATUS='REPLACE', ACTION='WRITE' )
        IsOpen = .True.
        WRITE( u, '(A)' ) '#      t*       bubble/L0     spike/L0'
    END IF

    WRITE( u, '(3ES13.5)' ) Time, ( Lowest_Crossing( 0 ) - Y0 ) / L0, &
                                  ( Lowest_Crossing( L0/2 ) - Y0 ) / L0
    FLUSH( u )

END SUBROUTINE Front_Diagnostics

!******************************************************************************

!== the lowest height at which C crosses 1/2 on the vertical line through the
!   cells that contain the finest-level position x (a cell edge: by symmetry
!   the cells on either side are alike). The leaf cells on that line are
!   spread to finest-level resolution (the interface itself is on the finest
!   level) and the crossing is interpolated linearly between cell centres.
FUNCTION Lowest_Crossing( x ) RESULT( y )
    USE Config_Mod, ONLY: nx, ny, dx_level
    USE Blocks_Mod, ONLY: Blk, Blk_C, nLeaf, LeafByLevel
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: x
    REAL(wp) :: y

    INTEGER :: m, n, dx, i, j, xb, yb, k

    DO m = 1, nLeaf
        n  = LeafByLevel(m)
        dx = dx_level( Blk(n)%Level )
        xb = Blk(n)%I * dx * nx
        IF( x < xb .OR. x >= xb + dx*nx ) CYCLE
        i  = (x - xb)/dx + 1
        yb = Blk(n)%J * dx * ny
        DO j = 1, ny
            column( yb + (j-1)*dx : yb + j*dx - 1 ) = Blk_C(i,j, n)
        END DO
    END DO

    y = 0.5_wp
    DO k = 0, UBOUND(column,1)
        IF( column(k) >= 0.5_wp ) THEN
            IF( k > 0 ) y = k - 0.5_wp + (0.5_wp - column(k-1)) / (column(k) - column(k-1))
            RETURN
        END IF
    END DO
    y = UBOUND(column,1) + 0.5_wp     ! no heavy fluid on the line

END FUNCTION Lowest_Crossing

!******************************************************************************

SUBROUTINE Close_Diagnostics
    IF( IsOpen ) CLOSE( u )
    IsOpen = .False.
END SUBROUTINE Close_Diagnostics

END MODULE Rayleigh_Taylor_Mod
