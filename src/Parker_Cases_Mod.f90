!================================================
!   Test cases of the two-species model (Model_Parker2009_Mod)
!
!   Parameters from the namelist &parker in the file named by AMR_CASE
!   (defaults below). Lattice units of the finest level; one time unit is
!   one time step. The domain is set in Config_Mod (default L0 x 3*L0,
!   periodic in x, walls at the bottom and top); the cases sit in its centre.
!
!   kind = 'droplet'  : a droplet of radius R, rich in species 1, in the
!                       species-2-rich phase, elongated by eps along x
!                       (eps = 0: static droplet, Laplace test; eps > 0:
!                       n = 2 shape oscillation). Written to droplet.dat:
!                       t, deformation (Ixx-Iyy)/(Ixx+Iyy), R_eq, p_in, p_out,
!                       mass of species 1, largest |u|.
!   kind = 'box'      : the whole domain filled with the mixture x1 = x_mix +
!                       noise (noise = 0.5, x_mix = 0.5: uniform in [0, 1]),
!                       e.g. in a periodic box (AMR_DOMAIN_HEIGHT=1, AMR_PERIODIC_Y=ON).
!                       Written to spinodal.dat like 'spinodal'.
!   kind = 'spinodal' : the middle L0 x L0 square filled with the mixture
!                       x1 = x_mix + noise; spinodal decomposition. Written
!                       to spinodal.dat: t, number of leaf cells, mass of
!                       species 1, largest |u|.
!
!   The coexisting mole fractions of the Margules model are computed at
!   start-up (Newton) and used for the initial phases.
!================================================
MODULE Parker_Cases_Mod
    USE Kinds_Mod, ONLY: wp, dp, i64
    USE Case_Mod, ONLY: AMR_Case
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Parker_Case

    !== namelist &parker
    CHARACTER(len=16) :: kind = 'droplet'
    REAL(wp) :: D_lu      = 0.2_wp      ! diffusivity
    REAL(wp) :: kappa     = 1._wp       ! surface-tension scaling
    REAL(wp) :: nu_lu     = 0.02_wp     ! kinematic viscosity
    REAL(wp) :: A_12      = 2.7627_wp   ! Margules parameters
    REAL(wp) :: A_21      = 2.4474_wp
    REAL(wp) :: c_total   = 1._wp       ! total concentration
    REAL(wp) :: Radius    = 20._wp      ! droplet radius
    REAL(wp) :: eps       = 0._wp       ! droplet elongation: semi-axes R(1+eps), R/(1+eps)
    REAL(wp) :: x_mix     = 0.5_wp      ! spinodal: mean mole fraction of species 1
    REAL(wp) :: noise     = 0.05_wp     ! spinodal: amplitude of the perturbation
    INTEGER  :: n_steps   = 12000       ! number of time steps
    !== mesh policy overrides (negative: keep the model's value)
    INTEGER  :: max_island = -1         ! largest coarse island refined (0: off)
    REAL(wp) :: refine_tol = -1         ! refinement / derefinement tolerances of max |grad ln gamma|
    REAL(wp) :: deref_tol  = -1

    REAL(wp) :: xa, xb                  ! coexisting mole fractions of species 1
    REAL(wp) :: X0, Y0                  ! centre of the droplet / of the square

    INTEGER :: u
    LOGICAL :: IsOpen = .FALSE.

CONTAINS

FUNCTION Parker_Case() RESULT( problem )
    TYPE(AMR_Case) :: problem

    problem%name          =  'two-species droplet / spinodal decomposition (namelist &parker)'
    problem%Setup         => Setup
    problem%Initial_State => Initial_State
    problem%Diagnostics   => Diagnostics
    problem%Finish        => Close_Diagnostics
    problem%Adjust_Model  => Adjust_Model

END FUNCTION Parker_Case

!******************************************************************************

SUBROUTINE Setup( t_end_out, time_per_step )
    USE Config_Mod, ONLY: L0, nBaseX, nBaseY
    USE Physics_Mod, ONLY: W, Sigma, taul, tauh, Rhol, Rhoh, M0, wc
    USE Model_Parker2009_Mod, ONLY: Diff, kappa_s, nu, A12, A21
    IMPLICIT NONE
    REAL(wp), INTENT(OUT) :: t_end_out, time_per_step

    NAMELIST /parker/ kind, D_lu, kappa, nu_lu, A_12, A_21, c_total, Radius, eps, x_mix, noise, n_steps, &
                      max_island, refine_tol, deref_tol

    CHARACTER(len=256) :: file
    INTEGER :: length, stat, un

    CALL GET_ENVIRONMENT_VARIABLE( 'AMR_CASE', file, length, stat )
    IF( stat == 0 .AND. length > 0 ) THEN
        OPEN( NEWUNIT=un, FILE=TRIM(file), STATUS='OLD', ACTION='READ', IOSTAT=stat )
        IF( stat /= 0 ) ERROR STOP 'AMR_CASE: cannot open the case file'
        READ( un, NML=parker, IOSTAT=stat )
        IF( stat /= 0 ) ERROR STOP 'AMR_CASE: error in namelist &parker'
        CLOSE( un )
        PRINT '(/2A)', 'Case file: ', TRIM(file)
    END IF
    IF( kind /= 'droplet' .AND. kind /= 'spinodal' .AND. kind /= 'box' ) &
        ERROR STOP 'namelist &parker: kind must be droplet, spinodal or box'

    Diff    = D_lu
    kappa_s = kappa
    nu      = nu_lu
    A12     = A_12
    A21     = A_21

    CALL Coexistence( xa, xb )

    X0 = 0.5_wp*nBaseX*L0
    Y0 = 0.5_wp*nBaseY*L0

    !== the solver's progress table: "Rho" is c_1 (Rhol = 0, Rhoh = 1), "Mass" its integral
    Rhol = 0;  Rhoh = 1
    W = 0;  Sigma = 0;  M0 = Diff;  wc = 1._wp/(3._wp*Diff + 0.5_wp)
    taul = 3._wp*nu;  tauh = taul

    PRINT '(/A,A)',        'kind:                    ', TRIM(kind)
    PRINT '(A,3F10.5)',    'D, kappa, nu:            ', Diff, kappa_s, nu
    PRINT '(A,2F10.5)',    'Margules A12, A21:       ', A12, A21
    PRINT '(A,2F10.5)',    'coexistence x1 (a, b):   ', xa, xb
    PRINT '(A,ES12.4)',    'kappa D^2:               ', kappa_s*Diff**2

    t_end_out     = n_steps
    time_per_step = 1

END SUBROUTINE Setup

!******************************************************************************

!== mesh policy from the namelist (negative values keep the model's defaults)
SUBROUTINE Adjust_Model( model )
    USE Model_Mod, ONLY: LBM_Model
    TYPE(LBM_Model), INTENT(INOUT) :: model

    IF( max_island >= 0 ) model%max_island_blocks = max_island
    IF( refine_tol >  0 ) model%refine_tol = refine_tol
    IF( deref_tol  >  0 ) model%deref_tol  = deref_tol
    PRINT '(A,2ES10.2,I4)', 'mesh: refine_tol, deref_tol, max_island_blocks:', &
        model%refine_tol, model%deref_tol, model%max_island_blocks

END SUBROUTINE Adjust_Model

!******************************************************************************

!== coexisting mole fractions: equal activities x_i gamma_i of both species
SUBROUTINE Coexistence( x_a, x_b )
    USE Model_Parker2009_Mod, ONLY: Ln_Gamma
    REAL(wp), INTENT(OUT) :: x_a, x_b

    REAL(wp) :: v(2), F(2), Fp(2), J(2,2), dv(2), h
    INTEGER  :: it, k

    v = [ 0.05_wp, 0.95_wp ]
    DO it = 1, 100
        F = Residual( v )
        h = 1.e-7_wp
        DO k = 1, 2
            dv = 0;  dv(k) = h
            Fp = Residual( v + dv )
            J(:,k) = ( Fp - F ) / h
        END DO
        dv = Solve2( J, -F )
        v  = v + dv
        IF( MAXVAL(ABS(dv)) < 1.e-13_wp ) EXIT
    END DO
    IF( v(1) >= v(2) .OR. v(1) <= 0 .OR. v(2) >= 1 ) ERROR STOP 'Coexistence: no two-phase solution (Margules A too small?)'
    x_a = v(1)
    x_b = v(2)

CONTAINS

    FUNCTION Residual( v ) RESULT( r )
        REAL(wp), INTENT(IN) :: v(2)
        REAL(wp) :: r(2), g1a, g2a, g1b, g2b
        CALL Ln_Gamma( v(1), g1a, g2a )
        CALL Ln_Gamma( v(2), g1b, g2b )
        r(1) = LOG(v(1)) + g1a - LOG(v(2)) - g1b
        r(2) = LOG(1._wp - v(1)) + g2a - LOG(1._wp - v(2)) - g2b
    END FUNCTION Residual

    FUNCTION Solve2( A, b ) RESULT( x )
        REAL(wp), INTENT(IN) :: A(2,2), b(2)
        REAL(wp) :: x(2), det
        det  = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        x(1) = (  A(2,2)*b(1) - A(1,2)*b(2) ) / det
        x(2) = ( -A(2,1)*b(1) + A(1,1)*b(2) ) / det
    END FUNCTION Solve2

END SUBROUTINE Coexistence

!******************************************************************************

!== a reproducible pseudo-random number in [-1, 1] of the integer position (i, j)
!   The multiplication by 1274126177 is meant modulo 2**64, but signed
!   overflow is not defined in Fortran. Only bits 0..46 of the product
!   reach the result (bits 0..30 after the 16-bit shift), so they are
!   formed exactly from 24- and 23-bit pieces of h, without overflow.
PURE FUNCTION Hash_Noise( i, j ) RESULT( r )
    INTEGER, INTENT(IN) :: i, j
    REAL(wp) :: r
    INTEGER(i64), PARAMETER :: mult = 1274126177_i64
    INTEGER(i64) :: h, lo, hi

    h  = IEOR( INT(i,i64)*73856093_i64, INT(j,i64)*19349663_i64 )
    h  = IEOR( h, ISHFT(h, -13) )
    lo = IBITS( h, 0, 24 )
    hi = IBITS( h, 24, 23 )
    h  = IAND( lo*mult + ISHFT( IBITS(hi*mult, 0, 23), 24 ), 140737488355327_i64 )   ! 2**47 - 1
    h  = IAND( IEOR( h, ISHFT(h, -16) ), 2147483647_i64 )
    r = 2._wp*REAL( h, wp )/2147483647._wp - 1._wp

END FUNCTION Hash_Noise

!******************************************************************************

!== mole fraction of species 1 at (x, y) (finest-level units)
FUNCTION X1_Initial( x, y ) RESULT( x1 )
    USE Config_Mod, ONLY: L0
    REAL(wp), INTENT(IN) :: x, y
    REAL(wp) :: x1

    REAL(wp) :: a, b, r, s

    IF( kind == 'droplet' ) THEN
        a  = Radius*(1._wp + eps)
        b  = Radius/(1._wp + eps)
        r  = SQRT( ((x - X0)/a)**2 + ((y - Y0)/b)**2 )
        s  = (r - 1._wp)*SQRT(a*b)                  ! approximate distance to the interface
        x1 = xa + (xb - xa)*0.5_wp*( 1._wp - TANH(0.5_wp*s) )
    ELSE IF( kind == 'box' ) THEN
        x1 = x_mix + noise*Hash_Noise( FLOOR(x), FLOOR(y) )
    ELSE
        !== the mixture in the square |x-X0|, |y-Y0| < L0/2 - 4 with a smooth edge,
        !   noise on the finest-level cells
        s  = MAX( ABS(x - X0), ABS(y - Y0) ) - (0.5_wp*L0 - 4._wp)
        a  = 0.5_wp*( 1._wp - TANH(0.5_wp*s) )
        x1 = xa + a*( x_mix + noise*Hash_Noise( FLOOR(x), FLOOR(y) ) - xa )
    END IF

END FUNCTION X1_Initial

!******************************************************************************

!== leaf block n at rest: C = c x1, C2 = c (1 - x1) (with ghost layers), P, U
SUBROUTINE Initial_State( n )
    USE Config_Mod, ONLY: nx, ny, dx_level
    USE Blocks_Mod, ONLY: Blk, Blk_C, Blk_C2, Blk_P, Blk_Ux, Blk_Uy
    INTEGER, INTENT(IN) :: n

    INTEGER  :: dx, i, j
    REAL(wp) :: x1

    dx = dx_level( Blk(n)%Level )
    DO j = -1, ny+2
    DO i = -1, nx+2
        x1 = X1_Initial( Blk(n)%I*dx*nx + dx*(i - 0.5_wp), Blk(n)%J*dx*ny + dx*(j - 0.5_wp) )
        Blk_C (i,j,n) = c_total*x1
        Blk_C2(i,j,n) = c_total*(1._wp - x1)
    END DO
    END DO
    Blk_Ux(:,:,n) = 0
    Blk_Uy(:,:,n) = 0
    Blk_P (:,:,n) = c_total/3._wp

END SUBROUTINE Initial_State

!******************************************************************************

SUBROUTINE Diagnostics( Time )
    USE Config_Mod, ONLY: nx, ny, dx_level
    USE Blocks_Mod, ONLY: Blk, nLeaf, LeafByLevel, Blk_C, Blk_C2, Blk_P, Blk_Ux, Blk_Uy
    REAL(wp), INTENT(IN) :: Time

    INTEGER  :: m, n, i, j, dx
    REAL(dp) :: A, Ixx, Iyy, Pin, Nin, Pout, Nout, Mass, w, x, y, r2, x1, area
    REAL(wp) :: Umax

    IF( .NOT. IsOpen ) THEN
        IF( kind == 'droplet' ) THEN
            OPEN( NEWUNIT=u, FILE='droplet.dat', STATUS='REPLACE', ACTION='WRITE' )
            WRITE( u, '(A)' ) '#      t        deformation      R_eq          p_in         p_out        mass_1       max|u|'
        ELSE
            OPEN( NEWUNIT=u, FILE='spinodal.dat', STATUS='REPLACE', ACTION='WRITE' )
            WRITE( u, '(A)' ) '#      t        leaf_cells     mass_1       max|u|'
        END IF
        IsOpen = .True.
    END IF

    A = 0;  Ixx = 0;  Iyy = 0;  Pin = 0;  Nin = 0;  Pout = 0;  Nout = 0;  Mass = 0;  Umax = 0
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED( Blk, Blk_C, Blk_C2, Blk_P, Blk_Ux, Blk_Uy, LeafByLevel, nLeaf, &
!$OMP   xa, xb, X0, Y0, Radius ) REDUCTION(+:A,Ixx,Iyy,Pin,Nin,Pout,Nout,Mass) REDUCTION(max:Umax)
    DO m = 1, nLeaf
        n  = LeafByLevel(m)
        dx = dx_level( Blk(n)%Level )
        area = REAL(dx,dp)**2
        DO j = 1, ny
        DO i = 1, nx
            x  = Blk(n)%I*dx*nx + dx*(i - 0.5_dp) - X0
            y  = Blk(n)%J*dx*ny + dx*(j - 0.5_dp) - Y0
            r2 = x*x + y*y
            x1 = Blk_C(i,j,n) / ( Blk_C(i,j,n) + Blk_C2(i,j,n) )
            !== droplet moments from r < 2.5 R only: the far field relaxes to the discrete
            !   coexistence value, slightly off xa, and would dominate in the tall domain
            w  = MIN( MAX( (x1 - xa)/(xb - xa), 0._dp ), 1._dp ) * area
            IF( r2 > (2.5_dp*Radius)**2 ) w = 0
            A   = A + w
            Ixx = Ixx + w*x*x
            Iyy = Iyy + w*y*y
            Mass = Mass + Blk_C(i,j,n)*area
            IF( r2 < (0.5_dp*Radius)**2 ) THEN
                Pin = Pin + Blk_P(i,j,n)*area;  Nin = Nin + area
            ELSE IF( r2 > (Radius + 12)**2 .AND. r2 < (2.5_dp*Radius)**2 ) THEN
                Pout = Pout + Blk_P(i,j,n)*area;  Nout = Nout + area
            END IF
            Umax = MAX( Umax, Blk_Ux(i,j,n)**2 + Blk_Uy(i,j,n)**2 )
        END DO
        END DO
    END DO
!$OMP END PARALLEL DO

    IF( kind == 'droplet' ) THEN
        WRITE( u, '(F10.1,6ES14.6)' ) Time, (Ixx - Iyy)/(Ixx + Iyy), SQRT(A/(4*ATAN(1._dp))), &
            Pin/MAX(Nin,1._dp), Pout/MAX(Nout,1._dp), Mass, SQRT(Umax)
    ELSE
        WRITE( u, '(F10.1,I12,2ES14.6)' ) Time, nLeaf*nx*ny, Mass, SQRT(Umax)
    END IF
    FLUSH( u )

END SUBROUTINE Diagnostics

!******************************************************************************

SUBROUTINE Close_Diagnostics
    IF( IsOpen ) CLOSE( u )
    IsOpen = .False.
END SUBROUTINE Close_Diagnostics

END MODULE Parker_Cases_Mod
