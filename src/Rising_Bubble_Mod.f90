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
!   Rising-bubble case: setup and time loop, for a given LBM model
!   (Model_Mod); the drivers main_*.f90 select the model.
!
!   Bubble of diameter D = L0/5 at (L0/2, L0/2) in an L0 x 3*L0 channel,
!   periodic sides, walls at the bottom and top (Fakhari, Geier & Lee,
!   J. Comput. Phys. 315 (2016), Sec. 5.1). The dimensionless numbers are
!   read from the namelist &rising_bubble in the file named by the
!   environment variable AMR_CASE (see cases/); without it the defaults
!   below are used. The rise of the bubble is written to bubble_rise.dat.
!================================================
MODULE Rising_Bubble_Mod
    USE Kinds_Mod, ONLY: wp, dp
    USE Mesh_Mod, ONLY: Initialize_Base_Block, Fix_Initial_Neighbors_3Blocks, Initial_Mesh, &
                        Determine_Min_Max_Levels, Print_Diagnostics
    USE Refine_Mod, ONLY: Refine_Derefine_Mesh, Block_Error
    USE LBM_Common_Mod, ONLY: Perfect_Shift_AMR, Streaming_LaxWendroff_AMR, Macroscopic_Values_h_AMR
    USE Model_Mod, ONLY: LBM_Model
    USE Ghost_C_Mod, ONLY: Fill_Ghost_Cells_C
    USE Ghost_hg_Mod, ONLY: Fill_Ghost_Cells_hg
    USE Output_Mod, ONLY: Select_Writer, Write_Output, Close_Output
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: AMR_LBM_Rising_Bubble

    !== refinement / derefinement tolerances of the interface criterion
    REAL(wp), PARAMETER :: refinTol = 0.002_wp
    REAL(wp), PARAMETER :: derefTol = 0.001_wp

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

CONTAINS

!== run the rising-bubble case with an LBM model and a refinement criterion
SUBROUTINE AMR_LBM_Rising_Bubble( model, criterion )
    USE Parameters_Mod, ONLY: minLevels, maxLevels, LevelMin, LevelMax
    IMPLICIT NONE
    TYPE(LBM_Model), INTENT(IN) :: model
    PROCEDURE(Block_Error)      :: criterion

    INTEGER    :: L
    INTEGER(8) :: t1, t2, Rate

    CALL SYSTEM_CLOCK( t1, Rate )

    LevelMin = minLevels
    LevelMax = minLevels

    CALL Initialize_Base_Block( 1 ) ! lower  part of the domain
    CALL Initialize_Base_Block( 2 ) ! middle part of the domain
    CALL Initialize_Base_Block( 3 ) ! upper  part of the domain
    CALL Fix_Initial_Neighbors_3Blocks

    CALL Initial_Mesh( LevelMin )

    CALL Initialize_AMR_Rising_Bubble( model )

    !== refine level by level; LevelMax tells the refinement routines the current depth
    DO L = minLevels, maxLevels-1
        LevelMax = L
        CALL Refine_Derefine_Mesh( criterion, model%Moments_g, refinTol, derefTol )
        CALL Initialize_AMR_Rising_Bubble( model )
    END DO

    CALL SYSTEM_CLOCK( t2 )

    CALL Main_Loop( model, criterion )

    PRINT '(A,F12.3)', ' t_initialization:', DBLE(t2-t1)/Rate

END

!******************************************************************************

!== bubble of diameter D = L0/5 at rest, Laplace pressure jump; the model
!   turns C, P, U into its distribution functions h, g
SUBROUTINE Initialize_AMR_Rising_Bubble( model )
    USE Parameters_Mod, ONLY: Blk, L0, nLeaf, LeafBlks, nx, ny
    USE D2Q9_Mod,   ONLY: W, D, Sigma
    IMPLICIT NONE
    TYPE(LBM_Model), INTENT(IN) :: model

    LOGICAL, SAVE :: IsFirstTime = .True.

    INTEGER :: dx, i, j, X0, Y0, n, m
    REAL(wp) :: R, Ri, Xc, Yc

    R  = L0/10._wp
    D  = 2._wp*R

    X0 = L0/2
    Y0 = L0/2

    IF( IsFirstTime ) THEN
        CALL Set_Rising_Bubble_Parameters
        IsFirstTime = .False.
    END IF

!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP   SHARED( Blk, LeafBlks, nLeaf, R, X0, Y0, W, Sigma, model )
    DO m = 1, nLeaf
        n = LeafBlks(m)

        dx = Blk(n)%dx

        !=============  first initialize the phase-field
        DO j = -1, ny+2
        DO i = -1, nx+2

            Xc = Blk(n)%I * dx * nx + dx * (i-0.5_wp)
            Yc = Blk(n)%J * dx * ny + dx * (j-0.5_wp)

            Ri = SQRT((Xc-X0)**2 + (Yc-Y0)**2)

            Blk(n)%C(i,j) = 0.5_wp + 0.5_wp*TANH(2._wp*(Ri-R)/W)

        END DO
        END DO
        !===========    now initialize the hydrodynamic properties
        DO j = 1, ny
        DO i = 1, nx

            Blk(n)%Ux(i,j) = 0
            Blk(n)%Uy(i,j) = 0

            Blk(n)%P (i,j) = -Sigma/R * Blk(n)%C(i,j)

        END DO
        END DO
        !=====

        CALL model%Initialize( n )

    END DO
!$OMP END PARALLEL DO

END

!******************************************************************************

SUBROUTINE Set_Rising_Bubble_Parameters
    USE Parameters_Mod, ONLY: L0
    USE D2Q9_Mod, ONLY: W, Sigma, Beta, kappa, taul, tauh, Rhol, Rhoh, Gy, Bo, Ar, D, M0, wc
    IMPLICIT NONE

    NAMELIST /rising_bubble/ Eo, Re, rho_ratio, mu_ratio, Gy_lu, Mobility, Pe, Width, Cn, t_end

    CHARACTER(len=256) :: file
    INTEGER  :: length, stat, u
    REAL(wp) :: mul

    CALL GET_ENVIRONMENT_VARIABLE( 'AMR_CASE', file, length, stat )
    IF( stat == 0 .AND. length > 0 ) THEN
        OPEN( NEWUNIT=u, FILE=TRIM(file), STATUS='OLD', ACTION='READ', IOSTAT=stat )
        IF( stat /= 0 ) ERROR STOP 'AMR_CASE: cannot open the case file'
        READ( u, NML=rising_bubble, IOSTAT=stat )
        IF( stat /= 0 ) ERROR STOP 'AMR_CASE: error in namelist &rising_bubble'
        CLOSE( u )
        PRINT '(/2A)', 'Case file: ', TRIM(file)
    END IF

    Rhoh = 1
    Rhol = Rhoh / rho_ratio

    Bo = Eo
    Ar = Re

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

    Sigma = Gy * (Rhoh-Rhol) * D**2 / Bo
    mu_h  = SQRT(Gy*Rhoh*(Rhoh-Rhol)*D**3) / Ar

    mul   = mu_h/mu_ratio

    tauh = 3._wp*(mu_h/Rhoh)
    taul = 3._wp*(mul/Rhol)

    Beta  = 12.0_wp * Sigma/W
    kappa =  1.5_wp * Sigma*W

    PRINT '(/A,4F12.3)', 'muh, mul, Re_h, Pe =', mu_h, mul, Ar, SQRT(Gy*D)*W/M0

END

!******************************************************************************

SUBROUTINE Main_Loop( model, criterion )
    USE Parameters_Mod, ONLY: Blk, T, nLeaf, LeafBlks, nx, ny
    USE D2Q9_Mod,   ONLY: W, M0, Bo, Ar, taul, tauh, wc, Sigma, Gy, D
    USE OMP_LIB, ONLY: OMP_GET_MAX_THREADS
    USE ISO_FORTRAN_ENV, ONLY: OUTPUT_UNIT
    IMPLICIT NONE
    TYPE(LBM_Model), INTENT(IN) :: model
    PROCEDURE(Block_Error)      :: criterion

    INTEGER :: n, m, Tf, nOutStep, nargs, ios
    CHARACTER(len=32) :: arg
    REAL(wp)    :: cfl
    LOGICAL    :: changed, moments_done, is_output, is_diagnostic
    INTEGER    :: nDiagStep
    INTEGER(8) :: t1, Rate
    INTEGER(8) :: t_collision, t_advection, t_Output, t_Ghost_C, t_Ghost_hg, t_refine
    INTEGER    :: nProgStep
    INTEGER(8) :: tp0, tpLast, cellUpdates

    !==================================
    Tf   = NINT(t_end * SQRT(D/Gy) /20)*20
    !== optional command-line override of the number of time steps
    nargs = COMMAND_ARGUMENT_COUNT()
    IF( nargs >= 1 ) THEN
        CALL GET_COMMAND_ARGUMENT( 1, arg )
        READ( arg, *, IOSTAT=ios ) Tf
        IF( ios /= 0 .OR. Tf < 1 ) STOP 'usage: amr_lbm [number_of_time_steps [output_interval_in_steps]]'
    END IF
    !== optional output interval in time steps (default: 20 outputs over the run);
    !   the initial and the final state are always written
    nOutStep = MAX( Tf/20, 1 )
    IF( nargs >= 2 ) THEN
        CALL GET_COMMAND_ARGUMENT( 2, arg )
        READ( arg, *, IOSTAT=ios ) nOutStep
        IF( ios /= 0 .OR. nOutStep < 1 ) STOP 'usage: amr_lbm [number_of_time_steps [output_interval_in_steps]]'
    END IF
    CALL Select_Writer
    !== the rise of the bubble is recorded 120 times per run and at every output
    nDiagStep = MAX( Tf/120, 1 )
    !== a progress line (speed, estimated time left) 20 times per run
    nProgStep = MAX( Tf/20, 1 )
    !==================================

    t_Output    = 0
    t_collision = 0
    t_advection = 0
    t_Ghost_C   = 0
    t_Ghost_hg  = 0
    t_refine    = 0

    !==========================================================
    PRINT '(/2A)', 'Model: ', TRIM(model%name)
    PRINT '(/A/)', '    Tf     Sigma    W      Eo    Re_h    M0    tauc   taul   tauh     Ref_Tols'
    PRINT '(I8,F8.4,F6.2,F7.1,F7.1,4F7.3,2F8.4)', Tf, Sigma, W, Bo, Ar, M0, 1/wc, taul, tauh, refinTol, derefTol
    PRINT '(/A/)','     T   Nblks  nLeaf  Unused  Rho_min    Rho_max   U_max      Mass     Levels'
    !==========================================================

    CALL SYSTEM_CLOCK( t1, Rate )
    tp0 = t1
    tpLast = t1
    cellUpdates = 0

    !== One time step, starting from h, g and C (interior) after streaming:
    !   1) ghost layers of C, as many as the model's collision reads; they also
    !      serve the refinement criterion and the moments of g
    !   2) mesh adaptation (from the second step on: the initialization adapted
    !      the mesh); if it changed, the ghost layers are filled again
    !   3) output of the current state and the rise of the bubble (pressure and
    !      velocity computed first)
    !   4) one parallel region: moments of g and collision in one block loop,
    !      ghost cells of h and g, streaming together with C = sum(h).
    !      The work of each stage is shared with orphaned !$OMP DO constructs
    !      (also inside the Fill_Ghost_Cells_* routines); their implicit
    !      barriers separate the stages. Timers are updated by the primary
    !      thread right after each barrier.
    !   The moments of g are not recomputed before the first collision (the
    !   initial P, U are used) nor twice on output steps.
    DO T = 0, Tf

        CALL Determine_Min_Max_Levels
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED( model )
        CALL Fill_Ghost_Cells_C( model%ghost_C )
!$OMP END PARALLEL
        CALL Lap( t_Ghost_C )

        IF( T > 0 ) THEN
            CALL Refine_Derefine_Mesh( criterion, model%Moments_g, refinTol, derefTol, changed )
            CALL Determine_Min_Max_Levels
            IF( changed ) THEN
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED( model )
                CALL Fill_Ghost_Cells_C( model%ghost_C )
!$OMP END PARALLEL
            END IF
            CALL Lap( t_refine )
        END IF

        moments_done  = ( T == 0 )
        is_output     = ( MOD(T,nOutStep) == 0 .OR. T == Tf )
        is_diagnostic = ( is_output .OR. MOD(T,nDiagStep) == 0 )
        IF( is_diagnostic ) THEN
            IF( .NOT. moments_done ) THEN
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED( LeafBlks, nLeaf, model )
                DO m = 1, nLeaf
                    CALL model%Moments_g( LeafBlks(m) )
                END DO
!$OMP END PARALLEL DO
                moments_done = .True.
            END IF
            IF( is_output ) THEN
                CALL Write_Output( T*SQRT(Gy/D) )
                CALL Print_Diagnostics
            END IF
            CALL Bubble_Diagnostics( T*SQRT(Gy/D) )
            CALL Lap( t_Output )
            IF( T==Tf ) THEN
                CALL Close_Output
                CALL Bubble_Diagnostics( -1._wp )       ! close the file
                EXIT
            END IF
        END IF

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED( Blk, LeafBlks, nLeaf, model, moments_done, &
!$OMP   t_collision, t_Ghost_hg, t_advection )

        !== pressure and velocity from g, then the collision, block by block
!$OMP DO SCHEDULE(STATIC)
        DO m = 1, nLeaf
            n = LeafBlks(m)
            IF( moments_done ) THEN
                CALL model%Collide( n )
            ELSE
                CALL model%Moments_and_Collide( n )
            END IF
        END DO
!$OMP END DO
!$OMP MASKED
        CALL Lap( t_collision )
!$OMP END MASKED

        CALL Fill_Ghost_Cells_hg
!$OMP MASKED
        CALL Lap( t_Ghost_hg )
!$OMP END MASKED

        !== streaming, and the phase field C = sum(h) while h of the block is in cache
!$OMP DO SCHEDULE(STATIC)
        DO m = 1, nLeaf
            n = LeafBlks(m)
            cfl = 1._wp/Blk(n)%dx
            IF( cfl > 0.9 ) THEN
                CALL Perfect_Shift_AMR( Blk(n)%h )
                CALL Perfect_Shift_AMR( Blk(n)%g )
            ELSE
                CALL Streaming_LaxWendroff_AMR( cfl, Blk(n)%h )
                CALL Streaming_LaxWendroff_AMR( cfl, Blk(n)%g )
            END IF
            CALL Macroscopic_Values_h_AMR( Blk(n)%h, Blk(n)%C )
        END DO
!$OMP END DO
!$OMP MASKED
        CALL Lap( t_advection )
!$OMP END MASKED

!$OMP END PARALLEL

        cellUpdates = cellUpdates + INT(nLeaf,8)*nx*ny
        IF( MOD(T+1,nProgStep) == 0 .AND. T+1 < Tf ) CALL Progress( T+1 )

    END DO

    PRINT '(/A)', '------------------------------'
    PRINT '(A,F12.3)',   ' t_Output:        ', DBLE(t_Output   )/Rate
    PRINT '(A,F12.3,A)', ' t_collision:     ', DBLE(t_collision)/Rate, '   (moments of g and collision)'
    PRINT '(A,F12.3,A)', ' t_advection:     ', DBLE(t_advection)/Rate, '   (streaming and C)'
    PRINT '(A,F12.3)',   ' t_Ghost_C :      ', DBLE(t_Ghost_C  )/Rate
    PRINT '(A,F12.3)',   ' t_Ghost_hg:      ', DBLE(t_Ghost_hg )/Rate
    PRINT '(A,F12.3)',   ' t_refine:        ', DBLE(t_refine   )/Rate

CONTAINS

    !== progress after step Tdone: the speed over the last interval in million
    !   lattice-cell updates per second (MLUPS; every leaf cell counts once per
    !   step whatever its level) and the time left estimated from the average
    !   time per step so far (it changes with the number of refined blocks)
    SUBROUTINE Progress( Tdone )
        INTEGER, INTENT(IN) :: Tdone
        INTEGER(8) :: tnow
        REAL(dp)   :: elapsed, mlups

        CALL SYSTEM_CLOCK( tnow )
        elapsed = REAL( tnow - tp0, dp ) / Rate
        mlups   = REAL( cellUpdates, dp ) / MAX( REAL(tnow-tpLast,dp)/Rate, 1.e-9_dp ) / 1.e6_dp
        PRINT '(A,I9,A,I9,A,F5.1,A,F8.3,A,I9,A,F8.2,A,F7.2,4A)', ' step', Tdone, ' /', Tf, &
            ' (', 100._dp*Tdone/Tf, '%)  t* =', Tdone*SQRT(Gy/D), '  leaf cells', nLeaf*nx*ny, &
            '  MLUPS', mlups, ' (', mlups/OMP_GET_MAX_THREADS(), ' per thread)  elapsed ', &
            Hms(elapsed), '  left ~', Hms( elapsed/Tdone*(Tf-Tdone) )
        FLUSH( OUTPUT_UNIT )
        tpLast = tnow
        cellUpdates = 0
    END SUBROUTINE Progress

    !== seconds as h:mm:ss
    FUNCTION Hms( sec ) RESULT( str )
        REAL(dp), INTENT(IN) :: sec
        CHARACTER(len=10) :: str
        INTEGER :: isec
        isec = NINT( sec )
        WRITE( str, '(I4,":",I2.2,":",I2.2)' ) isec/3600, MOD(isec/60,60), MOD(isec,60)
    END FUNCTION Hms

    !== add the time elapsed since the last call to the accumulator tacc
    SUBROUTINE Lap( tacc )
        INTEGER(8), INTENT(INOUT) :: tacc
        INTEGER(8) :: tnow
        CALL SYSTEM_CLOCK( tnow )
        tacc = tacc + (tnow - t1)
        t1 = tnow
    END SUBROUTINE Lap

END SUBROUTINE Main_Loop

!******************************************************************************

!== rise of the bubble: centroid y_c and rise velocity u_r of the light phase
!   (weight 1-C), and the rise Reynolds number Re_r = rho_h u_r D / mu_h,
!   appended to bubble_rise.dat; Time < 0 closes the file.
!   Cells with C >= 0.9 are left out: faint tails of C slightly below 1 in
!   coarse cells far from the bubble (wake, coarse-fine edges) carry a large
!   area and change whenever blocks are refined or derefined, which makes
!   y_c jump.
SUBROUTINE Bubble_Diagnostics( Time )
    USE Parameters_Mod, ONLY: Blk, nLeaf, LeafBlks, nx, ny
    USE D2Q9_Mod, ONLY: Rhoh, D
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: Time

    INTEGER, SAVE :: u
    LOGICAL, SAVE :: IsOpen = .False.
    INTEGER  :: m, n, i, j, dx
    REAL(dp) :: A, Ay, Av, wgt

    IF( Time < 0 ) THEN
        IF( IsOpen ) CLOSE( u )
        IsOpen = .False.
        RETURN
    END IF
    IF( .NOT. IsOpen ) THEN
        OPEN( NEWUNIT=u, FILE='bubble_rise.dat', STATUS='REPLACE', ACTION='WRITE' )
        IsOpen = .True.
        WRITE( u, '(A)' ) '#      t*           y_c          u_r          Re_r'
    END IF

    A  = 0
    Ay = 0
    Av = 0
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED( Blk, LeafBlks, nLeaf ) REDUCTION(+:A,Ay,Av)
    DO m = 1, nLeaf
        n  = LeafBlks(m)
        dx = Blk(n)%dx
        DO j = 1, ny
        DO i = 1, nx
            IF( Blk(n)%C(i,j) >= 0.9_wp ) CYCLE
            wgt = REAL( 1._wp - MAX( Blk(n)%C(i,j), 0._wp ), dp ) * dx**2
            A   = A  + wgt
            Ay  = Ay + wgt * ( Blk(n)%J*dx*ny + dx*(j-0.5_dp) )
            Av  = Av + wgt * Blk(n)%Uy(i,j)
        END DO
        END DO
    END DO
!$OMP END PARALLEL DO

    WRITE( u, '(4ES13.5)' ) Time, Ay/A, Av/A, Rhoh*(Av/A)*D/mu_h
    FLUSH( u )

END SUBROUTINE Bubble_Diagnostics

END MODULE Rising_Bubble_Mod
