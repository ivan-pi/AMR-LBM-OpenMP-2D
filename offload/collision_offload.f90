!================================================
!   OpenMP target-offload prototype of the collision kernel
!
!   Stand-alone test of the most expensive stage of the solver
!   (Collision_h_g_AMR, about half of the time step) on a GPU.
!   The block data are kept in a flat "block pool" layout
!   (last index = block), which is what the solver would need for
!   offloading: the BLOCK derived type with allocatable components
!   and pointer-linked neighbours cannot be mapped efficiently.
!
!   The kernel is run on the host (OpenMP threads) and on the device,
!   results are compared, and timings are reported for
!     - the kernel with data resident on the device, and
!     - the kernel plus copying all fields to/from the device
!       (the cost if only this stage were offloaded).
!
!   Build:   cmake -B build-ifx -DCMAKE_Fortran_COMPILER=ifx -DAMR_BUILD_OFFLOAD=ON
!            cmake --build build-ifx --target collision_offload
!            (block size and precision follow -DAMR_BLOCK_SIZE / -DAMR_PRECISION)
!   Run:     ./build-ifx/offload/collision_offload [nblocks] [nrepeat]
!
!   Note: the Arc A7xx GPUs have no FP64 hardware; double precision
!   needs the driver's FP64 emulation, e.g.
!      export OverrideDefaultFP64Settings=1 IGC_EnableDPEmulation=1
!================================================
MODULE Offload_Pool_Mod
    USE Kinds_Mod, ONLY: wp
    IMPLICIT NONE

    !== local copies of the D2Q9 vectors and weights of D2Q9_Mod (ifx 2026.1
    !   cannot link device code that indexes PARAMETER arrays of another module)
    !   and of the default mobility (a run-time variable in D2Q9_Mod)
    INTEGER,  PARAMETER :: ex(9) = [0, 1, 0,-1, 0, 1,-1,-1, 1]
    INTEGER,  PARAMETER :: ey(9) = [0, 0, 1, 0,-1, 1, 1,-1,-1]
    REAL(wp), PARAMETER :: Wa(9) = [16,4, 4, 4, 4, 1, 1, 1, 1] / 36._wp
    REAL(wp), PARAMETER :: M0 = 0.03_wp, wc = 1._wp/(0.5_wp + M0*3._wp)

    INTEGER, PARAMETER :: nx = AMR_BLOCK_SIZE, ny = AMR_BLOCK_SIZE

    ! physical parameters of the rising-bubble case (Set_Rising_Bubble_Parameters)
    TYPE Phys
        REAL(wp) :: W, Beta, kappa, taul, tauh, Rhol, Rhoh, Gy
    END TYPE Phys

CONTAINS

    !== collision of one cell (i,j) of block b: same math as Collision_h_g_AMR
    !   (same layout h(i,j,k), clamped relaxation rate, reduced MRT)
    SUBROUTINE Collide_Cell( i, j, dx, C, P, Ux, Uy, h, g, ph )
    !$omp declare target
        INTEGER,    INTENT(IN)    :: i, j, dx
        REAL(wp),   INTENT(IN)    :: C(-1:nx+2,-1:ny+2), P(nx,ny), Ux(nx,ny), Uy(nx,ny)
        REAL(wp),   INTENT(INOUT) :: h(0:nx+1,0:ny+1,9), g(0:nx+1,0:ny+1,9)
        TYPE(Phys), INTENT(IN)    :: ph

        INTEGER  :: k
        REAL(wp) :: eDc(9), eDm(9), Gam(9), heq(9), geq(9), hlp(9)
        REAL(wp) :: DcDx, DcDy, DmDx, DmDy, ni, nj, tmp, U, V, U2, eU, eN
        REAL(wp) :: D2C, mu, dRho3, Rho, Fb, x, s9, rdx, m8, m9

        rdx   = 1._wp / dx
        dRho3 = (ph%Rhoh - ph%Rhol)/3._wp

        D2C = ( C(i-1,j-1)+C(i+1,j-1)+C(i-1,j+1)+C(i+1,j+1) &
            +4*(C(i  ,j-1)+C(i-1,j  )+C(i+1,j  )+C(i  ,j+1)) - 20*C(i,j) )/6
        mu  = 4._wp*ph%Beta*C(i,j)*(C(i,j)-1._wp)*(C(i,j)-0.5_wp) - ph%kappa*D2C/dx**2

        DO k = 1, 9
            eDc(k) = ( C(i+ex(k),j+ey(k)) - C(i-ex(k),j-ey(k)) ) * 0.5_wp*rdx
            eDm(k) = ( - C(i-ex(k),j-ey(k)) - 3._wp*C(i,j) + 5._wp*C(i+ex(k),j+ey(k)) &
                       - C(i+2*ex(k),j+2*ey(k)) ) * 0.25_wp*rdx
        END DO
        DcDx = ( eDc(2)-eDc(4) )/3._wp + ( eDc(6)-eDc(7)-eDc(8)+eDc(9) )/12._wp
        DcDy = ( eDc(3)-eDc(5) )/3._wp + ( eDc(6)+eDc(7)-eDc(8)-eDc(9) )/12._wp
        DmDx = ( eDm(2)-eDm(4) )/3._wp + ( eDm(6)-eDm(7)-eDm(8)+eDm(9) )/12._wp
        DmDy = ( eDm(3)-eDm(5) )/3._wp + ( eDm(6)+eDm(7)-eDm(8)-eDm(9) )/12._wp

        tmp = SQRT( DcDx**2 + DcDy**2 + 1.e-32_wp )
        ni = DcDx / tmp
        nj = DcDy / tmp

        U = Ux(i,j); V = Uy(i,j); U2 = U*U + V*V
        tmp = ( 12*(1-C(i,j)) )/ph%W * M0
        DO k = 1, 9
            eU = ex(k)*U  + ey(k)*V
            eN = ex(k)*ni + ey(k)*nj
            Gam(k) = Wa(k) * ( 1._wp + eU*(3._wp + 4.5_wp*eU) - 1.5_wp*U2 )
            heq(k) = C(i,j) * ( Gam(k) + tmp * Wa(k) * eN )
        END DO

        h(i,j,:) = h(i,j,:) * (1._wp-wc) + heq(:)*wc

        Rho = ph%Rhol + C(i,j) * (ph%Rhoh - ph%Rhol)
        Fb  = (ph%Rhoh - Rho) * ph%Gy
        DO k = 1, 9
            hlp(k) = ( (Gam(k)-Wa(k))*dRho3 + Gam(k)*mu ) * ( eDc(k) - (U*DcDx + V*DcDy) ) &
                   + (ey(k)-V) * Gam(k) * Fb
            geq(k) = P(i,j)*Wa(k) + (Gam(k)-Wa(k))*Rho/3._wp - 0.5_wp*hlp(k)
            hlp(k) = ( (Gam(k)-Wa(k))*dRho3 + Gam(k)*mu ) * ( eDm(k) - (U*DmDx + V*DmDy) ) &
                   + (ey(k)-V) * Gam(k) * Fb
        END DO

        !== relaxation rate: 1/tau_fluid interpolated in the clamped C, tau = 1/x + 1/2
        x  = MIN( MAX( C(i,j), 0._wp ), 1._wp ) * (1._wp/ph%tauh - 1._wp/ph%taul) + 1._wp/ph%taul
        s9 = x / ( 1._wp + 0.5_wp*x )

        !== reduced MRT: only the stress moments m8, m9 relax with s9
        m8 = ( (g(i,j,2)-geq(2)) - (g(i,j,3)-geq(3)) + (g(i,j,4)-geq(4)) - (g(i,j,5)-geq(5)) ) * (1._wp-s9) * 0.25_wp
        m9 = ( (g(i,j,6)-geq(6)) - (g(i,j,7)-geq(7)) + (g(i,j,8)-geq(8)) - (g(i,j,9)-geq(9)) ) * (1._wp-s9) * 0.25_wp
        g(i,j,1) = geq(1)      + hlp(1)
        g(i,j,2) = geq(2) + m8 + hlp(2)
        g(i,j,3) = geq(3) - m8 + hlp(3)
        g(i,j,4) = geq(4) + m8 + hlp(4)
        g(i,j,5) = geq(5) - m8 + hlp(5)
        g(i,j,6) = geq(6) + m9 + hlp(6)
        g(i,j,7) = geq(7) - m9 + hlp(7)
        g(i,j,8) = geq(8) + m9 + hlp(8)
        g(i,j,9) = geq(9) - m9 + hlp(9)

    END SUBROUTINE Collide_Cell

END MODULE Offload_Pool_Mod

!******************************************************************************

PROGRAM Collision_Offload
    USE Kinds_Mod, ONLY: wp, dp
    USE Offload_Pool_Mod, ONLY: nx, ny, Wa, Phys, Collide_Cell
    USE OMP_LIB, ONLY: OMP_GET_WTIME, OMP_GET_MAX_THREADS, OMP_GET_NUM_DEVICES
    IMPLICIT NONE

    INTEGER  :: nb, nrep, b, i, j, k, r, nargs
    CHARACTER(len=32) :: arg
    INTEGER,  ALLOCATABLE :: dx(:)
    REAL(wp), ALLOCATABLE :: C(:,:,:), P(:,:,:), Ux(:,:,:), Uy(:,:,:)
    REAL(wp), ALLOCATABLE :: h(:,:,:,:), g(:,:,:,:), hd(:,:,:,:), gd(:,:,:,:)
    REAL(wp), ALLOCATABLE :: h0(:,:,:,:), g0(:,:,:,:)
    TYPE(Phys) :: ph
    REAL(wp) :: x, y, R0, muh, mul, Sigma, D
    REAL(dp) :: t0, t_host, t_dev, t_devx, err

    nb = 2000; nrep = 20
    nargs = COMMAND_ARGUMENT_COUNT()
    IF( nargs >= 1 ) THEN; CALL GET_COMMAND_ARGUMENT(1,arg); READ(arg,*) nb;   END IF
    IF( nargs >= 2 ) THEN; CALL GET_COMMAND_ARGUMENT(2,arg); READ(arg,*) nrep; END IF

    !== rising-bubble parameters (as in Set_Rising_Bubble_Parameters, L0 = 256)
    D  = 2*256/10._wp
    ph%Rhol = 0.001_wp; ph%Rhoh = 1; ph%Gy = 0.000005_wp; ph%W = 5
    Sigma = ph%Gy*(ph%Rhoh-ph%Rhol)*D**2/40
    muh   = SQRT(ph%Gy*ph%Rhoh*(ph%Rhoh-ph%Rhol)*D**3)/40; mul = muh/100
    ph%tauh = 3*(muh/ph%Rhoh); ph%taul = 3*(mul/ph%Rhol)
    ph%Beta = 12*Sigma/ph%W;  ph%kappa = 1.5_wp*Sigma*ph%W

    ALLOCATE( dx(nb), C(-1:nx+2,-1:ny+2,nb), P(nx,ny,nb), Ux(nx,ny,nb), Uy(nx,ny,nb) )
    ALLOCATE( h(0:nx+1,0:ny+1,9,nb), g(0:nx+1,0:ny+1,9,nb) )

    !== synthetic data: an interface crossing every block, small velocities
    R0 = 4
    DO b = 1, nb
        dx(b) = 2**MOD(b,4)
        DO j = -1, ny+2
        DO i = -1, nx+2
            x = i - 0.5_wp*nx + 0.3_wp*MOD(b,7); y = j - 0.5_wp*ny
            C(i,j,b) = 0.5_wp + 0.5_wp*TANH( 2*(SQRT(x*x+y*y)-R0)/ph%W )
        END DO
        END DO
        DO j = 1, ny
        DO i = 1, nx
            Ux(i,j,b) = 1.e-3_wp*SIN(0.1_wp*(i+b)); Uy(i,j,b) = 1.e-3_wp*COS(0.1_wp*(j+b))
            P(i,j,b)  = 1.e-4_wp*C(i,j,b)
        END DO
        END DO
        DO k = 1, 9
            h(:,:,k,b) = 0.1_wp*Wa(k)
            g(:,:,k,b) = 0.2_wp*Wa(k)
        END DO
    END DO
    hd = h; gd = g; h0 = h; g0 = g

    PRINT '(A,I0,A,I0,A,I0)', 'blocks: ', nb, '   cells: ', nb*nx*ny, '   repetitions: ', nrep
    PRINT '(A,I0)', 'real kind (bytes): ', STORAGE_SIZE(1._wp)/8
    PRINT '(A,I0,A,I0)', 'host threads: ', OMP_GET_MAX_THREADS(), '   devices: ', OMP_GET_NUM_DEVICES()

    !== host
    t0 = OMP_GET_WTIME()
    DO r = 1, nrep
!$omp parallel do collapse(3) schedule(static)
        DO b = 1, nb
        DO j = 1, ny
        DO i = 1, nx
            CALL Collide_Cell( i, j, dx(b), C(:,:,b), P(:,:,b), Ux(:,:,b), Uy(:,:,b), h(:,:,:,b), g(:,:,:,b), ph )
        END DO
        END DO
        END DO
    END DO
    t_host = (OMP_GET_WTIME() - t0)/nrep

    !== device, data resident
!$omp target enter data map(to: dx, C, P, Ux, Uy, hd, gd)
    ! warm-up (JIT compilation of the kernel), not timed; then reset the state
    CALL Device_Kernel
    hd = h0; gd = g0
!$omp target update to(hd, gd)
    t0 = OMP_GET_WTIME()
    DO r = 1, nrep
        CALL Device_Kernel
    END DO
    t_dev = (OMP_GET_WTIME() - t0)/nrep
!$omp target update from(hd, gd)

    err = MAX( MAXVAL(ABS(REAL(hd-h,dp))), MAXVAL(ABS(REAL(gd-g,dp))) )

    !== device, with all fields copied in and out every call
    t0 = OMP_GET_WTIME()
    DO r = 1, nrep
!$omp target update to(C, P, Ux, Uy, hd, gd)
        CALL Device_Kernel
!$omp target update from(hd, gd)
    END DO
    t_devx = (OMP_GET_WTIME() - t0)/nrep
!$omp target exit data map(delete: dx, C, P, Ux, Uy, hd, gd)

    PRINT '(/A,ES10.3)', 'max |device - host| after all repetitions: ', err
    PRINT '(A,F10.3,A)', 'host   kernel              : ', 1.e3_dp*t_host, ' ms'
    PRINT '(A,F10.3,A)', 'device kernel (resident)   : ', 1.e3_dp*t_dev,  ' ms'
    PRINT '(A,F10.3,A)', 'device kernel + transfers  : ', 1.e3_dp*t_devx, ' ms'
    PRINT '(A,F10.2,A)', 'host   throughput          : ', nb*nx*ny/t_host/1.e6_dp, ' Mcells/s'
    PRINT '(A,F10.2,A)', 'device throughput          : ', nb*nx*ny/t_dev /1.e6_dp, ' Mcells/s'

CONTAINS

    SUBROUTINE Device_Kernel
!$omp target teams distribute parallel do collapse(3)
        DO b = 1, nb
        DO j = 1, ny
        DO i = 1, nx
            CALL Collide_Cell( i, j, dx(b), C(:,:,b), P(:,:,b), Ux(:,:,b), Uy(:,:,b), hd(:,:,:,b), gd(:,:,:,b), ph )
        END DO
        END DO
        END DO
    END SUBROUTINE Device_Kernel

END PROGRAM Collision_Offload
