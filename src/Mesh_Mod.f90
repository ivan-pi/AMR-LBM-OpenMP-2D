!================================================
!   Block (quadtree) mesh: base blocks, children, leaf lists, diagnostics
!================================================
MODULE Mesh_Mod
    USE Kinds_Mod, ONLY: wp, dp
    USE Neighbors_Mod, ONLY: Find_Neighbors
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Initialize_Base_Block, Fix_Initial_Neighbors_3Blocks, Initial_Mesh, &
              Add_Child_Blocks, Array_of_Leaf_Blocks, Determine_Min_Max_Levels, &
              Print_Diagnostics

CONTAINS

SUBROUTINE Initialize_Base_Block( N )
    USE Parameters_Mod, ONLY: Blk, Nblks, L0, nx
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N

    INTEGER :: I

    Nblks = N

    Blk(N)%Level = 0
    Blk(N)%dx    = L0 / nx
    Blk(N)%I     = 0
    Blk(N)%J     = N-1  !N=1: the lower block

    ALLOCATE( Blk(N)%child(1:4) )
    Blk(N)%child(:)= -1         ! Children will be determined soon
    Blk(N)%Parent  = -1         ! This is the base block and has no parent

    !============ This is the numbering for the neighboring blocks
    !   6   4   8
    !   1       2
    !   5   3   7
    !============ This is the numbering for the neighboring blocks

    DO I = 1, 8
        ALLOCATE( Blk(N)%nghbr(I)%n )
        Blk(N)%nghbr(I)%n = -I  !External Boundaries
    END DO
    Blk(N)%nghbr(1)%n = N   !Periodic
    Blk(N)%nghbr(2)%n = N   !Periodic

    Blk(N)%IsLeaf  = .True.
    Blk(N)%hasLeaf = .False.
    Blk(N)%Refin   = .False.
    Blk(N)%Deref   = .False.

END SUBROUTINE Initialize_Base_Block

!******************************************************************************

SUBROUTINE Fix_Initial_Neighbors_3Blocks
    USE Parameters_Mod, ONLY: Blk
    IMPLICIT NONE

    Blk(1)%nghbr(4)%n = 2
    Blk(1)%nghbr(6)%n = 2
    Blk(1)%nghbr(8)%n = 2

    Blk(2)%nghbr(3)%n = 1
    Blk(2)%nghbr(4)%n = 3
    Blk(2)%nghbr(5)%n = 1
    Blk(2)%nghbr(6)%n = 3
    Blk(2)%nghbr(7)%n = 1
    Blk(2)%nghbr(8)%n = 3

    Blk(3)%nghbr(3)%n = 2
    Blk(3)%nghbr(5)%n = 2
    Blk(3)%nghbr(7)%n = 2

END

!******************************************************************************

SUBROUTINE Initial_Mesh( nLevels )
    USE Parameters_Mod, ONLY: Blk, Nblks
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nLevels

    INTEGER :: L, n

    DO L = 1, nLevels
        DO n = 1, Nblks     ! Nblks will change inside Sub:Add_Child_Blocks
            IF( Blk(n)%Level == L-1 ) THEN
                CALL Add_Child_Blocks( n, Nblks+1 )
            END IF
        END DO
    END DO

    CALL Array_of_Leaf_Blocks

END SUBROUTINE Initial_Mesh

!******************************************************************************

SUBROUTINE Add_Child_Blocks( par, FreeBlk )
    USE Parameters_Mod, ONLY: Blk, Nblks, MaxBlks, nx, ny
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: par, FreeBlk

    INTEGER :: N, grandpar, I, J

    N = FreeBlk - 1

    IF( N+4 > MaxBlks ) THEN
        CALL Print_Diagnostics
        PRINT '(/A/)', ' number of blocks exceeded MaxBlks! ... Program Terminated.'
        STOP
    END IF

    DO i = 1, 4
        Blk(par)%child(i) = N + i
    END DO
    Blk(par)%IsLeaf  = .False.
    Blk(par)%hasLeaf = .True.

    grandpar = Blk(par)%Parent
    IF( grandpar > 0 ) Blk(grandpar)%hasLeaf = .False.

    DO J = 0, 1
    DO I = 0, 1

        N = N + 1

        Blk(N)%Level = Blk(par)%Level + 1
        Blk(N)%dx    = Blk(par)%dx / 2

        Blk(N)%I = 2 * Blk(par)%I + I
        Blk(N)%J = 2 * Blk(par)%J + J

        Blk(N)%Parent = par

        ALLOCATE( Blk(N)%child(1:4) )
        Blk(N)%child(:) = N

        CALL Find_Neighbors( par, N, I, J )

        Blk(N)%IsLeaf  = .True.
        Blk(N)%hasLeaf = .False.
        Blk(N)%Refin   = .False.
        Blk(N)%Deref   = .False.

        ALLOCATE( Blk(N)%h(0:nx+1,0:ny+1,9) )
        ALLOCATE( Blk(N)%g(0:nx+1,0:ny+1,9) )
        ALLOCATE( Blk(N)%C(  -1:nx+2,-1:ny+2) )
        ALLOCATE( Blk(N)%P (nx,ny) )
        ALLOCATE( Blk(N)%Ux(nx,ny) )
        ALLOCATE( Blk(N)%Uy(nx,ny) )

        Blk(N)%h  = 0
        Blk(N)%g  = 0
        Blk(N)%C  = 0
        Blk(N)%P  = 0
        Blk(N)%Ux = 0
        Blk(N)%Uy = 0

    END DO
    END DO

    Nblks = MAX(Nblks,N)    ! because unallocated blocks might have been used

END SUBROUTINE Add_Child_Blocks

!******************************************************************************

SUBROUTINE Array_of_Leaf_Blocks
    USE Parameters_Mod, ONLY: Blk, Nblks, LeafBlks, nLeaf, LeafByLevel, LevelStart, maxLevels
    IMPLICIT NONE

    INTEGER :: n, i, L
    INTEGER :: tmpArray(Nblks), cnt(0:maxLevels+1)

    i = 0
    DO n = 1, Nblks
        IF( Blk(n)%IsLeaf ) THEN
            i = i + 1
            tmpArray(i) = n
        END IF
    END DO

    nLeaf = i
    IF( ALLOCATED(LeafBlks) ) DEALLOCATE( LeafBlks )
    ALLOCATE( LeafBlks(i) )

    LeafBlks(:) = tmpArray(1:i)

    !== counting sort of the leaves by level (stable: keeps block order per level)
    cnt = 0
    DO i = 1, nLeaf
        L = Blk(LeafBlks(i))%Level
        cnt(L) = cnt(L) + 1
    END DO
    LevelStart(0) = 1
    DO L = 0, maxLevels
        LevelStart(L+1) = LevelStart(L) + cnt(L)
    END DO

    IF( ALLOCATED(LeafByLevel) ) DEALLOCATE( LeafByLevel )
    ALLOCATE( LeafByLevel(nLeaf) )
    cnt(0:maxLevels) = LevelStart(0:maxLevels)
    DO i = 1, nLeaf
        n = LeafBlks(i)
        L = Blk(n)%Level
        LeafByLevel(cnt(L)) = n
        cnt(L) = cnt(L) + 1
    END DO

END

!******************************************************************************

SUBROUTINE Determine_Min_Max_Levels
    USE Parameters_Mod, ONLY: LevelMin, LevelMax, LevelStart, maxLevels
    IMPLICIT NONE

    INTEGER :: L

    LevelMin =  100
    LevelMax = -100
    DO L = 0, maxLevels
        IF( LevelStart(L+1) > LevelStart(L) ) THEN
            LevelMin = MIN( LevelMin, L )
            LevelMax = MAX( LevelMax, L )
        END IF
    END DO

END SUBROUTINE Determine_Min_Max_Levels

!******************************************************************************

!== one line of diagnostics: density range, maximum vertical velocity, total mass
SUBROUTINE Print_Diagnostics
    USE Parameters_Mod, ONLY: Blk, Nblks, LevelMin, LevelMax, nLeaf, LeafBlks, T, KK, nx, ny
    USE D2Q9_Mod,   ONLY: Rhol, Rhoh, Gy, D
    IMPLICIT NONE

    INTEGER :: m, n, dx, i, j
    REAL(wp) :: Rho, Rho_min, Rho_max, U_max
    REAL(dp) :: Mass    ! accumulated in double precision independent of wp

    Rho_min =  10
    Rho_max = -10
    Mass  = 0
    U_max = 0
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED( Blk, LeafBlks, nLeaf, Rhol, Rhoh ) &
!$OMP   REDUCTION(min:Rho_min) REDUCTION(max:Rho_max,U_max) REDUCTION(+:Mass)
    DO m = 1, nLeaf
        n = LeafBlks(m)

        dx = Blk(n)%dx
        DO j = 1, ny
        DO i = 1, nx
            Rho = Rhol + Blk(n)%C(i,j) * (Rhoh - Rhol)
            Rho_min = MIN( Rho_min, Rho )
            Rho_max = MAX( Rho_max, Rho )
            Mass    = Mass + REAL(Rho,dp) * dx**2
        END DO
        END DO

        U_max = MAX( U_max, MAXVAL(ABS(Blk(n)%Uy)) )
    END DO
!$OMP END PARALLEL DO
    WRITE(*,10) T*SQRT(Gy/D), Nblks, nLeaf, KK, Rho_min, Rho_max, U_max, Mass, LevelMin, LevelMax
10  FORMAT(F7.2,2I7,I6,2F11.5,F9.5,F12.1,2I4)

END SUBROUTINE Print_Diagnostics

END MODULE Mesh_Mod
