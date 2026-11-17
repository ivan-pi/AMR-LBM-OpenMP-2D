!================================================
!   Ghost cells of the phase field C (ng = 1 or 2 layers)
!================================================
MODULE Ghost_C_Mod
    USE Kinds_Mod, ONLY: wp
    USE Interpolations_Mod, ONLY: Biquadratic_from_Fine, Biquadratic_from_Coarse
    USE BC_Mod, ONLY: Symmetric_C_3, Symmetric_C_4
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Fill_Ghost_Cells_C

CONTAINS


SUBROUTINE Fill_Ghost_Cells_C( ng )
    !== fill ng ghost layers of C on all leaf blocks, level by level (coarse first)
    USE Blocks_Mod, ONLY: LevelMin, LevelMax, LeafByLevel, LevelStart, Blk_C, Neighbor_of
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ng

    INTEGER :: L, I, m, n, Nghbr(8)

    DO L = LevelMin, LevelMax
!$OMP DO schedule(static)
        DO m = LevelStart(L), LevelStart(L+1)-1
            n = LeafByLevel(m)

            DO I = 1, 8
                Nghbr(I) = Neighbor_of(n, I)
            END DO

            IF( Nghbr(1) > 0 ) CALL Fill_Left_C  ( n, L, Nghbr(1), ng )
            IF( Nghbr(2) > 0 ) CALL Fill_Right_C ( n, L, Nghbr(2), ng )
            IF( Nghbr(3) > 0 ) CALL Fill_Bottom_C( n, L, Nghbr(3), ng )
            IF( Nghbr(4) > 0 ) CALL Fill_Top_C   ( n, L, Nghbr(4), ng )

            IF( Nghbr(5) > 0 ) CALL Fill_Left_Bottom_C ( n, L, Nghbr(5), ng )
            IF( Nghbr(6) > 0 ) CALL Fill_Left_Top_C    ( n, L, Nghbr(6), ng )
            IF( Nghbr(7) > 0 ) CALL Fill_Right_Bottom_C( n, L, Nghbr(7), ng )
            IF( Nghbr(8) > 0 ) CALL Fill_Right_Top_C   ( n, L, Nghbr(8), ng )

            !   IF( Nghbr(1) < 0 ) apply BC ...
            !   IF( Nghbr(2) < 0 ) apply BC ...
            IF( Nghbr(3) < 0 ) CALL Symmetric_C_3( ng, Blk_C(:,:,n) )
            IF( Nghbr(4) < 0 ) CALL Symmetric_C_4( ng, Blk_C(:,:,n) )

        END DO
!$OMP END DO
    END DO

END

!******************************************************************************

SUBROUTINE Fill_Left_C( n, L, N1, ng )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk, Blk_C
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, L, N1, ng

    INTEGER :: L_nbr, child, i, j, Y
    REAL(wp) :: C(2+ng,2*ny), C2(0:2,0:ny/2+1), Cg(ng,ny)

    L_nbr = Blk(N1)%Level

    IF( L == L_nbr ) THEN

        IF( Blk(N1)%IsLeaf ) THEN
            DO j=1,ny
                Blk_C(1-ng:0,j, n) = Blk_C(nx+1-ng:nx,j, N1)
            END DO
        ELSE
        !== saving neighbor's children data to temporary arrays
            child = Blk(N1)%child(1)
            DO J = 1,    ny
            DO I = nx-1-ng, nx
                C(nx-I+1,J   ) = Blk_C(I,J, child+1)
                C(nx-I+1,J+ny) = Blk_C(I,J, child+3)
            END DO
            END DO
        !========================================================

            CALL Biquadratic_from_Fine( 2+ng,2*ny, ng,ny, C, Cg )
            Blk_C(0:1-ng:-1,1:ny, n) = Cg

        END IF

    ELSEIF( L == L_nbr + 1 ) THEN

        Y = (ny/2) * MOD(Blk(n)%J,2)

        DO j = 0, ny/2 + 1
            C2(2:0:-1,j) = Blk_C(nx-1:nx+1,j+Y, N1)
        END DO

        CALL Biquadratic_from_Coarse( 3,ny/2+2, ng,ny, C2, Cg )
        Blk_C(0:1-ng:-1,1:ny, n) = Cg

    ELSE
        PRINT '(A,4I5)', " Error in Fill_Left_C: ", n, N1, L, L_nbr
    END IF

END

!******************************************************************************

SUBROUTINE Fill_Right_C( n, L, N2, ng )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk, Blk_C
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, L, N2, ng

    INTEGER :: L_nbr, child, i, j, Y
    REAL(wp) :: C(2+ng,2*ny), C2(0:2,0:ny/2+1), Cg(ng,ny)

    L_nbr = Blk(N2)%Level

    IF( L == L_nbr ) THEN

        IF( Blk(N2)%IsLeaf ) THEN
            DO j=1,ny
                Blk_C(nx+1:nx+ng,j, n) = Blk_C(1:ng,j, N2)
            END DO
        ELSE
        !== saving neighbor's children data to temporary arrays
            child = Blk(N2)%child(1)
            DO J = 1, ny
            DO I = 1, 2+ng
                C(I,J   ) = Blk_C(I,J, child  )
                C(I,J+ny) = Blk_C(I,J, child+2)
            END DO
            END DO
        !==============================================

            CALL Biquadratic_from_Fine( 2+ng,2*ny, ng,ny, C, Cg )
            Blk_C(nx+1:nx+ng,1:ny, n) = Cg

        END IF

    ELSEIF( L == L_nbr + 1 ) THEN

        Y = (ny/2) * MOD(Blk(n)%J,2)

        DO j = 0, ny/2 + 1
            C2(0:2,j) = Blk_C(0:2,j+Y, N2)
        END DO

        CALL Biquadratic_from_Coarse( 3,ny/2+2, ng,ny, C2, Cg )
        Blk_C(nx+1:nx+ng,1:ny, n) = Cg

    ELSE
        PRINT '(A,4I5)', " Error in Fill_Right_C: ", n, N2, L, L_nbr
    END IF

END

!******************************************************************************

SUBROUTINE Fill_Bottom_C( n, L, N3, ng )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk, Blk_C
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, L, N3, ng

    INTEGER :: L_nbr, child, i, j, X
    REAL(wp) :: C(2*nx,2+ng), C2(0:nx/2+1,0:2), Cg(nx,ng)

    L_nbr = Blk(N3)%Level

    IF( L == L_nbr ) THEN

        IF( Blk(N3)%IsLeaf ) THEN
            DO i=1,nx
                Blk_C(i,1-ng:0, n) = Blk_C(i,ny+1-ng:ny, N3)
            END DO
        ELSE
        !== saving neighbor's children data to temporary arrays
            child = Blk(N3)%child(1)
            DO J = ny-1-ng, ny
            DO I = 1,    nx
                C(I   ,ny-J+1) = Blk_C(I,J, child+2)
                C(I+nx,ny-J+1) = Blk_C(I,J, child+3)
            END DO
            END DO
        !============================================

            CALL Biquadratic_from_Fine( 2*nx,2+ng, nx,ng, C, Cg )
            Blk_C(1:nx,0:1-ng:-1, n) = Cg

        END IF

    ELSEIF( L == L_nbr + 1 ) THEN

        X = (nx/2) * MOD(Blk(n)%I,2)

        DO i = 0, nx/2 + 1
            C2(i,2:0:-1) = Blk_C(i+X,ny-1:ny+1, N3)
        END DO

        CALL Biquadratic_from_Coarse( nx/2+2,3, nx,ng, C2, Cg )
        Blk_C(1:nx,0:1-ng:-1, n) = Cg

    ELSE
        PRINT '(A,4I5)', " Error in Fill_Bottom_C: ", n, N3, L, L_nbr
    END IF

END

!******************************************************************************

SUBROUTINE Fill_Top_C( n, L, N4, ng )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk, Blk_C
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, L, N4, ng

    INTEGER :: L_nbr, child, i, j, X
    REAL(wp) :: C(2*nx,2+ng), C2(0:nx/2+1,0:2), Cg(nx,ng)

    L_nbr = Blk(N4)%Level

    IF( L == L_nbr ) THEN

        IF( Blk(N4)%IsLeaf ) THEN
            DO i=1,nx
                Blk_C(i,ny+1:ny+ng, n) = Blk_C(i,1:ng, N4)
            END DO
        ELSE
        !== saving neighbor's children data to temporary arrays
            child = Blk(N4)%child(1)
            DO J = 1, 2+ng
            DO I = 1, nx
                C(I   ,J) = Blk_C(I,J, child  )
                C(I+nx,J) = Blk_C(I,J, child+1)
            END DO
            END DO
        !==============================================

            CALL Biquadratic_from_Fine( 2*nx,2+ng, nx,ng, C, Cg )
            Blk_C(1:nx,ny+1:ny+ng, n) = Cg

        END IF

    ELSEIF( L == L_nbr + 1 ) THEN

        X = (nx/2) * MOD(Blk(n)%I,2)

        DO i = 0, nx/2 + 1
            C2(i,0:2) = Blk_C(i+X,0:2, N4)
        END DO

        CALL Biquadratic_from_Coarse( nx/2+2,3, nx,ng, C2, Cg )
        Blk_C(1:nx,ny+1:ny+ng, n) = Cg

    ELSE
        PRINT '(A,4I5)', " Error in Fill_Top_C: ", n, N4, L, L_nbr
    END IF

END

!******************************************************************************

SUBROUTINE Fill_Left_Bottom_C( n, L, N5, ng )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk, Blk_C
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, L, N5, ng

    INTEGER :: L_nbr, child, i, j
    REAL(wp) :: C(2+ng,2+ng), C2(0:2,0:2), Cg(ng,ng)

    L_nbr = Blk(N5)%Level

    IF( L == L_nbr ) THEN

        IF( Blk(N5)%IsLeaf ) THEN
            Blk_C(1-ng:0,1-ng:0, n) = Blk_C(nx+1-ng:nx,ny+1-ng:ny, N5)
        ELSE
        !== saving neighbor's children data to temporary arrays
            child = Blk(N5)%child(1)
            DO J = ny-1-ng, ny
            DO I = nx-1-ng, nx
                C(nx-I+1,ny-J+1) = Blk_C(I,J, child+3)
            END DO
            END DO
        !========================================================

            CALL Biquadratic_from_Fine( 2+ng,2+ng, ng,ng, C, Cg )
            Blk_C(0:1-ng:-1,0:1-ng:-1, n) = Cg

        END IF

    ELSEIF( L == L_nbr + 1 ) THEN

        I = nx - (nx/2) * MOD(Blk(n)%I,2) - 1
        J = ny - (ny/2) * MOD(Blk(n)%J,2) - 1

        C2(2:0:-1,2:0:-1) = Blk_C(I:I+2,J:J+2, N5)

        CALL Biquadratic_from_Coarse( 3,3, ng,ng, C2, Cg )
        Blk_C(0:1-ng:-1,0:1-ng:-1, n) = Cg

    ELSE
        PRINT '(A,4I5)', " Error in Fill_Left_Bottom_C: ", n, N5, L, L_nbr
    END IF

END

!******************************************************************************

SUBROUTINE Fill_Left_Top_C( n, L, N6, ng )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk, Blk_C
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, L, N6, ng

    INTEGER :: L_nbr, child, i, j
    REAL(wp) :: C(2+ng,2+ng), C2(0:2,0:2), Cg(ng,ng)

    L_nbr = Blk(N6)%Level

    IF( L == L_nbr ) THEN

        IF( Blk(N6)%IsLeaf ) THEN
            Blk_C(1-ng:0,ny+1:ny+ng, n) = Blk_C(nx+1-ng:nx,1:ng, N6)
        ELSE
        !== saving neighbor's children data to temporary arrays
            child = Blk(N6)%child(1)
            DO J = 1,    2+ng
            DO I = nx-1-ng, nx
                C(nx-I+1,J) = Blk_C(I,J, child+1)
            END DO
            END DO
        !========================================================

            CALL Biquadratic_from_Fine( 2+ng,2+ng, ng,ng, C, Cg )
            Blk_C(0:1-ng:-1,ny+1:ny+ng, n) = Cg

        END IF

    ELSEIF( L == L_nbr + 1 ) THEN

        I = nx   - (nx/2) * MOD(Blk(n)%I,2) - 1
        J = ny/2 - (ny/2) * MOD(Blk(n)%J,2)

        C2(2:0:-1,0:2) = Blk_C(I:I+2,J:J+2, N6)

        CALL Biquadratic_from_Coarse( 3,3, ng,ng, C2, Cg )
        Blk_C(0:1-ng:-1,ny+1:ny+ng, n) = Cg

    ELSE
        PRINT '(A,4I5)', " Error in Fill_Left_Top_C: ", n, N6, L, L_nbr
    END IF

END

!******************************************************************************

SUBROUTINE Fill_Right_Bottom_C( n, L, N7, ng )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk, Blk_C
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, L, N7, ng

    INTEGER :: L_nbr, child, i, j
    REAL(wp) :: C(2+ng,2+ng), C2(0:2,0:2), Cg(ng,ng)

    L_nbr = Blk(N7)%Level

    IF( L == L_nbr ) THEN

        IF( Blk(N7)%IsLeaf ) THEN
            Blk_C(nx+1:nx+ng,1-ng:0, n) = Blk_C(1:ng,ny+1-ng:ny, N7)
        ELSE
        !== saving neighbor's children data to temporary arrays
            child = Blk(N7)%child(1)
            DO J = ny-1-ng, ny
            DO I = 1,    2+ng
                C(I,ny-J+1) = Blk_C(I,J, child+2)
            END DO
            END DO
        !========================================================

            CALL Biquadratic_from_Fine( 2+ng,2+ng, ng,ng, C, Cg )
            Blk_C(nx+1:nx+ng,0:1-ng:-1, n) = Cg

        END IF

    ELSEIF( L == L_nbr + 1 ) THEN

        I = nx/2 - (nx/2) * MOD(Blk(n)%I,2)
        J = ny   - (ny/2) * MOD(Blk(n)%J,2) - 1

        C2(0:2,2:0:-1) = Blk_C(I:I+2,J:J+2, N7)

        CALL Biquadratic_from_Coarse( 3,3, ng,ng, C2, Cg )
        Blk_C(nx+1:nx+ng,0:1-ng:-1, n) = Cg

    ELSE
        PRINT '(A,4I5)', " Error in Fill_Right_Bottom_C: ", n, N7, L, L_nbr
    END IF

END

!******************************************************************************

SUBROUTINE Fill_Right_Top_C( n, L, N8, ng )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk, Blk_C
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, L, N8, ng

    INTEGER :: L_nbr, child, i, j
    REAL(wp) :: C(2+ng,2+ng), C2(0:2,0:2), Cg(ng,ng)

    L_nbr = Blk(N8)%Level

    IF( L == L_nbr ) THEN

        IF( Blk(N8)%IsLeaf ) THEN
            Blk_C(nx+1:nx+ng,ny+1:ny+ng, n) = Blk_C(1:ng,1:ng, N8)
        ELSE
        !== saving neighbor's children data to temporary arrays
            child = Blk(N8)%child(1)
            DO J = 1, 2+ng
            DO I = 1, 2+ng
                C(I,J) = Blk_C(I,J, child)
            END DO
            END DO
        !========================================================

            CALL Biquadratic_from_Fine( 2+ng,2+ng, ng,ng, C, Cg )
            Blk_C(nx+1:nx+ng,ny+1:ny+ng, n) = Cg

        END IF

    ELSEIF( L == L_nbr + 1 ) THEN

        I = nx/2 - (nx/2) * MOD(Blk(n)%I,2)
        J = ny/2 - (ny/2) * MOD(Blk(n)%J,2)

        C2(0:2,0:2) = Blk_C(I:I+2,J:J+2, N8)

        CALL Biquadratic_from_Coarse( 3,3, ng,ng, C2, Cg )
        Blk_C(nx+1:nx+ng,ny+1:ny+ng, n) = Cg

    ELSE
        PRINT '(A,4I5)', " Error in Fill_Right_Top_C: ", n, N8, L, L_nbr
    END IF

END

END MODULE Ghost_C_Mod
