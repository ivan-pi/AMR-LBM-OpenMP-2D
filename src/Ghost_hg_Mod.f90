!================================================
!   Ghost cells of the distribution functions h and g
!================================================
MODULE Ghost_hg_Mod
    USE Kinds_Mod, ONLY: wp
    USE Interpolations_Mod, ONLY: Biquadratic_from_Fine_f, Biquadratic_from_Coarse_f
    USE BC_Mod, ONLY: Bounce_Back_3, Bounce_Back_4
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Fill_Ghost_Cells_hg

    !== Only the moving directions 2:9 get ghost values: the rest population
    !   (k=1) never streams. (Edge-specific subsets are not enough: the
    !   coarse-to-fine interpolation also reads the coarse neighbour's ghost
    !   rows and corners in the same direction.)
    !   Fine and coarse strips are gathered into small contiguous buffers,
    !   ordered from the interface outwards, interpolated, and stored.

CONTAINS


SUBROUTINE Fill_Ghost_Cells_hg
    USE Blocks_Mod, ONLY: Blk, LevelMin, LevelMax, LeafByLevel, LevelStart
    IMPLICIT NONE

    INTEGER :: L, I, m, n, Nghbr(8)

    DO L = LevelMin, LevelMax
!$OMP DO schedule(static)
        DO m = LevelStart(L), LevelStart(L+1)-1
            n = LeafByLevel(m)

            DO I = 1, 8
                Nghbr(I) = Blk(n)%nghbr(I)%n
            END DO

            IF( Nghbr(1) > 0 ) CALL Fill_Left_hg  ( n, L, Nghbr(1) )
            IF( Nghbr(2) > 0 ) CALL Fill_Right_hg ( n, L, Nghbr(2) )
            IF( Nghbr(3) > 0 ) CALL Fill_Bottom_hg( n, L, Nghbr(3) )
            IF( Nghbr(4) > 0 ) CALL Fill_Top_hg   ( n, L, Nghbr(4) )

            IF( Nghbr(5) > 0 ) CALL Fill_Left_Bottom_hg ( n, L, Nghbr(5) )
            IF( Nghbr(6) > 0 ) CALL Fill_Left_Top_hg    ( n, L, Nghbr(6) )
            IF( Nghbr(7) > 0 ) CALL Fill_Right_Bottom_hg( n, L, Nghbr(7) )
            IF( Nghbr(8) > 0 ) CALL Fill_Right_Top_hg   ( n, L, Nghbr(8) )

            !   IF( Nghbr(1) < 0 ) THEN
            !       apply BC ...
            !   END IF
            !   IF( Nghbr(2) < 0 ) THEN
            !       apply BC ...
            !   END IF
            IF( Nghbr(3) < 0 ) THEN
                CALL Bounce_Back_3( Blk(n)%h )
                CALL Bounce_Back_3( Blk(n)%g )
            END IF
            IF( Nghbr(4) < 0 ) THEN
                CALL Bounce_Back_4( Blk(n)%h )
                CALL Bounce_Back_4( Blk(n)%g )
            END IF

        END DO
!$OMP END DO
    END DO

END

!******************************************************************************

SUBROUTINE Fill_Left_hg( n, L, N1 )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk
    INTEGER, INTENT(IN) :: n, L, N1

    INTEGER :: child, Y
    REAL(wp) :: h(3,2*ny,2:9), H2(0:2,0:ny/2+1,2:9), bh(1,ny,2:9)
    REAL(wp) :: g(3,2*ny,2:9), G2(0:2,0:ny/2+1,2:9), bg(1,ny,2:9)

    IF( L == Blk(N1)%Level ) THEN

        IF( Blk(N1)%IsLeaf ) THEN
            Blk(n)%h(0,1:ny,2:9) = Blk(N1)%h(nx,1:ny,2:9)
            Blk(n)%g(0,1:ny,2:9) = Blk(N1)%g(nx,1:ny,2:9)
            RETURN
        END IF
        !== the neighbour is refined: restrict its children's strip
        child = Blk(N1)%child(1)
        h(3:1:-1,   1:  ny,:) = Blk(child+1)%h(nx-2:nx,1:ny,2:9)
        h(3:1:-1,ny+1:2*ny,:) = Blk(child+3)%h(nx-2:nx,1:ny,2:9)
        g(3:1:-1,   1:  ny,:) = Blk(child+1)%g(nx-2:nx,1:ny,2:9)
        g(3:1:-1,ny+1:2*ny,:) = Blk(child+3)%g(nx-2:nx,1:ny,2:9)
        CALL Biquadratic_from_Fine_f( 3,2*ny, 1,ny, h, bh, 8 )
        CALL Biquadratic_from_Fine_f( 3,2*ny, 1,ny, g, bg, 8 )

    ELSEIF( L == Blk(N1)%Level + 1 ) THEN

        !== the neighbour is coarser: interpolate its strip
        Y = (ny/2) * MOD(Blk(n)%J,2)
        H2(2:0:-1,:,:) = Blk(N1)%h(nx-1:nx+1,Y:Y+ny/2+1,2:9)
        G2(2:0:-1,:,:) = Blk(N1)%g(nx-1:nx+1,Y:Y+ny/2+1,2:9)
        CALL Biquadratic_from_Coarse_f( 3,ny/2+2, 1,ny, H2, bh, 8 )
        CALL Biquadratic_from_Coarse_f( 3,ny/2+2, 1,ny, G2, bg, 8 )

    ELSE
        PRINT '(A,4I5)', " Error in Fill_Left_hg: ", n, N1, L, Blk(N1)%Level
        RETURN
    END IF

    Blk(n)%h(0,1:ny,2:9) = bh(1,:,:)
    Blk(n)%g(0,1:ny,2:9) = bg(1,:,:)

END SUBROUTINE Fill_Left_hg

!******************************************************************************

SUBROUTINE Fill_Right_hg( n, L, N2 )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk
    INTEGER, INTENT(IN) :: n, L, N2

    INTEGER :: child, Y
    REAL(wp) :: h(3,2*ny,2:9), H2(0:2,0:ny/2+1,2:9), bh(1,ny,2:9)
    REAL(wp) :: g(3,2*ny,2:9), G2(0:2,0:ny/2+1,2:9), bg(1,ny,2:9)

    IF( L == Blk(N2)%Level ) THEN

        IF( Blk(N2)%IsLeaf ) THEN
            Blk(n)%h(nx+1,1:ny,2:9) = Blk(N2)%h(1,1:ny,2:9)
            Blk(n)%g(nx+1,1:ny,2:9) = Blk(N2)%g(1,1:ny,2:9)
            RETURN
        END IF
        !== the neighbour is refined: restrict its children's strip
        child = Blk(N2)%child(1)
        h(:,   1:  ny,:) = Blk(child  )%h(1:3,1:ny,2:9)
        h(:,ny+1:2*ny,:) = Blk(child+2)%h(1:3,1:ny,2:9)
        g(:,   1:  ny,:) = Blk(child  )%g(1:3,1:ny,2:9)
        g(:,ny+1:2*ny,:) = Blk(child+2)%g(1:3,1:ny,2:9)
        CALL Biquadratic_from_Fine_f( 3,2*ny, 1,ny, h, bh, 8 )
        CALL Biquadratic_from_Fine_f( 3,2*ny, 1,ny, g, bg, 8 )

    ELSEIF( L == Blk(N2)%Level + 1 ) THEN

        !== the neighbour is coarser: interpolate its strip
        Y = (ny/2) * MOD(Blk(n)%J,2)
        H2(:,:,:) = Blk(N2)%h(0:2,Y:Y+ny/2+1,2:9)
        G2(:,:,:) = Blk(N2)%g(0:2,Y:Y+ny/2+1,2:9)
        CALL Biquadratic_from_Coarse_f( 3,ny/2+2, 1,ny, H2, bh, 8 )
        CALL Biquadratic_from_Coarse_f( 3,ny/2+2, 1,ny, G2, bg, 8 )

    ELSE
        PRINT '(A,4I5)', " Error in Fill_Right_hg: ", n, N2, L, Blk(N2)%Level
        RETURN
    END IF

    Blk(n)%h(nx+1,1:ny,2:9) = bh(1,:,:)
    Blk(n)%g(nx+1,1:ny,2:9) = bg(1,:,:)

END SUBROUTINE Fill_Right_hg

!******************************************************************************

SUBROUTINE Fill_Bottom_hg( n, L, N3 )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk
    INTEGER, INTENT(IN) :: n, L, N3

    INTEGER :: child, X
    REAL(wp) :: h(2*nx,3,2:9), H2(0:nx/2+1,0:2,2:9), bh(nx,1,2:9)
    REAL(wp) :: g(2*nx,3,2:9), G2(0:nx/2+1,0:2,2:9), bg(nx,1,2:9)

    IF( L == Blk(N3)%Level ) THEN

        IF( Blk(N3)%IsLeaf ) THEN
            Blk(n)%h(1:nx,0,2:9) = Blk(N3)%h(1:nx,ny,2:9)
            Blk(n)%g(1:nx,0,2:9) = Blk(N3)%g(1:nx,ny,2:9)
            RETURN
        END IF
        !== the neighbour is refined: restrict its children's strip
        child = Blk(N3)%child(1)
        h(   1:  nx,3:1:-1,:) = Blk(child+2)%h(1:nx,ny-2:ny,2:9)
        h(nx+1:2*nx,3:1:-1,:) = Blk(child+3)%h(1:nx,ny-2:ny,2:9)
        g(   1:  nx,3:1:-1,:) = Blk(child+2)%g(1:nx,ny-2:ny,2:9)
        g(nx+1:2*nx,3:1:-1,:) = Blk(child+3)%g(1:nx,ny-2:ny,2:9)
        CALL Biquadratic_from_Fine_f( 2*nx,3, nx,1, h, bh, 8 )
        CALL Biquadratic_from_Fine_f( 2*nx,3, nx,1, g, bg, 8 )

    ELSEIF( L == Blk(N3)%Level + 1 ) THEN

        !== the neighbour is coarser: interpolate its strip
        X = (nx/2) * MOD(Blk(n)%I,2)
        H2(:,2:0:-1,:) = Blk(N3)%h(X:X+nx/2+1,ny-1:ny+1,2:9)
        G2(:,2:0:-1,:) = Blk(N3)%g(X:X+nx/2+1,ny-1:ny+1,2:9)
        CALL Biquadratic_from_Coarse_f( nx/2+2,3, nx,1, H2, bh, 8 )
        CALL Biquadratic_from_Coarse_f( nx/2+2,3, nx,1, G2, bg, 8 )

    ELSE
        PRINT '(A,4I5)', " Error in Fill_Bottom_hg: ", n, N3, L, Blk(N3)%Level
        RETURN
    END IF

    Blk(n)%h(1:nx,0,2:9) = bh(:,1,:)
    Blk(n)%g(1:nx,0,2:9) = bg(:,1,:)

END SUBROUTINE Fill_Bottom_hg

!******************************************************************************

SUBROUTINE Fill_Top_hg( n, L, N4 )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk
    INTEGER, INTENT(IN) :: n, L, N4

    INTEGER :: child, X
    REAL(wp) :: h(2*nx,3,2:9), H2(0:nx/2+1,0:2,2:9), bh(nx,1,2:9)
    REAL(wp) :: g(2*nx,3,2:9), G2(0:nx/2+1,0:2,2:9), bg(nx,1,2:9)

    IF( L == Blk(N4)%Level ) THEN

        IF( Blk(N4)%IsLeaf ) THEN
            Blk(n)%h(1:nx,ny+1,2:9) = Blk(N4)%h(1:nx,1,2:9)
            Blk(n)%g(1:nx,ny+1,2:9) = Blk(N4)%g(1:nx,1,2:9)
            RETURN
        END IF
        !== the neighbour is refined: restrict its children's strip
        child = Blk(N4)%child(1)
        h(   1:  nx,:,:) = Blk(child  )%h(1:nx,1:3,2:9)
        h(nx+1:2*nx,:,:) = Blk(child+1)%h(1:nx,1:3,2:9)
        g(   1:  nx,:,:) = Blk(child  )%g(1:nx,1:3,2:9)
        g(nx+1:2*nx,:,:) = Blk(child+1)%g(1:nx,1:3,2:9)
        CALL Biquadratic_from_Fine_f( 2*nx,3, nx,1, h, bh, 8 )
        CALL Biquadratic_from_Fine_f( 2*nx,3, nx,1, g, bg, 8 )

    ELSEIF( L == Blk(N4)%Level + 1 ) THEN

        !== the neighbour is coarser: interpolate its strip
        X = (nx/2) * MOD(Blk(n)%I,2)
        H2(:,:,:) = Blk(N4)%h(X:X+nx/2+1,0:2,2:9)
        G2(:,:,:) = Blk(N4)%g(X:X+nx/2+1,0:2,2:9)
        CALL Biquadratic_from_Coarse_f( nx/2+2,3, nx,1, H2, bh, 8 )
        CALL Biquadratic_from_Coarse_f( nx/2+2,3, nx,1, G2, bg, 8 )

    ELSE
        PRINT '(A,4I5)', " Error in Fill_Top_hg: ", n, N4, L, Blk(N4)%Level
        RETURN
    END IF

    Blk(n)%h(1:nx,ny+1,2:9) = bh(:,1,:)
    Blk(n)%g(1:nx,ny+1,2:9) = bg(:,1,:)

END SUBROUTINE Fill_Top_hg

!******************************************************************************

SUBROUTINE Fill_Left_Bottom_hg( n, L, N5 )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk
    INTEGER, INTENT(IN) :: n, L, N5

    INTEGER :: child, I, J
    REAL(wp) :: h(3,3,2:9), H2(0:2,0:2,2:9), bh(1,1,2:9)
    REAL(wp) :: g(3,3,2:9), G2(0:2,0:2,2:9), bg(1,1,2:9)

    IF( L == Blk(N5)%Level ) THEN

        IF( Blk(N5)%IsLeaf ) THEN
            Blk(n)%h(0,0,2:9) = Blk(N5)%h(nx,ny,2:9)
            Blk(n)%g(0,0,2:9) = Blk(N5)%g(nx,ny,2:9)
            RETURN
        END IF
        !== the neighbour is refined: restrict its children's strip
        child = Blk(N5)%child(1)
        h(3:1:-1,3:1:-1,:) = Blk(child+3)%h(nx-2:nx,ny-2:ny,2:9)
        g(3:1:-1,3:1:-1,:) = Blk(child+3)%g(nx-2:nx,ny-2:ny,2:9)
        CALL Biquadratic_from_Fine_f( 3,3, 1,1, h, bh, 8 )
        CALL Biquadratic_from_Fine_f( 3,3, 1,1, g, bg, 8 )

    ELSEIF( L == Blk(N5)%Level + 1 ) THEN

        !== the neighbour is coarser: interpolate its strip
        I = nx - (nx/2) * MOD(Blk(n)%I,2) - 1
        J = ny - (ny/2) * MOD(Blk(n)%J,2) - 1
        H2(2:0:-1,2:0:-1,:) = Blk(N5)%h(I:I+2,J:J+2,2:9)
        G2(2:0:-1,2:0:-1,:) = Blk(N5)%g(I:I+2,J:J+2,2:9)
        CALL Biquadratic_from_Coarse_f( 3,3, 1,1, H2, bh, 8 )
        CALL Biquadratic_from_Coarse_f( 3,3, 1,1, G2, bg, 8 )

    ELSE
        PRINT '(A,4I5)', " Error in Fill_Left_Bottom_hg: ", n, N5, L, Blk(N5)%Level
        RETURN
    END IF

    Blk(n)%h(0,0,2:9) = bh(1,1,:)
    Blk(n)%g(0,0,2:9) = bg(1,1,:)

END SUBROUTINE Fill_Left_Bottom_hg

!******************************************************************************

SUBROUTINE Fill_Left_Top_hg( n, L, N6 )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk
    INTEGER, INTENT(IN) :: n, L, N6

    INTEGER :: child, I, J
    REAL(wp) :: h(3,3,2:9), H2(0:2,0:2,2:9), bh(1,1,2:9)
    REAL(wp) :: g(3,3,2:9), G2(0:2,0:2,2:9), bg(1,1,2:9)

    IF( L == Blk(N6)%Level ) THEN

        IF( Blk(N6)%IsLeaf ) THEN
            Blk(n)%h(0,ny+1,2:9) = Blk(N6)%h(nx,1,2:9)
            Blk(n)%g(0,ny+1,2:9) = Blk(N6)%g(nx,1,2:9)
            RETURN
        END IF
        !== the neighbour is refined: restrict its children's strip
        child = Blk(N6)%child(1)
        h(3:1:-1,:,:) = Blk(child+1)%h(nx-2:nx,1:3,2:9)
        g(3:1:-1,:,:) = Blk(child+1)%g(nx-2:nx,1:3,2:9)
        CALL Biquadratic_from_Fine_f( 3,3, 1,1, h, bh, 8 )
        CALL Biquadratic_from_Fine_f( 3,3, 1,1, g, bg, 8 )

    ELSEIF( L == Blk(N6)%Level + 1 ) THEN

        !== the neighbour is coarser: interpolate its strip
        I = nx   - (nx/2) * MOD(Blk(n)%I,2) - 1
        J = ny/2 - (ny/2) * MOD(Blk(n)%J,2)
        H2(2:0:-1,:,:) = Blk(N6)%h(I:I+2,J:J+2,2:9)
        G2(2:0:-1,:,:) = Blk(N6)%g(I:I+2,J:J+2,2:9)
        CALL Biquadratic_from_Coarse_f( 3,3, 1,1, H2, bh, 8 )
        CALL Biquadratic_from_Coarse_f( 3,3, 1,1, G2, bg, 8 )

    ELSE
        PRINT '(A,4I5)', " Error in Fill_Left_Top_hg: ", n, N6, L, Blk(N6)%Level
        RETURN
    END IF

    Blk(n)%h(0,ny+1,2:9) = bh(1,1,:)
    Blk(n)%g(0,ny+1,2:9) = bg(1,1,:)

END SUBROUTINE Fill_Left_Top_hg

!******************************************************************************

SUBROUTINE Fill_Right_Bottom_hg( n, L, N7 )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk
    INTEGER, INTENT(IN) :: n, L, N7

    INTEGER :: child, I, J
    REAL(wp) :: h(3,3,2:9), H2(0:2,0:2,2:9), bh(1,1,2:9)
    REAL(wp) :: g(3,3,2:9), G2(0:2,0:2,2:9), bg(1,1,2:9)

    IF( L == Blk(N7)%Level ) THEN

        IF( Blk(N7)%IsLeaf ) THEN
            Blk(n)%h(nx+1,0,2:9) = Blk(N7)%h(1,ny,2:9)
            Blk(n)%g(nx+1,0,2:9) = Blk(N7)%g(1,ny,2:9)
            RETURN
        END IF
        !== the neighbour is refined: restrict its children's strip
        child = Blk(N7)%child(1)
        h(:,3:1:-1,:) = Blk(child+2)%h(1:3,ny-2:ny,2:9)
        g(:,3:1:-1,:) = Blk(child+2)%g(1:3,ny-2:ny,2:9)
        CALL Biquadratic_from_Fine_f( 3,3, 1,1, h, bh, 8 )
        CALL Biquadratic_from_Fine_f( 3,3, 1,1, g, bg, 8 )

    ELSEIF( L == Blk(N7)%Level + 1 ) THEN

        !== the neighbour is coarser: interpolate its strip
        I = nx/2 - (nx/2) * MOD(Blk(n)%I,2)
        J = ny   - (ny/2) * MOD(Blk(n)%J,2) - 1
        H2(:,2:0:-1,:) = Blk(N7)%h(I:I+2,J:J+2,2:9)
        G2(:,2:0:-1,:) = Blk(N7)%g(I:I+2,J:J+2,2:9)
        CALL Biquadratic_from_Coarse_f( 3,3, 1,1, H2, bh, 8 )
        CALL Biquadratic_from_Coarse_f( 3,3, 1,1, G2, bg, 8 )

    ELSE
        PRINT '(A,4I5)', " Error in Fill_Right_Bottom_hg: ", n, N7, L, Blk(N7)%Level
        RETURN
    END IF

    Blk(n)%h(nx+1,0,2:9) = bh(1,1,:)
    Blk(n)%g(nx+1,0,2:9) = bg(1,1,:)

END SUBROUTINE Fill_Right_Bottom_hg

!******************************************************************************

SUBROUTINE Fill_Right_Top_hg( n, L, N8 )
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk
    INTEGER, INTENT(IN) :: n, L, N8

    INTEGER :: child, I, J
    REAL(wp) :: h(3,3,2:9), H2(0:2,0:2,2:9), bh(1,1,2:9)
    REAL(wp) :: g(3,3,2:9), G2(0:2,0:2,2:9), bg(1,1,2:9)

    IF( L == Blk(N8)%Level ) THEN

        IF( Blk(N8)%IsLeaf ) THEN
            Blk(n)%h(nx+1,ny+1,2:9) = Blk(N8)%h(1,1,2:9)
            Blk(n)%g(nx+1,ny+1,2:9) = Blk(N8)%g(1,1,2:9)
            RETURN
        END IF
        !== the neighbour is refined: restrict its children's strip
        child = Blk(N8)%child(1)
        h(:,:,:) = Blk(child)%h(1:3,1:3,2:9)
        g(:,:,:) = Blk(child)%g(1:3,1:3,2:9)
        CALL Biquadratic_from_Fine_f( 3,3, 1,1, h, bh, 8 )
        CALL Biquadratic_from_Fine_f( 3,3, 1,1, g, bg, 8 )

    ELSEIF( L == Blk(N8)%Level + 1 ) THEN

        !== the neighbour is coarser: interpolate its strip
        I = nx/2 - (nx/2) * MOD(Blk(n)%I,2)
        J = ny/2 - (ny/2) * MOD(Blk(n)%J,2)
        H2(:,:,:) = Blk(N8)%h(I:I+2,J:J+2,2:9)
        G2(:,:,:) = Blk(N8)%g(I:I+2,J:J+2,2:9)
        CALL Biquadratic_from_Coarse_f( 3,3, 1,1, H2, bh, 8 )
        CALL Biquadratic_from_Coarse_f( 3,3, 1,1, G2, bg, 8 )

    ELSE
        PRINT '(A,4I5)', " Error in Fill_Right_Top_hg: ", n, N8, L, Blk(N8)%Level
        RETURN
    END IF

    Blk(n)%h(nx+1,ny+1,2:9) = bh(1,1,:)
    Blk(n)%g(nx+1,ny+1,2:9) = bg(1,1,:)

END SUBROUTINE Fill_Right_Top_hg

END MODULE Ghost_hg_Mod
