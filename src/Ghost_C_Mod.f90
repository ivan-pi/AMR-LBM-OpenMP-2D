!================================================
!   Ghost cells of the phase field C (ng = 1 or 2 layers)
!
!   As for h and g (Ghost_hg_Mod, where the strip layout is described): a
!   copy from a same-level leaf, the restriction of the facing children of a
!   refined neighbour, or the interpolation of a coarser one, with one
!   routine per strip shape. The strips across an edge are 2+ng fine or 3
!   coarse cells deep.
!================================================
MODULE Ghost_C_Mod
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: nx, ny
    USE Blocks_Mod, ONLY: Blk, Blk_C, di, dj, facing
    USE Interpolations_Mod, ONLY: Biquadratic_from_Fine, Biquadratic_from_Coarse
    USE BC_Mod, ONLY: Symmetric_C_3, Symmetric_C_4
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Fill_Ghost_Cells_C

CONTAINS


!== fill ng ghost layers of C on all leaf blocks, level by level (coarse
!   first: an interpolation from a coarser neighbour reads its ghost cells)
SUBROUTINE Fill_Ghost_Cells_C( ng )
    USE Blocks_Mod, ONLY: LevelMin, LevelMax, LeafByLevel, LevelStart
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ng

    INTEGER :: L, k, m, n

    DO L = LevelMin, LevelMax
!$OMP DO schedule(static)
        DO m = LevelStart(L), LevelStart(L+1)-1
            n = LeafByLevel(m)

            DO k = 1, 2
                IF( Blk(n)%nbr(k) > 0 ) CALL Fill_Side_X( n, Blk(n)%nbr(k), di(k), ng )
            END DO
            DO k = 3, 4
                IF( Blk(n)%nbr(k) > 0 ) CALL Fill_Side_Y( n, Blk(n)%nbr(k), dj(k), ng )
            END DO
            DO k = 5, 8
                IF( Blk(n)%nbr(k) > 0 ) CALL Fill_Corner( n, k, ng )
            END DO

            IF( Blk(n)%nbr(3) < 0 ) CALL Symmetric_C_3( ng, Blk_C(:,:,n) )
            IF( Blk(n)%nbr(4) < 0 ) CALL Symmetric_C_4( ng, Blk_C(:,:,n) )

        END DO
!$OMP END DO
    END DO

END SUBROUTINE Fill_Ghost_Cells_C

!******************************************************************************

!== the ng ghost columns of block n towards its left (s = -1) or right (s = +1)
!   neighbour nb. Copies are element loops: as sections of the same pool
!   array they would go through a temporary (no proof that they don't overlap).
SUBROUTINE Fill_Side_X( n, nb, s, ng )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nb, s, ng

    INTEGER  :: ig, i0, c, l, j, Y
    REAL(wp) :: Cf(2+ng,2*ny), Cc(3,ny/2+2), Cg(ng,ny)

    ig = MERGE( 0, nx+1, s < 0 )

    IF( Blk(nb)%Level < Blk(n)%Level ) THEN
        i0 = MERGE( nx+1, 0, s < 0 )
        Y  = (ny/2) * MOD(Blk(n)%J,2)
        DO l = 0, 2
            Cc(l+1,:) = Blk_C(i0+s*l,Y:Y+ny/2+1, nb)
        END DO
        CALL Biquadratic_from_Coarse( 3,ny/2+2, ng,ny, Cc, Cg )
    ELSE
        i0 = MERGE( nx, 1, s < 0 )
        IF( Blk(nb)%IsLeaf ) THEN
            DO j = 1, ny
            DO l = 0, ng-1
                Blk_C(ig+s*l,j, n) = Blk_C(i0+s*l,j, nb)
            END DO
            END DO
            RETURN
        END IF
        c = Blk(nb)%child + MERGE( 1, 0, s < 0 )      ! the facing children c, c+2
        DO l = 0, 1+ng
            Cf(l+1,   1:  ny) = Blk_C(i0+s*l,1:ny, c  )
            Cf(l+1,ny+1:2*ny) = Blk_C(i0+s*l,1:ny, c+2)
        END DO
        CALL Biquadratic_from_Fine( 2+ng,2*ny, ng,ny, Cf, Cg )
    END IF

    DO j = 1, ny
    DO l = 0, ng-1
        Blk_C(ig+s*l,j, n) = Cg(l+1,j)
    END DO
    END DO

END SUBROUTINE Fill_Side_X

!******************************************************************************

!== the ng ghost rows of block n towards its bottom (s = -1) or top (s = +1)
!   neighbour nb
SUBROUTINE Fill_Side_Y( n, nb, s, ng )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nb, s, ng

    INTEGER  :: jg, j0, c, l, i, X
    REAL(wp) :: Cf(2*nx,2+ng), Cc(nx/2+2,3), Cg(nx,ng)

    jg = MERGE( 0, ny+1, s < 0 )

    IF( Blk(nb)%Level < Blk(n)%Level ) THEN
        j0 = MERGE( ny+1, 0, s < 0 )
        X  = (nx/2) * MOD(Blk(n)%I,2)
        DO l = 0, 2
            Cc(:,l+1) = Blk_C(X:X+nx/2+1,j0+s*l, nb)
        END DO
        CALL Biquadratic_from_Coarse( nx/2+2,3, nx,ng, Cc, Cg )
    ELSE
        j0 = MERGE( ny, 1, s < 0 )
        IF( Blk(nb)%IsLeaf ) THEN
            DO l = 0, ng-1
            DO i = 1, nx
                Blk_C(i,jg+s*l, n) = Blk_C(i,j0+s*l, nb)
            END DO
            END DO
            RETURN
        END IF
        c = Blk(nb)%child + MERGE( 2, 0, s < 0 )      ! the facing children c, c+1
        DO l = 0, 1+ng
            Cf(   1:  nx,l+1) = Blk_C(1:nx,j0+s*l, c  )
            Cf(nx+1:2*nx,l+1) = Blk_C(1:nx,j0+s*l, c+1)
        END DO
        CALL Biquadratic_from_Fine( 2*nx,2+ng, nx,ng, Cf, Cg )
    END IF

    DO l = 0, ng-1
    DO i = 1, nx
        Blk_C(i,jg+s*l, n) = Cg(i,l+1)
    END DO
    END DO

END SUBROUTINE Fill_Side_Y

!******************************************************************************

!== the ng x ng ghost corner of block n towards its neighbour in direction k = 5..8
SUBROUTINE Fill_Corner( n, k, ng )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, k, ng

    INTEGER  :: nb, sx, sy, ig, jg, i0, j0, px, py, c, l, m
    REAL(wp) :: Cf(2+ng,2+ng), Cc(3,3), Cg(ng,ng)

    nb = Blk(n)%nbr(k)
    sx = di(k)
    sy = dj(k)
    ig = MERGE( 0, nx+1, sx < 0 )
    jg = MERGE( 0, ny+1, sy < 0 )

    IF( Blk(nb)%Level < Blk(n)%Level ) THEN
        px = MOD(Blk(n)%I,2)
        py = MOD(Blk(n)%J,2)
        i0 = MERGE( nx+1-(nx/2)*px, (nx/2)*(1-px), sx < 0 )
        j0 = MERGE( ny+1-(ny/2)*py, (ny/2)*(1-py), sy < 0 )
        DO m = 0, 2
        DO l = 0, 2
            Cc(l+1,m+1) = Blk_C(i0+sx*l,j0+sy*m, nb)
        END DO
        END DO
        CALL Biquadratic_from_Coarse( 3,3, ng,ng, Cc, Cg )
    ELSE
        i0 = MERGE( nx, 1, sx < 0 )
        j0 = MERGE( ny, 1, sy < 0 )
        c  = nb
        IF( .NOT. Blk(nb)%IsLeaf ) c = Blk(nb)%child + facing(1,k)     ! the facing child
        DO m = 0, MERGE( ng-1, 1+ng, c == nb )
        DO l = 0, MERGE( ng-1, 1+ng, c == nb )
            Cf(l+1,m+1) = Blk_C(i0+sx*l,j0+sy*m, c)
        END DO
        END DO
        IF( c == nb ) THEN
            Cg = Cf(1:ng,1:ng)
        ELSE
            CALL Biquadratic_from_Fine( 2+ng,2+ng, ng,ng, Cf, Cg )
        END IF
    END IF

    DO m = 0, ng-1
    DO l = 0, ng-1
        Blk_C(ig+sx*l,jg+sy*m, n) = Cg(l+1,m+1)
    END DO
    END DO

END SUBROUTINE Fill_Corner

END MODULE Ghost_C_Mod
