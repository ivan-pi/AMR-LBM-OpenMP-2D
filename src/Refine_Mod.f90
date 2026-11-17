!================================================
!   Flagging, refinement and derefinement of blocks
!================================================
MODULE Refine_Mod
    USE Kinds_Mod, ONLY: wp
    USE Mesh_Mod, ONLY: Add_Child_Blocks, Array_of_Leaf_Blocks
    USE Prolong_Restrict_Mod, ONLY: Prolong_Bilinear, Restrict_Bilinear
    USE Model_Mod, ONLY: Block_Kernel, Block_Error
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Refine_Derefine_Mesh

CONTAINS


!== adapt the mesh: criterion decides where to refine / derefine; Moments_g
!   (the model's kernel) completes the pressure and velocity of new blocks;
!   changed tells whether any block was refined or derefined
SUBROUTINE Refine_Derefine_Mesh( criterion, Moments_g, refinTol, derefTol, changed )
    USE Blocks_Mod, ONLY: Blk, Nblks
    IMPLICIT NONE
    PROCEDURE(Block_Error)  :: criterion
    PROCEDURE(Block_Kernel) :: Moments_g
    REAL(wp), INTENT(IN) :: refinTol, derefTol
    LOGICAL, INTENT(OUT), OPTIONAL :: changed

    LOGICAL :: mesh_changed

    CALL Flag_for_refinement  ( criterion, refinTol )
    CALL Flag_for_derefinement( criterion, derefTol )

    mesh_changed = ANY( Blk(1:Nblks)%Refin ) .OR. ANY( Blk(1:Nblks)%Deref )
    IF( PRESENT(changed) ) changed = mesh_changed
    IF( .NOT. mesh_changed ) RETURN

    CALL Derefine
    CALL Refine( Moments_g )

    !== leaf lists and neighbour table of the new mesh
    CALL Array_of_Leaf_Blocks

END SUBROUTINE Refine_Derefine_Mesh

!******************************************************************************

SUBROUTINE Flag_for_refinement( criterion, refinTol )
    USE Config_Mod, ONLY: minLevels, maxLevels
    USE Blocks_Mod, ONLY: Blk, Nblks
    IMPLICIT NONE
    PROCEDURE(Block_Error) :: criterion
    REAL(wp), INTENT(IN) :: refinTol

    INTEGER :: n

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED( Blk, Nblks, refinTol )
    DO n = 1, Nblks
        Blk(n)%Refin = .False.

        IF( Blk(n)%IsLeaf ) THEN
            IF( Blk(n)%Level < minLevels ) THEN
                Blk(n)%Refin = .True.
                CYCLE
            ELSEIF( Blk(n)%Level < maxLevels ) THEN

                IF( criterion(n) > refinTol ) THEN
                    Blk(n)%Refin = .True.
                    CYCLE
                END IF

            END IF
        END IF

    END DO
!$OMP END PARALLEL DO

    DO n = 1, Nblks
        IF( Blk(n)%Refin ) THEN
            CALL Enforce_1_Level_Jump( n )
        END IF
    END DO

END SUBROUTINE Flag_for_refinement

!******************************************************************************

RECURSIVE SUBROUTINE Enforce_1_Level_Jump( n )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n

    INTEGER :: I

    DO I = 1, 8
        CALL Check_Level_Jump( n, I )
    END DO

END

!******************************************************************************

RECURSIVE SUBROUTINE Check_Level_Jump( n, I )
    USE Blocks_Mod, ONLY: Blk
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n, I

    INTEGER :: b

    b = Blk(n)%nbr(I)

    IF( b > 0 ) THEN
        IF( .NOT. Blk(b)%Refin .AND. Blk(b)%Level == Blk(n)%Level-1 ) THEN
            IF( Blk(b)%IsLeaf ) THEN
                Blk(b)%Refin = .True.
                CALL Enforce_1_Level_Jump( b )
            END IF
        END IF
    END IF

END SUBROUTINE Check_Level_Jump

!******************************************************************************

SUBROUTINE Flag_for_derefinement( criterion, derefTol )
    USE Config_Mod, ONLY: minLevels
    USE Blocks_Mod, ONLY: Blk, Nblks, LevelMin, LevelMax
    IMPLICIT NONE
    PROCEDURE(Block_Error) :: criterion
    REAL(wp), INTENT(IN) :: derefTol

    INTEGER :: n, L, c1, c4, c

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED( Blk, Nblks, derefTol )
    DO n = 1, Nblks
        Blk(n)%Deref = .False.

        IF( Blk(n)%hasLeaf ) THEN

            IF( Blk(n)%Level < minLevels ) CYCLE

            c1 = Blk(n)%child(1)
            c4 = c1 + 3

            IF( ANY(Blk(c1:c4)%Refin) ) CYCLE

        !== derefinement based on the children
            Blk(n)%Deref = .True.

            DO c = c1, c4
                IF( criterion(c) > derefTol ) THEN
                    Blk(n)%Deref = .False.
                    EXIT
                END IF
            END DO
        !== derefinement based on the children

        END IF

    END DO
!$OMP END PARALLEL DO

    DO L = LevelMax-1, LevelMin-1, -1
        DO n = 1, Nblks
            IF( Blk(n)%Level == L ) THEN
                IF( Blk(n)%Deref ) THEN
                    CALL enforce_1_Level_jump_deref_2D( n )
                END IF
            END IF
        END DO
    END DO

END SUBROUTINE Flag_for_derefinement

!******************************************************************************

SUBROUTINE enforce_1_Level_jump_deref_2D( n )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n

    CALL Check_Level_Jump_Deref_2D( n,1, 2, 4 ) ! Left
    CALL Check_Level_Jump_Deref_2D( n,2, 1, 3 ) ! Right
    CALL Check_Level_Jump_Deref_2D( n,3, 3, 4 ) ! Bottom
    CALL Check_Level_Jump_Deref_2D( n,4, 1, 2 ) ! Top

    CALL Check_Level_Jump_Deref_2D( n,5, 4, 4 ) ! SW (Left -Bottom)
    CALL Check_Level_Jump_Deref_2D( n,6, 2, 2 ) ! NW (Left -Top)
    CALL Check_Level_Jump_Deref_2D( n,7, 3, 3 ) ! SE (Right-Bottom)
    CALL Check_Level_Jump_Deref_2D( n,8, 1, 1 ) ! NE (Right-Top)

CONTAINS

    SUBROUTINE Check_Level_Jump_Deref_2D( n, Neighbor, c1, c2 )
        USE Blocks_Mod, ONLY: Blk
        INTEGER, INTENT(IN) :: n, Neighbor, c1, c2

        INTEGER :: b, c, I

        IF( .NOT. Blk(n)%Deref ) RETURN

        b = Blk(n)%nbr(Neighbor)

        IF( b > 0 ) THEN
            IF( .NOT. Blk(b)%IsLeaf ) THEN
                c = Blk(b)%child(c1)
                DO I = 1, 2
                    IF( Blk(c)%Refin .OR. (.Not.Blk(c)%IsLeaf .AND. .Not.Blk(c)%Deref) ) THEN
                        Blk(n)%Deref = .False.
                        RETURN
                    END IF
                    c = Blk(b)%child(c2)
                END DO
            END IF
        END IF

    END SUBROUTINE Check_Level_Jump_Deref_2D

END SUBROUTINE enforce_1_Level_jump_deref_2D


!******************************************************************************

SUBROUTINE Derefine
    USE Blocks_Mod, ONLY: Blk, Nblks
    IMPLICIT NONE

    INTEGER :: n

    DO n = 1, Nblks
        IF( Blk(n)%Deref ) THEN
            CALL Restrict_Bilinear( n )
            CALL Remove_Children( n )
        END IF
    END DO

CONTAINS

    SUBROUTINE Remove_Children( par )
        USE Blocks_Mod, ONLY: Blk, Release_Siblings
        INTEGER, INTENT(IN) :: par

        INTEGER :: c1, c4, grandpa, s1, s4

        c1 = Blk(par)%child(1)
        c4 = c1 + 3
        Blk(c1:c4)%Level  = -1
        Blk(c1:c4)%IsLeaf = .False.

        Blk(par)%child(:)= -1
        Blk(par)%IsLeaf  = .True.
        Blk(par)%hasLeaf = .False.

        grandpa = Blk(par)%Parent
        IF( grandpa > 0 ) THEN
            s1 = Blk(grandpa)%child(1)  ! first sibling
            s4 = Blk(grandpa)%child(4)  ! last  sibling
            IF( ALL(Blk(s1:s4)%IsLeaf) ) THEN
                Blk(grandpa)%hasLeaf = .True.
            END IF
        END IF

        CALL Release_Siblings( c1 )

    END SUBROUTINE Remove_Children

END SUBROUTINE Derefine


!******************************************************************************

SUBROUTINE Refine( Moments_g )
    USE Blocks_Mod, ONLY: Blk, Nblks, LevelMin, LevelMax, Take_Free_Siblings
    IMPLICIT NONE
    PROCEDURE(Block_Kernel) :: Moments_g

    INTEGER :: n, L, c1

    DO L = LevelMin, LevelMax
        DO n = 1, Nblks
            IF( Blk(n)%Level == L )THEN
                IF( Blk(n)%Refin ) THEN

                    CALL Take_Free_Siblings( c1 )
                    IF( c1 < 0 ) c1 = Nblks+1
                    CALL Add_Child_Blocks( n, c1 )

                    CALL Prolong_Bilinear( n, c1, Moments_g )

                END IF
            END IF
        END DO
    END DO

END SUBROUTINE Refine

END MODULE Refine_Mod
