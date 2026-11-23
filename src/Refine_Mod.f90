!================================================
!   Flagging, refinement and derefinement of blocks
!================================================
MODULE Refine_Mod
    USE Kinds_Mod, ONLY: wp
    USE Blocks_Mod, ONLY: Block_t, Blk, Nblks, LevelMin, LevelMax, facing
    USE Mesh_Mod, ONLY: Add_Child_Blocks, Array_of_Leaf_Blocks
    USE Prolong_Restrict_Mod, ONLY: Prolong_Bilinear, Restrict_Bilinear
    USE Model_Mod, ONLY: LBM_Model, Block_Error
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Refine_Derefine_Mesh

CONTAINS


!== adapt the mesh with the model's refinement criterion and tolerances; its
!   Moments_g completes the pressure and velocity of new blocks. changed
!   tells whether any block was refined or derefined.
SUBROUTINE Refine_Derefine_Mesh( model, changed )
    IMPLICIT NONE
    TYPE(LBM_Model), INTENT(IN) :: model
    LOGICAL, INTENT(OUT), OPTIONAL :: changed

    LOGICAL :: mesh_changed

    CALL Flag_for_refinement  ( model%Refinement_Error, model%refine_tol, model%max_island_blocks )
    CALL Flag_for_derefinement( model%Refinement_Error, model%deref_tol, model%max_island_blocks )

    mesh_changed = ANY( Blk(1:Nblks)%Refin ) .OR. ANY( Blk(1:Nblks)%Deref )
    IF( PRESENT(changed) ) changed = mesh_changed
    IF( .NOT. mesh_changed ) RETURN

    CALL Derefine
    CALL Refine( model )

    !== leaf lists and neighbour table of the new mesh
    CALL Array_of_Leaf_Blocks

END SUBROUTINE Refine_Derefine_Mesh

!******************************************************************************

!== leaves below the finest level whose error exceeds refinTol, plus the
!   coarser neighbours needed to keep the 2:1 balance. With max_island > 0,
!   also small groups of leaves enclosed by finer blocks (e.g. the interior of
!   a large droplet), so that no coarse/fine edges are left inside a phase.
SUBROUTINE Flag_for_refinement( criterion, refinTol, max_island )
    USE Config_Mod, ONLY: maxLevels
    IMPLICIT NONE
    PROCEDURE(Block_Error) :: criterion
    REAL(wp), INTENT(IN) :: refinTol
    INTEGER,  INTENT(IN) :: max_island

    INTEGER :: n

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED( Blk, Nblks, refinTol )
    DO n = 1, Nblks
        Blk(n)%Refin = Blk(n)%IsLeaf .AND. Blk(n)%Level < maxLevels
        IF( Blk(n)%Refin ) Blk(n)%Refin = criterion(n) > refinTol
    END DO
!$OMP END PARALLEL DO

    IF( max_island > 0 ) CALL Refine_Coarse_Islands( max_island, after_deref = .FALSE. )

    DO n = 1, Nblks
        IF( Blk(n)%Refin ) CALL Enforce_1_Level_Jump( n )
    END DO

END SUBROUTINE Flag_for_refinement

!******************************************************************************

!== coarse islands: connected groups (through the four sides) of at most
!   max_island leaves on one level below the finest that touch only finer
!   blocks, e.g. the
!   interior of a large droplet whose gradient-based criterion is zero (a group
!   that reaches a wall or a coarser block is not an island). Their
!   leaves are flagged for refinement. after_deref: the mesh after the planned
!   derefinement, where parents flagged Deref are leaves; an island made of
!   them keeps its children instead.
SUBROUTINE Refine_Coarse_Islands( max_island, after_deref )
    USE Config_Mod, ONLY: maxLevels
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: max_island
    LOGICAL, INTENT(IN) :: after_deref

    INTEGER, ALLOCATABLE :: comp(:), stack(:), members(:)
    INTEGER :: n, m, b, k, top, nm, ic
    LOGICAL :: enclosed

    ALLOCATE( comp(Nblks), stack(Nblks), members(Nblks) )
    comp = 0
    ic   = 0
    DO n = 1, Nblks
        IF( comp(n) /= 0 .OR. .NOT. Coarse_Leaf( n ) ) CYCLE
        ic = ic + 1
        comp(n) = ic
        top = 1;  stack(1) = n
        nm  = 0
        enclosed = .TRUE.
        DO WHILE( top > 0 )
            m = stack(top);  top = top - 1
            nm = nm + 1;  members(nm) = m
            DO k = 1, 4
                b = Blk(m)%nbr(k)
                IF( b <= 0 ) THEN
                    enclosed = .FALSE.                                  ! a wall: the far field, not an island
                ELSE IF( Blk(b)%Level < Blk(m)%Level ) THEN
                    enclosed = .FALSE.                                  ! a coarser neighbour
                ELSE IF( Leaf_After( b ) ) THEN
                    IF( comp(b) == 0 ) THEN                             ! same-level leaf: same island
                        comp(b) = ic
                        top = top + 1;  stack(top) = b
                    END IF
                END IF                                                  ! else finer
            END DO
        END DO
        IF( .NOT. enclosed .OR. nm > max_island ) CYCLE
        DO k = 1, nm
            m = members(k)
            IF( Blk(m)%Deref ) THEN
                Blk(m)%Deref = .FALSE.
            ELSE
                Blk(m)%Refin = .TRUE.
            END IF
        END DO
    END DO

CONTAINS

    !== n is a leaf (after the planned derefinement, if after_deref) that is not flagged for refinement
    LOGICAL FUNCTION Leaf_After( n )
        INTEGER, INTENT(IN) :: n
        Leaf_After = ( Blk(n)%IsLeaf .AND. .NOT. Blk(n)%Refin ) .OR. ( after_deref .AND. Blk(n)%Deref )
    END FUNCTION Leaf_After

    LOGICAL FUNCTION Coarse_Leaf( n )
        INTEGER, INTENT(IN) :: n
        Coarse_Leaf = Blk(n)%Level >= 0 .AND. Blk(n)%Level < maxLevels .AND. Leaf_After( n )
    END FUNCTION Coarse_Leaf

END SUBROUTINE Refine_Coarse_Islands

!******************************************************************************

!== block n will be refined: so must every coarser neighbour (a coarser
!   neighbour is a leaf, else its child at n's level would be the neighbour)
RECURSIVE SUBROUTINE Enforce_1_Level_Jump( n )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    INTEGER :: k, b

    DO k = 1, 8
        b = Blk(n)%nbr(k)
        IF( b <= 0 ) CYCLE
        IF( .NOT. Blk(b)%Refin .AND. Blk(b)%Level == Blk(n)%Level-1 ) THEN
            Blk(b)%Refin = .TRUE.
            CALL Enforce_1_Level_Jump( b )
        END IF
    END DO

END SUBROUTINE Enforce_1_Level_Jump

!******************************************************************************

!== parents of four leaves, none flagged for refinement, all below derefTol;
!   kept where removing the children would break the 2:1 balance
SUBROUTINE Flag_for_derefinement( criterion, derefTol, max_island )
    USE Config_Mod, ONLY: minLevels
    IMPLICIT NONE
    PROCEDURE(Block_Error) :: criterion
    REAL(wp), INTENT(IN) :: derefTol
    INTEGER,  INTENT(IN) :: max_island

    INTEGER :: n, L, c, q

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED( Blk, Nblks, derefTol )
    DO n = 1, Nblks
        Blk(n)%Deref = .FALSE.
        c = Blk(n)%child
        IF( c < 0 .OR. Blk(n)%Level < minLevels ) CYCLE      ! a leaf, or an unused slot
        IF( .NOT. ALL(Blk(c:c+3)%IsLeaf) .OR. ANY(Blk(c:c+3)%Refin) ) CYCLE

        Blk(n)%Deref = .TRUE.
        DO q = c, c+3
            IF( criterion(q) > derefTol ) THEN
                Blk(n)%Deref = .FALSE.
                EXIT
            END IF
        END DO
    END DO
!$OMP END PARALLEL DO

    !== parents that would form a small island of leaves enclosed by finer
    !   blocks keep their children (they would be refined again at the next step)
    IF( max_island > 0 ) CALL Refine_Coarse_Islands( max_island, after_deref = .TRUE. )

    !== finest first: a block kept refined may keep its coarser neighbours refined
    DO L = LevelMax-1, LevelMin-1, -1
        DO n = 1, Nblks
            IF( Blk(n)%Level == L .AND. Blk(n)%Deref ) CALL Keep_Level_Jump( n )
        END DO
    END DO

CONTAINS

    !== n stays refined if a neighbour's child touching it will be refined
    !   or stays refined
    SUBROUTINE Keep_Level_Jump( n )
        INTEGER, INTENT(IN) :: n

        INTEGER :: k, b, i, c

        DO k = 1, 8
            b = Blk(n)%nbr(k)
            IF( b <= 0 ) CYCLE
            IF( Blk(b)%IsLeaf ) CYCLE
            DO i = 1, 2
                c = Blk(b)%child + facing(i,k)
                IF( Blk(c)%Refin .OR. (.NOT.Blk(c)%IsLeaf .AND. .NOT.Blk(c)%Deref) ) THEN
                    Blk(n)%Deref = .FALSE.
                    RETURN
                END IF
            END DO
        END DO

    END SUBROUTINE Keep_Level_Jump

END SUBROUTINE Flag_for_derefinement

!******************************************************************************

!== restrict the flagged groups of siblings to their parent and queue their slots
SUBROUTINE Derefine
    USE Blocks_Mod, ONLY: Release_Siblings
    IMPLICIT NONE

    INTEGER :: n, c

    DO n = 1, Nblks
        IF( .NOT. Blk(n)%Deref ) CYCLE
        CALL Restrict_Bilinear( n )
        c = Blk(n)%child
        Blk(c:c+3)    = Block_t()       ! unused slots
        Blk(n)%child  = -1
        Blk(n)%IsLeaf = .TRUE.
        CALL Release_Siblings( c )
    END DO

END SUBROUTINE Derefine

!******************************************************************************

!== children of the flagged leaves, coarse levels first, in free slots
SUBROUTINE Refine( model )
    USE Blocks_Mod, ONLY: Take_Free_Siblings
    IMPLICIT NONE
    TYPE(LBM_Model), INTENT(IN) :: model

    INTEGER :: n, L, c1, c

    DO L = LevelMin, LevelMax
        DO n = 1, Nblks
            IF( Blk(n)%Level /= L .OR. .NOT. Blk(n)%Refin ) CYCLE

            CALL Take_Free_Siblings( c1 )
            IF( c1 < 0 ) c1 = Nblks+1
            CALL Add_Child_Blocks( n, c1 )

            CALL Prolong_Bilinear( n, c1 )
            DO c = c1, c1+3
                CALL model%Moments_g( c )
            END DO
        END DO
    END DO

END SUBROUTINE Refine

END MODULE Refine_Mod
