!================================================
!   Block storage of the quadtree mesh
!
!   Metadata and field data are stored apart. Blk(1:Nblks) holds the
!   metadata of every block, leaf or not, as plain integers and logicals
!   (no pointers). The fields of block n are the last-index slices of
!   one pool array per field, e.g. Blk_h(:,:,:,n) or Blk_C(:,:,n), so a
!   block index is all a kernel needs.
!
!   The pool is allocated once with a capacity of MaxBlks slots, the
!   number of blocks of the fully refined quadtree, which bounds the
!   number of live blocks. Only slots that get used are touched, and
!   the OS commits a page of the pool when it is first written. A smaller
!   start capacity can be set with the environment variable
!   AMR_BLOCK_CAPACITY. Reserve_Blocks then grows the pool (and the
!   metadata) whenever a refinement would overflow it. The block indices
!   stay valid because everything refers to blocks by index.
!
!   Slots of derefined blocks are kept in a FIFO queue of sibling groups
!   (four consecutive slots, nUnused slots in total) and reused by the
!   next refinement. The leaves are listed in LeafBlks and, sorted by
!   level, in LeafByLevel.
!================================================
MODULE Blocks_Mod
    USE Kinds_Mod, ONLY: wp
    USE Config_Mod, ONLY: nx, ny, maxLevels, MaxBlks
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Block_t, Blk, Blk_h, Blk_g, Blk_C, Blk_P, Blk_Ux, Blk_Uy, &
              Nblks, nLeaf, LeafBlks, LeafByLevel, LevelStart, LevelMin, LevelMax, &
              nUnused, &
              Init_Block_Pool, Reserve_Blocks, Release_Siblings, Take_Free_Siblings

    !== metadata of a block
    !   nbr(k) is the neighbour in direction k (numbering below): the block
    !   at the neighbouring position on the same level, else the coarser
    !   block covering it, or the boundary code -k beyond the bottom and top
    !   walls. It is rebuilt from the block positions (Level, I, J) after
    !   every mesh change (Rebuild_Neighbor_Table in Mesh_Mod).
    !       6   4   8
    !       1       2
    !       5   3   7
    TYPE :: Block_t
        INTEGER :: Level  = -1
        INTEGER :: dx     =  0      ! cell size in finest-level lattice units
        INTEGER :: I = 0, J = 0     ! block position on its level
        INTEGER :: Parent = -1
        INTEGER :: child(4)     = -1    ! first child: child(1), siblings follow; -1 for a leaf
        INTEGER :: nbr(8)       =  0
        LOGICAL :: IsLeaf  = .FALSE.
        LOGICAL :: hasLeaf = .FALSE.
        LOGICAL :: Refin   = .FALSE.
        LOGICAL :: Deref   = .FALSE.
    END TYPE Block_t

    TYPE(Block_t), ALLOCATABLE :: Blk(:)

    !== block fields, last index = block
    REAL(wp), ALLOCATABLE :: Blk_h(:,:,:,:), Blk_g(:,:,:,:)     ! (0:nx+1,0:ny+1,9,:)
    REAL(wp), ALLOCATABLE :: Blk_C(:,:,:)                       ! (-1:nx+2,-1:ny+2,:)
    REAL(wp), ALLOCATABLE :: Blk_P(:,:,:), Blk_Ux(:,:,:), Blk_Uy(:,:,:)   ! (nx,ny,:)

    INTEGER :: Nblks = 0, nLeaf = 0

    !== coarsest and finest level that currently has leaves
    INTEGER :: LevelMin, LevelMax

    INTEGER, ALLOCATABLE :: LeafBlks(:)

    !== leaf blocks sorted by refinement level, in Z-order within a level
    !   (rebuilt in Array_of_Leaf_Blocks): the leaves on level L are
    !   LeafByLevel( LevelStart(L) : LevelStart(L+1)-1 )
    INTEGER, ALLOCATABLE :: LeafByLevel(:)
    INTEGER :: LevelStart(0:maxLevels+1)

    !== FIFO ring of the first slots of free sibling groups
    INTEGER :: free_first((MaxBlks+1)/4)
    INTEGER :: free_head = 1, free_len = 0
    INTEGER :: nUnused = 0

CONTAINS

!******************************************************************************

!== allocate the metadata and the field pool
SUBROUTINE Init_Block_Pool
    IMPLICIT NONE

    CHARACTER(len=32) :: str
    INTEGER :: length, stat, capacity

    capacity = MaxBlks
    CALL GET_ENVIRONMENT_VARIABLE( 'AMR_BLOCK_CAPACITY', str, length, stat )
    IF( stat == 0 .AND. length > 0 ) THEN
        READ( str, *, IOSTAT=stat ) capacity
        IF( stat /= 0 ) ERROR STOP 'AMR_BLOCK_CAPACITY: not an integer'
        capacity = MIN( MAX( capacity, 3 ), MaxBlks )
    END IF

    CALL Allocate_Pool( capacity )
    Nblks     = 0
    free_head = 1
    free_len  = 0
    nUnused   = 0

    PRINT '(/A,I0,A,I0,A,I0,A)', 'Block pool: ', capacity, ' of at most ', MaxBlks, &
        ' blocks (', CEILING( capacity * Slot_Bytes() / 2._wp**20 ), ' MiB reserved)'

END SUBROUTINE Init_Block_Pool

!******************************************************************************

!== make room for blocks 1:nmax; grows the pool by doubling, keeping 1:Nblks
SUBROUTINE Reserve_Blocks( nmax )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nmax

    INTEGER :: old

    old = SIZE( Blk )
    IF( nmax <= old ) RETURN
    IF( nmax > MaxBlks ) THEN
        PRINT '(/A,2I9)', ' blocks needed, MaxBlks:', nmax, MaxBlks
        ERROR STOP 'more blocks than in the fully refined quadtree'
    END IF

    CALL Allocate_Pool( MIN( MAX( nmax, 2*old ), MaxBlks ) )
    PRINT '(A,I0,A,I0,A)', ' block pool grown from ', old, ' to ', SIZE(Blk), ' blocks'

END SUBROUTINE Reserve_Blocks

!******************************************************************************

!== (re)allocate everything with the given capacity, keeping blocks 1:Nblks
SUBROUTINE Allocate_Pool( capacity )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: capacity

    TYPE(Block_t), ALLOCATABLE :: meta(:)
    REAL(wp), ALLOCATABLE :: f4(:,:,:,:), f3(:,:,:)

    ALLOCATE( meta(capacity) )
    IF( ALLOCATED(Blk) ) meta(1:Nblks) = Blk(1:Nblks)
    CALL MOVE_ALLOC( meta, Blk )

    !== no initialization: a slot is zeroed when it is assigned to a block
    ALLOCATE( f4(0:nx+1,0:ny+1,9,capacity) )
    IF( ALLOCATED(Blk_h) ) f4(:,:,:,1:Nblks) = Blk_h(:,:,:,1:Nblks)
    CALL MOVE_ALLOC( f4, Blk_h )

    ALLOCATE( f4(0:nx+1,0:ny+1,9,capacity) )
    IF( ALLOCATED(Blk_g) ) f4(:,:,:,1:Nblks) = Blk_g(:,:,:,1:Nblks)
    CALL MOVE_ALLOC( f4, Blk_g )

    ALLOCATE( f3(-1:nx+2,-1:ny+2,capacity) )
    IF( ALLOCATED(Blk_C) ) f3(:,:,1:Nblks) = Blk_C(:,:,1:Nblks)
    CALL MOVE_ALLOC( f3, Blk_C )

    ALLOCATE( f3(nx,ny,capacity) )
    IF( ALLOCATED(Blk_P) ) f3(:,:,1:Nblks) = Blk_P(:,:,1:Nblks)
    CALL MOVE_ALLOC( f3, Blk_P )

    ALLOCATE( f3(nx,ny,capacity) )
    IF( ALLOCATED(Blk_Ux) ) f3(:,:,1:Nblks) = Blk_Ux(:,:,1:Nblks)
    CALL MOVE_ALLOC( f3, Blk_Ux )

    ALLOCATE( f3(nx,ny,capacity) )
    IF( ALLOCATED(Blk_Uy) ) f3(:,:,1:Nblks) = Blk_Uy(:,:,1:Nblks)
    CALL MOVE_ALLOC( f3, Blk_Uy )

END SUBROUTINE Allocate_Pool

!******************************************************************************

!== bytes of field data per block slot
PURE FUNCTION Slot_Bytes() RESULT( bytes )
    IMPLICIT NONE
    REAL(wp) :: bytes

    bytes = REAL( STORAGE_SIZE(1._wp)/8, wp ) * &
            ( 2*9*(nx+2)*(ny+2) + (nx+4)*(ny+4) + 3*nx*ny )

END FUNCTION Slot_Bytes

!******************************************************************************

!== queue the slots c1:c1+3 of removed siblings for reuse
SUBROUTINE Release_Siblings( c1 )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: c1

    IF( free_len == SIZE(free_first) ) ERROR STOP 'free-slot queue overflow'
    free_first( 1 + MOD( free_head-1 + free_len, SIZE(free_first) ) ) = c1
    free_len = free_len + 1
    nUnused  = nUnused + 4

END SUBROUTINE Release_Siblings

!******************************************************************************

!== the oldest free sibling group c1:c1+3, or c1 = -1 if there is none
SUBROUTINE Take_Free_Siblings( c1 )
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: c1

    c1 = -1
    IF( free_len == 0 ) RETURN

    c1 = free_first( free_head )
    free_head = 1 + MOD( free_head, SIZE(free_first) )
    free_len  = free_len - 1
    nUnused   = nUnused - 4

END SUBROUTINE Take_Free_Siblings

END MODULE Blocks_Mod
