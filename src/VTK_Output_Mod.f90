!================================================
!
!   VTK (XML, binary-appended) output of the AMR solution and mesh
!
!   For every output time two VTK UnstructuredGrid files are written
!   into the directory vtk/ :
!
!     Rising_Bubble_2D_NNNN.vtu  all leaf cells as quads, with the cell
!                                data C, P, Velocity, Level and Block
!     AMR_Blocks_NNNN.vtu        one quad per leaf block (the AMR block
!                                structure), with Level, Block, dx and
!                                Parent
!
!   plus the ParaView time-series files Rising_Bubble_2D.pvd and
!   AMR_Blocks.pvd (rewritten after every output, so that they stay valid
!   if the run is interrupted). Coordinates are scaled by L0 as in the
!   TecPlot output.
!
!   If the program is built with the fVTKHDF library (-DAMR_WITH_FVTKHDF=ON,
!   https://github.com/nncarlson/fvtkhdf), the same two data sets can instead
!   be written as time series into two single VTKHDF files,
!   Rising_Bubble_2D.vtkhdf and AMR_Blocks.vtkhdf (mesh changes with time).
!
!   The writer used in a run is selected in Output_Mod (AMR_OUTPUT).
!
!================================================
MODULE VTK_Output_Mod
    USE Kinds_Mod, ONLY: wp
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT8, INT32, INT64, REAL64
#ifdef USE_FVTKHDF
    USE vtkhdf_ug_file_type, ONLY: vtkhdf_ug_file, vtkhdf_cell_data_handle, VTKHDF_DYNAMIC_MESH
#endif
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Write_VTU_Step
#ifdef USE_FVTKHDF
    PUBLIC :: Write_VTKHDF_Step, Close_VTKHDF
#endif

    !== an unstructured grid of quads with cell data (0-based connectivity)
    TYPE Mesh_Data
        REAL(wp),          ALLOCATABLE :: Pts(:,:)      ! (3, nPts)
        INTEGER(INT32),    ALLOCATABLE :: Conn(:)       ! 0-based point ids
        INTEGER(INT32),    ALLOCATABLE :: Offs(:)       ! end offset of each cell in Conn
        INTEGER(INT8),     ALLOCATABLE :: Types(:)
        CHARACTER(len=16), ALLOCATABLE :: RNames(:), INames(:)
        REAL(wp),          ALLOCATABLE :: RData(:,:)    ! (nCells, nR) scalar real fields
        CHARACTER(len=16)              :: VName = ''    ! optional 3-component real field
        REAL(wp),          ALLOCATABLE :: VData(:,:)    ! (3, nCells)
        INTEGER(INT32),    ALLOCATABLE :: IData(:,:)    ! (nCells, nI) scalar integer fields
    END TYPE Mesh_Data

    CHARACTER(len=*), PARAMETER :: dir = 'vtk'

    INTEGER, PARAMETER :: MaxSteps = 10000
    INTEGER  :: nSteps = 0
    REAL(wp) :: Times(MaxSteps)

    INTEGER(INT8), PARAMETER :: VTK_QUAD = 9_INT8

#ifdef USE_FVTKHDF
    !== one time-series file per data set (cells, blocks)
    TYPE VTKHDF_Series
        TYPE(vtkhdf_ug_file) :: file
        LOGICAL :: IsOpen = .False.
        TYPE(vtkhdf_cell_data_handle), ALLOCATABLE :: R(:), I(:)
        TYPE(vtkhdf_cell_data_handle) :: V
    END TYPE VTKHDF_Series
    TYPE(VTKHDF_Series), SAVE :: H5Cells, H5Blocks
#endif

CONTAINS

!******************************************************************************

    !== VTK XML series: the leaf cells and the block structure at the current
    !   time as two .vtu files in vtk/, plus the two .pvd time-series files
    SUBROUTINE Write_VTU_Step( Time )
        REAL(wp), INTENT(IN) :: Time    ! physical time stored with the data

        LOGICAL, SAVE :: IsFirstCall = .True.
        CHARACTER(len=4) :: Id
        TYPE(Mesh_Data)  :: Cells, Blocks

        IF( IsFirstCall ) CALL EXECUTE_COMMAND_LINE( 'mkdir -p '//dir )
        IsFirstCall = .False.

        IF( nSteps == MaxSteps ) THEN
            PRINT '(A)', ' Write_VTU_Step: too many output steps, skipped.'
            RETURN
        END IF
        nSteps = nSteps + 1
        Times(nSteps) = Time
        WRITE( Id, '(I4.4)' ) nSteps-1

        CALL Build_Cells ( Cells  )
        CALL Build_Blocks( Blocks )

        CALL Write_VTU( dir//'/Rising_Bubble_2D_'//Id//'.vtu', Cells  )
        CALL Write_VTU( dir//'/AMR_Blocks_'//Id//'.vtu',       Blocks )

        CALL Write_PVD( 'Rising_Bubble_2D' )
        CALL Write_PVD( 'AMR_Blocks' )

    END SUBROUTINE Write_VTU_Step

!******************************************************************************

    !== all leaf cells (nx*ny quads per leaf block)
    SUBROUTINE Build_Cells( md )
        USE Config_Mod, ONLY: nx, ny, L0
        USE Blocks_Mod, ONLY: Blk, nLeaf, LeafBlks
        TYPE(Mesh_Data), INTENT(OUT) :: md

            INTEGER :: nPts, nCells, m, b, i, j, k, e, p0, dx

            nPts   = nLeaf * (nx+1) * (ny+1)
            nCells = nLeaf *  nx    *  ny

            ALLOCATE( md%Pts(3,nPts), md%RData(nCells,2), md%VData(3,nCells), md%IData(nCells,2) )
            ALLOCATE( md%Conn(4*nCells), md%Offs(nCells), md%Types(nCells) )
            md%RNames = [character(len=16) :: 'C', 'P']
            md%VName  = 'Velocity'
            md%INames = [character(len=16) :: 'Level', 'Block']

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED( Blk, LeafBlks, nLeaf, md )
            DO m = 1, nLeaf
                b  = LeafBlks(m)
                dx = Blk(b)%dx
                p0 = (m-1) * (nx+1) * (ny+1)    ! points before this block (0-based)

                DO j = 1, ny+1
                DO i = 1, nx+1
                    k = p0 + (j-1)*(nx+1) + i
                    md%Pts(1,k) = REAL( Blk(b)%I * dx * nx + dx * (i-1), wp ) / L0
                    md%Pts(2,k) = REAL( Blk(b)%J * dx * ny + dx * (j-1), wp ) / L0
                    md%Pts(3,k) = 0
                END DO
                END DO

                DO j = 1, ny
                DO i = 1, nx
                    e = (m-1)*nx*ny + (j-1)*nx + i
                    k = p0 + (j-1)*(nx+1) + i - 1   ! 0-based lower-left point
                    md%Conn(4*e-3:4*e) = [ k, k+1, k+1+nx+1, k+nx+1 ]
                    md%Offs(e)    = 4*e
                    md%Types(e)   = VTK_QUAD
                    md%RData(e,1) = Blk(b)%C (i,j)  ! C
                    md%RData(e,2) = Blk(b)%P (i,j)  ! P
                    md%VData(:,e) = [ Blk(b)%Ux(i,j), Blk(b)%Uy(i,j), 0.0_wp ]  ! Velocity
                    md%IData(e,1) = Blk(b)%Level    ! Level
                    md%IData(e,2) = b               ! Block
                END DO
                END DO
            END DO
!$OMP END PARALLEL DO

    END SUBROUTINE Build_Cells

!******************************************************************************

    !== the block structure: one quad per leaf block
    SUBROUTINE Build_Blocks( md )
        USE Config_Mod, ONLY: nx, ny, L0
        USE Blocks_Mod, ONLY: Blk, nLeaf, LeafBlks
        TYPE(Mesh_Data), INTENT(OUT) :: md

            INTEGER :: m, b, dx, k

            ALLOCATE( md%Pts(3,4*nLeaf), md%RData(nLeaf,1), md%VData(3,0), md%IData(nLeaf,3) )
            ALLOCATE( md%Conn(4*nLeaf), md%Offs(nLeaf), md%Types(nLeaf) )
            md%RNames = [character(len=16) :: 'dx']
            md%INames = [character(len=16) :: 'Level', 'Block', 'Parent']

            ASSOCIATE( Pts => md%Pts, Conn => md%Conn, Offs => md%Offs, Types => md%Types, &
                       BSize => md%RData, IData => md%IData )

            DO m = 1, nLeaf
                b  = LeafBlks(m)
                dx = Blk(b)%dx
                k  = 4*(m-1)
                Pts(:,k+1) = [ REAL( Blk(b)%I   *dx*nx, wp), REAL( Blk(b)%J   *dx*ny, wp), 0.0_wp ] / L0
                Pts(:,k+2) = [ REAL((Blk(b)%I+1)*dx*nx, wp), REAL( Blk(b)%J   *dx*ny, wp), 0.0_wp ] / L0
                Pts(:,k+3) = [ REAL((Blk(b)%I+1)*dx*nx, wp), REAL((Blk(b)%J+1)*dx*ny, wp), 0.0_wp ] / L0
                Pts(:,k+4) = [ REAL( Blk(b)%I   *dx*nx, wp), REAL((Blk(b)%J+1)*dx*ny, wp), 0.0_wp ] / L0
                Conn(k+1:k+4) = [ k, k+1, k+2, k+3 ]
                Offs(m)    = 4*m
                Types(m)   = VTK_QUAD
                IData(m,:) = [ Blk(b)%Level, b, Blk(b)%Parent ]
                BSize(m,1) = REAL(dx,wp) / L0
            END DO
            END ASSOCIATE

    END SUBROUTINE Build_Blocks

!******************************************************************************
!   Generic writer of an UnstructuredGrid with scalar real cell fields,
!   one optional 3-component real cell field and scalar integer cell fields.
!   All data go to an appended raw-binary section (UInt64 size headers).

    SUBROUTINE Write_VTU( FileName, md )
        CHARACTER(len=*), INTENT(IN) :: FileName
        TYPE(Mesh_Data),  INTENT(IN) :: md

        CHARACTER(len=*), PARAMETER :: LF = ACHAR(10)
        CHARACTER(len=7) :: RType
        INTEGER(INT64)   :: Off, nb
        INTEGER :: u, k, nPts, nCells, rb

        ASSOCIATE( Pts => md%Pts, Conn => md%Conn, Offs => md%Offs, Types => md%Types, &
                   RNames => md%RNames, RData => md%RData, VName => md%VName, &
                   VData => md%VData, INames => md%INames, IData => md%IData )

        nPts   = SIZE(Pts,2)
        nCells = SIZE(Types)
        rb     = STORAGE_SIZE(Pts) / 8
        IF( rb == 4 ) THEN
            RType = 'Float32'
        ELSE
            RType = 'Float64'
        END IF

        OPEN( NEWUNIT=u, FILE=FileName, ACCESS='STREAM', FORM='UNFORMATTED', STATUS='REPLACE' )

        WRITE(u) '<?xml version="1.0"?>'//LF
        WRITE(u) '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="'//Byte_Order()// &
                 '" header_type="UInt64">'//LF
        WRITE(u) '<UnstructuredGrid>'//LF
        WRITE(u) '<Piece NumberOfPoints="'//Str(nPts)//'" NumberOfCells="'//Str(nCells)//'">'//LF

        Off = 0
        WRITE(u) '<CellData>'//LF
        DO k = 1, SIZE(RNames)
            CALL Header( TRIM(RNames(k)), RType, 1, INT(nCells,INT64)*rb )
        END DO
        IF( LEN_TRIM(VName) > 0 ) CALL Header( TRIM(VName), RType, 3, 3_INT64*nCells*rb )
        DO k = 1, SIZE(INames)
            CALL Header( TRIM(INames(k)), 'Int32', 1, 4_INT64*nCells )
        END DO
        WRITE(u) '</CellData>'//LF

        WRITE(u) '<Points>'//LF
        CALL Header( 'Points', RType, 3, 3_INT64*nPts*rb )
        WRITE(u) '</Points>'//LF

        WRITE(u) '<Cells>'//LF
        CALL Header( 'connectivity', 'Int32', 1, 4_INT64*SIZE(Conn) )
        CALL Header( 'offsets',      'Int32', 1, 4_INT64*nCells )
        CALL Header( 'types',        'UInt8', 1, 1_INT64*nCells )
        WRITE(u) '</Cells>'//LF

        WRITE(u) '</Piece>'//LF
        WRITE(u) '</UnstructuredGrid>'//LF
        WRITE(u) '<AppendedData encoding="raw">'//LF//'_'

        !== the binary data, in the same order as the headers above
        DO k = 1, SIZE(RNames)
            nb = INT(nCells,INT64)*rb ; WRITE(u) nb, RData(:,k)
        END DO
        IF( LEN_TRIM(VName) > 0 ) THEN
            nb = 3_INT64*nCells*rb ; WRITE(u) nb, VData
        END IF
        DO k = 1, SIZE(INames)
            nb = 4_INT64*nCells ; WRITE(u) nb, IData(:,k)
        END DO
        nb = 3_INT64*nPts*rb      ; WRITE(u) nb, Pts
        nb = 4_INT64*SIZE(Conn)   ; WRITE(u) nb, Conn
        nb = 4_INT64*nCells       ; WRITE(u) nb, Offs
        nb = 1_INT64*nCells       ; WRITE(u) nb, Types

        WRITE(u) LF//'</AppendedData>'//LF//'</VTKFile>'//LF
        CLOSE(u)

        END ASSOCIATE

    CONTAINS

        SUBROUTINE Header( Name, DType, nComp, nBytes )
            CHARACTER(len=*), INTENT(IN) :: Name, DType
            INTEGER,          INTENT(IN) :: nComp
            INTEGER(INT64),   INTENT(IN) :: nBytes

            WRITE(u) '  <DataArray type="'//TRIM(DType)//'" Name="'//Name// &
                     '" NumberOfComponents="'//Str(nComp)// &
                     '" format="appended" offset="'//Str(Off)//'"/>'//LF
            Off = Off + 8 + nBytes      ! 8 bytes for the UInt64 size header

        END SUBROUTINE Header

        FUNCTION Byte_Order() RESULT( order )
            CHARACTER(len=:), ALLOCATABLE :: order
            INTEGER(INT8) :: bytes(4)

            bytes = TRANSFER( 1_INT32, bytes )
            IF( bytes(1) == 1 ) THEN
                order = 'LittleEndian'
            ELSE
                order = 'BigEndian'
            END IF

        END FUNCTION Byte_Order

        FUNCTION Str( i ) RESULT( s )
            CLASS(*), INTENT(IN) :: i
            CHARACTER(len=:), ALLOCATABLE :: s
            CHARACTER(len=24) :: buf

            SELECT TYPE( i )
            TYPE IS (INTEGER)
                WRITE( buf, '(I0)' ) i
            TYPE IS (INTEGER(INT64))
                WRITE( buf, '(I0)' ) i
            END SELECT
            s = TRIM(buf)

        END FUNCTION Str

    END SUBROUTINE Write_VTU

!******************************************************************************

    SUBROUTINE Write_PVD( Base )
        CHARACTER(len=*), INTENT(IN) :: Base

        INTEGER :: u, k
        CHARACTER(len=4)  :: Id
        CHARACTER(len=24) :: tstr

        OPEN( NEWUNIT=u, FILE=Base//'.pvd', STATUS='REPLACE', ACTION='WRITE' )
        WRITE(u,'(A)') '<?xml version="1.0"?>'
        WRITE(u,'(A)') '<VTKFile type="Collection" version="0.1">'
        WRITE(u,'(A)') '<Collection>'
        DO k = 1, nSteps
            WRITE( Id,   '(I4.4)'    ) k-1
            WRITE( tstr, '(ES16.8E3)') Times(k)
            WRITE(u,'(A)') '  <DataSet timestep="'//TRIM(ADJUSTL(tstr))//'" part="0" file="'// &
                           dir//'/'//Base//'_'//Id//'.vtu"/>'
        END DO
        WRITE(u,'(A)') '</Collection>'
        WRITE(u,'(A)') '</VTKFile>'
        CLOSE(u)

    END SUBROUTINE Write_PVD

#ifdef USE_FVTKHDF
!******************************************************************************
!   VTKHDF output through the fVTKHDF library: an UnstructuredGrid whose
!   mesh changes in time (VTKHDF_DYNAMIC_MESH), one file per data set.

    SUBROUTINE Write_VTKHDF_Step( Time )
        REAL(wp), INTENT(IN) :: Time    ! physical time stored with the data

        TYPE(Mesh_Data) :: Cells, Blocks

        CALL Build_Cells ( Cells  )
        CALL Build_Blocks( Blocks )

        CALL Write_Series( H5Cells,  'Rising_Bubble_2D.vtkhdf', Cells  )
        CALL Write_Series( H5Blocks, 'AMR_Blocks.vtkhdf',       Blocks )

    CONTAINS

        SUBROUTINE Write_Series( S, FileName, md )
            TYPE(VTKHDF_Series), INTENT(INOUT) :: S
            CHARACTER(len=*),    INTENT(IN)    :: FileName
            TYPE(Mesh_Data),     INTENT(IN)    :: md

            INTEGER :: stat, k
            CHARACTER(len=:), ALLOCATABLE :: errmsg

            IF( .NOT. S%IsOpen ) THEN
                CALL S%file%create( FileName, stat, errmsg, mode=VTKHDF_DYNAMIC_MESH )
                IF( stat /= 0 ) ERROR STOP 'VTKHDF: '//errmsg
                ALLOCATE( S%R(SIZE(md%RNames)), S%I(SIZE(md%INames)) )
                DO k = 1, SIZE(md%RNames)
                    S%R(k) = S%file%register_temporal_cell_data( TRIM(md%RNames(k)), 0.0_wp )
                END DO
                IF( LEN_TRIM(md%VName) > 0 ) &
                    S%V = S%file%register_temporal_cell_data( TRIM(md%VName), [0.0_wp, 0.0_wp, 0.0_wp] )
                DO k = 1, SIZE(md%INames)
                    S%I(k) = S%file%register_temporal_cell_data( TRIM(md%INames(k)), 0_INT32 )
                END DO
                S%IsOpen = .True.
            END IF

            CALL S%file%start_time_step( REAL(Time,REAL64) )
            !== fVTKHDF uses 1-based connectivity and start offsets (ncells+1)
            CALL S%file%write_mesh( md%Pts, md%Conn+1, [1_INT32, md%Offs+1], md%Types )
            DO k = 1, SIZE(md%RNames)
                CALL S%file%write_cell_data( S%R(k), md%RData(:,k) )
            END DO
            IF( LEN_TRIM(md%VName) > 0 ) CALL S%file%write_cell_data( S%V, md%VData )
            DO k = 1, SIZE(md%INames)
                CALL S%file%write_cell_data( S%I(k), md%IData(:,k) )
            END DO
            CALL S%file%finalize_time_step
            CALL S%file%flush       ! the file is complete and readable after every step

        END SUBROUTINE Write_Series

    END SUBROUTINE Write_VTKHDF_Step

!******************************************************************************

    SUBROUTINE Close_VTKHDF
        IF( H5Cells%IsOpen  ) CALL H5Cells%file%close
        IF( H5Blocks%IsOpen ) CALL H5Blocks%file%close
        H5Cells%IsOpen  = .False.
        H5Blocks%IsOpen = .False.
    END SUBROUTINE Close_VTKHDF
#endif


END MODULE VTK_Output_Mod
