!================================================
!   amrviz with fields in GPU memory (Fortran, OpenMP target offload)
!
!   The Fortran counterpart of sycl_demo.cpp: a quadtree of 8x8-cell blocks
!   refined around a moving ring on the host, two fields with one ghost
!   layer in a device-resident pool, and an offloaded kernel that packs the
!   leaf-block interiors into one array per field. The packed fields reach
!   the viewer in the two ways amrviz supports:
!     field 0: omp_target_memcpy straight into amrviz_field_buffer(0)
!     field 1: amrviz_field_packed(1, device address) through a copy hook
!
!   Build with viewer/examples/CMakeLists.txt (ifx), run ./offload_demo [frames].
!================================================
MODULE Demo_Copy_Mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_PTR, C_SIZE_T
    USE OMP_LIB, ONLY: omp_target_memcpy, omp_get_initial_device, omp_get_default_device
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: Copy_From_Device

CONTAINS

    !== copy hook for amrviz_field_packed: device address -> host buffer
    FUNCTION Copy_From_Device( dst, src, bytes, user ) RESULT( stat ) BIND(C)
        TYPE(C_PTR),       VALUE :: dst, src, user
        INTEGER(C_SIZE_T), VALUE :: bytes
        INTEGER(C_INT) :: stat
        stat = omp_target_memcpy( dst, src, bytes, 0_C_SIZE_T, 0_C_SIZE_T, &
                                  omp_get_initial_device(), omp_get_default_device() )
    END FUNCTION Copy_From_Device

END MODULE Demo_Copy_Mod

!******************************************************************************

PROGRAM Offload_Demo
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_LONG, C_DOUBLE, C_SIZE_T, &
                                           C_LOC, C_FUNLOC, C_F_POINTER, C_NULL_PTR, C_NULL_CHAR
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: sp => REAL32
    USE OMP_LIB, ONLY: omp_target_memcpy, omp_get_initial_device, omp_get_default_device, omp_get_num_devices
    USE amrviz, ONLY: amrviz_open, amrviz_set_field_name, amrviz_is_open, amrviz_begin_snapshot, &
                      amrviz_blocks_buffer, amrviz_field_buffer, amrviz_set_copy, amrviz_field_packed, &
                      amrviz_end_snapshot, amrviz_wait_if_paused, amrviz_close
    USE Demo_Copy_Mod, ONLY: Copy_From_Device
    IMPLICIT NONE

    INTEGER, PARAMETER :: NX = 8, NY = 8, L0 = 256, LMAX = 5, NFIELDS = 2, CAP = 4096

    REAL(sp), ALLOCATABLE         :: pool(:,:,:,:)          ! (0:NX+1, 0:NY+1, CAP, NFIELDS), device
    REAL(sp), ALLOCATABLE, TARGET :: pack(:,:,:,:)          ! (NX, NY, CAP, NFIELDS), device
    REAL(sp) :: leaves(4,CAP)                               ! x0, y0, dx, level
    REAL(sp), POINTER :: meta(:,:)
    REAL(sp) :: t, cx, cy, R, x, y, rr
    INTEGER  :: frames, frame, n, b, i, j, f, nargs, stat
    CHARACTER(len=16) :: arg

    frames = 2000
    nargs = COMMAND_ARGUMENT_COUNT()
    IF( nargs >= 1 ) THEN
        CALL GET_COMMAND_ARGUMENT( 1, arg )
        READ( arg, * ) frames
    END IF
    PRINT '(A,I0)', 'offload devices: ', omp_get_num_devices()

    ALLOCATE( pool(0:NX+1, 0:NY+1, CAP, NFIELDS), pack(NX, NY, CAP, NFIELDS) )
!$omp target enter data map(alloc: pool, pack, leaves)

    IF( amrviz_open( 800_C_INT, 800_C_INT, NX, NY, REAL(L0,C_DOUBLE), REAL(L0,C_DOUBLE), NFIELDS ) /= 0 ) STOP 1
    CALL amrviz_set_field_name( 0_C_INT, 'C'    // C_NULL_CHAR )
    CALL amrviz_set_field_name( 1_C_INT, 'wave' // C_NULL_CHAR )
    CALL amrviz_set_copy( C_FUNLOC(Copy_From_Device), C_NULL_PTR )

    DO frame = 0, frames-1
        IF( amrviz_is_open() == 0 ) EXIT
        t  = 0.01_sp * frame
        cx = 128 + 60 * COS(t)
        cy = 128 + 60 * SIN(t)
        R  = 40
        CALL Build_Mesh( cx, cy, R, leaves, n )
!$omp target update to(leaves(:,1:n))

        !== "solver": both fields in the pool, ghost layers included
!$omp target teams distribute parallel do collapse(3) private(x, y, rr)
        DO b = 1, n
            DO j = 0, NY+1
                DO i = 0, NX+1
                    x  = leaves(1,b) + (i - 0.5_sp) * leaves(3,b)
                    y  = leaves(2,b) + (j - 0.5_sp) * leaves(3,b)
                    rr = SQRT( (x-cx)**2 + (y-cy)**2 )
                    pool(i,j,b,1) = 0.5_sp * (1 + TANH( (rr - R) / 3 ))
                    pool(i,j,b,2) = SIN( 0.1_sp*x - 3*t ) * COS( 0.07_sp*y )
                END DO
            END DO
        END DO

        !== pack the interiors of the leaf blocks, field by field
!$omp target teams distribute parallel do collapse(4)
        DO f = 1, NFIELDS
            DO b = 1, n
                DO j = 1, NY
                    DO i = 1, NX
                        pack(i,j,b,f) = pool(i,j,b,f)
                    END DO
                END DO
            END DO
        END DO

        CALL amrviz_begin_snapshot( REAL(t,C_DOUBLE), INT(frame,C_LONG) )
        CALL C_F_POINTER( amrviz_blocks_buffer( INT(n,C_INT) ), meta, [4, n] )
        meta = leaves(:,1:n)
!$omp target data use_device_addr(pack)
        stat = omp_target_memcpy( amrviz_field_buffer(0_C_INT), C_LOC(pack(1,1,1,1)), &
                                  INT(n,C_SIZE_T)*NX*NY*4, 0_C_SIZE_T, 0_C_SIZE_T, &
                                  omp_get_initial_device(), omp_get_default_device() )
        CALL amrviz_field_packed( 1_C_INT, C_LOC(pack(1,1,1,2)) )
!$omp end target data
        IF( stat /= 0 ) STOP 'omp_target_memcpy failed'
        CALL amrviz_end_snapshot
        CALL amrviz_wait_if_paused
    END DO

    CALL amrviz_close( 1_C_INT )
!$omp target exit data map(delete: pool, pack, leaves)

CONTAINS

    !== leaf blocks of a quadtree refined to LMAX where a block touches the ring
    SUBROUTINE Build_Mesh( cx, cy, R, leaves, n )
        REAL(sp), INTENT(IN)  :: cx, cy, R
        REAL(sp), INTENT(OUT) :: leaves(:,:)
        INTEGER,  INTENT(OUT) :: n
        REAL(sp) :: stack(4,4*LMAX+4), blk(4), w, px, py, dmin, dmax, h
        INTEGER  :: top, k

        n = 0
        top = 1
        stack(:,1) = [ 0._sp, 0._sp, REAL(L0/NX,sp), 0._sp ]
        DO WHILE( top > 0 )
            blk = stack(:,top)
            top = top - 1
            w  = blk(3) * NX
            px = MIN( MAX(cx, blk(1)), blk(1) + w )
            py = MIN( MAX(cy, blk(2)), blk(2) + w )
            dmin = HYPOT( px - cx, py - cy )
            dmax = HYPOT( MAX( ABS(blk(1)-cx), ABS(blk(1)+w-cx) ), MAX( ABS(blk(2)-cy), ABS(blk(2)+w-cy) ) )
            IF( blk(4) < LMAX .AND. ( (dmin <= R+6 .AND. dmax >= R-6) .OR. blk(4) < 1 ) ) THEN
                h = 0.5_sp * blk(3)
                DO k = 0, 3
                    top = top + 1
                    stack(:,top) = [ blk(1) + MOD(k,2)*h*NX, blk(2) + (k/2)*h*NX, h, blk(4) + 1 ]
                END DO
            ELSE IF( n < SIZE(leaves,2) ) THEN
                n = n + 1
                leaves(:,n) = blk
            END IF
        END DO
    END SUBROUTINE Build_Mesh

END PROGRAM Offload_Demo
