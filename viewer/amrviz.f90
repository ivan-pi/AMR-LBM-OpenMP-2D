!================================================
!   Fortran interface of the amrviz C library (amrviz.h)
!
!   Fields are numbered from 0 as in C. amrviz_block_field is generic over
!   real32 and real64 block arrays: cell (i,j), i,j from 0, is element
!   offset + i + j*ld of a. The packed calls return C pointers to the
!   viewer's host buffers, e.g. as the destination of omp_target_memcpy
!   from a device array, or for C_F_POINTER.
!================================================
MODULE amrviz
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_LONG, C_DOUBLE, C_FLOAT, C_CHAR, C_PTR, C_FUNPTR
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: amrviz_open, amrviz_set_field_name, amrviz_is_open, amrviz_begin_snapshot, &
              amrviz_add_block, amrviz_block_field, amrviz_blocks_buffer, amrviz_field_buffer, &
              amrviz_set_copy, amrviz_field_packed, &
              amrviz_end_snapshot, amrviz_wait_if_paused, amrviz_close

    INTERFACE
        FUNCTION amrviz_open( win_w, win_h, nx, ny, dom_w, dom_h, nfields ) RESULT( stat ) &
                BIND(C, name='amrviz_open')
            IMPORT :: C_INT, C_DOUBLE
            INTEGER(C_INT), VALUE :: win_w, win_h, nx, ny, nfields
            REAL(C_DOUBLE), VALUE :: dom_w, dom_h
            INTEGER(C_INT) :: stat
        END FUNCTION amrviz_open

        SUBROUTINE amrviz_set_field_name( field, name ) BIND(C, name='amrviz_set_field_name')
            IMPORT :: C_INT, C_CHAR
            INTEGER(C_INT), VALUE :: field
            CHARACTER(kind=C_CHAR), INTENT(IN) :: name(*)
        END SUBROUTINE amrviz_set_field_name

        FUNCTION amrviz_is_open() RESULT( is_open ) BIND(C, name='amrviz_is_open')
            IMPORT :: C_INT
            INTEGER(C_INT) :: is_open
        END FUNCTION amrviz_is_open

        SUBROUTINE amrviz_begin_snapshot( time, step ) BIND(C, name='amrviz_begin_snapshot')
            IMPORT :: C_DOUBLE, C_LONG
            REAL(C_DOUBLE),  VALUE :: time
            INTEGER(C_LONG), VALUE :: step
        END SUBROUTINE amrviz_begin_snapshot

        FUNCTION amrviz_add_block( level, x0, y0, dx ) RESULT( slot ) BIND(C, name='amrviz_add_block')
            IMPORT :: C_INT, C_DOUBLE
            INTEGER(C_INT), VALUE :: level
            REAL(C_DOUBLE), VALUE :: x0, y0, dx
            INTEGER(C_INT) :: slot
        END FUNCTION amrviz_add_block

        !== packed snapshots: the block table (4 floats per block: x0, y0, dx,
        !   level) and the field arrays (nblocks*nx*ny floats) to fill
        FUNCTION amrviz_blocks_buffer( nblocks ) RESULT( meta ) BIND(C, name='amrviz_blocks_buffer')
            IMPORT :: C_INT, C_PTR
            INTEGER(C_INT), VALUE :: nblocks
            TYPE(C_PTR) :: meta
        END FUNCTION amrviz_blocks_buffer

        FUNCTION amrviz_field_buffer( field ) RESULT( buf ) BIND(C, name='amrviz_field_buffer')
            IMPORT :: C_INT, C_PTR
            INTEGER(C_INT), VALUE :: field
            TYPE(C_PTR) :: buf
        END FUNCTION amrviz_field_buffer

        !== copy hook (a BIND(C) function: int copy(void *dst, const void *src,
        !   size_t bytes, void *user)) and a field from any pointer it can read
        SUBROUTINE amrviz_set_copy( copy, user ) BIND(C, name='amrviz_set_copy')
            IMPORT :: C_FUNPTR, C_PTR
            TYPE(C_FUNPTR), VALUE :: copy
            TYPE(C_PTR),    VALUE :: user
        END SUBROUTINE amrviz_set_copy

        SUBROUTINE amrviz_field_packed( field, src ) BIND(C, name='amrviz_field_packed')
            IMPORT :: C_INT, C_PTR
            INTEGER(C_INT), VALUE :: field
            TYPE(C_PTR),    VALUE :: src
        END SUBROUTINE amrviz_field_packed

        SUBROUTINE amrviz_end_snapshot() BIND(C, name='amrviz_end_snapshot')
        END SUBROUTINE amrviz_end_snapshot

        SUBROUTINE amrviz_wait_if_paused() BIND(C, name='amrviz_wait_if_paused')
        END SUBROUTINE amrviz_wait_if_paused

        SUBROUTINE amrviz_close( wait_for_user ) BIND(C, name='amrviz_close')
            IMPORT :: C_INT
            INTEGER(C_INT), VALUE :: wait_for_user
        END SUBROUTINE amrviz_close
    END INTERFACE

    INTERFACE amrviz_block_field
        SUBROUTINE amrviz_block_field_d( slot, field, a, ld, offset ) BIND(C, name='amrviz_block_field_d')
            IMPORT :: C_INT, C_DOUBLE
            INTEGER(C_INT), VALUE :: slot, field, ld, offset
            REAL(C_DOUBLE), INTENT(IN) :: a(ld,*)
        END SUBROUTINE amrviz_block_field_d
        SUBROUTINE amrviz_block_field_f( slot, field, a, ld, offset ) BIND(C, name='amrviz_block_field_f')
            IMPORT :: C_INT, C_FLOAT
            INTEGER(C_INT), VALUE :: slot, field, ld, offset
            REAL(C_FLOAT), INTENT(IN) :: a(ld,*)
        END SUBROUTINE amrviz_block_field_f
    END INTERFACE amrviz_block_field

END MODULE amrviz
