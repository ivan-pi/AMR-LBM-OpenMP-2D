!================================================
!   TecPlot ASCII output
!================================================
MODULE TecPlot_Mod
    USE Kinds_Mod, ONLY: wp
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: TecPlot_CellCenter

CONTAINS


SUBROUTINE TecPlot_CellCenter( Time )
    USE Parameters_Mod, ONLY: Blk, nLeaf, LeafBlks, nx, ny, L0
    USE D2Q9_Mod,   ONLY: Sigma, D
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: Time    ! physical time written into the zone header

    INTEGER :: I, J, dx, b, m
    INTEGER :: Nn, Ne, N, E
    INTEGER, ALLOCATABLE :: connect(:,:), X(:), Y(:)
    REAL(wp), ALLOCATABLE :: u(:), v(:), c(:), p(:)

    CHARACTER(Len= 8):: Id
    CHARACTER(Len=32):: FileName
    INTEGER, SAVE    :: counter = -1

    counter = counter + 1

    IF( counter < 100 ) THEN
        WRITE( Id,'(I2.2)' ) counter
    ELSE
        WRITE( Id,'(I0)' ) counter
    END IF

    FileName = 'Rising_Bubble_2D_' // TRIM(Id) // '.dat'

    OPEN (2, file = FileName)
    WRITE(2,*) 'Variables = X, Y, Ux, Uy, C, P'

    Nn = nLeaf * (nx+1) * (ny+1)
    Ne = nLeaf *  nx    *  ny

    WRITE(2,20) 'Zone N =', Nn, ', E =', Ne, ', ZONETYPE=FEquadrilateral, DATAPACKING=Block'
    WRITE(2,* ) 'VARLOCATION=([3,4,6]=CELLCENTERED), T = "', Time, '"'
20  FORMAT(2(A,I6),A)

    !----   allocate the connectivity matrix for writing output
    ALLOCATE( connect(4,Ne) )
    ALLOCATE( X(Nn), Y(Nn), u(Ne), v(Ne), c(Nn), p(Ne) )

    !----   construct the connectivity matrix
    N = 1
    E = 1
    DO m = 1, nLeaf
        b = LeafBlks(m)

        dx = Blk(b)%dx

        DO J = 1, ny+1
        DO I = 1, nx+1

            X(N) = Blk(b)%I * dx * nx + dx * (I-1)
            Y(N) = Blk(b)%J * dx * ny + dx * (J-1)

            c(N) = 0.25_wp * ( Blk(b)%C(i-1,j-1) + Blk(b)%C(i,j-1) &
                            + Blk(b)%C(i-1,j  ) + Blk(b)%C(i,j  )) - 0.5_wp

            IF( I <= nx .AND. J <= ny ) THEN
                u(E) = Blk(b)%Ux(I,J) !/U0
                v(E) = Blk(b)%Uy(I,J) !/U0

                p(E) = Blk(b)%P(I,J) * 0.5_wp*D/Sigma

                connect(1,E) = N
                connect(2,E) = N + 1
                connect(3,E) = N + 1 + nx+1
                connect(4,E) = N     + nx+1

                E = E + 1
            END IF

            N = N + 1

        END DO
        END DO
    END DO

    WRITE(2,'(9F10.3)') ( X(i) / REAL(L0,wp), i=1,Nn )
    WRITE(2,'(9F10.3)') ( Y(i) / REAL(L0,wp), i=1,Nn )
    WRITE(2,'(8E14.6)') ( u(i), i=1,Ne )
    WRITE(2,'(8E14.6)') ( v(i), i=1,Ne )
    WRITE(2,'(9E14.6)') ( c(i), i=1,Nn )
    WRITE(2,'(8E14.6)') ( p(i), i=1,Ne )
    WRITE(2,'(8I7)'   ) ( connect(1:4,i), i=1,Ne )

END

END MODULE TecPlot_Mod
