
SUBROUTINE Fill_Ghost_Cells_hg
	USE parameters, ONLY: Blk, LevelMin, LevelMax, nLeaf, LeafBlks, nx, ny
	IMPLICIT NONE

	INTEGER	:: L, I, m, n, Nghbr(8)

	DO L = LevelMin, LevelMax
!$OMP PARALLEL DO private( n, Nghbr )
		DO m = 1, nLeaf
			n = LeafBlks(m)

			IF( L == Blk(n)%Level ) THEN

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

			!	IF( Nghbr(1) < 0 ) THEN
			!		apply BC ...
			!	END IF
			!	IF( Nghbr(2) < 0 ) THEN
			!		apply BC ...
			!	END IF
				IF( Nghbr(3) < 0 ) THEN
					CALL Bounce_Back_3( nx, ny, Blk(n)%h )
					CALL Bounce_Back_3( nx, ny, Blk(n)%g )
				END IF
				IF( Nghbr(4) < 0 ) THEN
					CALL Bounce_Back_4( nx, ny, Blk(n)%h )
					CALL Bounce_Back_4( nx, ny, Blk(n)%g )
				END IF

			END IF
		END DO
!$OMP END PARALLEL DO
	END DO

END

!******************************************************************************

SUBROUTINE Fill_Left_hg( n, L, N1 )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, L, N1

	INTEGER :: L_nbr, child, i, j, Y
	REAL(8) :: h(9, 3,2*ny), H2(9, 0:2,0:ny/2+1)
	REAL(8) :: g(9, 3,2*ny), G2(9, 0:2,0:ny/2+1)

	L_nbr = Blk(N1)%Level

	IF( L == L_nbr ) THEN

		IF( Blk(N1)%IsLeaf ) THEN
			DO j=1,ny
				Blk(n)%h(:, 0,j) = Blk(N1)%h(:, nx,j)
				Blk(n)%g(:, 0,j) = Blk(N1)%g(:, nx,j)
			END DO
		ELSE
		!== saving neighbor's children data to temporary arrays
			child = Blk(N1)%child(1)
			DO J = 1,    ny
			DO I = nx-2, nx
				h(:, nx-I+1,J   ) = Blk(child+1)%h(:, I,J)
				h(:, nx-I+1,J+ny) = Blk(child+3)%h(:, I,J)
				g(:, nx-I+1,J   ) = Blk(child+1)%g(:, I,J)
				g(:, nx-I+1,J+ny) = Blk(child+3)%g(:, I,J)
			END DO
			END DO
		!========================================================

			CALL Biquadratic_from_Fine_f( 3,2*ny, 1,ny, h, Blk(n)%h(:, 0,1:ny) )
			CALL Biquadratic_from_Fine_f( 3,2*ny, 1,ny, g, Blk(n)%g(:, 0,1:ny) )

		END IF

	ELSEIF( L == L_nbr + 1 ) THEN

		Y = (ny/2) * MOD(Blk(n)%J,2)

		DO j = 0, ny/2 + 1
			H2(:, 2:0:-1,j) = Blk(N1)%h(:, nx-1:nx+1,j+Y)
			G2(:, 2:0:-1,j) = Blk(N1)%g(:, nx-1:nx+1,j+Y)
		END DO

		CALL Biquadratic_from_Coarse_f( 3,ny/2+2, 1,ny, H2, Blk(n)%h(:, 0,1:ny) )
		CALL Biquadratic_from_Coarse_f( 3,ny/2+2, 1,ny, G2, Blk(n)%g(:, 0,1:ny) )

	ELSE
		PRINT '(A,4I5)', " Error in SUB_Fill_Left_hg: ", n, N1, L, L_nbr
	END IF

END

!******************************************************************************

SUBROUTINE Fill_Right_hg( n, L, N2 )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, L, N2

	INTEGER :: L_nbr, child, i, j, Y
	REAL(8) :: h(9, 3,2*ny), H2(9, 0:2,0:ny/2+1)
	REAL(8) :: g(9, 3,2*ny), G2(9, 0:2,0:ny/2+1)

	L_nbr = Blk(N2)%Level

	IF( L == L_nbr ) THEN

		IF( Blk(N2)%IsLeaf ) THEN
			DO j=1,ny
				Blk(n)%h(:, nx+1,j) = Blk(N2)%h(:, 1,j)
				Blk(n)%g(:, nx+1,j) = Blk(N2)%g(:, 1,j)
			END DO
		ELSE
		!== saving neighbor's children data to temporary arrays
			child = Blk(N2)%child(1)
			DO J = 1, ny
			DO I = 1, 3
				h(:, I,J   ) = Blk(child  )%h(:, I,J)
				h(:, I,J+ny) = Blk(child+2)%h(:, I,J)
				g(:, I,J   ) = Blk(child  )%g(:, I,J)
				g(:, I,J+ny) = Blk(child+2)%g(:, I,J)
			END DO
			END DO
		!==============================================

			CALL Biquadratic_from_Fine_f( 3,2*ny, 1,ny, h, Blk(n)%h(:, nx+1,1:ny) )
			CALL Biquadratic_from_Fine_f( 3,2*ny, 1,ny, g, Blk(n)%g(:, nx+1,1:ny) )

		END IF

	ELSEIF( L == L_nbr + 1 ) THEN

		Y = (ny/2) * MOD(Blk(n)%J,2)

		DO j = 0, ny/2 + 1
			H2(:, 0:2,j) = Blk(N2)%h(:, 0:2,j+Y)
			G2(:, 0:2,j) = Blk(N2)%g(:, 0:2,j+Y)
		END DO

		CALL Biquadratic_from_Coarse_f( 3,ny/2+2, 1,ny, H2, Blk(n)%h(:, nx+1,1:ny) )
		CALL Biquadratic_from_Coarse_f( 3,ny/2+2, 1,ny, G2, Blk(n)%g(:, nx+1,1:ny) )

	ELSE
		PRINT '(A,4I5)', " Error in SUB_Fill_Right_hg: ", n, N2, L, L_nbr
	END IF

END

!******************************************************************************

SUBROUTINE Fill_Bottom_hg( n, L, N3 )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, L, N3

	INTEGER :: L_nbr, child, i, j, X
	REAL(8) :: h(9, 2*nx,3), H2(9, 0:nx/2+1,0:2)
	REAL(8) :: g(9, 2*nx,3), G2(9, 0:nx/2+1,0:2)

	L_nbr = Blk(N3)%Level

	IF( L == L_nbr ) THEN

		IF( Blk(N3)%IsLeaf ) THEN
			DO i=1,nx
				Blk(n)%h(:, i,0) = Blk(N3)%h(:, i,ny)
				Blk(n)%g(:, i,0) = Blk(N3)%g(:, i,ny)
			END DO
		ELSE
		!== saving neighbor's children data to temporary arrays
			child = Blk(N3)%child(1)
			DO J = ny-2, ny
			DO I = 1,    nx
				h(:, I   ,ny-J+1) = Blk(child+2)%h(:, I,J)
				h(:, I+nx,ny-J+1) = Blk(child+3)%h(:, I,J)
				g(:, I   ,ny-J+1) = Blk(child+2)%g(:, I,J)
				g(:, I+nx,ny-J+1) = Blk(child+3)%g(:, I,J)
			END DO
			END DO
		!============================================

			CALL Biquadratic_from_Fine_f( 2*nx,3, nx,1, h, Blk(n)%h(:, 1:nx,0) )
			CALL Biquadratic_from_Fine_f( 2*nx,3, nx,1, g, Blk(n)%g(:, 1:nx,0) )

		END IF

	ELSEIF( L == L_nbr + 1 ) THEN

		X = (nx/2) * MOD(Blk(n)%I,2)

		DO i = 0, nx/2 + 1
			H2(:, i,2:0:-1) = Blk(N3)%h(:, i+X,ny-1:ny+1)
			G2(:, i,2:0:-1) = Blk(N3)%g(:, i+X,ny-1:ny+1)
		END DO

		CALL Biquadratic_from_Coarse_f( nx/2+2,3, nx,1, H2, Blk(n)%h(:, 1:nx,0) )
		CALL Biquadratic_from_Coarse_f( nx/2+2,3, nx,1, G2, Blk(n)%g(:, 1:nx,0) )

	ELSE
		PRINT '(A,4I5)', " Error in SUB_Fill_Bottom_hg: ", n, N3, L, L_nbr
	END IF

END

!******************************************************************************

SUBROUTINE Fill_Top_hg( n, L, N4 )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, L, N4

	INTEGER :: L_nbr, child, i, j, X
	REAL(8) :: h(9, 2*nx,3), H2(9, 0:nx/2+1,0:2)
	REAL(8) :: g(9, 2*nx,3), G2(9, 0:nx/2+1,0:2)

	L_nbr = Blk(N4)%Level

	IF( L == L_nbr ) THEN

		IF( Blk(N4)%IsLeaf ) THEN
			DO i=1,nx
				Blk(n)%h(:, i,ny+1) = Blk(N4)%h(:, i,1)
				Blk(n)%g(:, i,ny+1) = Blk(N4)%g(:, i,1)
			END DO
		ELSE
		!== saving neighbor's children data to temporary arrays
			child = Blk(N4)%child(1)
			DO J = 1, 3
			DO I = 1, nx
				h(:, I   ,J) = Blk(child  )%h(:, I,J)
				h(:, I+nx,J) = Blk(child+1)%h(:, I,J)
				g(:, I   ,J) = Blk(child  )%g(:, I,J)
				g(:, I+nx,J) = Blk(child+1)%g(:, I,J)
			END DO
			END DO
		!==============================================

			CALL Biquadratic_from_Fine_f( 2*nx,3, nx,1, h, Blk(n)%h(:, 1:nx,ny+1) )
			CALL Biquadratic_from_Fine_f( 2*nx,3, nx,1, g, Blk(n)%g(:, 1:nx,ny+1) )

		END IF

	ELSEIF( L == L_nbr + 1 ) THEN

		X = (nx/2) * MOD(Blk(n)%I,2)

		DO i = 0, nx/2 + 1
			H2(:, i,0:2) = Blk(N4)%h(:, i+X,0:2)
			G2(:, i,0:2) = Blk(N4)%g(:, i+X,0:2)
		END DO

		CALL Biquadratic_from_Coarse_f( nx/2+2,3, nx,1, H2, Blk(n)%h(:, 1:nx,ny+1) )
		CALL Biquadratic_from_Coarse_f( nx/2+2,3, nx,1, G2, Blk(n)%g(:, 1:nx,ny+1) )

	ELSE
		PRINT '(A,4I5)', " Error in SUB_Fill_Top_hg: ", n, N4, L, L_nbr
	END IF

END

!******************************************************************************

SUBROUTINE Fill_Left_Bottom_hg( n, L, N5 )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, L, N5

	INTEGER :: L_nbr, child, i, j
	REAL(8) :: h(9, 3,3), H2(9, 0:2,0:2)
	REAL(8) :: g(9, 3,3), G2(9, 0:2,0:2)

	L_nbr = Blk(N5)%Level

	IF( L == L_nbr ) THEN

		IF( Blk(N5)%IsLeaf ) THEN
			Blk(n)%h(:, 0,0) = Blk(N5)%h(:, nx,ny)
			Blk(n)%g(:, 0,0) = Blk(N5)%g(:, nx,ny)
		ELSE
		!== saving neighbor's children data to temporary arrays
			child = Blk(N5)%child(1)
			DO J = ny-2, ny
			DO I = nx-2, nx
				h(:, nx-I+1,ny-J+1) = Blk(child+3)%h(:, I,J)
				g(:, nx-I+1,ny-J+1) = Blk(child+3)%g(:, I,J)
			END DO
			END DO
		!========================================================

			CALL Biquadratic_from_Fine_f( 3,3, 1,1, h, Blk(n)%h(:, 0,0) )
			CALL Biquadratic_from_Fine_f( 3,3, 1,1, g, Blk(n)%g(:, 0,0) )

		END IF

	ELSEIF( L == L_nbr + 1 ) THEN

		I = nx - (nx/2) * MOD(Blk(n)%I,2) - 1
		J = ny - (ny/2) * MOD(Blk(n)%J,2) - 1

		H2(:, 2:0:-1,2:0:-1) = Blk(N5)%h(:, I:I+2,J:J+2)
		G2(:, 2:0:-1,2:0:-1) = Blk(N5)%g(:, I:I+2,J:J+2)

		CALL Biquadratic_from_Coarse_f( 3,3, 1,1, H2, Blk(n)%h(:, 0,0) )
		CALL Biquadratic_from_Coarse_f( 3,3, 1,1, G2, Blk(n)%g(:, 0,0) )

	ELSE
		PRINT '(A,4I5)', " Error in Fill_Left_Bottom_hg: ", n, N5, L, L_nbr
	END IF

END

!******************************************************************************

SUBROUTINE Fill_Left_Top_hg( n, L, N6 )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, L, N6

	INTEGER :: L_nbr, child, i, j
	REAL(8) :: h(9, 3,3), H2(9, 0:2,0:2)
	REAL(8) :: g(9, 3,3), G2(9, 0:2,0:2)

	L_nbr = Blk(N6)%Level

	IF( L == L_nbr ) THEN

		IF( Blk(N6)%IsLeaf ) THEN
			Blk(n)%h(:, 0,ny+1) = Blk(N6)%h(:, nx,1)
			Blk(n)%g(:, 0,ny+1) = Blk(N6)%g(:, nx,1)
		ELSE
		!== saving neighbor's children data to temporary arrays
			child = Blk(N6)%child(1)
			DO J = 1,    3
			DO I = nx-2, nx
				h(:, nx-I+1,J) = Blk(child+1)%h(:, I,J)
				g(:, nx-I+1,J) = Blk(child+1)%g(:, I,J)
			END DO
			END DO
		!========================================================

			CALL Biquadratic_from_Fine_f( 3,3, 1,1, h, Blk(n)%h(:, 0,ny+1) )
			CALL Biquadratic_from_Fine_f( 3,3, 1,1, g, Blk(n)%g(:, 0,ny+1) )

		END IF

	ELSEIF( L == L_nbr + 1 ) THEN

		I = nx   - (nx/2) * MOD(Blk(n)%I,2) - 1
		J = ny/2 - (ny/2) * MOD(Blk(n)%J,2)

		H2(:, 2:0:-1,0:2) = Blk(N6)%h(:, I:I+2,J:J+2)
		G2(:, 2:0:-1,0:2) = Blk(N6)%g(:, I:I+2,J:J+2)

		CALL Biquadratic_from_Coarse_f( 3,3, 1,1, H2, Blk(n)%h(:, 0,ny+1) )
		CALL Biquadratic_from_Coarse_f( 3,3, 1,1, G2, Blk(n)%g(:, 0,ny+1) )

	ELSE
		PRINT '(A,4I5)', " Error in Fill_Left_Top_hg: ", n, N6, L, L_nbr
	END IF

END

!******************************************************************************

SUBROUTINE Fill_Right_Bottom_hg( n, L, N7 )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, L, N7

	INTEGER :: L_nbr, child, i, j
	REAL(8) :: h(9, 3,3), H2(9, 0:2,0:2)
	REAL(8) :: g(9, 3,3), G2(9, 0:2,0:2)

	L_nbr = Blk(N7)%Level

	IF( L == L_nbr ) THEN

		IF( Blk(N7)%IsLeaf ) THEN
			Blk(n)%h(:, nx+1,0) = Blk(N7)%h(:, 1,ny)
			Blk(n)%g(:, nx+1,0) = Blk(N7)%g(:, 1,ny)
		ELSE
		!== saving neighbor's children data to temporary arrays
			child = Blk(N7)%child(1)
			DO J = ny-2, ny
			DO I = 1,    3
				h(:, I,ny-J+1) = Blk(child+2)%h(:, I,J)
				g(:, I,ny-J+1) = Blk(child+2)%g(:, I,J)
			END DO
			END DO
		!========================================================

			CALL Biquadratic_from_Fine_f( 3,3, 1,1, h, Blk(n)%h(:, nx+1,0) )
			CALL Biquadratic_from_Fine_f( 3,3, 1,1, g, Blk(n)%g(:, nx+1,0) )

		END IF

	ELSEIF( L == L_nbr + 1 ) THEN

		I = nx/2 - (nx/2) * MOD(Blk(n)%I,2)
		J = ny   - (ny/2) * MOD(Blk(n)%J,2) - 1

		H2(:, 0:2,2:0:-1) = Blk(N7)%h(:, I:I+2,J:J+2)
		G2(:, 0:2,2:0:-1) = Blk(N7)%g(:, I:I+2,J:J+2)

		CALL Biquadratic_from_Coarse_f( 3,3, 1,1, H2, Blk(n)%h(:, nx+1,0) )
		CALL Biquadratic_from_Coarse_f( 3,3, 1,1, G2, Blk(n)%g(:, nx+1,0) )

	ELSE
		PRINT '(A,4I5)', " Error in Fill_Right_Bottom_hg: ", n, N7, L, L_nbr
	END IF

END

!******************************************************************************

SUBROUTINE Fill_Right_Top_hg( n, L, N8 )
	USE parameters, ONLY: Blk, nx, ny
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n, L, N8

	INTEGER :: L_nbr, child, i, j
	REAL(8) :: h(9, 3,3), H2(9, 0:2,0:2)
	REAL(8) :: g(9, 3,3), G2(9, 0:2,0:2)

	L_nbr = Blk(N8)%Level

	IF( L == L_nbr ) THEN

		IF( Blk(N8)%IsLeaf ) THEN
			Blk(n)%h(:, nx+1,ny+1) = Blk(N8)%h(:, 1,1)
			Blk(n)%g(:, nx+1,ny+1) = Blk(N8)%g(:, 1,1)
		ELSE
		!== saving neighbor's children data to temporary arrays
			child = Blk(N8)%child(1)
			DO J = 1, 3
			DO I = 1, 3
				h(:, I,J) = Blk(child)%h(:, I,J)
				g(:, I,J) = Blk(child)%g(:, I,J)
			END DO
			END DO
		!========================================================

			CALL Biquadratic_from_Fine_f( 3,3, 1,1, h, Blk(n)%h(:, nx+1,ny+1) )
			CALL Biquadratic_from_Fine_f( 3,3, 1,1, g, Blk(n)%g(:, nx+1,ny+1) )

		END IF

	ELSEIF( L == L_nbr + 1 ) THEN

		I = nx/2 - (nx/2) * MOD(Blk(n)%I,2)
		J = ny/2 - (ny/2) * MOD(Blk(n)%J,2)

		H2(:, 0:2,0:2) = Blk(N8)%h(:, I:I+2,J:J+2)
		G2(:, 0:2,0:2) = Blk(N8)%g(:, I:I+2,J:J+2)

		CALL Biquadratic_from_Coarse_f( 3,3, 1,1, H2, Blk(n)%h(:, nx+1,ny+1) )
		CALL Biquadratic_from_Coarse_f( 3,3, 1,1, G2, Blk(n)%g(:, nx+1,ny+1) )

	ELSE
		PRINT '(A,4I5)', " Error in Fill_Right_Top_hg: ", n, N8, L, L_nbr
	END IF

END
