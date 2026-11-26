    REAL(wp), INTENT(IN)    :: inv_dx
    REAL(wp), INTENT(IN)    :: C(-1:nx+2,-1:ny+2), mu(1,1)
    REAL(wp), INTENT(OUT)   :: P (nx,ny), Ux(nx,ny), Uy(nx,ny)
    REAL(wp), INTENT(INOUT) :: h(0:nx+1,0:ny+1,9), g(0:nx+1,0:ny+1,9)
