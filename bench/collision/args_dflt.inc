    REAL(wp), INTENT(IN)    :: inv_dx
    REAL(wp), INTENT(IN)    :: C(nx+4,ny+4), mu(nx,ny)
    REAL(wp), INTENT(IN)    :: P (nx,ny), Ux(nx,ny), Uy(nx,ny)
    REAL(wp), INTENT(INOUT) :: h(nx+2,ny+2,9), g(nx+2,ny+2,9)
