    REAL(wp), INTENT(IN)    :: inv_dx
    REAL(wp), INTENT(IN)    :: C(-1:ldp-2,-1:ny+2), mu(1,1)
    REAL(wp), INTENT(IN)    :: P (nx,ny), Ux(nx,ny), Uy(nx,ny)
    REAL(wp), INTENT(INOUT) :: h(0:ldp-1,0:ny+1,9), g(0:ldp-1,0:ny+1,9)
