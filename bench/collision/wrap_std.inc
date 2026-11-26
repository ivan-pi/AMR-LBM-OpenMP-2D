SUBROUTINE VNAME( n )
    INTEGER, INTENT(IN) :: n
    REAL(wp) :: mu(nx,ny), inv_dx
    inv_dx = inv_dx_level( Blk(n)%Level )
    CALL Chemical_Potential_AMR( inv_dx, Blk_C(:,:,n), mu )
    CALL KNAME( inv_dx, Blk_C(:,:,n), mu, Blk_P(:,:,n), &
                Blk_Ux(:,:,n), Blk_Uy(:,:,n), Blk_h(:,:,:,n), Blk_g(:,:,:,n) )
END SUBROUTINE VNAME
