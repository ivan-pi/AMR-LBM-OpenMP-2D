        g2 = DcDx*DcDx + DcDy*DcDy
        rn = MERGE( 1._wp, 0._wp, g2 > grad2_min ) / SQRT( MAX( g2, grad2_min ) )
        ni = DcDx*rn
        nj = DcDy*rn
