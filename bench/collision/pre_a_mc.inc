        g2 = DcDx*DcDx + DcDy*DcDy
        T_GX = DcDx
        T_GY = DcDy
        T_RN = MERGE( 1._wp, 0._wp, g2 > grad2_min ) / SQRT( MAX( g2, grad2_min ) )
        t = MIN( MAX( C0, 0._wp ), 1._wp ) * dInvTau + invTaul
        T_KM = 0.25_wp*( 2._wp - t )/( 2._wp + t )
        T_RR = 1._wp/( Rhol + C0*dRho )
