        g2  = DcDx*DcDx + DcDy*DcDy
        Rho = Rhol + C0*dRho
        t   = MIN( MAX( C0, 0._wp ), 1._wp ) * dInvTau + invTaul
        km  = 0.125_wp + 1.e-3_wp*t
        rRho = 1._wp + 1.e-3_wp*Rho
        rn   = MERGE( 1.e-3_wp + g2, 0._wp, g2 > grad2_min )
        ni = DcDx*rn
        nj = DcDy*rn
