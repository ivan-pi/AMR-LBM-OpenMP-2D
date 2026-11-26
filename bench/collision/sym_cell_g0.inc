        vxx   = vx*vx
        vyy   = vy*vy
        base0 = 1._wp - 1.5_wp*( vxx + vyy )
        uDc   = vx*DcDx + vy*DcDy
        D0    = 0.5_wp*uDc - ( vx*DmDx + vy*DmDy )
        rho3  = Rho*c13
        hF    = (Rhoh - Rho)*hGy
        amu   = dRho3 + xmu
        Pr    = Pc - rho3
        r0    = rho3 - hF*vy
        rp    = r0 + hF
        rm    = r0 - hF
        Q     = 2._wp*( rho3 + hF*vy ) + amu*uDc
        t  = MIN( MAX( C0, 0._wp ), 1._wp ) * dInvTau + invTaul
        km = 0.25_wp*( 2._wp - t )/( 2._wp + t )
        m8 = km*( ( (G_(2) + G_(4)) - (G_(3) + G_(5)) ) &
                - ( ws45*(vxx - vyy)*Q + ws3*( 2._wp*vy*hF - amu*( vx*eC2 - vy*eC3 ) ) ) )
        m9 = km*( ( (G_(6) + G_(8)) - (G_(7) + G_(9)) ) &
                - ( wd18*vx*vy*Q - wd12*hF*vx - wd3*amu*( (vx + vy)*eC6 - (vy - vx)*eC7 ) ) )
