        xmu = c4B*C0*(C0-1._wp)*(C0-0.5_wp) - cK*( ( C_(-1,-1) + C_(+1,-1) + C_(-1,+1) + C_(+1,+1) ) &
            + 4._wp*( C_(+0,-1) + C_(-1,+0) + C_(+1,+0) + C_(+0,+1) ) - 20._wp*C0 )
        S0 = G_(1) + ( (G_(2) + G_(3)) + (G_(4) + G_(5)) ) + ( (G_(6) + G_(7)) + (G_(8) + G_(9)) )
        Sx = ( G_(2) - G_(4) ) + ( (G_(6) - G_(7)) + (G_(9) - G_(8)) )
        Sy = ( G_(3) - G_(5) ) + ( (G_(6) + G_(7)) - (G_(8) + G_(9)) )
        Rho  = Rhol + C0*dRho
        rRho = 1._wp/Rho
        vx = ( 3._wp*Sx + 0.5_wp*( xmu*DcDx ) )*rRho
        vy = ( 3._wp*Sy + 0.5_wp*( xmu*DcDy + (Rhoh - Rho)*Gy ) )*rRho
        Pc = S0 + dRho6*( vx*DcDx + vy*DcDy )
        UX_ = vx
        UY_ = vy
        P_  = Pc
