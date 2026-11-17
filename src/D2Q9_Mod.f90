
MODULE D2Q9_Mod
    USE Kinds_Mod, ONLY: wp

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: ex, ey, Wa, M0, wc, W, Beta, kappa, taul, tauh, Rhol, Rhoh, Sigma, Gy, D, &
              Ar, Bo

    INTEGER,PARAMETER :: ex(9) = [0, 1, 0,-1, 0, 1,-1,-1, 1]
    INTEGER,PARAMETER :: ey(9) = [0, 0, 1, 0,-1, 1, 1,-1,-1]
    REAL(wp),PARAMETER :: Wa(9) = [16,4, 4, 4, 4, 1, 1, 1, 1] / 36._wp

    !== set by the case (Set_Rising_Bubble_Parameters):
    !   M0 mobility, wc = 1/(tau_phi + 1/2) of the phase field, W interface width,
    !   Sigma surface tension (Beta, kappa), taul/tauh = 3 mu/rho of the light and
    !   heavy fluid, Rhol/Rhoh densities, Gy gravity, D bubble diameter,
    !   Bo Eotvos number, Ar gravity Reynolds number Re_h
    REAL(wp) :: M0 = 0.03_wp
    REAL(wp) :: wc = 1._wp/(0.5_wp + 0.03_wp*3._wp)

    REAL(wp) :: W, Sigma, Beta, kappa, taul, tauh, Rhol, Rhoh
    REAL(wp) :: Gy, D, Ar, Bo

END MODULE D2Q9_Mod
