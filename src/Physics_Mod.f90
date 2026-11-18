!================================================
!   Physical parameters of the two-phase flow (lattice units), set by the
!   case (e.g. Rising_Bubble_Mod) before the initialization:
!   M0 mobility, wc = 1/(tau_phi + 1/2) of the phase field, W interface width,
!   Sigma surface tension (Beta, kappa), taul/tauh = 3 mu/rho of the light and
!   heavy fluid, Rhol/Rhoh densities, Gy gravity (along y), D characteristic
!   length (bubble diameter)
!================================================
MODULE Physics_Mod
    USE Kinds_Mod, ONLY: wp
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: M0, wc, W, Beta, kappa, taul, tauh, Rhol, Rhoh, Sigma, Gy, D

    REAL(wp) :: M0 = 0.03_wp
    REAL(wp) :: wc = 1._wp/(0.5_wp + 0.03_wp*3._wp)
    REAL(wp) :: W, Sigma, Beta, kappa, taul, tauh, Rhol, Rhoh
    REAL(wp) :: Gy, D

END MODULE Physics_Mod
