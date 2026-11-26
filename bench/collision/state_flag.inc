        IF( with_moments ) THEN
#include "state_moments_fusemu.inc"
        ELSE
#include "state_fusemu.inc"
        END IF
