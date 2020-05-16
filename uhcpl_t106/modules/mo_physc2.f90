MODULE mo_physc2

  USE mo_parameters

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module *mo_physc2* constants to communicate between the main program
  !                    and the physical subroutines (except radiation ones).
  !
  ! ----------------------------------------------------------------

  REAL :: clam            !  *asymptotic mixing length for momentum.
  REAL :: ckap            !  *karman constant.
  REAL :: cb              !  *stability parameter near neutrality.
  REAL :: cc              !  *stability parameter for unstable cases.
  REAL :: cd              !  *stability parameter for stable cases.
  REAL :: cchar           !  *charnock constant.
  REAL :: crgcg           !  *heat capacity x density of the soil.
  REAL :: cdif            !  *thermal diffusivity of soil.
  REAL :: cwcap           !  *soil water content at field capacity (=wsmax)
  REAL :: cwcrit          !  *soil water content above which evaporation
                          !   occurs at potential rate
  REAL :: cwpwp           !  *soil water content at wilting point
  REAL :: clice           !  *thermal conductivity factor of snow/ice
  REAL :: cgh2o           !  *heat capacity of water
  REAL :: cqwevap         !  *inverse of (cwcrit-cwpwp)
  REAL :: cqcon           !  *soil hydraulic conductivity
  REAL :: cqdif           !  *soil hydraulic diffusivity
  REAL :: cqwsbcr         !  *inverse of critical wetness for bare ground
  REAL :: cd1             !  *thickness of 1st ground level.
  REAL :: cd2             !  *thickness of 2nd ground level.
  ! ----------------------- ADDED FOR AMIP2 -----------------------------
  REAL :: cqsncr          !  *inverse of equivalent water height for
                          !   snow cover calculations
  REAL :: cdel(5)         !  *thickness of soil layers
  REAL :: cmid(5)         !  *depth of mids of soil layers
  REAL :: csncri          !  *critical snow depth for soil computations
  ! ---------------------------------------------------------------------
  REAL :: cwlmax          !  *maximum moisture content of
                          !   the skin reservoir
  REAL :: cvdifts         !   factor for timestep weighting
                          !   in *rhs* of *vdiff* and *scv*.
  REAL, ALLOCATABLE :: cevapcu(:) !   *evaporation coefficient for kuo0.
  REAL :: ctfreez         !   temperature at which sea
                          !   starts freezing/melting
  REAL :: cz0ice          !   roughness over sea-ice
  REAL :: corvari         !  *minimum sub-grid scale orography variance
                          !   to allow surface runoff
  REAL :: corvars         !  *sub-grid scale orography variance
                          !   calibration in surface runoff

END MODULE mo_physc2
