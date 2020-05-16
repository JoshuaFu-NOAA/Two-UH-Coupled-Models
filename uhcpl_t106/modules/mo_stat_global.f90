MODULE mo_stat_global

  USE mo_parameters

  IMPLICIT NONE

  ! ---------------------------------------------------------------
  !
  ! module *mo_stat_global* - arrays and parameters for statistics.
  !
  ! ---------------------------------------------------------------

  LOGICAL :: ldiad    !  *true when dynamical diagnostics are computed
                      !  and printed.
  LOGICAL :: lverdia  !  *true when dynamical diagnostics are computed
                      !  and printed for each vertical level.

  REAL, ALLOCATABLE :: voh(:) !  horizontal rms of vorticity.
  REAL, ALLOCATABLE :: dh(:)  !  horizontal rms of divergence.
  REAL, ALLOCATABLE :: qh(:)  !  horizontal average of humidity.
  REAL, ALLOCATABLE :: th(:)  !  horizontal average of temperature.
  REAL, ALLOCATABLE :: xh(:)  !  horizontal average of extra variable

  REAL :: gvo         !  global rms of vorticity.
  REAL :: gd          !  global rms of divergence.
  REAL :: gq          !  global average of humidity.
  REAL :: gt          !  global average of temperature.
  REAL :: gx          !  global average of extra variable
  REAL :: gps         !  mean surface pressure.
  REAL :: gke         !  global kinetic energy.
  REAL :: gpe         !  global potential energy.
  REAL :: gte         !  global total energy.
  REAL :: glq         !  global latent heat energy.
  REAL :: gtpe        !  *gte+glq*.
  REAL :: gqm         !  equivalent water content.
  REAL :: gts         !  land top layer energy.
  REAL :: gtd         !  land deep layer energy.
  REAL :: gws         !  land top layer water content.
  REAL :: gwd         !  land deep layer water content.
  REAL :: gsn         !  land snow equivalent depth.

END MODULE mo_stat_global
