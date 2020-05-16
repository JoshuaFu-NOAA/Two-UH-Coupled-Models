MODULE mo_hdiff

  USE mo_parameters

  IMPLICIT NONE

  ! ---------------------------------------------------------------
  !
  ! module *mo_hdiff* - coefficients for horizontal diffusion.
  !
  ! ---------------------------------------------------------------

  REAL    :: dampth
  REAL    :: difvo           !  coefficient for vorticity.
  REAL    :: difd            !  coefficient for divergence.
  REAL    :: dift            !  coefficient for temperature.
  REAL, ALLOCATABLE    :: diftcor(:) !  correction profile for temperature
  REAL    :: cdrag           !  drag coefficient in seconds
                             !  or -(drag coefficent) in days
  REAL    :: enstdif         !  *factor by which stratospheric
                             !  horizontal diffusion is increased from one
                             !  level to next level above.
  INTEGER :: nlvstd1         !  *last (uppermost) layer at which
                             !  stratospheric horizontal diffusion is
                             !  enhanced.
  INTEGER :: nlvstd2         !  *first (lowest) layer at which
                             !  stratospheric horizontal diffusion is
                             !  enhanced.
  LOGICAL :: ldiahdf         !  .true. for statistics of horizontal diffusion
  LOGICAL :: ldrag           !  .true. to apply drag to upper 2 levels

END MODULE mo_hdiff
