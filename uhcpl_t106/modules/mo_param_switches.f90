MODULE mo_param_switches

  IMPLICIT NONE

  !
  ! ----------------------------------------------------------------
  !
  ! module *mo_param_switches* switches related to the parameterisations of
  !                            diabatic processes except for radiation.
  !
  ! ----------------------------------------------------------------

  LOGICAL :: lphys    !   *true for parameterisation of diabatic processes.
  LOGICAL :: lvdiff   !   *true for vertical diffusion.
  LOGICAL :: lcond    !   *true for large scale condensation scheme.
  LOGICAL :: lsurf    !   *true for surface exchanges.
  LOGICAL :: lconv    !   *true to allow convection
  LOGICAL :: lgwdrag  !   *true for gravity wave drag scheme
  LOGICAL :: lice     !   *true for sea-ice temperature calculation

  INTEGER :: ndiapfr  !   *frequency of physics budgets.
                      !   if(nradfr.gt.0):rad every *nradfr* time step.
                      !   if(nradfr.lt.0):rad every *-nradfr* hours.


END MODULE mo_param_switches
