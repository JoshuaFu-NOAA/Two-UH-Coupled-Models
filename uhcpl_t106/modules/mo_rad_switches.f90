MODULE mo_rad_switches

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module *mo_rad_switches* switches related to the radiation scheme.
  !
  ! ----------------------------------------------------------------

  LOGICAL :: lrad     !    *true for radiation.
  LOGICAL :: ldiur    !    *true for diurnal cycle.
  LOGICAL :: lsolc    !    *true for solar clear sky diagnostic
  LOGICAL :: laer     !    *true for aerosols
  LOGICAL :: lcfc     !    *true for cfc's
  LOGICAL :: lgadsrh  !    *true for rel.hum. dependency
                      !     of *gads* aerosol optical parameters
  REAL    :: co2fac

  INTEGER :: nmonth
  INTEGER :: nradfr   !    *frequency of full radiation computations.
                      !    if(nradfr.gt.0):rad every *nradfr* time step.
                      !    if(nradfr.lt.0):rad every *-nradfr* hours.
  INTEGER :: nradpfr  !    *print frequency for radiation statistics
                      !    (in number of radiation steps).
  INTEGER :: nradpla  !    *print radiation statistics every *nradpla*
                      !    latitude line.

END MODULE mo_rad_switches
