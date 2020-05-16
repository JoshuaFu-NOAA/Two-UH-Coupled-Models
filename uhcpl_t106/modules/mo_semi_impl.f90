MODULE mo_semi_impl

  USE mo_parameters

  IMPLICIT NONE

  ! ---------------------------------------------------------------
  !
  ! module *mo_semi_impl* - quantities needed for the semi-implicit scheme
  !
  ! ---------------------------------------------------------------

  REAL :: betadt       !  =0 : explicit scheme(for *d*,*t*,*alps*).
                       !  =1.:semi implicit scheme.
  REAL :: betazq       !  =0 : explicit scheme(for *v0*,*q*).
                       !  =1.:semi implicit scheme.
  REAL :: apr          !  *reference surface pressure for semi-implicit scheme.
  REAL :: tr           !  *reference temperature for semi-implicit scheme.
  REAL :: ulat         !  *latitudeof maximum !u!+!v! (real winds).
  REAL :: uvmax        !  *max(!u!+!v!) (real winds).
  REAL :: ulm          !  *linearization wind profile on the latitude
                       !  and at the level of maximum !u!+!v!.
  REAL :: vcrit        !  *critical velocity above which
                       !  horizontal diffusion is enhanced for
                       !  t63 with dt=20min.
  REAL :: vcheck       !  threshold value for check of high windspeed
  REAL :: hdamp        !  *damping factor for strong
                       !  stratospheric damping.
  REAL, ALLOCATABLE :: vmax(:)

  INTEGER :: nulev     !  *level of maximum !u!+!v! (real winds).

END MODULE mo_semi_impl
