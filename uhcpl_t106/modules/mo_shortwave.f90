MODULE mo_shortwave

  IMPLICIT NONE

  !      -----------------------------------------------------------------
  !
  ! *    coefficients for the shortwave radiation subroutines
  !
  !      -----------------------------------------------------------------

  !  solar fraction in spectral intervals
  REAL, PARAMETER :: sun_data(2) = (/ 0.441676 ,  0.558324 /)

  REAL :: apad(2,3,7)  !  pade approximants numerator
  REAL :: bpad(2,3,7)  !  pade approximants denominator
  REAL :: d(2,3)       !  transmission limit for infinite absorber amount
  REAL :: rayc(2,6)    !  rayleigh scattering coefficients
  REAL :: rpdh1        !  1 + exponent pressure dependence h2o
  REAL :: rpdu1        !  1 + exponent pressure dependence uniformly mixed gases
  REAL :: rpnh         !  reference pressure factor for h2o
  REAL :: rpnu         !  reference pressure factor for uniformly mixed gases
  REAL :: rswce        !  e-type, h2o continuum absorption coefficient over 0.68-4 mu
  REAL :: rswcp        !  p-type, h2o continuum absorption coefficient over 0.68-4 mu
  REAL :: rtdh2o       !  exponent temperature dependence h2o
  REAL :: rtdumg       !  exponent temperature dependence uniformly mixed gases
  REAL :: rth2o        !  reference temperature h2o
  REAL :: rtumg        !  reference temperature uniformly mixed gases

!
  INTEGER, PRIVATE :: I, J, K
!
      DATA RTDH2O,RTDUMG           / 0.40  , 0.375 /
      DATA RTH2O ,RTUMG            / 240., 240.  /
      DATA RSWCE ,RSWCP            / 0.000 , 0.000 /
!
      DATA (D(1,K),K = 1,3) / 0.00, 0.00, 0.00 /
!* DERIVED FROM HITRAN APRIL 1991
!       H2O:  Pref=300 hPa, Tref=240K, Pdep=0.8
!       O3 :  unchanged
!
      DATA ((APAD(1,I,J),I=1,3),J=1,7) /                   &
    &  0.912418292E+05, 0.000000000E-00, 0.925887084E-04,  &
    &  0.723613782E+05, 0.000000000E-00, 0.129353723E-01,  &
    &  0.596037057E+04, 0.000000000E-00, 0.800821928E+00,  &
    &  0.000000000E-00, 0.000000000E-00, 0.242715973E+02,  &
    &  0.000000000E-00, 0.000000000E-00, 0.878331486E+02,  &
    &  0.000000000E-00, 0.000000000E-00, 0.191559725E+02,  &
    &  0.000000000E-00, 0.000000000E-00, 0.000000000E+00 /
!
      DATA ((BPAD(1,I,J),I=1,3),J=1,7) /                   &
    &  0.912418292E+05, 0.000000000E-00, 0.925887084E-04,  &
    &  0.724555318E+05, 0.000000000E-00, 0.131812683E-01,  &
    &  0.602593328E+04, 0.000000000E-00, 0.812706117E+00,  &
    &  0.100000000E+01, 0.000000000E-00, 0.249863591E+02,  &
    &  0.000000000E-00, 0.000000000E-00, 0.931071925E+02,  &
    &  0.000000000E-00, 0.000000000E-00, 0.252233437E+02,  &
    &  0.000000000E-00, 0.000000000E-00, 0.100000000E+01 /
!
!
      DATA (RAYC(1,K),K=1,6) /                 &
    &  .428937E-01, .890743E+00,-.288555E+01,  &
    &  .522744E+01,-.469173E+01, .161645E+01/
!
!
      DATA (D(2,K),K=1,3) / 0.000000000, 0.000000000, 0.800000000 /
!
!* INTERVAL 2:  0.68 - 4.00 MICRONS
!* DERIVED FROM HITRAN APRIL 1991
!       H2O:  Pref=300 hPa, Tref=240K, Pdep=0.80
!       UMG:  Pref=300 hPa, Tref=240K, Pdep=0.75 (CO2+O2+CH4+N2O+CO)
!       O3 :  unchanged
!
      DATA ((APAD(2,I,J),I=1,3),J=1,7) /                   &
    &  0.376655383E-08, 0.739646016E-08, 0.410177786E+03,  &
    &  0.978576773E-04, 0.131849595E-03, 0.672595424E+02,  &
    &  0.387714006E+00, 0.437772681E+00, 0.000000000E-00,  &
    &  0.118461660E+03, 0.151345118E+03, 0.000000000E-00,  &
    &  0.119079797E+04, 0.233628890E+04, 0.000000000E-00,  &
    &  0.293353397E+03, 0.797219934E+03, 0.000000000E-00,  &
    &  0.000000000E+00, 0.000000000E+00, 0.000000000E+00 /
!
      DATA ((BPAD(2,I,J),I=1,3),J=1,7) /                   &
    &  0.376655383E-08, 0.739646016E-08, 0.410177786E+03,  &
    &  0.979023421E-04, 0.131861712E-03, 0.731185438E+02,  &
    &  0.388611139E+00, 0.437949001E+00, 0.100000000E+01,  &
    &  0.120291383E+03, 0.151692730E+03, 0.000000000E+00,  &
    &  0.130531005E+04, 0.237071130E+04, 0.000000000E+00,  &
    &  0.415049409E+03, 0.867914360E+03, 0.000000000E+00,  &
    &  0.100000000E+01, 0.100000000E+01, 0.000000000E+00 /
!
!
      DATA (RAYC(2,K),K=1,6) /                 &
    &  .697200E-02, .173297E-01,-.850903E-01,  &
    &  .248261E+00,-.302031E+00, .129662E+00/
!


END MODULE mo_shortwave
