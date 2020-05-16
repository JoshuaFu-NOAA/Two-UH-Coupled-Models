!+ for solar orbital parameters.
!+ $Id: orbit.f90,v 1.4 1998/10/28 12:31:38 m214003 Exp $

SUBROUTINE orbit(pclock,pytime,pdisse,pzen1,pzen2,pzen3,prae)

  ! Description:
  !
  ! Computes the solar constant.
  ! 
  ! Method:
  !
  ! This routine computes the solar constant normalised by its
  ! annual mean and three orbital parameters depending on the time
  ! of the day as well as of the year (both in radians). Time scale
  ! origin is 1900/01/01 00 gmt. Also returned is a constant for the
  ! effect of the earth's curvature on the cosine of the solar zenith
  ! angle.
  !
  ! Staightforward. Intermediate variables are the solar
  ! declination and the equation of time.
  !
  ! *orbit* is called from *physc* at the first latitude row.
  ! There are seven dummy arguments: 
  !     *pclock* is the time of the day
  !     *pytime* is the time of the year (both in radians).
  !     *pdisse* is the ratio of the solar constant to its annual mean
  !     *pzen1*, *pzen2* and *pzen3* are zenithal parameters.
  !     *prae*   is the ratio of the height of the equivalent atmosphere 
  !              to the radius of the earth.
  !
  ! Authors:
  !
  ! J. F. Geleyn, ECMWF, June 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  IMPLICIT NONE

  !  Scalar arguments 
  REAL :: pclock, pdisse, prae, pytime, pzen1, pzen2, pzen3

  !  Local scalars: 
  REAL :: zc1yt, zc2yt, zclock, zdecli, zdisse, zeqtim, zrae, zs1yt, zs2yt, &
&      zytime, zzen1, zzen2, zzen3

  !  Local arrays: 
  REAL :: zcdec(5), zcdis(5), zceqt(5)

  !  Intrinsic functions 
  INTRINSIC COS, SIN

  !  Data statements 
  ! *zcdis*, *zcdec* and *zceqt* are arrays for a second order
  ! Fourier development for respectively: Solar constant, solar
  ! declination and equation of time. *zrae* is the value for *prae*.
  DATA zcdis/ + 1.000110, + 0.034221, + 0.001280, + 0.000719, + 0.000077/
  DATA zcdec/ + 0.006918, -0.399912, + 0.070257, -0.006758, + 0.000907/
  DATA zceqt/ + 0.000075, + 0.001868, -0.032077, -0.014615, -0.040849/
  DATA zrae/ + .1277E-02/


  !  Executable statements 

!-- 1. Preliminary setting

  zclock = pclock
  zytime = pytime

!-- 2. Computations

  zc1yt = COS(zytime)
  zs1yt = SIN(zytime)
  zc2yt = zc1yt**2 - zs1yt**2
  zs2yt = 2.*zs1yt*zc1yt
  zdisse = zcdis(1) + zcdis(2)*zc1yt + zcdis(3)*zs1yt + zcdis(4)*zc2yt + &
&      zcdis(5)*zs2yt
  zdecli = zcdec(1) + zcdec(2)*zc1yt + zcdec(3)*zs1yt + zcdec(4)*zc2yt + &
&      zcdec(5)*zs2yt
  zeqtim = zceqt(1) + zceqt(2)*zc1yt + zceqt(3)*zs1yt + zceqt(4)*zc2yt + &
&      zceqt(5)*zs2yt
  zzen1 = SIN(zdecli)
  zzen2 = COS(zdecli)*COS(zclock+zeqtim)

  zzen3 = COS(zdecli)*SIN(zclock+zeqtim)

  pdisse = zdisse
  pzen1 = zzen1
  pzen2 = zzen2
  pzen3 = zzen3
  prae = zrae

  RETURN
END SUBROUTINE orbit
