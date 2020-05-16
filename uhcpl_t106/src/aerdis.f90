!+ parameters for the vertical distributions of aerosols.
!+ $Id: aerdis.f90,v 1.5 1999/04/26 07:47:08 m214003 Exp $

SUBROUTINE aerdis(petah,pvdaes,pvdael,pvdaeu,pvdaed,klevp1,ptrbga,pvobga, &
&      pstbga,paeops,paeopl,paeopu,paeopd,ptrpt,paeadk,paeadm)

  ! Description:
  !
  ! parameters for the vertical distributions of aerosols.
  !
  ! Method:
  !
  ! This routine computes the values *pvdaen* (*n=*s,*l,*u or *d
  ! for sea,land,urban or desert) of a surface-normalised vertical
  ! distribution of aerosols' optical dephts from the argument *petah*
  ! (vertical coordinate) at *klevp1* levels. It also sets values for
  ! non-geographically weighted total optical depths (at 0.55 e-06
  ! wave-length) *paeopn* for the same four types and similear optical
  ! dephts divided by pressure for background well-mixed aerosols
  ! of three types *pmnbga* (*mn*=*tr*,*vo* or *st* for tropospheric,
  ! volcanic (stratospheric ashes) or stratospheric (sulfuric type)).
  ! It finally set values for the power to be applied to a temperature
  ! ratio smaller than one in order to obtain an idex one in the
  ! stratosphere and zero in the troposphere with a relatively smooth
  ! transition (*ptrpt*), as well as for adsorption coefficients for
  ! water to the three type of tropospheric aerosols (*paeadk*) with
  ! a minimum value (in the whole atmosphere) for the sum of the
  ! products of *paeadk* by the optical depths divided by presure
  ! thickness: *paeadm*.
  !
  ! *aerdis* is called from *physc*.
  ! there are sixteen dummy arguments:
  ! 
  !   *petah*  is the vertical coordinate.
  !   *pvdaen* (*n=*s,*l,*u or*d) are the normalised vertical distributions.
  !   *klevp1* is the number of levels.
  !   *pmnbga* (*mn*=*tr*,*vo* or *st*) 
  !            are the background optical depths divided by pressure.
  !   *paeopn* (*n=*s,*l,*u or *d) are the total optical dephts 
  !            for the vertically varying aerosols.
  !   *ptrpt*  is the temperature exponent for the stratospheric definition.
  !   *paeadk* (1,2,3) and
  !   *paeadm* are the constants for the definition of the quantity
  !   of water vapour that will be adsorbed to the dry aerosols to form
  !   moist aerosols.
  !
  ! straightforward, equivalent heigths are given in meters (8434
  ! for the atmosphere) and tropospheric and stratospheric pressure
  ! boundary values are set at 101325 and 19330 *pascal.
  !
  ! Authors:
  !
  ! J. F. Geleyn, ECMWF, November 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  IMPLICIT NONE

  !  Scalar arguments 
  REAL :: paeadm, paeopd, paeopl, paeops, paeopu, pstbga, ptrbga, ptrpt, &
&      pvobga
  INTEGER :: klevp1

  !  Array arguments 
  REAL :: paeadk(3), petah(klevp1), pvdaed(klevp1), pvdael(klevp1), &
&      pvdaes(klevp1), pvdaeu(klevp1)

  !  Local scalars: 
  REAL :: zhsd, zhsl, zhss, zhsu
  INTEGER :: ilevp1, jk

  !  Intrinsic functions 
  INTRINSIC MAX


  !  Executable statements 

!-- 1. Preliminary setting

  ilevp1 = klevp1

!-- 2. Computations

  zhss = MAX(1.,8434./1000.)
  zhsl = MAX(1.,8434./1000.)
  zhsu = MAX(1.,8434./1000.)
  zhsd = MAX(1.,8434./3000.)
  pvdaes(1) = 0.
  pvdael(1) = 0.
  pvdaeu(1) = 0.
  pvdaed(1) = 0.
  IF (ABS(petah(1)) > 0.) THEN
    pvdaes(1) = petah(1)**zhss
    pvdael(1) = petah(1)**zhsl
    pvdaeu(1) = petah(1)**zhsu
    pvdaed(1) = petah(1)**zhsd
  END IF
  DO jk = 2, ilevp1
    pvdaes(jk) = petah(jk)**zhss
    pvdael(jk) = petah(jk)**zhsl
    pvdaeu(jk) = petah(jk)**zhsu
    pvdaed(jk) = petah(jk)**zhsd
  END DO
  ptrbga = 0.03/(101325.-19330.)
  pvobga = 0.007/19330.
  ! PVOBGA=1.E-45
  pstbga = 0.045/19330.
  ! PSTBGA=1.E-45
  paeops = 0.05
  paeopl = 0.2
  paeopu = 0.1
  paeopd = 1.9
  ptrpt = 30.
  paeadk(1) = + .3876E-03
  paeadk(2) = + .6693E-02
  paeadk(3) = + .8563E-03

  paeadm = 2.6E-10

  RETURN
END SUBROUTINE aerdis
