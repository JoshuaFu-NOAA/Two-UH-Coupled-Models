!+ organizes the longwave calculations
!+ $Id: lw.f90,v 1.5 1998/12/11 13:23:33 m214003 Exp $

SUBROUTINE lw(kdlon,kflev,kewaer,kaerh,kaer,kcfc,pcco2,pcldlw,pdp,ppmb,pqof, &
&      ptl,paer,ptave,pwv,pflux,pfluc,pcfcabs)

  ! Description:
  !
  ! Organizes the longwave calculations.
  !
  ! Method:
  !
  ! Depending on kmode, computes longwave fluxes and/or
  ! radiances
  !
  ! *lw* is called from *radlsw*
  !
  ! Explicit arguments :
  ! pcco2  :                       ; concentration in CO2 [Pa/Pa]
  ! pqof   : (kdlon,kflev)         ; concentration in ozone [Pa/Pa]
  ! ptave  : (kdlon,kflev)         ; temperature
  ! ptl    : (kdlon,kflev+1)       ; half level temperature
  ! ppmb   : (kdlon,kflev+1)       ; half level pressure
  ! pwv    : (kdlon,kflev)         ; specific humidity pa/pa
  ! pcldlw : (kdlon,kflev)         ; cloud fractional cover
  ! paer   : (kdlon,kflev,5+kewaer); aerosol optical thickness (1,,5)
  !                                  Tanre et al., 1984
  !                                  aerosol mass mixing ratio [kg/kg]
  !                                  (6,.,5+kewaer) computed in echam4
  ! ==== outputs ===
  ! if kmode = 0, 1, 2
  ! pflux(kdlon,2,kflev)           ; radiative fluxes :
  !                      1  ==>  upward   flux total
  !                      2  ==>  downward flux total
  ! pfluc(kdlon,2,kflev)           ; radiative fluxes clear sky:
  !                      1  ==>  upward   flux total
  !                      2  ==>  downward flux total
  !
  ! 1. Computes the pressure and temperature weighted amounts of
  !    absorbers.
  ! 2. Computes the planck functions on the interfaces and the
  !    gradient of planck functions in the layers.
  ! 3. Performs the vertical integration distinguishing the con-
  !    tributions of the adjacent and distant layers and those from the
  !    boundaries.
  ! 4. Computes the clear-sky downward and upward emissivities.
  ! 5. Introduces the effects of the clouds on the fluxes.
  !
  ! Reference:
  ! See radiation's part of the model's documentation and
  ! ECMWF research department documentation of the "in core model"
  !
  ! Authors:
  !
  ! J.-J. Morcrette, ECMWF, July 1989, original source
  ! R. Van Dorland, KNMI, May 1995, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_longwave
  USE mo_radiation

  IMPLICIT NONE

  !  Scalar arguments 
  REAL :: pcco2
  INTEGER :: kaer, kcfc, kdlon, kewaer, kflev

  !  Array arguments 
  REAL :: paer(kdlon,kflev,5+kewaer), pcfcabs(4), pcldlw(kdlon,kflev), &
&      pdp(kdlon,kflev), pfluc(kdlon,2,kflev+1), pflux(kdlon,2,kflev+1), &
&      ppmb(kdlon,kflev+1), pqof(kdlon,kflev), ptave(kdlon,kflev), &
&      ptl(kdlon,kflev+1), pwv(kdlon,kflev)
  INTEGER :: kaerh(kdlon,kflev)

  !  Local arrays: 
  REAL :: zabcu(kdlon,nua,3*kflev+1), zb(kdlon,nsint,kflev+1), &
&      zbint(kdlon,kflev+1), zbsui(kdlon), zbtop(kdlon,nsint), &
&      zcntrb(kdlon,kflev+1,kflev+1), zdbsl(kdlon,nsint,kflev*2), &
&      zga(kdlon,8,2,kflev), zgasur(kdlon,8,2), zgatop(kdlon,8,2), &
&      zgb(kdlon,8,2,kflev), zgbsur(kdlon,8,2), zgbtop(kdlon,8,2), &
&      zgc(kdlon,5,2,kflev), zgcsur(kdlon,5,2), zgctop(kdlon,5,2), &
&      zgd(kdlon,5,2,kflev), zgdsur(kdlon,5,2), zgdtop(kdlon,5,2)

  !  External subroutines 
  EXTERNAL lwb, lwc, lwu, lwv


  !  Executable statements 

!-- 1. Initialization

!-- 1.1 Computes absorber amounts

  CALL lwu(kdlon,kflev,kewaer,kaerh,paer,pcco2,pdp,ppmb,pqof,ptave,pwv,zabcu, &
&      pcfcabs)

!-- 2. Computes planck functions

  CALL lwb(kdlon,kflev,ptave,ptl,zb,zbint,zbsui,zbtop,zdbsl,zga,zgb,zgasur, &
&      zgbsur,zgatop,zgbtop,zgc,zgd,zgcsur,zgdsur,zgctop,zgdtop)

!-- 3. Performs the vertical integration

  CALL lwv(kdlon,kflev,nuaer,ntraer,kaer,kcfc,zabcu,zbint,zbtop,zdbsl,zga, &
&      zgb,zgasur,zgbsur,zgatop,zgbtop,zgc,zgd,zgcsur,zgdsur,zgctop,zgdtop, &
&      zcntrb,pfluc)

!-- 4. Introduces the effects of clouds

  CALL lwc(kdlon,kflev,zbint,zbsui,pcldlw,zcntrb,pfluc,pflux)

  RETURN
END SUBROUTINE lw
