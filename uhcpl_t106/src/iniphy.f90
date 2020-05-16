!+ initialises physical constants of uncertain value.
!+ $Id: iniphy.f90,v 1.6 2000/09/19 07:57:50 m214003 Exp $

SUBROUTINE iniphy

  ! Description:
  !
  ! Initialises physical constants of uncertain value.
  !
  ! Method:
  !
  ! This routine sets the values for the physical constants used
  ! in the parameterization routines (except for the radiation
  ! black-box) whenever these values are not well enough known to
  ! forbid any tuning or whenever they are subject to an arbitrary
  ! choice of the modeller. These constants will be in *comph2*.
  !
  ! *iniphy* is called from *physc*.
  !
  ! Authors:
  !
  ! J. F. Geleyn, ECMWF, December 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_parameters
  USE mo_control
  USE mo_constants
  USE mo_physc2
  USE mo_diagnostics
  USE mo_hyb
  USE mo_vegetation
  USE mo_cumulus_flux, ONLY: cuparam

  IMPLICIT NONE

  !  Local scalars: 
  INTEGER :: jk

  !  Intrinsic functions 
  INTRINSIC SQRT


  !  Executable statements 

!-- 1. Setting of single variables

!-- 1.1 Constants for *vdiff*

  clam = 300.
  ckap = 0.4
  cb = 5.
  cc = 5.
  cd = 5.

  IF (nn==21 .AND. .NOT. lamip2) THEN
    cchar = 0.018
  ELSE
    cchar = 0.032
  END IF

  cvdifts = 1.5

!-- 1.2 Constants for *vdiff*, *clsst* and *atmice*

  IF (lamip2) THEN
    ctfreez = 271.38
  ELSE
    ctfreez = 273.16 - 1.79
  END IF
  cz0ice = 0.001

!-- 1.3 Constant for *vdiff* and *surf*

  IF (lamip2) THEN
    cqsncr = 1./0.006
  ELSE
    cqsncr = 1./0.015
  END IF

!-- 1.4 Constant for massflux convection scheme

  CALL cuparam

  DO jk = 1, nlev
    cevapcu(jk) = 1.93E-6*261.*SQRT(1.E3/(38.3*0.293)*SQRT(ceta(jk)))*0.5/g
  END DO

!-- 1.5 Constants for *surf*

  IF (lamip2) THEN

    !  thickness of soil layers

    cdel(1) = 0.065
    cdel(2) = 0.254
    cdel(3) = 0.913
    cdel(4) = 2.902
    cdel(5) = 5.700

    !  depth of mids of soil layers

    cmid(1) = cdel(1)*0.5
    cmid(2) = cdel(1)+cdel(2)*0.5
    cmid(3) = cdel(1)+cdel(2)+cdel(3)*0.5
    cmid(4) = cdel(1)+cdel(2)+cdel(3)+cdel(4)*0.5
    cmid(5) = cdel(1)+cdel(2)+cdel(3)+cdel(4)+cdel(5)*0.5

    csncri = 0.  ! critical snow depth for soil computations
  ELSE
    cd1 = 0.07
    cd2 = 0.42
  END IF
  clice = 0.3
  cqcon = 1.E-10
  cqdif = 1.E-7

!-- 1.6 Constants for diagnostics in *surf*

  IF (lamip2) THEN
    cdiawd = 0.
  ELSE
    cdiawd = cd2/cd1
  END IF

!-- 1.7 Constants for *vdiff* and *surf*

  cgh2o  = 4.18E6
  cwlmax = 2.E-4

  cvroots = .50
  cvrootd = .50
  cvrootc = 0.
  cvrad   = 0.55
  IF (lamip2) THEN
    cvinter = 0.5
  ELSE
    cvinter = 1.
  END IF
  corvari = 50.**2
  corvars = 750.**2

  cva = 5000.
  cvb = 10.
  cvc = 100.
  cvk = .9
  cvbc  = cvb*cvc
  cvkc  = cvk*cvc
  cvabc = (cva+cvbc)/cvc

  RETURN
END SUBROUTINE iniphy
