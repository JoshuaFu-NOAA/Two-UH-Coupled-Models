MODULE mo_cumulus_flux

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module *mo_cumulus_flux* - parameters for cumulus massflux scheme
  !
  ! ----------------------------------------------------------------

  REAL :: entrpen      !    entrainment rate for penetrative convection
  REAL :: entrscv      !    entrainment rate for shallow convection
  REAL :: entrmid      !    entrainment rate for midlevel convection
  REAL :: entrdd       !    entrainment rate for cumulus downdrafts
  REAL :: cmfctop      !    relat. cloud massflux at level above nonbuoyanc
  REAL :: cmfcmax      !    maximum massflux value allowed for
  REAL :: cmfcmin      !    minimum massflux value (for safety)
  REAL :: cmfdeps      !    fractional massflux for downdrafts at lfs
  REAL :: rhcdd        !    relative saturation in downdrafts
  REAL :: cprcon       !    coefficients for determining conversion
                       !    from cloud water to rain
  LOGICAL :: lmfpen    !    true if penetrative convection is switched on
  LOGICAL :: lmfscv    !    true if shallow     convection is switched on
  LOGICAL :: lmfmid    !    true if midlevel    convection is switched on
  LOGICAL :: lmfdd     !    true if cumulus downdraft      is switched on
  LOGICAL :: lmfdudv   !    true if cumulus friction       is switched on

CONTAINS

SUBROUTINE cuparam

  ! Description:
  !
  ! Defines disposable parameters for massflux scheme
  !
  ! Method:
  !
  ! This routine is called from *iniphy*
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, February 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A. Rhodin, MPI, Jan 1999, subroutine cuparam -> module mo_cumulus_flux
  ! 
  ! for more details see file AUTHORS
  ! 

  USE mo_control, ONLY: lamip2

  IMPLICIT NONE


  !  Executable Statements 

!-- 1. Specify parameters for massflux-scheme

!fu++  entrpen =1.0E-4 ! Average entrainment rate for penetrative convection
       entrpen =1.0E-4 ! Average entrainment rate for penetrative convection

!fu++  entrscv =3.0E-4 ! Average entrainment rate for shallow convection
!       entrscv =5.0E-4 ! Average entrainment rate for shallow convection
!fu++ test on Sep/16/2008
       entrscv =1.0E-3 ! enhanced shallow convection entrainment
!fu++ test on Aug/17/2010
!	entrscv =0.0    !turn off
!
!fu++
  entrmid = 1.0E-4 ! Average entrainment rate for midlevel convection

  entrdd  = 2.0E-4 ! Average entrainment rate for downdrafts

  IF (lamip2) THEN
    cmfctop = 0.33 ! Relative cloud massflux at level above nonbuoyancy level
  ELSE
    cmfctop = 0.1
  END IF

  cmfcmax = 1.0    ! Maximum massflux value allowed for updrafts etc

  cmfcmin = 1.E-10 ! Minimum massflux value (for safety)

  cmfdeps = 0.3    ! Fractional massflux for downdrafts at lfs

  cprcon  = 6.E-4  ! Coefficients for determining conversion from cloud water

  ! Next value is relative saturation in downdrafrs
  ! but is no longer used ( formulation implies saturation)

  rhcdd = 1.

  RETURN
END SUBROUTINE cuparam

END MODULE mo_cumulus_flux
