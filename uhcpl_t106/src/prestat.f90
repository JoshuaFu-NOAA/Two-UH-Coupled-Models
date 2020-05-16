!+ set up logical switches and blank statistics at beginning of time
!  step.
!+ $Id: prestat.f90,v 1.8 1999/07/21 12:56:22 m214089 Exp $

SUBROUTINE prestat

  ! Description:
  !
  ! Set up logical switches and blank statistics at beginning of time
  ! step.
  !
  ! Method:
  !
  ! *prestat* is called from *nnsc1*.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, June 1983, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !


  USE mo_parameters
  USE mo_control
  USE mo_tmp_buffer
  USE mo_diagnostics
  USE mo_param_switches
  USE mo_forecast_switches
  USE mo_constants
  USE mo_start_dataset
  USE mo_stat_global
  USE mo_stat_zonal

  IMPLICIT NONE

  !  Local scalars: 
  REAL :: zalf, zals, zalsmv, zalv, znorm
  INTEGER :: idiapfr

  !  Intrinsic functions 
  INTRINSIC MOD


  !  Executable statements 

!-- 1. Set up logical switches

  ldiap = .FALSE.
  ldiad = .FALSE.

  IF (ndiapfr /= 0) THEN
    idiapfr = ndiapfr
    IF (MOD(nstep,idiapfr) == 0) ldiap = .TRUE.
  END IF

  lverdia = .FALSE.
  IF (ndiadfr /= 0) THEN
    IF (MOD(nstep,ndiadfr) == 0) ldiad = .TRUE.

    IF (ndiavfr /= 0) THEN
      IF (MOD(nstep,ndiavfr) == 0) lverdia = .TRUE.
    END IF
  END IF

!-- 2. Blank statistics at beginning of time step

  IF (ldiad) THEN
    voh(1:nlev)   = 0.
    dh(1:nlev)    = 0.
    qh(1:nlev)    = 0.
    xh(1:nlev)    = 0.
    th(1:nlev)    = 0.

    delph(1:nlev) = 0.

    gvo  = 0.
    gd   = 0.
    gq   = 0.
    gx   = 0.
    gt   = 0.
    gps  = 0.
    gke  = 0.
    gpe  = 0.
    glq  = 0.
    gte  = 0.
    gtpe = 0.
  END IF

!-- 3. Prepare budgets for the physics

  IF (ldiap) THEN

    ! 3.1    Blank diagnostics.

    gts = 0.
    gtd = 0.
    gws = 0.
    gwd = 0.
    gsn = 0.

    ! 3.2    Set up some constants.

    znorm  = 1000.
    zalv   = alv*rhoh2o
    zals   = als*rhoh2o
    zalf   = alf*rhoh2o
    zalsmv = zals - zalv

    ! 3.3    Preserve diagnostics before the physics.

    IF (nstep>nstart) THEN
      dia(1)  = dsrad0
      dia(2)  = dtrad0
      dia(3)  = -dsrads
      dia(4)  = dssrad
      dia(5)  = -dtrads
      dia(6)  = dstrad
      dia(7)  = dvdis
      dia(8)  = -dvdis
      dia(9)  = -dhfs
      dia(10) = (zalv*(dsevw+dsevdw+dsevcw)+zals*dsevi) + dshfl
      dia(11) = -devap*znorm
      dia(12) = -dia(10) + dshfl
      dia(13) = -dsevi*znorm
      dia(14) = -dsevw*znorm
      dia(15) = dgwdis
      dia(16) = -dgwdis
      dia(17) = dcvfr
      dia(18) = -dcvfr
      dia(19) = -znorm*(dcvmoi+dcvgr+dcvgs)
      dia(20) = dcvmoi*znorm
      dia(21) = 0.
      dia(22) = zalv*dcvgr
      dia(23) = 0.
      dia(24) = zals*dcvgs
      dia(25) = 0.
      dia(26) = -alf*dcvms
      dia(27) = -zalv*dcver
      dia(28) = dcver*znorm
      dia(29) = -zals*dcves
      dia(30) = dcves*znorm
      dia(31) = dscvs*znorm
      dia(32) = dscvr*znorm
      dia(33) = zalv*dlsgr
      dia(34) = -dlsgr*znorm
      dia(35) = zals*dlsgs
      dia(36) = -dlsgs*znorm
      dia(37) = -zalsmv*dlsms
      dia(38) = -zalv*dlser
      dia(39) = dlser*znorm
      dia(40) = -zals*dlses
      dia(41) = dlses*znorm
      dia(42) = dslss*znorm
      dia(43) = dslsr*znorm
      dia(44) = dsdtfl
      dia(45) = ddctfl - dsdtfl
      dia(46) = dsdwfl*znorm
      dia(47) = (ddcwfl-dsdwfl)*znorm
      dia(48) = dstsml
      dia(49) = -dssnmt*znorm
      dia(50) = dssnmt*znorm
      dia(51) = -dsros*znorm
      dia(52) = -dsrod*znorm
      dia(53) = -dadcon
      dia(54) = dadcon
      dia(55) = -dsevdw*znorm
      dia(56) = dstdml

      ! 3.4  Reset tendencies

      dsrad0 = 0.0
      dtrad0 = 0.0
      dsrads = 0.0
      dtrads = 0.0
      dvdis = 0.0
      dhfs = 0.0
      devap = 0.0
      dcvfr = 0.0
      dcvqac = 0.0
      dcvmoi = 0.0
      dcvgr = 0.0
      dcvgs = 0.0
      dcvms = 0.0
      dcver = 0.0
      dcves = 0.0
      dlsgr = 0.0
      dlsgs = 0.0
      dlsms = 0.0
      dlser = 0.0
      dlses = 0.0
      dssrad = 0.0
      dstrad = 0.0
      dshfl = 0.0
      dsdtfl = 0.0
      dslsr = 0.0
      dslss = 0.0
      dscvr = 0.0
      dscvs = 0.0
      dsevw = 0.0
      dsevi = 0.0
      dsdwfl = 0.0
      dssnmt = 0.0
      ddctfl = 0.0
      ddcwfl = 0.0
      dsros = 0.0
      dsrod = 0.0
      dadcon = 0.0
      dgwdis = 0.0
      dsevdw = 0.0
      dsevcw = 0.0
      dstsml = 0.0
      dstdml = 0.0

    END IF

  END IF

END SUBROUTINE prestat
