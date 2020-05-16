!+ preset some constants for radiation before first line of latitude.

SUBROUTINE prerad

  ! Description:
  !
  ! This subroutine computes some parameters for the radiation
  ! scheme which remain unchanged during scan over latitude lines.
  !
  ! Method:
  !
  ! *prerad* is called from *nnsc1*.
  !
  ! Externals:
  !   *orbit*   compute orbital parameters.
  !   *ozone*   compute parameters for ozone.
  !   *aerosol* compute parameters for aerosols.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, June 1983, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_tmp_buffer
  USE mo_parameters
  USE mo_control
  USE mo_rad_switches
  USE mo_start_dataset
  USE mo_physc1
  USE mo_constants
  USE mo_rad1
  USE mo_year

  IMPLICIT NONE

  !  Local scalars: 
  REAL :: zclock, zytime, zzytime
  INTEGER :: jk, jn
  INTEGER :: iddd, iy, nobase
  LOGICAL :: lodia
  REAL :: zvemar0, zve00, zvebas, zyearve, zvetim

  !  External subroutines 
  EXTERNAL aerosol, orbit, ozone

  !  Intrinsic functions 
  INTRINSIC MOD


  !  Executable statements 

  ! Compute orbital,ozone and aerosol parameters

  IF (ly365) THEN
    CALL cd2dy(ncbase, iddd, iy)
    nobase = iddd
  ELSE
    nobase = ncbase
  ENDIF

  IF (lamip2) THEN
    nobase = nobase + 25202
    !   DETERMINE DAY OF VERNAL EQUINOX IN MARCH 1900 
    !   ASSUMING YEARL=365.2422

    zvemar0 = 20.41-0.0078*(1900.-1987.)+0.25*MOD(1900.,4.)

    !   DETERMINE DAY IN 1900 

    zve00 = zvemar0+31.+28.-1.

    !   REMAINDER OF THE YEAR

    zvebas = yearl-zve00
  END IF

  IF (nmonth/=0) THEN
    zzytime = (nmonth*30-15)/yearl*2.*api
  END IF

  zclock = MOD((ntbase*ntimst+nstep*dtime)/dayl,1.)*2*api
  IF (lamip2) THEN
    zyearve = (nobase-1.)+(ntbase+nstep*dtime)/dayl+zvebas
    zvetim = MOD(zyearve/yearl,1.)*2.*api
  ELSE
    zytime = MOD(((nobase-1)+(ntbase*ntimst+nstep*dtime)/dayl)/yearl,1.)*2*api
  ENDIF

  IF (nmonth/=0) THEN
    zytime = zzytime
  END IF

  IF (lamip2) THEN
    CALL orbit2(zclock,zvetim,cdisse,czen1,czen2,czen3,crae)
  ELSE
    CALL orbit(zclock,zytime,cdisse,czen1,czen2,czen3,crae)
  ENDIF

  IF (MOD(nstep,nradfr) == 0) THEN
    CALL ozone(zytime,cozqc,cozqs,cozhc,cozhs)
    CALL aerosol(caesc,caess,caelc,caels,caeuc,caeus,caedc,caeds)

    ! Compute orbital parameters for full radiation time step.

    zclock = MOD((ntbase*ntimst+nstep*dtime+0.5*nradfr*dtime)/dayl,1.)*2*api
    IF (lamip2) THEN
      zyearve= (nobase-1.)+(ntbase+nstep*dtime+0.5*nradfr*dtime)/dayl+zvebas
      zvetim = MOD(zyearve/yearl,1.)*2.*api
    ELSE
      zytime = MOD(((nobase-1)+(ntbase*ntimst+nstep*dtime+ &
               0.5*nradfr*dtime)/dayl)/yearl,1.)*2*api
    ENDIF

    IF (nmonth/=0) THEN
      zytime = zzytime
    END IF

    IF (lamip2) THEN
      CALL orbit2(zclock,zvetim,cdissem,czen1m,czen2m,czen3m,crae)
    ELSE
      CALL orbit(zclock,zytime,cdissem,czen1m,czen2m,czen3m,crae)
    ENDIF
  END IF

  ! Blank diagnostics

  IF ( .NOT. lrad) RETURN

  lodia = .FALSE.
  IF (nradpfr/=0) lodia = MOD(nstep,nradpfr) == 0
  IF (lodia) THEN
    DO jn = 1, 4
      DO jk = 1, nlevp1
        diag(jk,jn) = 0.
      END DO
    END DO
    DO jn = 1, 7
      dia1(jn) = 0.
    END DO
  END IF

  RETURN
END SUBROUTINE prerad
