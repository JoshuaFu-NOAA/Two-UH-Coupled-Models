!+ calculate seaice skin-temperature

MODULE mo_skintem

CONTAINS

!OCL NOALIAS

SUBROUTINE skintem (klp2, klevp1, klon,                                       &
                    albedom, emterm, seaice, siced, slmm, teff, teffm, tsm1m, &
                    dhft, srfl, thfl, tskin, tskinm, tskin1, tskin1m)

  !!! This dummy argument is not used in the subprogram.(name:tskin1m) !!!

  ! Description:
  !
  ! Calculate ice skin-temperature as a prognostic variable
  !
  ! Method:
  !
  ! *skintem* called from physc
  ! *physc* called gpc
  !
  ! Surface fluxes are linearised around surface-temperature.
  !
  !     q(tskin)=q(tsm1m)+dqdt*(tskin-tsm1m)
  !
  ! The heat-balance equation: q(tskin)*dt=cp*(tskin-tskinm)
  ! leeds to:
  !
  !     qc(tsm1m,tskinm)=tskin*qd(tsm1m)  ==> tskin=qc/qd
  !     qc=constant flux-terms only dependet on tsm1m and tskinm
  !     qd=dirivations of q
  !
  ! Authors:
  !
  ! F. Lunkeit, MI, April 1991, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A, Rhodin, MPI, Jan 1999, argument list added
  !
  ! for more details see file AUTHORS
  !

  USE mo_control,        ONLY: lmidatm, twodt, lamip2
  USE mo_param_switches, ONLY: lice
  USE mo_constants,      ONLY: stbo, tmelt
  USE mo_physc2,         ONLY: ctfreez

  USE mo_radint,         ONLY: cemiss

  IMPLICIT NONE

  ! Arguments

  INTEGER ,INTENT(in)  :: klp2   ! first  (longitudinal) bound of arrays
  INTEGER ,INTENT(in)  :: klevp1 ! second (vertical)     bound of emterm
  INTEGER ,INTENT(in)  :: klon   ! first  (longitudinal) loop index bound

  REAL    ,INTENT(in)  :: albedom (klp2)
  REAL    ,INTENT(in)  :: emterm  (klp2,klevp1)
  REAL    ,INTENT(in)  :: seaice  (klp2)
  REAL    ,INTENT(in)  :: siced   (klp2)
  REAL    ,INTENT(in)  :: slmm    (klp2)
  REAL    ,INTENT(out) :: teff    (klp2)
  REAL    ,INTENT(in)  :: teffm   (klp2)
  REAL    ,INTENT(in)  :: tsm1m   (klp2)

  REAL    ,INTENT(in)  :: dhft    (klp2)
  REAL    ,INTENT(in)  :: srfl    (klp2)
  REAL    ,INTENT(in)  :: thfl    (klp2)

  REAL    ,INTENT(out) :: tskin   (klp2) ! auxil1
  REAL    ,INTENT(in)  :: tskinm  (klp2) ! auxil1m
  REAL    ,INTENT(out) :: tskin1  (klp2) ! auxil2
  REAL    ,INTENT(in)  :: tskin1m (klp2) ! auxil2m

  !  Local scalars: 
  REAL :: zalb, zalbmax, zalbmin, zalpha, zcpcon, zcpdt, zcpice, zdalb, &
&      zdiagt, zdice, zdicefl, zdqice, zdqidt, zdsflx, zdthfl, zdtrfl, &
&      zicefl, zrhoice, zsflx, zsig, zsofl, ztalb, zthfl, ztmelt, ztmst, &
&      ztrfl, ztsea, ztsnmelt
  INTEGER :: jl

  !  Local arrays: 
  REAL :: zemterm(klon)

  !  Intrinsic functions 
  INTRINSIC MIN


  !  Executable statements 

!-- 1. Set up constants

  zalbmin  = 0.5
  zalbmax  = 0.75

  IF (lamip2) THEN
    ztmelt = ctfreez
    ztalb  = ctfreez - 10.
  ELSE
    ztmelt = tmelt
    ztalb  = ztmelt - 10.
  END IF

  ztsnmelt = ztmelt
  zalpha   = 2.0
  zcpice   = 2090.
  zrhoice  = 1000.
  zsig     = stbo
  ztmst    = twodt
  zdiagt   = 0.5*twodt
  zdice    = 0.10

  zcpcon = zrhoice*zcpice*zdice
  zcpdt  = zcpcon/zdiagt

!-- 2. Compute new skin-temperature

  IF (lice) THEN

    zemterm(:) = emterm(:klon,klevp1)

    DO jl = 1, klon

      IF (slmm(jl)<0.5 .AND. siced(jl)>zdice) THEN

        ztsea = ctfreez

        ! Calculate ice albedo

        IF (tskinm(jl)>=ztsnmelt) THEN
          zalb = zalbmin
        ELSE IF (tskinm(jl)<ztalb) THEN
          zalb = zalbmax
        ELSE
          zdalb = (zalbmax-zalbmin)/(ztsnmelt-ztalb)
          zalb  = zalbmin + zdalb*(ztsnmelt-tskinm(jl))
        END IF

        ! Constant flux terms

        IF (lmidatm) THEN
          ztrfl = zemterm(jl) + 4.*cemiss*zsig*tsm1m(jl)**4
        ELSE
          ztrfl = zemterm(jl)*zsig*tsm1m(jl)**4 + 4.*cemiss*zsig*tsm1m(jl)**4
        ENDIF
        zsofl  = (1.-zalb)*srfl(jl)/(1.-albedom(jl))
        zthfl  = thfl(jl) - dhft(jl)*tsm1m(jl)
        zicefl = zalpha*ztsea/siced(jl)
        zdqice = zcpdt*tskinm(jl)
        zsflx  = -ztrfl - zthfl - zsofl - zicefl - zdqice

        ! Dirivations

        zdtrfl  = -4.*cemiss*zsig*tsm1m(jl)**3
        zdthfl  = dhft(jl)
        zdicefl = -zalpha/siced(jl)
        zdqidt  = -zcpdt
        zdsflx  = zdtrfl + zdthfl + zdicefl + zdqidt

        tskin(jl) = zsflx/zdsflx

        tskin(jl) = MIN(tskin(jl),ztsea)

      ELSE
        ztsea = tsm1m(jl)
        tskin(jl) = tsm1m(jl)
      END IF

      ! Accumulate skin-temperature

      tskin1(jl) = seaice(jl)*tskin(jl) + (1.-seaice(jl))*ztsea
      teff(jl)   = teffm(jl) + zdiagt*tskin1(jl)

    END DO

    ! Necessary computations if subroutine is bypassed
  ELSE
    DO jl = 1, klon
      tskin1(jl) = tsm1m(jl)
      teff(jl)   = teffm(jl) + zdiagt*tskin1(jl)
    END DO
  END IF

  RETURN
END SUBROUTINE skintem

END MODULE mo_skintem
