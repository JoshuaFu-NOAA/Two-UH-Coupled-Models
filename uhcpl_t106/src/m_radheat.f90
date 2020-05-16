!+ computes temperature changes due to radiation.

!MODULE m_radheat

!CONTAINS

SUBROUTINE radheat(kidia,kfdia,klon,klp2,klof,ktdia,klev,klevp1,paphm1,papm1,              &
                   pemtefm,pemterm,pqm1,ptm1,ptrsofm,ptrsolm,paclcvm,palbedom,palbm,       &
                   psclf0m,psclfsm,psrad0m,psrad0um,psradsm,psradsum,psraf0m,psrafsm,      &
                   ptclf0m,ptclfsm,ptrad0m,ptradsm,ptradsum,ptraf0m,ptrafsm,ptsm1m,pemtef, &
                   pemter,ptrsof,ptrsol,paclcv,palbedo,palb,psclf0,psclfs,psrad0,psrad0u,  &
                   psrads,psradsu,psraf0,psrafs,psrfl,ptclf0,ptclfs,ptrad0,ptrads,ptradsu, &
                   ptraf0,ptrafs,pte)

  ! This dummy argument is not used in the subprogram.(name:ktdia)

  ! Description:
  !
  ! Computes temperature changes due to radiation.
  !
  ! Method:
  !
  ! This routine computes the tendencies of the atmosphere's
  ! temperature due to the effects of long wave and short wave
  ! radiation. The computation is done on the t-1 time level using
  ! values of atmospheric transmisivities and emissivities that have
  ! been stored at the last full radiation time step. The surface
  ! solar flux later to be used in the soil process calculations is
  ! also stored.
  !
  ! *radheat* is called from *physc*.
  !
  ! A call to subroutine *solang* gives fields of solar zenith
  ! angles and relative day length from which an effective solar
  ! influx is computed. The results are of course different depending
  ! on the switch on or off of the diurnal cycle. Product of solar
  ! influx by transmissivities leads to solar fluxes. Then the
  ! temperatures are interpolated/extrapolated to the layer boundaries
  ! (at the bottom one takes the surface temperature) and a product by
  ! emissivities of sigma*t**4 gives thermal fluxes. The two fluxes
  ! are added and divergences computed to give heating rates.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, June 1995, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, IfM, September 2005, ocean coupling
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_control,           ONLY: lmidatm, nrow, twodt, ltdiag, lcouple,&
                                  ldsst      !fu++
  USE mo_gaussgrid,         ONLY: budw, sqcst, twomu
  USE mo_rad_switches,      ONLY: lrad
  USE mo_physc1,            ONLY: cdisse, czen1, czen2, czen3
  USE mo_constants,         ONLY: cpd, g, solc, stbo, vtmpc2
  USE mo_diagnostics_zonal, ONLY: dsrad0z, dsradsz, dtrad0z, dtradsz
!  USE m_solang,             ONLY: solang ! module procedure
  USE mo_radint,            ONLY: cemiss
  USE mo_diag_tendency,     ONLY: pdiga
  USE mo_couple,            ONLY: trads_local, srads_local
!fu++
  USE mo_dsst,              ONLY: srads_diurnal

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN)  :: kfdia, kidia, klev, klevp1, klon, klp2, ktdia
  INTEGER, INTENT (IN)  :: klof ! longitude offset to global field

  !  Array arguments with intent(In):
  REAL, INTENT (IN)  :: paclcvm(klp2), palbedom(klp2), palbm(klp2),             &
                        paphm1(klp2,klevp1), papm1(klp2,klev), pemtefm(klp2,2), &
                        pemterm(klp2,klevp1), pqm1(klp2,klev), psclf0m(klp2),   &
                        psclfsm(klp2), psrad0m(klp2), psrad0um(klp2),           &
                        psradsm(klp2), psradsum(klp2), psraf0m(klp2),           &
                        psrafsm(klp2), ptclf0m(klp2), ptclfsm(klp2),            &
                        ptm1(klp2,klev), ptrad0m(klp2), ptradsm(klp2),          &
                        ptradsum(klp2), ptraf0m(klp2), ptrafsm(klp2),           &
                        ptrsofm(klp2,2), ptrsolm(klp2,klevp1), ptsm1m(klp2)

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: pte(klp2,klev)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: pemtef(klp2,2), pemter(klp2,klevp1), ptrsof(klp2,2),      &
                        ptrsol(klp2,klevp1), paclcv(klp2), palbedo(klp2),         &
                        palb(klp2), psclf0(klp2), psclfs(klp2), psrad0(klp2),     &
                        psrad0u(klp2), psrads(klp2), psradsu(klp2), psrafs(klp2), &
                        psrfl(klp2), ptclf0(klp2), ptclfs(klp2), ptrad0(klp2),    &
                        ptrads(klp2), ptradsu(klp2), ptraf0(klp2), ptrafs(klp2),  &
                        psraf0(klp2)

  !  Local scalars: 
  REAL :: zbud, zcons2, zcons3, zdiagt, zdtdt, zflb, zflt, zqemiss, &
          zsr0u, zsrsu, ztim1, ztim2, ztim3, ztrsu
  REAL :: zflts, zfltt, zffact ! for tendency diagnostics
  INTEGER :: irow, jk, jl, jrow
  !  Local arrays: 
  REAL :: zamu0(klp2), zflxs(klon,klevp1), zflxt(klon,klevp1), zi0(klp2), &
          zrdayl(klp2), zti(klon,klevp1)

  !  Intrinsic functions 
  INTRINSIC SUM


  !  Executable Statements 

  zdiagt = 0.5*twodt
  zcons2 = cdisse*solc
  zcons3 = g/cpd
  zqemiss = (1.-cemiss)/cemiss

  irow = nrow(1)
  jrow = nrow(2)

  IF (lrad) THEN

!-- 2. SOLAR ANGLE COMPUTATIONS  

!-- 2.1 INTRODUCE THE LATITUDE DEPENDENCY

    ztim1 = czen1*.5*twomu(irow)
    ztim2 = -czen2*sqcst(irow)

!-- 2.2 CALL TO *SOLANG*

    ztim3 = czen3*sqcst(irow)

    CALL solang(klon,klp2,klof,ztim1,ztim2,ztim3,zamu0,zrdayl)

!-- 2.3 CREATE THE SOLAR INFLUX

    DO jl = kidia, kfdia
      zi0(jl) = zcons2*zamu0(jl)*zrdayl(jl)
    END DO

!-- 3. TEMPERATURES AT LAYERS' BOUDARIES

!-- 3.1 INTERPOLATION PROPER

    DO jk = 2, klev
      DO jl = kidia, kfdia
        zti(jl,jk) = (ptm1(jl,jk-1)*papm1(jl,jk-1)*(papm1(jl,jk)-paphm1(jl,jk)) &
                   +ptm1(jl,jk)*papm1(jl,jk)*(paphm1(jl,jk)-papm1(jl,jk-1)))    &
                   /(paphm1(jl,jk)*(papm1(jl,jk)-papm1(jl,jk-1)))
      END DO

!-- 3.2 SURFACE AND TOP OF ATMOSPHERE TEMPERATURE

    END DO
    DO jl = kidia, kfdia
      zti(jl,klevp1) = ptsm1m(jl)
      zti(jl,1) = ptm1(jl,1) - papm1(jl,1)*(ptm1(jl,1)-zti(jl,2))/ &
                 (papm1(jl,1)-paphm1(jl,2))
    END DO

!-- 4. UPDATE FLUXES AND COMPUTE HEATING RATES

    DO jl = kidia, kfdia
      zflxs(jl,1) = zi0(jl)*ptrsolm(jl,1)
      IF (lmidatm) THEN
        zflxt(jl,1) = pemterm(jl,1)
      ELSE
        zflxt(jl,1) = stbo*zti(jl,1)**4*pemterm(jl,1)
      ENDIF
    END DO

    DO jk = 1, klev
      DO jl = kidia, kfdia
        zflt = zflxs(jl,jk) + zflxt(jl,jk)
        zflxs(jl,jk+1) = zi0(jl)*ptrsolm(jl,jk+1)
        IF (lmidatm) THEN
          zflxt(jl,jk+1) = pemterm(jl,jk+1)
        ELSE
          zflxt(jl,jk+1) = stbo*zti(jl,jk+1)**4*pemterm(jl,jk+1)
        ENDIF
        zflb = zflxs(jl,jk+1) + zflxt(jl,jk+1)
        zdtdt = -zcons3*(zflb-zflt)/((paphm1(jl,jk+1)-paphm1(jl,jk)) &
              *(1.+vtmpc2*pqm1(jl,jk)))
        pte(jl,jk) = pte(jl,jk) + zdtdt
        IF (ltdiag) THEN
           ! tendency diagnostics
           zflts = zflxs(jl,jk+1)-zflxs(jl,jk)
           zfltt = zflxt(jl,jk+1)-zflxt(jl,jk)
           zffact = - zcons3/((paphm1(jl,jk+1)-paphm1(jl,jk)) &
                  * (1.+vtmpc2*pqm1(jl,jk)) )
           pdiga(jl,jk,24,jrow) = pdiga(jl,jk,24,jrow) + zfltt*zffact
           pdiga(jl,jk,25,jrow) = pdiga(jl,jk,25,jrow) + zflts*zffact
        END IF
      END DO
    END DO

!-- 5. DIAGNOSTICS

!-- 5.1 TOP AND SURFACE FLUXES

    IF (lmidatm) THEN
       DO jl = kidia, kfdia
          psrad0(jl) = psrad0m(jl) + zdiagt*zflxs(jl,1)
          ptrad0(jl) = ptrad0m(jl) + zdiagt*zflxt(jl,1)
          zsr0u = zflxs(jl,1) - zi0(jl)
          psrad0u(jl) = psrad0um(jl) + zdiagt*zsr0u

          psrads(jl) = psradsm(jl) + zdiagt*zflxs(jl,klevp1)
          ptrads(jl) = ptradsm(jl) + zdiagt*zflxt(jl,klevp1)
!ik ocean coupling
          IF (lcouple) THEN
             srads_local(jl,jrow) = srads_local(jl,jrow) + zdiagt*zflxs(jl,klevp1)
             trads_local(jl,jrow) = trads_local(jl,jrow) + zdiagt*zflxt(jl,klevp1)
          END IF
!fu++ SST diurnal cycle 
          IF (ldsst) THEN
             srads_diurnal(jl,jrow) = zflxs(jl,klevp1)
          END IF

          zsrsu = -zflxs(jl,klevp1)*(1./(1.-palbedom(jl))-1.)
          ztrsu = -(stbo*zti(jl,klevp1)**4+pemterm(jl,klevp1)*zqemiss)

          psradsu(jl) = psradsum(jl) + zdiagt*zsrsu
          ptradsu(jl) = ptradsum(jl) + zdiagt*ztrsu

          psraf0(jl) = psraf0m(jl) + zdiagt*zi0(jl)*ptrsofm(jl,1)
          ptraf0(jl) = ptraf0m(jl) + zdiagt*pemtefm(jl,1)
          psrafs(jl) = psrafsm(jl) + zdiagt*zi0(jl)*ptrsofm(jl,2)
          ptrafs(jl) = ptrafsm(jl) + zdiagt*pemtefm(jl,2)

          psclf0(jl) = psrad0(jl) - psraf0(jl)
          ptclf0(jl) = ptrad0(jl) - ptraf0(jl)
          psclfs(jl) = psrads(jl) - psrafs(jl)
          ptclfs(jl) = ptrads(jl) - ptrafs(jl)

          psrfl(jl) = zflxs(jl,klevp1)
       END DO
    ELSE
! this OCL is need for fjsamp !?
!OCL SCALAR
       DO jl = kidia, kfdia
          psrad0(jl) = psrad0m(jl) + zdiagt*zflxs(jl,1)
          ptrad0(jl) = ptrad0m(jl) + zdiagt*zflxt(jl,1)
          zsr0u = zflxs(jl,1) - zi0(jl)
          psrad0u(jl) = psrad0um(jl) + zdiagt*zsr0u

          psrads(jl) = psradsm(jl) + zdiagt*zflxs(jl,klevp1)
          ptrads(jl) = ptradsm(jl) + zdiagt*zflxt(jl,klevp1)

!ik ocean coupling
          IF (lcouple) THEN
             srads_local(jl,jrow) = srads_local(jl,jrow) + zdiagt*zflxs(jl,klevp1)
             trads_local(jl,jrow) = trads_local(jl,jrow) + zdiagt*zflxt(jl,klevp1)
          END IF
!fu++ SST diurnal cycle 
          IF (ldsst) THEN
             srads_diurnal(jl,jrow) = zflxs(jl,klevp1)
          END IF

          zsrsu = -zflxs(jl,klevp1)*(1./(1.-palbedom(jl))-1.)
          ztrsu = -(1.+pemterm(jl,klevp1)*zqemiss)*stbo*zti(jl,klevp1)**4

          psradsu(jl) = psradsum(jl) + zdiagt*zsrsu
          ptradsu(jl) = ptradsum(jl) + zdiagt*ztrsu

          psraf0(jl) = psraf0m(jl) + zdiagt*zi0(jl)*ptrsofm(jl,1)
          ptraf0(jl) = ptraf0m(jl) + zdiagt*stbo*zti(jl,1)**4*pemtefm(jl,1)
          psrafs(jl) = psrafsm(jl) + zdiagt*zi0(jl)*ptrsofm(jl,2)
          ptrafs(jl) = ptrafsm(jl) + zdiagt*stbo*zti(jl,klevp1)**4*pemtefm(jl,2)

          psclf0(jl) = psrad0(jl) - psraf0(jl)
          ptclf0(jl) = ptrad0(jl) - ptraf0(jl)
          psclfs(jl) = psrads(jl) - psrafs(jl)
          ptclfs(jl) = ptrads(jl) - ptrafs(jl)

          psrfl(jl) = zflxs(jl,klevp1)
       END DO
    ENDIF

!-- 5.2 ZONAL BUDGETS

    zbud = budw(irow)
    dsrad0z(irow) = zdiagt*zbud*SUM(zflxs(1:klon,1))
    dtrad0z(irow) = zdiagt*zbud*SUM(zflxt(1:klon,1))
    dsradsz(irow) = zdiagt*zbud*SUM(zflxs(1:klon,klevp1))
    dtradsz(irow) = zdiagt*zbud*SUM(zflxt(1:klon,klevp1))

!-- 6. SAVE VARIABLES FOR NEXT TIMESTEP

    DO jk = 1, klevp1
      DO jl = kidia, kfdia
        pemter(jl,jk) = pemterm(jl,jk)
        ptrsol(jl,jk) = ptrsolm(jl,jk)
      END DO
    END DO

    DO jk = 1, 2
      DO jl = kidia, kfdia
        pemtef(jl,jk) = pemtefm(jl,jk)
        ptrsof(jl,jk) = ptrsofm(jl,jk)
      END DO
    END DO

    DO jl = kidia, kfdia
      palb(jl) = palbm(jl)
      paclcv(jl) = paclcvm(jl)
      palbedo(jl) = palbedom(jl)
    END DO

!-- 7. NECESSARY COMPUTATIONS IF SUBROUTINE IS BY-PASSED

  ELSE
    DO jl = kidia, kfdia
      psrads(jl) = psradsm(jl)
      ptrads(jl) = ptradsm(jl)
      psrad0(jl) = psrad0m(jl)
      ptrad0(jl) = ptrad0m(jl)
      psrad0u(jl) = psrad0um(jl)
      psradsu(jl) = psradsum(jl)
      ptradsu(jl) = ptradsum(jl)
      psrafs(jl) = psrafsm(jl)
      ptrafs(jl) = ptrafsm(jl)
      psraf0(jl) = psraf0m(jl)
      ptraf0(jl) = ptraf0m(jl)
      ptclf0(jl) = ptclf0m(jl)
      ptclfs(jl) = ptclfsm(jl)
      psclf0(jl) = psclf0m(jl)
      psclfs(jl) = psclfsm(jl)
      palbedo(jl) = palbedom(jl)
      palb(jl) = palbm(jl)
      psrfl(jl) = 0.
    END DO
    ptrsol(:,:) = ptrsolm(:,:)
    pemter(:,:) = pemterm(:,:)
    dsrad0z(irow) = 0.
    dtrad0z(irow) = 0.
    dsradsz(irow) = 0.
    dtradsz(irow) = 0.
  END IF

  RETURN
END SUBROUTINE radheat

!END MODULE m_radheat
