!+ updates land values of temperature, moisture and snow.

!MODULE m_surf

!CONTAINS

  SUBROUTINE surf (                    &
       klon, klp2, klevp1,             &
       wsmx , wsmxm,                   &
       rgcgn, rgcgnm,                  &
       sodif, sodifm,                  &
       vlt, vltm,                      &
       loland, loglac, slmm ,          &! land-sea mask, glacier mask
       ts, tsm, tsm1, tsm1m,           &! surface temperature
       ws, wsm, wsm1, wsm1m,           &! soil water amount
       sn, snm, snm1, snm1m,           &! depth of the snow pack
       wl, wlm, wlm1, wlm1m,           &! skin reservoir water content
       td, tdm, tdm1, tdm1m,           &! temperature, layer 4
       drain, drainm,                  &! drainage
       srfl, thfl, qhfl, xhfl,         &! -- surface fluxes  --
       rsfc, ssfc, rsfl, ssfl,         &
       ahfli,                          &
       dhft,                           &! sens. heat flux with resp. to surf. temp.
       dhfqw,                          &! moist. flux with respect to skin reservoir
       dhfqs,                          &! moist. flux over snow with resp. to snow d
       evap, evapm,                    &! -- evaporation ---
       vgrat, cvs, wlmx,               &! -- vegetation  --
       emterm,                         &! -- emissivities ---
       dsnac, dsnacm,                  &! -- remaining elements in 'surf' --
       runoff, runoffm,                &
       snmel, snmelm,                  &
       tdcl , tdclm , tdclm1, tdclm1m, &! temperature, layer 5
       td3  , td3m  , td3m1 , td3m1m , &! temperature, layer 1
       td4  , td4m  , td4m1 , td4m1m , &! temperature, layer 2
       td5  , td5m  , td5m1 , td5m1m , &! temperature, layer 3
       tslin, tslinm,                  &
       tsmax , tsmaxm,                 &
       tsmin, tsminm,                  &
       tsn  , tsnm  , tsnm1 , tsnm1m , &! temperature in the middle of the snow pack
       tsurf, tsurfm, varor            )

    ! Description:
    !
    ! Updates land values of temperature, moisture and snow.
    !
    ! Method:
    !
    ! This routine updates the land values of surface temperature,
    ! deep temperature, skin soil water,
    ! surface soil moisture (expressed as a water
    ! content), deep soil moisture (in comparable numbers) and snow
    ! depth (in water equivalent). This is done via a forward time step
    ! damped with some implicit linear considerations: As if all fluxes
    ! that explicitely depend on the variable had only a linear
    ! variation around the t-1 value. For consistency with the
    ! atmospheric treatment a time filter is applied on all five surface
    ! pronostic variables. However only the forward half of the time
    ! filter is used to avoid keeping indefinitely traces of snow for
    ! example. Climatic temperature and moisture, in a third, deeper
    ! layer, are used as lower boundary conditions.
    !
    ! Straightforward once the definition of the constants is
    ! understood. For this refer to documentation. For the time filter
    ! see corresponding part of the documentation of the adiabatic code.
    !
    ! The routine takes its input from the long-term storage:
    ! ts,td,wl,ws,wd,sn at t-1,
    ! Surface fluxes computed in other parts of
    ! the physic (the long-wave radiative flux is recomputed here from
    ! the emissivity), climatic t and w and land-sea mask. It returns
    ! its output to the same space: Same variables at t+1 and filtered
    ! values of the same variables at t.
    !
    ! *surf* is called from *physc*.
    !
    ! Reference:
    ! See soil processes' part of the model's documentation for
    ! details about the mathematics of this routine.
    !
    ! Modification:
    ! 1. Scheme modified for climate studies by using a five layer
    !    scheme as proposed by warrilow et al.(1986)
    !    soil temperatures are td3,td4,td5,td,tdcl,
    ! 2. Extra temperature variable for snow layer tsn if snow is
    !    deeper than zsncri
    ! 3. Soil hydrology scheme based on catchment considerations
    !
    ! Authors:
    !
    ! J. F. Geleyn, ECMWF, June 1982, original source
    ! C. B. Blondin, ECMWF, December 1986, changed
    ! L. Dumenil, MPI, May 1988, changed
    ! J.-P. Schulz, MPI, 1997, implementation of implicit coupling between 
    !                          land surface and atmosphere
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A. Rhodin, MPI, January 1999, parameters passed by argument list
    ! 
    ! for more details see file AUTHORS
    !


    !-------------
    ! Modules used
    !-------------
    USE mo_control,        ONLY: &
         eps     ,&! time filtering coefficient.
         lmidatm ,&! .true. for middle atmosphere model version
         nrow    ,&! urrent latitude line. (one entry per task).
         twodt   ,&! 2.*dtime.
         dtime   ,&! timestep in seconds
         nlev    ,&
         lamip2
    USE mo_gaussgrid,      ONLY: &
         aslm    ,&! number of land points on each latitude line.
         budw      ! weights for global budgets *budw=gw/klon*
    USE mo_param_switches, ONLY: &
         lsurf     ! true for surface exchanges
    USE mo_physc2,         ONLY: &
         cd1     ,&! thickness of 1st ground level
         cd2     ,&! thickness of 2nd ground level.
         cgh2o   ,&! heat capacity of water
         clice   ,&! thermal conductivity factor of snow/ice
         cqcon   ,&! soil hydraulic conductivity
         cqdif   ,&! soil hydraulic diffusivity
         cvdifts ,&! factor for timestep weighting in *rhs* of *vdiff* and *scv*.
         cwlmax  ,&! maximum moisture content of the skin reservoir
         cqsncr  ,&! inverse of equivalent water height when snow
         ! is considered to cover completely the ground
         ! in the box.
    csncri  ,&! critical snow depth for soil computations ADDED for IMPL
         cdel      ! thickness of soil layers
    USE mo_constants,      ONLY: &
         alf     ,&! latent heat for fusion.
         rhoh2o  ,&! density of liquid water.
         stbo    ,&! stephan boltzmann constant.
         tmelt     ! temperature of fusion of ice.
    USE mo_start_dataset,  ONLY: &
         nstart  ,&! time step for start/restart.
         nstep     ! current time step.
    USE mo_vegetation,     ONLY: &
         cvinter ,&! efficency of interception of precipitation as rain.
         cvrootc ,&! percentage of roots in the 3rd soil layer
         cvrootd ,&! percentage of roots in the 2nd soil layer
         cvroots   ! percentage of roots in the 1st soil layer
    USE mo_diagnostics_zonal, ONLY: & 
         dscvrz  ,&! zonal means of quantities related
         dscvsz  ,&!   to global physical diagnostics
         dsdtflz ,&!
         dseviz  ,&!
         dsevwz  ,&!
         dshflz  ,&!
         dslsrz  ,&!
         dslssz  ,&!
         dsrosz  ,&!
         dssnmtz ,&!
         dssradz ,&!
         dstradz   !
    USE mo_soil_impl
    USE mo_radint, ONLY: &
         cemiss    ! surface longwave emissivity

    !--------------------------------------------
    !  Variables now accessed via parameter list:
    !--------------------------------------------
    USE mo_decomposition, ONLY: dc => local_decomposition

    IMPLICIT NONE
    !-----------
    ! Parameters
    !-----------
    INTEGER ,INTENT(in) ::               &
         klon, klp2, klevp1
    REAL ,INTENT(out) ::                 &
         wsmx (klp2),                       &! field capacity   (just copy wsmxm)
         rgcgn(klp2),                       &! heat capacity    (just copy rgcgnm)
         sodif(klp2),                       &! soil diffusivity (just copy sodifm)
         drain(klp2),                       &! drainage
         vlt  (klp2)                         ! leaf area index
    REAL ,INTENT(in) ::                  &
         wsmxm (klp2),                      &! field capacity of soil
         rgcgnm(klp2),                      &! heat capacity of soil
         sodifm(klp2),                      &! soil diffusivity
         drainm(klp2),                      &! drainage
         vltm  (klp2),                      &! leaf area index
         slmm (klp2)                         ! some weight for diagnostics ?
    LOGICAL ,INTENT(in) ::               &
         loland(klp2), loglac(klp2)          ! land-sea mask, glacier mask
    REAL ,INTENT(out) ::                 &!==== surf temp and moist fields ====
         ts(klp2), tsm1(klp2),              &! surface temperature
         ws(klp2), wsm1(klp2),              &! soil wetness
         sn(klp2), snm1(klp2),              &! depth of the snow pack
         wl(klp2), wlm1(klp2),              &! skin reservoir content
         td(klp2), tdm1(klp2)                ! soil temperature layer 4
    REAL ,INTENT(in) ::                  &
         tsm(klp2), tsm1m(klp2),            &! surface temperature
         wsm(klp2), wsm1m(klp2),            &! soil wetness
         snm(klp2), snm1m(klp2),            &! depth of the snow pack
         wlm(klp2), wlm1m(klp2),            &! skin reservoir content
         tdm(klp2), tdm1m(klp2)              ! soil temperature layer 4
    REAL ,INTENT(in) ::                  &!======== surface fluxes ========
         srfl(klp2),                        &
         thfl(klp2),                        &
         xhfl(klp2),                        &
         ssfc(klp2),                        &! surf.flux convective snow
         ssfl(klp2),                        &! surf.flux largescale snow
         ahfli(klp2)                         ! sea ice heat flux
    REAL ,INTENT(inout) ::               &
         qhfl(klp2),                        &! 
         rsfl(klp2),                        &! largescale rain water
         rsfc(klp2)                          ! convective rain water
    REAL ,INTENT(in) ::                  &!======== flux derivatives ========
         dhft(klp2),                        &
         dhfqs(klp2)
    REAL ,INTENT(in) ,TARGET ::          &
         dhfqw(klp2)
    REAL ,INTENT(in) ::                  &
         evap(klp2), evapm(klp2),           &! evaporation
         vgrat(klp2),                       &! vegetation ratio 
         cvs(klp2),                         &!
         wlmx(klp2),                        &!
         emterm(klp2,klevp1)                         ! emissivities
    REAL, INTENT (in) ::                 &!==== remaining elements in 'surf' ====
         dsnacm(klp2),                      &
         runoffm(klp2),                     &! runoff accumulated
         snmelm(klp2),                      &! snow melt
         tdclm (klp2), tdclm1m(klp2),       &! temperature, layer 5
         td3m  (klp2), td3m1m (klp2),       &! temperature, layer 1
         td4m  (klp2), td4m1m (klp2),       &! temperature, layer 2
         td5m  (klp2), td5m1m (klp2),       &! temperature, layer 3
         tslinm(klp2),                      &! land: residual head budget
         ! sea ice: conductive heat flux
    tsmaxm(klp2),                      &! maximum surface temperature
         tsminm(klp2),                      &! minimum surface temperature
         tsnm  (klp2), tsnm1m (klp2),       &! temperature in the snow pack
         tsurfm(klp2),                      &! surface temperature (mean)
         varor (klp2)                        ! orographic variance (for surf.runoff)
    REAL ,INTENT(out) ::                 &
         dsnac(klp2),                       &
         runoff(klp2),                      &! runoff accumulated
         snmel(klp2),                       &! snow melt
         tdcl (klp2), tdclm1(klp2),         &! temperature, layer 5
         td3  (klp2), td3m1 (klp2),         &! temperature, layer 1
         td4  (klp2), td4m1 (klp2),         &! temperature, layer 2
         td5  (klp2), td5m1 (klp2),         &! temperature, layer 3
         tslin(klp2),                       &
         tsmax(klp2),                       &
         tsmin(klp2),                       &
         tsn  (klp2), tsnm1 (klp2),         &! temperature in the snow pack
         tsurf(klp2)

    !  Local scalars: 

    REAL :: zb1, zbase, zbm, zbud, zbws, zconb2, zcons11, zcons12, zcons14, &
         &      zcons15, zcons16, zcons17, zcons18, zcons19, zcons20, zcons21, zcons5, &
         &      zconss10, zconss9, zconw1, zconw2, zconw3, zconw4, zcvseps, zd1, zd2, &
         &      zdeff, zdhfqw, zdiags, zdiagt, zdiagw, zdifiz, zdisc, zdrexp, zdrmax, &
         &      zdrmin, zepext, zepprcp, zfac, zgh2o, zgice, ziprcp, zlice, &
         &      zlyeps, zlysic, zmprcp, zorvari, zorvars, zpmax, zpsfr, zqcon, zqdif, &
         &      zqhflw, zqsncr, zrici, zroeff, zsa, zsb, zsdtfl, zshfl, zsig0, zsig1, &
         &      zsig10, zsig11, zsig2, zsig3, zsig4, zsig5, zsig6, zsn, zsncoe, zsncri, &
         &      zsnlam, zsnmlt, zsros, zssnmt, zssrad, zstrad, ztmst, ztpfac1, &
         &      ztpfac2, ztpfac3, ztprcp, zvinter, zvrootc, zvrootd, zvroots, zwdtr, &
         &      zwlmax, zwmax, zwptr, zwslim, zwsup, zzdp, zzqsnfl

    INTEGER :: nglat     ! Number of gaussian latitudes.
    INTEGER :: irow, jrow, ivpos, jk, jkk, jl, nsl, jkg
    LOGICAL :: lo

    !  Local arrays: 

    REAL :: tem1(klp2), tem2(klp2), zair(klp2), zcfh(klp2,5), zconb1(klp2), &
         &      zdair(klp2), zdel(5), zdifi(klp2), zdrain(klp2), zdsn(klp2), &
         &      zeb(klp2,5), zemterm(klp2), zinfil(klp2), zlac(klp2), zper(klp2), &
         &      zpsfl(klp2), zqfld(klp2), zros(klp2), zsdif(klp2,5), zsndp(klp2), &
         &      zsnfl(klp2), zsnmel(klp2), zsufl(klp2), zt(klp2), ztcoe(klp2), &
         &      ztsdfl(klp2), zvol(klp2), zwl(klp2)

    !  Implicit soil scheme 
    REAL :: zcapav(klon)     ! ground volumetric heat capacity [J/m**3/K]
    REAL :: ztsvar(klon)     ! variance of surface temperature

    REAL, TARGET :: ztrfl(klp2 )
    REAL, POINTER :: zprfl(:), zqrfld(:)

#ifdef DRAIN
    ! Pointer for additional G3X fields
    REAL, POINTER :: drain(:), drainm(:)
#endif

    !  Intrinsic functions 
    INTRINSIC ABS, DOT_PRODUCT, MAX, MERGE, MIN, SQRT


    !  Executable statements 

    !  Local array bounds
    nglat = dc%nglat

    ! Physical constants.

    ! *zrgcg* is the product of density by heat capacity fot the
    ! soil, *zdif* is the thermal diffusivity of the soil, *zd1* and
    ! *zd2* are the thicknesses of the two ground layers (s,d) and *ztq*
    ! is the ratio of diffusivities for heat and moisture. *zqsncr* is
    ! is the inverse of a critical value for snow depth (see *vdiff*).

    zsncoe = 634500.
    ! J/(m**3 K)

    zsnlam = 0.22
    ! W/m/K

    IF (lamip2) THEN
       zsncri = csncri  ! m water equivalent critical snow height
    ELSE
       zsncri = 0.025
    END IF

    zrici  = 2.09E+06
    zdifiz = 12.E-07
    IF (.NOT. lamip2) THEN
       zd1 = cd1
       zd2 = cd2
    END IF
    zqcon   = cqcon
    zqdif   = cqdif
    zqsncr  = cqsncr
    zvroots = cvroots
    zvrootd = cvrootd
    zvrootc = cvrootc
    zgh2o   = cgh2o
    zlice   = clice
    zgice   = 0.5*zgh2o
    zwlmax  = cwlmax
    zvinter = cvinter
    zorvari = 100.
    IF (dc%nlat>=64) THEN
       zorvars = 1000.
    ELSE
       zorvars = 1500.
    END IF

    ! Security parameters

    zepext  = 1.E-20
    zcvseps = 1.E-20
    zepprcp = 1.E-20

    ! Computational constants.

    ztmst = twodt
    IF (nstep==nstart) ztmst = 0.5*twodt
    zdiagt = 0.5*twodt
    zdiagw = zdiagt/rhoh2o
    zdiags = zdiagt/ztmst

    zcons5 = ztmst/rhoh2o
    IF (.NOT. lamip2) THEN
       zcons11 = zd1/zd2
       zcons14 = zd2/zd1
       zcons15 = 1./(1.+zcons14)
       zcons16 = 1./(0.5*(zd1+zd2))
       zcons17 = 1./(0.5*zd2*(zd1+zd2))
       zcons18 = 1./(0.5*zd1*(zd1+zd2))
       zcons19 = 1./(0.5*zd1*zd1)
       zcons20 = 1./(zd2*zd2)
       zcons21 = 1./zd2
    END IF

    zpsfr  = 1.
    zdrmin = 0.001/(3600.*1000.)
    zdrmax = 0.1/(3600.*1000.)
    zdrexp = 1.5
    zpmax  = 1.

    zbase  = 0.

    irow = nrow(1)
    jrow = nrow(2)

    DO jl = 1, klon
       wsmx(jl)  = wsmxm(jl)
       rgcgn(jl) = rgcgnm(jl)
       sodif(jl) = sodifm(jl)
       vlt(jl)   = vltm(jl)
    END DO

    zprfl  => ztrfl

    zqrfld => dhfqw

    IF (lsurf) THEN

       !-- 1. If no land point by-pass computations

       IF (aslm(irow)<=0.5) THEN
          DO jl = 1, klon
             ts(jl)     = tsm1m(jl)
             td(jl)     = tdm1m(jl)
             zwl(jl)    = 0.
             ws(jl)     = 0.
             sn(jl)     = 0.
             tsm1(jl)   = tsm(jl)
             tdm1(jl)   = tdm(jl)
             wlm1(jl)   = 0.
             wsm1(jl)   = 0.
             snm1(jl)   = 0.
             tsn(jl)    = tsnm1m(jl)
             td3(jl)    = td3m1m(jl)
             td4(jl)    = td4m1m(jl)
             td5(jl)    = td5m1m(jl)
             tdcl(jl)   = tdclm1m(jl)
             tsnm1(jl)  = tsnm(jl)
             td3m1(jl)  = td3m(jl)
             td4m1(jl)  = td4m(jl)
             td5m1(jl)  = td5m(jl)
             tdclm1(jl) = tdclm(jl)
             tsurf(jl)  = tsurfm(jl) + zdiagt*ts(jl)
             !tsvari(jl) = tsvarim(jl)
             runoff(jl) = 0.
             drain(jl) = 0.
             snmel(jl) = 0.
             dsnac(jl) = 0.
             tslin(jl) = tslinm(jl)
          END DO

          dssradz(irow) = 0.
          dstradz(irow) = 0.
          dshflz(irow)  = 0.

          dsdtflz(irow) = 0.
          dslsrz(irow)  = 0.
          dslssz(irow)  = 0.
          dscvrz(irow)  = 0.
          dscvsz(irow)  = 0.
          dsevwz(irow)  = 0.
          dseviz(irow)  = 0.

          dssnmtz(irow) = 0.

          dsrosz(irow)  = 0.

       ELSE

          !-- 2.     Surface temperature

          !-- 2.1.1. Input for energy balance, snow temperature, snow melt

          zemterm(:) = emterm(:,klevp1)
          DO jl = 1, klon
             IF (lmidatm) THEN
                ztrfl(jl) = zemterm(jl)
             ELSE
                ztrfl(jl) = (zemterm(jl)*stbo*tsm1m(jl)**3)*tsm1m(jl)
             ENDIF
             zair(jl) = -(srfl(jl)+ztrfl(jl)+thfl(jl))
             zdair(jl) = 4.*cemiss*stbo*tsm1m(jl)**3 - dhft(jl)
          END DO
          DO jl = 1, klon
             IF (loglac(jl)) THEN
                zzqsnfl  = (evap(jl)-evapm(jl))/zdiagw
                zlac(jl) = rsfl(jl) + rsfc(jl) + ssfl(jl) + ssfc(jl) + zzqsnfl
             ELSE
                zlac(jl) = 0.
             END IF
          END DO

          zbud   = budw(irow)
          zssrad = zdiagt*zbud*DOT_PRODUCT(slmm(1:klon),srfl(1:klon))
          zstrad = zdiagt*zbud*DOT_PRODUCT(slmm(1:klon),ztrfl(1:klon))
          zshfl  = zdiagt*zbud*DOT_PRODUCT(slmm(1:klon),thfl(1:klon))

          IF (lamip2) THEN
             ! Time filter for surface temperature ts

             IF (nstep /= nstart) THEN
                DO jl = 1,klon
                   tsm1(jl) = tsm(jl)+eps*(tsm1m(jl)-2.*tsm(jl)+ts(jl))
                END DO
             ELSE
                DO jl = 1,klon
                   tsm1(jl) = tsm(jl)
                END DO
             END IF
          END IF

!DIR$ IVDEP
!OCL NOVREC
          DO jl = 1, klon
             zsnmel(jl) = 0.
             IF (.NOT. lamip2) THEN
                IF ((snm1m(jl)>0.) .AND. (loland(jl)) .AND. ( .NOT. loglac(jl))) THEN

                   zsndp(jl) = (1000./300.)*snm1m(jl)

                   IF ((snm1m(jl)<zsncri) .AND. (loland(jl))) THEN
                      IF (ABS(tsnm1m(jl)-td3m1m(jl))>0.2) THEN
                         tem1(jl) = -(zair(jl)+zsncoe*zsndp(jl)*(td3m1m(jl)-tsnm1m(jl))/ &
                              &                    ztmst)
                         tem2(jl) = -(zdair(jl)+zsncoe*zsndp(jl)/ztmst)
                      ELSE
                         tem1(jl) = -zair(jl)
                         tem2(jl) = -zdair(jl)
                      END IF
                   ELSE

                      zsufl(jl) = 2.*zsnlam*(td3m1m(jl)-tsnm1m(jl))/zsndp(jl)
                      zdsn(jl)  = -1./(zsncoe*zsndp(jl))*(zdair(jl)-2.*zsnlam/zsndp(jl))
                      zsnfl(jl) = -zair(jl) + zsufl(jl)

                      ! Snow melt

                      IF ((tsnm1m(jl)<tmelt) .OR. ((tsnm1m(jl)>= &
                           &                    tmelt) .AND. (zsnfl(jl)<=0.))) THEN
                         tsn(jl) = tsnm1m(jl) + (ztmst*zsnfl(jl)/(zsndp(jl)*zsncoe)/(1.- &
                              &                    0.5*ztmst*zdsn(jl)))
                         tem1(jl) = -zsufl(jl)
                         tem2(jl) = -2.*zsnlam/zsndp(jl)

                         ! Correction if snow surface temp would exceed tmelt

                         ts(jl) = 2.*tsn(jl) - td3m1m(jl)
                         IF (ts(jl)>=tmelt) THEN
                            zt(jl) = 0.5*(tmelt+td3m1m(jl))
                            tem1(jl) = tem1(jl) + zsncoe*zsndp(jl)*(tsn(jl)-zt(jl))/ztmst
                            tem2(jl) = -0.5*zsncoe*zsndp(jl)/ztmst
                            tsn(jl) = zt(jl)
                            ts(jl) = tmelt
                         END IF

                      ELSE

                         IF (td3m1m(jl)>=tmelt) THEN
                            zsnmel(jl) = MIN(snm1m(jl),ztmst*zsnfl(jl)/(rhoh2o*alf))
                            tem1(jl)   = -MIN((snm1m(jl)-ztmst*zsnfl(jl)/(rhoh2o*alf))*alf* &
                                 &                      rhoh2o/ztmst,0.)
                         ELSE
                            tem1(jl)   = zsnfl(jl)
                            zsnmel(jl) = 0.
                         END IF
                         tem2(jl) = -zdair(jl) + 2.*zsnlam/zsndp(jl)
                         tsn(jl)  = tmelt
                         ts(jl)   = tmelt
                      END IF

                   END IF

                ELSE

                   ! No snow on the ground

                   tem1(jl) = -zair(jl)
                   tem2(jl) = -zdair(jl)

                END IF
             END IF
          END DO

          !-- 2.1.3 Solution of soil diffusion eq.

          ! Set constants for ts scheme

          IF (lamip2) THEN
             DO jkg = 1,5
                zdel(jkg) = cdel(jkg)
             END DO
          ELSE
             zdel(1) = 0.065
             zdel(2) = 0.254
             zdel(3) = 0.913
             zdel(4) = 2.902
             zdel(5) = 5.700
          END IF

          zconb2 = 1./(0.5*(zdel(1)+zdel(2)))

          IF (.NOT.lamip2) THEN
             ztpfac1 = cvdifts
             ztpfac2 = 1./ztpfac1
             ztpfac3 = 1. - ztpfac2
          END IF

          DO jl = 1, klon
             zconb1(jl) = 1./MERGE(zrici,rgcgnm(jl),loglac(jl))
             ! NEW FOR IMPLICIT SOIL SCHEME
             IF (lamip2) zcapav(jl) = MERGE(zrici,rgcgnm(jl),loglac(jl))
             zdifi(jl) = MERGE(zdifiz,sodifm(jl),loglac(jl))
          END DO

          IF (lamip2) THEN
             DO jl = 1, klon
                tslin(jl) = tslinm(jl)  ! CHANGED FOR IMPLICIT SOIL SCHEME
             END DO
          ELSE
             DO jl = 1, klon
                IF (loland(jl)) THEN
                   zsa = ztmst*zconb1(jl)/zdel(1)
                   zsb = ztmst*zconb2*zdifi(jl)/zdel(1)
                   zdeff = -zair(jl) - tem1(jl)/(1.-0.5*(zsa*MIN(zsb/zsa,tem2(jl))-zsb))
                ELSE
                   zdeff = 0.
                END IF
                IF ( .NOT. loglac(jl) .AND. snm1m(jl)>zsncri .AND. loland(jl)) THEN
                   zdeff = -zair(jl) + zair(jl)/(1.-0.5*ztmst*zdsn(jl))
                END IF
                tslin(jl) = tslinm(jl) + zdiagt*zdeff
             END DO
          END IF

          IF (.NOT.lamip2) THEN

             ! Setting of right hand side

             ivpos = 1310
             DO jl = 1, klon
                IF (loland(jl)) THEN
                   zsdif(jl,1) = ztpfac2*td3m1m(jl)
                   zsdif(jl,2) = ztpfac2*td4m1m(jl)
                   zsdif(jl,3) = ztpfac2*td5m1m(jl)
                   zsdif(jl,4) = ztpfac2*tdm1m(jl)
                   zsdif(jl,5) = ztpfac2*tdclm1m(jl)
                END IF
             END DO

             ! Set a(k), b=1+a+c, c(k)=a(k-1)

             ivpos = 1320

             ! Number of vertical levels

             nsl = 5

             DO jk = 1, nsl - 1
                DO jl = 1, klon
                   IF (loland(jl)) THEN
                      zcfh(jl,jk) = ztpfac1*2.*ztmst*zdifi(jl)/(zdel(jk)+zdel(jk+1))
                   END IF
                END DO
             END DO

             ! Top layer elimination

             ivpos = 1330
             DO jl = 1, klon
                IF (loland(jl)) THEN
                   ztcoe(jl) = zcfh(jl,1)
                   zzdp = 1./zdel(1)
                   zdisc = 1./(1.+zcfh(jl,1)*zzdp)
                   zeb(jl,1) = zdisc*(zcfh(jl,1)*zzdp)
                   zsa = ztmst*zconb1(jl)*zzdp
                   zsb = ztmst*zconb2*zdifi(jl)*zzdp
                   zsdif(jl,1) = zdisc*(zsdif(jl,1)+tem1(jl)*zsa/(1.-0.5*(zsa*MIN(zsb/ &
                        &                zsa,tem2(jl))-zsb)))
                END IF
             END DO

             ! Elimination for middle layers

             ivpos = 1340
             DO jk = 2, nsl - 1
                DO jl = 1, klon
                   IF (loland(jl)) THEN
                      zzdp = 1./zdel(jk)
                      zfac = ztcoe(jl)*zzdp
                      ztcoe(jl) = zcfh(jl,jk)
                      zdisc = 1./(1.+zfac*(1.-zeb(jl,jk-1))+zcfh(jl,jk)*zzdp)
                      zeb(jl,jk) = zdisc*(zcfh(jl,jk)*zzdp)
                      zsdif(jl,jk) = zdisc*(zsdif(jl,jk)+zfac*zsdif(jl,jk-1))
                   END IF
                END DO
             END DO

             ! Bottom layer elimination

             ivpos = 1350
             DO jl = 1, klon
                IF (loland(jl)) THEN
                   zzdp  = 1./zdel(nsl)
                   zfac  = ztcoe(jl)*zzdp
                   zdisc = 1./(1.+zfac*(1.-zeb(jl,nsl-1)))
                   zsdif(jl,nsl) = zdisc*(zsdif(jl,nsl)+zfac*zsdif(jl,nsl-1))
                END IF
             END DO

             ! Back-substitution

             ivpos = 1360
             DO jkk = 1, nsl - 1
                jk = nsl - jkk
                DO jl = 1, klon
                   IF (loland(jl)) THEN
                      zsdif(jl,jk) = zsdif(jl,jk) + zeb(jl,jk)*zsdif(jl,jk+1)
                   END IF
                END DO
             END DO
          END IF

          IF (lamip2) THEN
             ! Compute soil heat quantities

             CALL soil2ech(jrow, klon, dtime, ts, zdifi, zcapav, snm1m)
          END IF

          ! Return to ts values

          ivpos = 1370
          DO jl = 1, klon
             IF (loland(jl)) THEN
                IF (lamip2) THEN
                   td3(jl)  = ptn(jl,jrow,1) 
                   td4(jl)  = ptn(jl,jrow,2)
                   td5(jl)  = ptn(jl,jrow,3)
                   td(jl)   = ptn(jl,jrow,4)
                   tdcl(jl) = ptn(jl,jrow,5)
                ELSE
                   td3(jl)  = zsdif(jl,1) + ztpfac3*td3m1m(jl)
                   td4(jl)  = zsdif(jl,2) + ztpfac3*td4m1m(jl)
                   td5(jl)  = zsdif(jl,3) + ztpfac3*td5m1m(jl)
                   td(jl)   = zsdif(jl,4) + ztpfac3*tdm1m(jl)
                   tdcl(jl) = zsdif(jl,5) + ztpfac3*tdclm1m(jl)
                END IF

                IF (snm1m(jl) < zsncri .OR. loglac(jl)) THEN
                   ts(jl)  = td3(jl)
                   tsn(jl) = td3(jl)
                END IF
             ELSE
                tsn(jl)  = tsnm1m(jl)
                ts(jl)   = tsm1m(jl)
                td3(jl)  = td3m1m(jl)
                td4(jl)  = td4m1m(jl)
                td5(jl)  = td5m1m(jl)
                td(jl)   = tdm1m(jl)
                tdcl(jl) = tdclm1m(jl)
             END IF
          END DO

          !-- 2.2 Moisture changes

          ! Evaporation of the skin reservoir

          DO jl = 1, klon
             IF (lamip2) THEN
                ! No evaporation over gaciers !
                ! Now possible because of changes in VDIFF !
                qhfl(jl) = qhfl(jl)-cvs(jl)*dhfqs(jl)
             END IF
             lo = dhfqw(jl) < 0. .AND. loland(jl)
             zdhfqw   = zcons5*dhfqw(jl)
             zqhflw   = zdhfqw*wlm1m(jl)
             zwl(jl)  = wlm1m(jl) + MERGE(zqhflw,0.,lo)
             zwl(jl)  = MAX(0.,zwl(jl))
             zqhflw   = (zwl(jl)-wlm1m(jl))/zcons5
             qhfl(jl) = qhfl(jl) - zqhflw
          END DO

          ! Collection of dew by the skin reservoir

          DO jl = 1, klon
             zqhflw = zcons5*qhfl(jl)
             lo = zqhflw > 0. .AND. loland(jl)
             zwl(jl)  = zwl(jl) + MERGE(zqhflw,0.,lo)
             qhfl(jl) = MERGE(MAX(0.,zwl(jl)-wlmx(jl))/zcons5,qhfl(jl),lo)
             zwl(jl)  = MIN(zwl(jl),wlmx(jl))
          END DO

!OCL NOALIAS
          DO jl = 1, klon
             zqrfld(jl) = (zwl(jl)-wlm1m(jl))/zcons5
          END DO

          zsig0  = zdiagw*zbud*DOT_PRODUCT(slmm(1:klon),zqrfld(1:klon))
          zsig10 = zdiagw*zbud*DOT_PRODUCT(slmm(1:klon),rsfl(1:klon))
          zsig11 = zdiagw*zbud*DOT_PRODUCT(slmm(1:klon),rsfc(1:klon))

          ! Interception of precipitation by the vegetation

          ! Large scale precipitation

          DO jl = 1, klon
             ztprcp = rsfl(jl)
             IF (lamip2) THEN
                ziprcp = MERGE(ztprcp*zvinter,0.,loland(jl))
             ELSE
                ziprcp = MERGE(ztprcp*vgrat(jl)*zvinter,0.,loland(jl))
             END IF
             zmprcp   = MIN(ziprcp,(wlmx(jl)-zwl(jl))/zcons5)
             zwl(jl)  = zwl(jl) + zmprcp*zcons5
             rsfl(jl) = rsfl(jl) - zmprcp
          END DO

          ! Convective precipitation

          DO jl = 1, klon
             ztprcp = rsfc(jl)
             IF (lamip2) THEN
                ziprcp = MERGE(ztprcp*zvinter,0.,loland(jl))
             ELSE
                ziprcp = MERGE(ztprcp*vgrat(jl)*zvinter,0.,loland(jl))
             END IF
             zmprcp   = MIN(ziprcp,(wlmx(jl)-zwl(jl))/zcons5)*zpsfr
             zwl(jl)  = zwl(jl) + zmprcp*zcons5
             rsfc(jl) = rsfc(jl) - zmprcp
          END DO

!OCL NOALIAS
          DO jl = 1, klon
             lo = qhfl(jl) >= 0.
             zqfld(jl) = qhfl(jl)
             zcons12 = 1.
             zprfl(jl) = rsfl(jl) + rsfc(jl) + zqfld(jl)*MERGE(1.,zcons12,lo)
             IF (lamip2) THEN
                zpsfl(jl) = ssfl(jl) + ssfc(jl) + cvs(jl)*(dhfqs(jl)+xhfl(jl))
             ELSE
                zpsfl(jl) = ssfl(jl) + ssfc(jl) + cvs(jl)*dhfqs(jl)
             END IF
             zros(jl)   = 0.
             zinfil(jl) = 0.
             zper(jl)   = 0.
             zdrain(jl) = 0.
             zvol(jl)   = 0.
             ztsdfl(jl) = sodifm(jl)*rgcgnm(jl)*zdel(1)/ &
                  &            (0.5*zdel(1)*(zdel(1)+zdel(2)))*(td4m1m(jl)-td3m1m(jl))
          END DO

          zsdtfl = zdiagt*zbud*DOT_PRODUCT(slmm(1:klon),ztsdfl(1:klon))
          zsig1  = zdiagw*zbud*DOT_PRODUCT(slmm(1:klon),zprfl(1:klon))
          zsig2  = zdiagw*zbud*DOT_PRODUCT(slmm(1:klon),zpsfl(1:klon))
          zsig3  = zdiagw*zbud*DOT_PRODUCT(slmm(1:klon),rsfl(1:klon))
          zsig4  = zdiagw*zbud*DOT_PRODUCT(slmm(1:klon),ssfl(1:klon))
          zsig5  = zdiagw*zbud*DOT_PRODUCT(slmm(1:klon),rsfc(1:klon))
          zsig6  = zdiagw*zbud*DOT_PRODUCT(slmm(1:klon),ssfc(1:klon))

          !-- 2.3 Snow changes

          DO jl = 1, klon
             IF (loland(jl)) THEN
                IF (loglac(jl)) THEN
                   sn(jl) = snm1m(jl)
                ELSE
                   sn(jl) = snm1m(jl) + zcons5*zpsfl(jl)
                   sn(jl) = MAX(sn(jl),0.)
                END IF
             ELSE
                sn(jl) = 0.
             END IF
          END DO

          !-- 3. Snow melt and corrections (including run off)

          !-- 3.1 Snow melt and snow correction

          IF (lamip2) THEN
            DO jl = 1, klon
              zconss9 = zcons5*zdel(1)/(alf*ztmst*zconb1(jl))
              zconss10 = 1./zconss9
              IF ((sn(jl) > 0.) .AND. (loland(jl))) THEN
                IF ((sn(jl) >= zsncri).AND.(.NOT.loglac(jl))) THEN
                  IF (ts(jl) > tmelt) THEN
                    IF (td3(jl) <= tmelt) THEN
                      pfluxgrd(jl,jrow) = pfluxgrd(jl,jrow) &
                           + pcapcal(jl,jrow)*(ts(jl)-tmelt)/dtime
                      zsnmel(jl) = 0.
                      ts(jl)  = tmelt
                      tsn(jl) = td3(jl)
                    ELSE
                      zsnmel(jl) = MIN(sn(jl),pcapcal(jl,jrow)/alf &
                           * (td3(jl)-tmelt)/rhoh2o)
                      td3(jl) = td3(jl)-zsnmel(jl)*alf/pcapcal(jl,jrow)*rhoh2o
                      ts(jl)  = td3(jl)
                      tsn(jl) = td3(jl)
                    END IF
                  ELSE
                    tsn(jl) = td3(jl)
                  END IF
                END IF

                IF (td3(jl) > tmelt .AND. (sn(jl) < zsncri .OR. loglac(jl))) THEN
                  zsnmel(jl) = MIN(sn(jl),zconss9*(td3(jl)-tmelt))
                  td3(jl) = td3(jl) - zsnmel(jl)*zconss10
                  tsn(jl) = td3(jl)
                  ts(jl)  = td3(jl)
                END IF
                IF (loglac(jl)) THEN
                  zsnmel(jl) = 0.
                END IF
                sn(jl) = sn(jl) - zsnmel(jl)
              END IF
              ! ---------------------------------------------------------------------------------

              IF ((sn(jl)<zsncri) .AND. loland(jl)) THEN
                ts(jl)  = td3(jl)
                tsn(jl) = td3(jl)
              END IF

              ! Allow for some water to remain on the ground or on the leaves

              zsnmlt  = MAX(0.,zsnmel(jl))
              zwl(jl) = zwl(jl) + zvinter*zsnmlt
              zsnmlt  = (1.-zvinter)*zsnmlt + MAX(0.,zwl(jl)-wlmx(jl))
              zwl(jl) = MIN(wlmx(jl),zwl(jl))
              IF (loland(jl) .AND. .NOT. loglac(jl)) THEN
                zprfl(jl) = zsnmlt*rhoh2o/ztmst
              ELSE
                zprfl(jl) = 0.
              END IF
            END DO
          ELSE
            DO jl = 1, klon
              zconss9 = zcons5*zdel(1)/(alf*ztmst*zconb1(jl))
              zconss10 = 1./zconss9
              IF ((sn(jl) > 0.) .AND. (loland(jl))) THEN
                IF (td3(jl) > tmelt .AND. (sn(jl) < zsncri .OR. loglac(jl))) THEN
                  zsnmel(jl) = MIN(sn(jl),zconss9*(td3(jl)-tmelt))
                  td3(jl) = td3(jl) - zsnmel(jl)*zconss10
                  tsn(jl) = td3(jl)
                  ts(jl)  = td3(jl)
                END IF
                IF (loglac(jl)) THEN
                  zsnmel(jl) = 0.
                END IF
                sn(jl) = sn(jl) - zsnmel(jl)
              END IF
              ! ---------------------------------------------------------------------------------

              IF ((sn(jl)<zsncri) .AND. loland(jl)) THEN
                ts(jl)  = td3(jl)
                tsn(jl) = td3(jl)
              END IF

              ! Allow for some water to remain on the ground or on the leaves

              zsnmlt  = MAX(0.,zsnmel(jl))
              zwl(jl) = zwl(jl) + zvinter*zsnmlt
              zsnmlt  = (1.-zvinter)*zsnmlt + MAX(0.,zwl(jl)-wlmx(jl))
              zwl(jl) = MIN(wlmx(jl),zwl(jl))
              IF (loland(jl) .AND. .NOT. loglac(jl)) THEN
                zprfl(jl) = zsnmlt*rhoh2o/ztmst
              ELSE
                zprfl(jl) = 0.
              END IF
            END DO
          END IF

          ! Summation of budgets.

          zssnmt = zdiags*zbud*DOT_PRODUCT(slmm(1:klon),zsnmel(1:klon))

          !-- 3.2 Compute soil hydrology

!DIR$ IVDEP
!OCL NOVREC
          DO jl = 1, klon

             ! Account for surface runoff due to sloping terrain

             zwmax  = wsmxm(jl)
             zwptr  = 0.90*zwmax
             zwdtr  = 0.90*zwmax
             zconw2 = zwmax - zwdtr
             zconw3 = zdrmax - zdrmin
             zconw4 = zwmax - zwptr
             zroeff = MAX(0.,SQRT(varor(jl))-zorvari)/(SQRT(varor(jl))+zorvars)
             zbws   = MAX(MIN(zroeff,0.5),0.01)
             zb1    = 1. + zbws
             zbm    = 1./zb1
             zconw1 = zwmax*zb1

             IF (loland(jl)) THEN

                lo = qhfl(jl) >= 0.
                zprfl(jl) = zprfl(jl) + rsfl(jl) + rsfc(jl) + &
                     &              qhfl(jl)*MERGE(1.,0.,lo)

                ! Compute infiltration from rainfall and runoff

                zlyeps = 0.
                IF (td3(jl)<tmelt) THEN
                   zros(jl) = zprfl(jl)*zcons5
                   ws(jl) = wsm1m(jl)
                ELSE
                   IF (zprfl(jl)>0.) THEN
                      IF (wsm1m(jl)>zwmax) THEN
                         zlyeps = wsm1m(jl) - zwmax
                      ELSE
                         zlyeps = 0.
                      END IF
                      zlysic = (wsm1m(jl)-zlyeps)/zwmax
                      zlysic = MIN(zlysic,1.)
                      zvol(jl) = (1.-zlysic)**zbm - ztmst*zprfl(jl)/(rhoh2o*zconw1)
                      zros(jl) = zprfl(jl)*ztmst/rhoh2o - (zwmax-wsm1m(jl))
                      IF (zvol(jl)>0.) THEN
                         zros(jl) = zros(jl) + zwmax*zvol(jl)**zb1
                      END IF
                      zros(jl)   = MAX(zros(jl),0.)
                      zinfil(jl) = zprfl(jl)*ztmst/rhoh2o - zros(jl)
                   END IF
                   ws(jl) = wsm1m(jl) + zinfil(jl)
                END IF

                IF (.NOT.lamip2 .AND. loglac(jl)) THEN
                   zros(jl) = 0.
                   ws(jl)   = wsm1m(jl)
                END IF

                ! Subtract evaporation

                IF (lamip2) THEN
                   IF (.NOT. loglac(jl)) THEN
                      IF (qhfl(jl) < 0.) THEN
                         ws(jl) = ws(jl)+qhfl(jl)*zcons5
                      END IF
                      ws(jl)   = ws(jl)+(1.-cvs(jl))*xhfl(jl)*zcons5
                   ELSE
                      zros(jl) = 0.
                      ws(jl)   = 0.
                   END IF
                ELSE
                   IF (qhfl(jl) < 0.) THEN
                      ws(jl) = ws(jl) + qhfl(jl)*ztmst/rhoh2o
                   END IF
                   ws(jl) = ws(jl) + xhfl(jl)*zcons5
                END IF

                IF (ws(jl)<0.) THEN
                   ws(jl) = 0.
                ELSE

                   ! Soil water depletion

                   IF ( .NOT. loglac(jl) .AND. td3(jl)>tmelt) THEN

                      zdrain(jl) = zdrmin*ws(jl)/zwmax
                      IF (ws(jl)>zwdtr) THEN
                         zdrain(jl) = zdrain(jl) + zconw3*((ws(jl)-zwdtr)/zconw2)** &
                              &                    zdrexp
                      END IF
                      zdrain(jl) = zdrain(jl)*ztmst

                      zwslim = 0.05*wsmxm(jl)
                      zdrain(jl) = MIN(zdrain(jl),MAX(ws(jl)-zwslim,0.))
                      ws(jl) = ws(jl) - zdrain(jl)

                   ELSE
                      zdrain(jl) = 0.
                   END IF

                   ! Soil water depletion : percolation

                   !  No percolation 

                   zper(jl) = 0.

                END IF

                zwsup    = MAX(ws(jl)-wsmxm(jl),0.)
                ws(jl)   = ws(jl) - zwsup
!                zros(jl) = zros(jl) + zdrain(jl) + zper(jl) + zwsup
                zros(jl) = zros(jl) + zper(jl) + zwsup
             ELSE
                IF (lamip2) THEN
                   ws(jl) = 0.
                ELSE
                   ws(jl) = wsm1m(jl)
                END IF
             END IF
          END DO

          IF (lamip2) THEN

             ! Compute variance of surface temperature

             DO jl = 1,klon
                ztsvar(jl) = (ts(jl)-tsm1(jl))*(ts(jl)-tsm1(jl))/(dtime*dtime)
             END DO
          END IF


          !-- 3.3 Accumulate surface parameters for diagnostics

          DO jl = 1, klon
             tsurf(jl) = tsurfm(jl) + zdiagt*ts(jl)
             !tsvari(jl) = tsvarim(jl) + zdiagt*ztsvar(jl) ! IMPLICIT SOIL SCHEME
             runoff(jl) = runoffm(jl) + zdiags*zros(jl)
             drain(jl) = drainm(jl) + zdiags*zdrain(jl)
             snmel(jl) = snmelm(jl) + zdiags*zsnmel(jl)
             zsn = (sn(jl)-snm1m(jl))
             dsnac(jl)  = dsnacm(jl) + zdiags*zsn + zdiagw*zlac(jl)
          END DO

          zsros = zdiags*zbud*DOT_PRODUCT(slmm(1:klon),zros(1:klon))

          dssradz(irow) = zssrad
          dstradz(irow) = zstrad
          dshflz(irow)  = zshfl

          dsdtflz(irow) = zsdtfl
          dslsrz(irow)  = zsig10
          dslssz(irow)  = zsig4
          dscvrz(irow)  = zsig11
          dscvsz(irow)  = zsig6
          dsevwz(irow)  = (zsig3+zsig5-zsig0-zsig1)
          dseviz(irow)  = (zsig4+zsig6-zsig2)

          dssnmtz(irow) = zssnmt

          dsrosz(irow)  = zsros

          !-- 4. Time filter (Temperature values are copied because of the full
          !                   implicit scheme)

          IF (nstep/=nstart) THEN
             DO jl = 1, klon
                ! Replaced by implicit soil scheme
                IF (lamip2) THEN
                   tdm1(jl) = tdm(jl)
                ELSE
                   tsm1(jl) = tsm(jl) + eps*(tsm1m(jl)-2.*tsm(jl)+ts(jl))
                   tdm1(jl) = tdm(jl) + eps*(tdm1m(jl)-2.*tdm(jl)+td(jl))
                END IF
                wlm1(jl) = wlm(jl) + eps*(wlm1m(jl)-2.*wlm(jl)+zwl(jl))
                wsm1(jl) = wsm(jl) + eps*(wsm1m(jl)-2.*wsm(jl)+ws(jl))
                snm1(jl) = snm(jl) + eps*(snm1m(jl)-2.*snm(jl)+sn(jl))
                ! Replaced by implicit soil scheme
                IF (lamip2) THEN
                   tsnm1(jl)  = tsnm(jl)
                   td3m1(jl)  = td3m(jl)
                   td4m1(jl)  = td4m(jl)
                   td5m1(jl)  = td5m(jl)
                   tdclm1(jl) = tdclm(jl)
                ELSE
                   tsnm1(jl)  = tsnm(jl)  + eps*(tsnm1m(jl) -2.*tsnm(jl) +tsn(jl))
                   td3m1(jl)  = td3m(jl)  + eps*(td3m1m(jl) -2.*td3m(jl) +td3(jl))
                   td4m1(jl)  = td4m(jl)  + eps*(td4m1m(jl) -2.*td4m(jl) +td4(jl))
                   td5m1(jl)  = td5m(jl)  + eps*(td5m1m(jl) -2.*td5m(jl) +td5(jl))
                   tdclm1(jl) = tdclm(jl) + eps*(tdclm1m(jl)-2.*tdclm(jl)+tdcl(jl))
                END IF
             END DO
          ELSE
             DO jl = 1, klon
                ! REMOVED FOR IMPLICIT SOIL SCHEME
                IF(.NOT. lamip2)  tsm1(jl) = tsm(jl)
                tdm1(jl)   = tdm(jl)
                wlm1(jl)   = wlm(jl)
                wsm1(jl)   = wsm(jl)
                snm1(jl)   = snm(jl)
                tsnm1(jl)  = tsnm(jl)
                td3m1(jl)  = td3m(jl)
                td4m1(jl)  = td4m(jl)
                td5m1(jl)  = td5m(jl)
                tdclm1(jl) = tdclm(jl)
             END DO
          END IF
       END IF

       !-- 4.1 Surface maximum and minimum temperature

       DO jl = 1, klon
          tsmax(jl) = MAX(ts(jl),tsmaxm(jl))
          tsmin(jl) = MIN(ts(jl),tsminm(jl))
       END DO

       !-- 5. Necessary computations if subroutine is by-passed

    ELSE
       DO jl = 1, klon
          tsm1(jl)   = tsm(jl)
          tdm1(jl)   = tdm(jl)
          wlm1(jl)   = wlm(jl)
          wsm1(jl)   = wsm(jl)
          snm1(jl)   = snm(jl)
          ts(jl)     = tsm1m(jl)
          td(jl)     = tdm1m(jl)
          zwl(jl)    = wlm1m(jl)
          ws(jl)     = wsm1m(jl)
          sn(jl)     = snm1m(jl)
          tslin(jl)  = tslinm(jl)
          runoff(jl) = runoffm(jl)
          drain(jl)  = drainm(jl)
          tsurf(jl)  = tsurfm(jl) + zdiagt*ts(jl)
          tsn(jl)    = tsnm1m(jl)
          td3(jl)    = td3m1m(jl)
          td4(jl)    = td4m1m(jl)
          td5(jl)    = td5m1m(jl)
          tdcl(jl)   = tdclm1m(jl)
          tsnm1(jl)  = tsnm(jl)
          td3m1(jl)  = td3m(jl)
          td4m1(jl)  = td4m(jl)
          td5m1(jl)  = td5m(jl)
          tdclm1(jl) = tdclm(jl)
          IF (lamip2) THEN
             !drain(jl)  = drainm(jl)
             !tsvari(jl) = tsvarim(jl)
             tsmax(jl) = -99.
             tsmin(jl) = 999.
             snmel(jl) =   0.
          END IF
       END DO
       dssradz(irow) = 0.
       dstradz(irow) = 0.
       dshflz(irow)  = 0.
       
       dsdtflz(irow) = 0.
       dslsrz(irow)  = 0.
       dslssz(irow)  = 0.
       dscvrz(irow)  = 0.
       dscvsz(irow)  = 0.
       dsevwz(irow)  = 0.
       dseviz(irow)  = 0.

       dssnmtz(irow) = 0.

       dsrosz(irow)  = 0.
    END IF

    ! Add sea ice heat flux to tslin

    DO jl = 1, klon
       wl(jl)    = zwl(jl)
       tslin(jl) = tslin(jl) + zdiagt*ahfli(jl)
    END DO

    RETURN
  END SUBROUTINE surf

!END MODULE m_surf
