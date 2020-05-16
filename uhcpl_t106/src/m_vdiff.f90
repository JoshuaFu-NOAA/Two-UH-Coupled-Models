
!#define DEBUG

!MODULE m_vdiff
!
!CONTAINS

  SUBROUTINE vdiff(kidia,kfdia,klon,klp2,ktdia,klev,klevm1,klevp1,ktrac,pxtm1,             &
                   paclcm,paphm1,papm1,pgeom1,pqm1,ptkem,ptkem1m,ptm1,pum1,pvm1,pxm1,pemterm,ptvm1, &
                   laland,pahflm,pahfsm,paz0m,pdew2m,pevapm,pforestm,pseaice,psnm1m,psrfl, &
                   ptemp2m,ptsm1m,pt2maxm,pt2minm,pustar3m,pustrm,pu10m,pvdism,pvstrm,     &
                   pv10m,pwimaxm,pwind10m,pwlm1m,pwsm1m,pwsmxm,pvltm,ptke,ptkem1,ktropo,   &
                   pahfl,pahfs,paz0,pcvs,pcvw,pdew2,pdhfqs,pdhfqw,pdhft,pevap,pts,pqhfl,prsfl, &
                   ptemp2,pthfl,pt2max,pt2min,pustar3,pustr,pu10,pvdis,pvstr,pv10,pwimax,  &
                   pwlmx,pwind10,pxhfl,zqsnew,zqsnfl,zqlfl,zqvfl,zqgfl,pxtte,pvol,pvom,    &
                   pqte,ptte,pxte,pvgrat)

    ! Description:
    !
    ! Does the vertical exchange of u,v,t,q by turbulence.
    !
    ! Method:
    !
    ! This routine computes the physical tendencies of the four
    ! prognostic variables u,v,t and q due to the vertical exchange by
    ! turbulent (= non-moist convective) processes. These tendencies are
    ! obtained as the difference between the result of an implicit
    ! time-step starting from values at t-1 and these t-1 values. All
    ! the diagnostic computations (exchange coefficients, .) are done
    ! from the t-1 values. As a by-product the roughness length over sea
    ! is updated accordingly to the *Charnock formula. Heat and moisture
    ! surface fluxes and their derivatives against ts, ws and wl
    ! (the latter will be later weighted with the snow factor in
    ! *vdiff*), later to be used for soil processes treatment, are also
    ! computed as well as a stability value to be used as a diagnostic
    ! of the depth of the well mixed layer in convective computations.
    !
    ! *vdiff* is called from *physc*.
    !
    ! First an auxialiary variable cp(q)t+gz is created on which
    ! the vertical diffusion process will work like on u,v and q. Then
    ! along the vertical and at the surface, exchange coefficients (with
    ! the dimension of a pressure thickness) are computed for momentum
    ! and for heat (sensible plus latent). The letters m and h are used
    ! to distinguish them. the diffusioncoefficents depend on the
    ! turbulent kinetic energy (tke) calculated by an additional
    ! prognostic equation, which considers advektion of tke.
    ! 
    ! In the second part of the routine the implicit linear
    ! systems for u,v first and t,q second are solved by a *Gaussian
    ! elimination back-substitution method. For t and q the lower
    ! boundary condition depends on the surface state.
    ! For tke the lower boundary condition depends on the square of
    ! the frictional velocity.
    ! 
    ! Over land, two different regimes of evaporation prevail:
    ! A stomatal resistance dependent one over the vegetated part
    ! and a soil relative humidity dependent one over the
    ! bare soil part of the grid mesh.
    ! Potential evaporation takes place over the sea, the snow
    ! covered part and the liquid water covered part of the
    ! grid mesh as well as in case of dew deposition.
    ! 
    ! Finally one returns to the variable temperature to compute
    ! its tendency and the later is modified by the dissipation's effect
    ! (one assumes no storage in the turbulent kinetic energy range) and
    ! the effect of moisture diffusion on cp. z0 is updated and the
    ! surface fluxes of t and q and their derivatives are prepared and
    ! stored like the difference between the implicitely obtained
    ! cp(q)t+gz and cp(q)t at the surface.
    ! 
    ! Reference:
    ! See vertical diffusion's part of the model's documentation
    ! for details about the mathematics of this routine.
    !
    ! Authors:
    !
    ! J. F. Geleyn, ECMWF, in 1982, original source
    ! C. B. Blondin, ECMWF, in 1986, changed
    ! J. Feichter, MI, in 1991, changed
    ! S. Brinkop, MPI, in 1992, changed
    ! U. Schlese, DKRZ, February 1993, changed
    ! M. Claussen, MPI, in 1993, changed
    ! E. Roeckner, MPI, in 1994, changed
    ! J.-P. Schulz, MPI, in 1997, Implementation of the implicit coupling
    !                             between land surface and atmosphere.
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A. Rhodin, MPI, December 1998, lookup table removed
    !                 January 1999,  subroutine put into module
    ! I. Kirchner, IfM, September 2005, coupling with ocean
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_exception,        ONLY:finish
    USE mo_control,          ONLY:dtime   ,&!time step (in seconds)
                                  eps     ,&!time filtering coefficien
                                  nrow    ,&!current latit. line(one entry/task)
                                  twodt   ,&!2.*dtime
                                  nlev    ,&!number of levels
                                  nlon    ,&!
                                  lcouple ,&!ocean coupling
                                  ldsst   ,&!turn on SST diurnal-cycle
                                  lamip2
    USE mo_gaussgrid,        ONLY:budw    ,&!weights for global budgets
                                  coriol    !coriolis parameter:2*omega*mu
    USE mo_diagnostics_zonal,ONLY:devapz  ,&
                                  dhfsz   ,&
                                  dvdisz
    USE mo_param_switches,   ONLY:lvdiff    !true for vertical diffusion
    USE mo_physc2,           ONLY:cb      ,&!stability parameter near neutrality
                                  cc      ,&!stability parameter for unst. cases
                                  cchar   ,&!charnock constant
                                  cd      ,&!stability parameter for stable cases
                                  ckap    ,&!karman constant
                                  clam    ,&!asymptotic mixing length for moment.
                                  cqsncr  ,&!inv. of equiv. water height whensnow
                                  cvdifts ,&!factor for timestep weighting
                                  cwlmax  ,&!max. moist. content of skin reserv.
                                  cz0ice    !roughness over sea-ice
    USE mo_constants,        ONLY:als     ,&!latent heat for sublimation
                                  alv     ,&!latent heat for vaporisation
                                  api     ,&!2.*arcsin(1.)
                                  c2es    ,&!constants used for computation 
                                  c3ies   ,&!  of saturation mixing ratio 
                                  c3les   ,&!  over liquid water(*c_les*) or
                                  c4ies   ,&!  ice(*c_ies*).
                                  c4les   ,&!
                                  c5ies   ,&!
                                  c5les   ,&!
                                  cpd     ,&!pec. heat at const. press. (dry air)
                                  g       ,&!gravity acceleration
                                  rd      ,&!gas constant for dry air
                                  rhoh2o  ,&!density of liquid water
                                  tmelt   ,&!temperature of fusion of ice
                                  vtmpc1  ,&!vtmpc1=rv/rd-1
                                  vtmpc2  ,&!vtmpc2=cpv/cpd-1
                                  stbo
    USE mo_start_dataset,    ONLY:nstart  ,&!time step for start/restart
                                  nstep     !current time step
    USE mo_vegetation,       ONLY:cva     ,&!const. to define stomat. resistance
                                  cvb     ,&!const. to define stomat. resistance
                                  cvc     ,&!minimum stomatal resistance
                                  cvabc   ,&!(cva+cvbc)/cvc
                                  cvbc    ,&!cvb*cvc
                                  cvk     ,&!?
                                  cvkc    ,&!cvk*cvc
                                  cvrad     !frac. of net s.w rad. contr.to p.a.r
    USE mo_tracer,           ONLY: lemis, ntrac, xtemiss, lxtvdiff
#ifndef NOLOOKUP
    USE mo_convect_tables,   ONLY: tlucua, jptlucu1, jptlucu2, &
                                   lookuperror, lookupoverflow
#endif
    USE mo_doctor,           ONLY: nout
    USE mo_radint,           ONLY: cemiss
    USE mo_soil_impl,        ONLY: pfluxgrd, pcapcal, evapotr4, surftemp

    USE mo_couple,           ONLY: wind10_local, vstr_local, ustr_local,&
         ahfs_local, ahfl_local, accu_time
    !fu++
    USE mo_dsst,             ONLY: wind10_diurnal

    IMPLICIT NONE

    !  Scalar arguments:
    INTEGER, INTENT (IN) :: &
         kidia, &
         kfdia, &
         klon,  &
         klp2,  &
         ktdia, &
         klev,  &
         klevm1,&
         klevp1,&
         ktrac
    !  Array arguments:
    REAL, INTENT (IN) ::           &
         pxtm1   (klon,klev,ktrac),&!tracer variables (t-dt)
         paclcm  (klp2,klev)      ,&! cloud cover (old value)
         paphm1  (klp2,klevp1)    ,&! half level pressure (t-dt)
         papm1   (klp2,klev)      ,&! full level pressure (t-dt)
         pgeom1  (klp2,klev)      ,&! geopotential above surface (t-dt)
         pqm1    (klp2,klev)      ,&! humidity (t-dt)
         ptm1    (klp2,klev)      ,&! temperature (t-dt)
         pum1    (klp2,klev)      ,&! zonal wind (t-dt)
         pvm1    (klp2,klev)      ,&! meridional wind (t-dt)
         pxm1    (klp2,klev)        ! cloud water (t-dt)

    ! INCLUDED FOR THE IMPLICIT SOIL SCHEME 
    REAL, INTENT (IN) ::        &
         pemterm (klp2,klevp1)    ,&
         ptvm1   (klp2,klev)         ! virtual temperature at t-dt

    LOGICAL, INTENT (IN) :: &
         laland  (klp2)         ! land-sea flag
    REAL, INTENT (IN) ::    &
         pahflm  (klp2)       ,&! surface latent heat flux (old value)
         pahfsm  (klp2)       ,&! surface sensible heat flux (old value)
         paz0m   (klp2)       ,&! roughness length (old value)
         pdew2m  (klp2)       ,&! dew point temperature at 2 meter(accum.,old value)
         pevapm  (klp2)       ,&! surface evaporation (accumulated, old value)
         pforestm(klp2)       ,&! ?
         pseaice (klp2)       ,&! sea ice cover (new value)
         psnm1m  (klp2)       ,&! snow depth (t-dt)
         psrfl   (klp2)       ,&! net solar radiative flux at the surface
         ptemp2m (klp2)       ,&! temperature at 2 meter (accumulated, old value)
         ptsm1m  (klp2)       ,&! surface temperature (t-dt)
         pt2maxm (klp2)       ,&! max. temp. at 2m between output intervals (old v.)
         pt2minm (klp2)       ,&! min. temp. at 2m between output intervals (old v.)
         pustar3m(klp2)       ,&! tke for ocean mixed layer (accumulated, old value)
         pustrm  (klp2)       ,&! u-stress (accumulated, old value)
         pu10m   (klp2)       ,&! u-wind at 10 meter (accumulated, old value)
         pvdism  (klp2)       ,&! boundary layer dissipation(accumulated, old value)
         pvstrm  (klp2)       ,&! v-stress (accumulated, old value)
         pv10m   (klp2)       ,&! v-wind at 10 meter (accumulated, old value)
         pwimaxm (klp2)       ,&! max windspeed at 10m. betw. outp. interv. (old v.)
         pwind10m(klp2)       ,&! wind speed at 10 meter (accumulated, old value)
         pwlm1m  (klp2)       ,&! skin reservoir content (t-dt)
         pwsmxm  (klp2)       ,&! field capacity of soil
         pvltm   (klp2)       ,&! leaf area index
         pvgrat  (klp2)         ! vegetation ratio
    REAL, INTENT (INOUT) ::    &
         pxtte  (klon,klev,ktrac),&! tendencies of tracer variables
         pvol   (klp2,klev)      ,&! tendency of meridional wind
         pvom   (klp2,klev)      ,&! tendency of zonal wind
         pqte   (klp2,klev)      ,&! tendency of humidity
         ptte   (klp2,klev)      ,&! tendency of temperature
         pxte   (klp2,klev)      ,&! tendency of cloud water
         pwsm1m (klp2)           ,&! surface soil wetness (t-dt)
         ptkem1m(klp2,klev)      ,&! turbulent kinetic energy (t-dt)
         ptkem  (klp2,klev)        ! turbulent kinetic energy
    REAL, INTENT (OUT) :: &
         ptke   (klp2,klev) ,&! turbulent kinetic energy (t+dt)
         ptkem1 (klp2,klev) ,&! turbulent kinetic energy (filtered)
         pahfl  (klp2)      ,&! surface latent heat flux (new value)
         pahfs  (klp2)      ,&! surface sensible heat flux (new value)
         paz0   (klp2)      ,&! roughness length (new value)
         pcvs   (klp2)      ,&! snow cover fraction
         pcvw   (klp2)      ,&! wet skin fraction
         pdew2  (klp2)      ,&! dew point temperature at 2m (accumulated, new value)
         pdhfqs (klp2)      ,&! deriv. of moist. flux over snow with resp. to snow d
         pdhfqw (klp2)      ,&! deriv. of moist. flux with respect to skin reservoir
         pdhft  (klp2)      ,&! deriv. of sens. heat flux with resp. to surf. temp.
         pevap  (klp2)      ,&! surface evaporation (accumulated, new value)
         pqhfl  (klp2)      ,&! moisture flux at the surface
         prsfl  (klp2)      ,&! large scale rain flux at the surface
         ptemp2 (klp2)      ,&! temperature at 2 meter (accumulated, new value)
         pthfl  (klp2)      ,&! total heat flux at the surface
         pt2max (klp2)      ,&! max temp. at 2m between output intervals (new value)
         pt2min (klp2)      ,&! min temp. at 2m between output intervals (new value)
         pustar3(klp2)      ,&! tke for ocean mixed layer (accumulated, new value)
         pustr  (klp2)      ,&! u-stress (accumulated, new value)
         pu10   (klp2)      ,&! u-wind at 10 meter (accumulated, new value)
         pvdis  (klp2)      ,&! boundary layer dissipation (accumulated, new value)
         pvstr  (klp2)      ,&! v-stress (accumulated, new value)
         pv10   (klp2)      ,&! v-wind at 10 meter (accumulated, new value)
         pwimax (klp2)      ,&! max windspeed at 10m betw. outp. interv. (new value)
         pwlmx  (klp2)      ,&! maximum skin reservoir contnet
         pwind10(klp2)      ,&! wind speed at 10 meter (accumulated, new value)
         pxhfl  (klp2)        ! liquid water flux at the surface

    REAL :: zqhfl  (klp2)

    ! INCLUDED FOR THE IMPLICIT SOIL SCHEME 
    REAL, INTENT (OUT) ::  pts    (klp2)

    ! This dummy argument is not used in the subprogram.(name:ktropo)
    INTEGER, INTENT (INOUT) :: ktropo(klon) ! tropopause index

    !  Local scalars: 
    REAL :: z0h, z1dgam, z2geomf, zabcs, zalf, zalh2, zalo, zaloh, zalvs,           &
            zaph2m, zb, zbet, zbud, zbuoy, zc, zc3ies, zc3les, zc4ies, zc4les, zca, &
            zcbn, zcbs, zcbu, zcdn2m, zcdnr, zcfm2m, zchar, zchneu, zcoeff, zcons,  &
            zcons10, zcons11, zcons12, zcons13, zcons14, zcons15, zcons16, zcons17, &
            zcons18, zcons2, zcons23, zcons25, zcons3, zcons5, zcons6, zcons8,      &
            zcons9, zconvs, zcor, zcpd, zcvm3, zcvm4, zcvm5, zd, zda1, zdiagt,      &
            zdiagw, zdisc, zdisl, zdisq, zdisxt, zdivv, zdivv1, zdqdt, zdqtot, zds, &
            zdtdt, zdthv, zdu2, zdus1, zdus2, zdxmdt, zdxtdt, zdz,                  &
            zepcor, zepdu2, zepevap, zephum, zeps, zepsec, zepsr, zepsw, zepz0o,    &
            zepzzo, zes, zfac, zfox, zfrac, zfreec, zfux, zgam, zghabl, zh1, zh2,   &
            zh2m, zhexp, zhtq, zhuv, zkap, zkappa, zktest, zlam, zln1, zln2, zm1,   &
            zm2, zm4, zmix, zmonob, zmult1, zmult2, zmult3, zmult4, zmult5, zplmax, &
            zplmin, zq2m, zqddif, zqdp, zqlwi1, zqlwi2, zqmitte, zqnlev, zqs1,      &
            zqs2, zqsmit, zqtmit, zqwevap, zrat, zrd, zrdrv, zred, zrh2m,           &
            zrhos, zrsi, zrvrd, zsdep1, zsdep2, zsh, zshear, zshn, zsm, zsmn,       &
            zsoil, zspeed, zsrfl, zstabf, zt2, ztau, ztaux, ztauy, zteldif,         &
            ztemitte, ztest, ztkemin, ztkesq, ztmelt, ztmit, ztmst, ztnlev,         &
            ztpfac1, ztpfac2, ztpfac3, ztpfac4, zu10, zustf, zusus1, zv10, zva,     &
            zvabc, zvb, zvbc, zvc, zvirmitte, zvk, zvkc, zvklt, zvrad, zvxmklt,     &
            zvxpklt, zwcrit, zwlmax, zwpwp, zwslev, zwstf, zwstop, zxhfl, zxnlev,   &
            zz2geo, zzb, zzcpts, zzqs, zzzh, zzzlam, zzzm
    REAL :: zepshr, zcons29, zcons30, zztvm

    INTEGER :: irow, jrow, itop, itopp1, jk, jl, jt
#ifndef NOLOOKUP
    INTEGER :: it
#endif
    LOGICAL :: lo, lo1

    !  Local arrays: 
    REAL :: z1mxtm1(klon), zbh(klon), zbm(klon), zbn(klon), zcair(klon),            &
            zccover(klon,klevm1), zcdn(klon), zcdum(klon,klev), zcfh(klon,klev),    &
            zcfm(klon,klev), zcfnc(klon), zcfnch(klon), zch(klon), zchn(klon),      &
            zcptgz(klon,klev), zcpts(klon), zcr(klon), zcsat(klon),                 &
            zdis(klon,klev), zdqs(klon), zebsh(klon,klev), zebsm(klon,klev),        &
            zedif(klon,klev), zfaxe(klon,klev), zfaxen(klon,klevm1), zhdyn(klon),   &
            zhh(klon,klevm1), zhsoil(klon), zhum(klon), zlteta1(klon,klev),         &
            zlwcmit(klon,klevm1), zqdif(klon,klev), zqmit(klon,klevm1), zqs(klon),  &
            zqss(klon,klev), zqssm(klon,klevm1), zri(klon), zricls(klon),           &
            zscf(klon), ztcoe(klon), ztdif(klon,klev), ztemit(klon,klevm1),         &
            ztess(klon), zteta1(klon,klev), ztkevn(klon,klev),                      &
            ztmitte(klon,klevm1), ztvir1(klon,klev), ztvirmit(klon,klevm1),         &
            ztvs(klon), zucf(klon), zudif(klon,klev), zustar(klon),                 &
            zvdif(klon,klev), zvidis(klon), zwet(klon), zwlmxi(klon), zwst(klon),   &
            zxdif(klon,klev), zxtdif(klon,klev,ktrac), zxtems(klon,ktrac)

    !
    ! ---------------------- INCLUDED FOR THE IMPLICIT SOIL SCHEME --------------------
    !
    !     The following variables are needed for the subroutines surftemp
    !     and evapotr4.

    REAL :: zetn(klon)       ! richtmyer-morton-coefficients
    REAL :: zftn(klon)       ! for dry static energy.
    REAL :: zeqn(klon)       ! richtmyer-morton-coefficients
    REAL :: zfqn(klon)       ! for moisture.
    REAL :: zcsurf(klon)     ! surface heat capacity
    REAL :: zgrhfx(klon)     ! ground heat flux

    REAL :: zcpq(klon)       ! spec. heat cap. of air as used in vdiff
    REAL :: zsold(klon)      ! old surface dry static energy
    REAL :: zqnatm(klon)     ! old atm. spec. humidity for old temp. ptm1
    REAL :: zqsold(klon)     ! surf. sat. spec. humidity for old temp.
    REAL :: zdqsol(klon)     ! derivative of zqsold

    REAL :: znetr(klon)      ! surface net radiation (old)
    REAL :: zcdrag(klon)     ! drag coefficient for heat and moisture

    ! tilde-tilde values of:
    REAL :: zsstt(klon)      ! surface dry static energy
    REAL :: tstiltil(klon)   ! surface temperature
    REAL :: zqstt(klon)      ! surf. sat. spec. humidity
    REAL :: zsatt(klon)      ! atm. dry static energy
    REAL :: zqatt(klon)      ! atm. spec. humidity

    REAL :: zsnew(klon)      ! new surface dry static energy
    REAL :: zqsnew(klon)     ! new surf. sat. spec. humidity

    REAL :: zeveff(klon)     ! evaporation efficiency (transpiration)
    REAL :: zrelhum(klon)    ! relative surface air humidity
    ! (bare soil evaporation)
    REAL :: zqsnfl(klon)     ! sublimation
    REAL :: zqlfl(klon)      ! interception loss
    REAL :: zqvfl(klon)      ! transpiration
    REAL :: zqgfl(klon)      ! bare soil evaporation

    ! ------------------------------------------------------------------------------

    INTEGER :: ihpbl(klon), ihpblc(klon), ihpbld(klon)

    !  Intrinsic functions 
#ifdef ECLIB
    REAL :: EXPHF,SQRTHF,ALOGHF
!DIR$ VFUNCTION EXPHF,SQRTHF
!DIR$ VFUNCTION ALOGHF
#define EXP(x)  EXPHF(x)
#define LOG(x)  ALOGHF(x)
#define SQRT  SQRTHF
#else
    INTRINSIC EXP, LOG, SQRT
#endif
    INTRINSIC ABS, COS, MAX, MERGE, MIN, SUM


    !  Executable statements 

#ifndef NOLOOKUP
    lookupoverflow = .FALSE.
#endif

    irow = nrow(1)
    jrow = nrow(2)

    ! Physical constants:
    !
    ! *zqwssat* are the inverses of critical values for
    ! soil water and snow depth that are used in the computation of the
    ! evapotranspiration's efficiency.

    zlam    = clam     ! asymptotic mixing length for momentum exchange
    zkap    = ckap     ! von karman constant
    zb      = cb       ! constants for the formulae about stability dependency
    zc      = cc       !   respectively near the neutral case, in the unstable
    zd      = cd       !   case and in the stable case
    zchar   = cchar    ! constant of the *charnock formula
    zva     = cva
    zvb     = cvb
    zvc     = cvc
    zvbc    = cvbc
    zvk     = cvk
    zvkc    = cvkc
    zvabc   = cvabc
    zvrad   = cvrad
    zwlmax  = cwlmax
    ztkemin = 1.E-4
    ztmelt  = tmelt
    zrvrd   = vtmpc1 + 1.
    zrdrv   = 1./zrvrd

    zcpd   = cpd
    zrd    = rd
    zkappa = zrd/zcpd
    zc3les = c3les
    zc3ies = c3ies
    zc4les = c4les
    zc4ies = c4ies

    ! Parameters for boundary layer diagnostics

    zhuv   = 10.*g
    zhtq   = 2.*g
    zephum = 5.E-2
    zrhos  = rhoh2o*1.025

    ! Security parameters.

    zepdu2 = 0.1     ! minimum squared wind increment in Ri calculation
    IF (lamip2) zepshr = 1.E-5
    zepzzo = 1.5E-05 ! minimum roughness length
    zepz0o = 2.
    zepcor = 5.E-05
    zepsw  = 1.E-3    ! minimum relative humidity of the ground
    zepsr  = 1.E-10   ! minimum value for the radiation in the
    !   visible part of the spectrum used to compute the
    !   canopy resistance.
    zepevap = 1.E-10 ! minimum atmospheric demand
    zepsec  = 1.E-2   ! minimum value for the drag coefficient

    ! Computational constants.

    ztmst = twodt
    IF (nstep == nstart) ztmst = 0.5*twodt
    zdiagt = 0.5*twodt
    zdiagw = zdiagt/rhoh2o

    ztpfac1 = cvdifts
    ztpfac2 = 1./ztpfac1
    ztpfac3 = 1. - ztpfac2
    ztpfac4 = 1. + ztpfac3

    zzzlam  = 30.
    zcons2  = 0.5*zkap/g
    zcons3  = zlam
    zcons5  = 3.*zb*zc*g**2
    zcons6  = 1./3.
    zcons8  = 2.*zb
    zcons9  = 3.*zb
    zcons10 = 1./cpd
    zcons11 = 3.*zb*zc
    zcons12 = ztpfac1*ztmst*g/rd
    zcons13 = 1./ztmst
    zcons14 = zchar*rd/(g**2*ztmst)
    zcons15 = 1./(g*ztmst)
    zcons16 = cpd*vtmpc2
    zcons18 = ztpfac1*ztmst*g**2
    zcons17 = 1./zkap**2
    zcons25 = zcons2/zcons3
    IF (lamip2) THEN
       zcons29 = ztpfac1*ztmst
       zcons30 = 1./(ztpfac1*ztmst*g)
    END IF

    zplmax = 0.75
    zplmin = 0.35
    zchneu = .3
    IF (lamip2) THEN
       zfreec = 0.001
    ELSE 
       zfreec = 0.0016  ! Constant for free convection
    END IF
    zgam = 1.25      ! exponent for the interpolation between free convection
    z1dgam = 1./zgam !   and neutral over sea

    ! Neutral stability functions (Mellor/Yamada, 1982)

    zh1  = 2.22
    zh2  = 0.22
    zm1  = 1.24
    zm2  = 2.37
    zm4  = 3.69
    zshn = zh1*zh2*SQRT(2.)
    zsmn = zshn*zm1*zm2/zm4
    IF (lamip2) THEN
       zda1  = 1./zsmn**3
       zustf = 1./zsmn**2
       zwstf = 0.2
    ELSE
       zda1  = 15.0
       zustf = 3.75
       zwstf = 0.2
    END IF

    itop = 1
    itopp1 = itop + 1

    !-- 1. New thermodynamic variable and boundary conditions

    !-- 1.1 Replace t by cp(q)*t+gz in the atmosphere

    DO jk = ktdia, klev
       DO jl = kidia, kfdia
          zcptgz(jl,jk) = pgeom1(jl,jk) + ptm1(jl,jk)*cpd*(1.+vtmpc2*pqm1(jl,jk))
          zteta1(jl,jk) = ptm1(jl,jk)*(100000./papm1(jl,jk))**zkappa
          ztvir1(jl,jk) = zteta1(jl,jk)*(1.+vtmpc1*pqm1(jl,jk)-pxm1(jl,jk))

          lo = ptm1(jl,jk) >= ztmelt
          zfaxe(jl,jk) = MERGE(alv,als,lo)
          zbet = zfaxe(jl,jk)/zcpd
          zusus1 = zbet*zteta1(jl,jk)/ptm1(jl,jk)*pxm1(jl,jk)
          zlteta1(jl,jk) = zteta1(jl,jk) - zusus1

#ifndef NOLOOKUP
          it = INT(ptm1(jl,jk)*1000.)
          IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
          IF (it<jptlucu1 .OR. it>jptlucu2) THEN
             WRITE(*,*) 'vdiff  : 1 it=',it,jl,jk
          END IF
#endif
          it  = MAX(MIN(it,jptlucu2),jptlucu1)
          zes = tlucua(it)/papm1(jl,jk)
#else
          lo  = ptm1(jl,jk) > tmelt
          zes = c2es*EXP(MERGE(c3les,c3ies,lo)*(ptm1(jl,jk)-tmelt) &
               / (ptm1(jl,jk)-MERGE(c4les,c4ies,lo)))/papm1(jl,jk)
#endif
          zes  = MIN(zes,0.5)
          zcor = 1./(1.-vtmpc1*zes)
          zqss(jl,jk) = zes*zcor
       END DO
    END DO
#ifndef NOLOOKUP
    IF (lookupoverflow) CALL lookuperror ('vdiff (1)   ')
#endif

    IF (lamip2) THEN
       DO jl = kidia,kfdia
          zqnatm(jl) = pqm1(jl,klev)
       END DO
    END IF

    DO jk = ktdia, klevm1
       DO jl = kidia, kfdia
          zhh(jl,jk) = (pgeom1(jl,jk)-pgeom1(jl,jk+1))
          zsdep1 = (paphm1(jl,jk)-paphm1(jl,jk+1))/(paphm1(jl,jk)-paphm1(jl,jk+2))
          zsdep2 = (paphm1(jl,jk+1)-paphm1(jl,jk+2))/(paphm1(jl,jk)-paphm1(jl,jk+2))
          zqssm(jl,jk)    = zsdep1*zqss(jl,jk)   + zsdep2*zqss(jl,jk+1)
          ztmitte(jl,jk)  = zsdep1*ptm1(jl,jk)   + zsdep2*ptm1(jl,jk+1)
          ztvirmit(jl,jk) = zsdep1*ztvir1(jl,jk) + zsdep2*ztvir1(jl,jk+1)
          zfaxen(jl,jk)   = zsdep1*zfaxe(jl,jk)  + zsdep2*zfaxe(jl,jk+1)
          zlwcmit(jl,jk)  = zsdep1*pxm1(jl,jk)   + zsdep2*pxm1(jl,jk+1)
          zqmit(jl,jk)    = zsdep1*pqm1(jl,jk)   + zsdep2*pqm1(jl,jk+1)
          ztemit(jl,jk)   = zsdep1*zteta1(jl,jk) + zsdep2*zteta1(jl,jk+1)
          zccover(jl,jk)  = paclcm(jl,jk)*zsdep1 + paclcm(jl,jk+1)*zsdep2
       END DO
    END DO

    ! Compute fractional surface coverages

    DO jl = kidia, kfdia
       IF (lamip2) THEN
          pcvs(jl) = 1. - EXP(-psnm1m(jl)*cqsncr)
       ELSE
          pcvs(jl) = MIN(1.,psnm1m(jl)*cqsncr)
       END IF
       pwlmx(jl)  = zwlmax*((1.-pvgrat(jl))+pvgrat(jl)*pvltm(jl))
       zwlmxi(jl) = 1./pwlmx(jl)
       pcvw(jl)   = pwlm1m(jl)*zwlmxi(jl)
    END DO

    !-- 1.2 Saturation parameters,
    !       relative humidity over the bare land part
    !       and virtual temperature at the surface.

    DO jl = kidia, kfdia
       lo = (ptsm1m(jl)-tmelt) > 0.
       zcvm3 = MERGE(c3les,c3ies,lo)
       zcvm4 = MERGE(c4les,c4ies,lo)
       zcvm5 = MERGE(c5les,c5ies,lo)
#ifndef NOLOOKUP
       it = INT(ptsm1m(jl)*1000.)
       IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
       IF (it<jptlucu1 .OR. it>jptlucu2) THEN
          WRITE(*,*) 'vdiff  : 2 it=',it,jl
       END IF
#endif
       it  = MAX(MIN(it,jptlucu2),jptlucu1)
       zes = tlucua(it)/paphm1(jl,klevp1)
#else
       zes = c2es*EXP(zcvm3*(ptsm1m(jl)-tmelt) &
            / (ptsm1m(jl)-zcvm4))/paphm1(jl,klevp1)
#endif
       zcor = 1./(1.-vtmpc1*zes)
       zqs(jl) = zes*zcor
       zdqs(jl) = zqs(jl)*zcvm5*zcor*(1./(ptsm1m(jl)-zcvm4))**2
       pwsm1m(jl) = MIN(pwsm1m(jl),pwsmxm(jl))
       zwstop = MIN(0.1,pwsmxm(jl))
       zwslev = pwsmxm(jl) - zwstop
       IF (pwsm1m(jl)>zwslev) THEN
          zhum(jl) = 0.5*(1.-COS((pwsm1m(jl)-zwslev)*api/zwstop))
       ELSE
          zhum(jl) = 0.
       END IF
       zhsoil(jl) = pcvs(jl) + (1.-pcvs(jl))*(pcvw(jl)+(1.-pcvw(jl))*zhum(jl))
       zhsoil(jl) = MERGE(zhsoil(jl),1.,laland(jl))
       lo = pqm1(jl,klev) > zqs(jl)
       zhsoil(jl) = MERGE(1.,zhsoil(jl),lo)
       ztess(jl)  = ptsm1m(jl)*(1.E5/paphm1(jl,klevp1))**zkappa
       ztvs(jl)   = ztess(jl)*(1.+vtmpc1*zhsoil(jl)*zqs(jl))
    END DO
#ifndef NOLOOKUP
    IF (lookupoverflow) CALL lookuperror ('vdiff (2)   ')
#endif

    !-- 1.3 Definition of the stomatal resistance

    DO jl = kidia, kfdia
       zwcrit   = zplmax*pwsmxm(jl)
       zwpwp    = zplmin*pwsmxm(jl)
       zqwevap  = 1./(zwcrit-zwpwp)
       zsoil    = MAX(zepsw,MIN(1.,(pwsm1m(jl)-zwpwp)*zqwevap))
       zsrfl    = MAX(zepsr,psrfl(jl)*zvrad)
       zabcs    = (zva+zvbc)/(zvc*zsrfl)
       zvklt    = zvk*pvltm(jl)
       zvxpklt  = EXP(zvklt)
       zvxmklt  = EXP(-zvklt)
       zln1     = LOG((zabcs*zvxpklt+1.)/(zabcs+1.))
       zln2     = LOG((zabcs+zvxmklt)/(zabcs+1.))
       zrsi     = (zvb*zln1/zvabc-zln2)/zvkc
       zwet(jl) = 1./(zrsi*zsoil)
       lo       = pqm1(jl,klev) > zqs(jl)
       zwet(jl) = MERGE(0.,zwet(jl),lo)
    END DO

    IF (lvdiff) THEN

       !-- 2. Computation of the exchange coefficients

       ! The surface layer is now computed before the other levels

       !-- 2.1 Computation of basic quantities: Wind shear,
       !       Richardson number,squared mixing lengths, unstable
       !       and stable case common factors and neutral case
       !       common part of the drag coefficients.

       DO jl = kidia, kfdia
          zdu2      = MAX(zepdu2,pum1(jl,klev)**2+pvm1(jl,klev)**2)
          zqmitte   = (pqm1(jl,klev)+zqs(jl)*zhsoil(jl))/2.
          zqtmit    = pxm1(jl,klev)*0.5 + zqmitte
          ztmit     = (ptm1(jl,klev)+ptsm1m(jl))/2.
          ztemitte  = (zteta1(jl,klev)+ztess(jl))/2.
          zvirmitte = (ztvir1(jl,klev)+ztvs(jl))/2.
          zqsmit    = (zqss(jl,klev)+zqs(jl))/2.
          zqlwi1    = pqm1(jl,klev) + pxm1(jl,klev)
          zqlwi2    = zqs(jl)*zhsoil(jl)
          zfux      = zfaxe(jl,klev)/(zcpd*ztmit)
          zfox      = zfaxe(jl,klev)/(zrd*ztmit)
          zmult1    = 1. + vtmpc1*zqtmit
          zmult2    = zfux*zmult1 - zrvrd
          zmult3    = zrdrv*zfox*zqsmit/(1.+zrdrv*zfox*zfux*zqsmit)
          zmult5    = zmult1 - zmult2*zmult3
          zmult4    = zfux*zmult5 - 1.

          zdus1   = paclcm(jl,klev)*zmult5 + (1.-paclcm(jl,klev))*zmult1
          zdus2   = paclcm(jl,klev)*zmult4 + (1.-paclcm(jl,klev))*vtmpc1
          zteldif = zlteta1(jl,klev) - ztess(jl)
          zqddif  = zqlwi1 - zqlwi2
          zbuoy   = zdus1*zteldif + zdus2*ztemitte*zqddif

          zri(jl)    = pgeom1(jl,klev)*zbuoy/(zvirmitte*zdu2)
          zricls(jl) = zri(jl)
          zcdn(jl)   = (zkap/LOG(1.+pgeom1(jl,klev)/(g*paz0m(jl))))**2
          z0h        = paz0m(jl)*EXP(2.-86.276*paz0m(jl)**0.375)
          zalo       = LOG(1.+pgeom1(jl,klev)/(g*paz0m(jl)))
          zaloh      = LOG(1.+pgeom1(jl,klev)/(g*z0h))
          zchn(jl)   = zkap**2/(zalo*zaloh)
          zucf(jl)   = 1./(1.+zcons11*zcdn(jl)*SQRT(ABS(zri(jl))* &
                          (1.+pgeom1(jl,klev)/(g*paz0m(jl)))))
          IF (lamip2) THEN
             zscf(jl) = SQRT(1.+ABS(zri(jl)))
          ELSE
             zscf(jl) = SQRT(1.+zd*ABS(zri(jl)))
          END IF
          zcons      = zcons12*paphm1(jl,klevp1)/(ptm1(jl,klev)* &
                      (1.+vtmpc1*pqm1(jl,klev)-pxm1(jl,klev)))
          zcfnc(jl)  = zcons*SQRT(zdu2)*zcdn(jl)
          zcfnch(jl) = zcons*SQRT(zdu2)*zchn(jl)
          zdthv      = MAX(0.,(ztvs(jl)-ztvir1(jl,klev)))
          zwst(jl)   = zdthv*SQRT(zdu2)/zvirmitte
          zcr(jl)    = (zfreec/(zchn(jl)*SQRT(zdu2)))*ABS(zbuoy)**zcons6

          ! Canopy resistance

          zwet(jl) = zwet(jl)/zcons
       END DO

       !-- 2.2 Dimensionless heat transfer coefficients multiplied
       !       by pressure thicknesses for momentum and heat exchange

       DO jl = kidia, kfdia
          IF (zri(jl)>=0.) THEN
             zcfm(jl,klev) = zcfnc(jl)/(1.+zcons8*zri(jl)/zscf(jl))
             IF (laland(jl)) THEN
                IF (lamip2) THEN
                   zcfh(jl,klev) = zcfnc(jl)/(1.+zcons8*zri(jl)*zscf(jl))
                ELSE
                   zcfh(jl,klev) = zcfnc(jl)/(1.+zcons9*zri(jl)*zscf(jl))
                END IF
                zch(jl)       = zcfh(jl,klev)/zcfnc(jl)*zcdn(jl)
             ELSE
                IF (lamip2) THEN
                   zcfh(jl,klev) = zcfnch(jl)/(1.+zcons8*zri(jl)*zscf(jl))
                ELSE
                   zcfh(jl,klev) = zcfnch(jl)/(1.+zcons9*zri(jl)*zscf(jl))
                END IF
                zch(jl)       = zcfh(jl,klev)/zcfnch(jl)*zchn(jl)
             END IF
          ELSE
             zcfm(jl,klev) = zcfnc(jl)*(1.-zcons8*zri(jl)*zucf(jl))
             IF (laland(jl)) THEN
                zcfh(jl,klev) = zcfnc(jl)*(1.-zcons9*zri(jl)*zucf(jl))
                zch(jl)       = zcfh(jl,klev)/zcfnc(jl)*zcdn(jl)
             ELSE

                ! Special free convection limit over sea

                zcfh(jl,klev) = zcfnch(jl)*((1.+zcr(jl)**zgam)**z1dgam)
                zch(jl)       = zcfh(jl,klev)/zcfnch(jl)*zchn(jl)
             END IF
          END IF

          zcdum(jl,klev) = zcfm(jl,klev)

          ! Interpolationfunctions for diagnostics

          zbn(jl) = zkap/SQRT(zcdn(jl))
          zbm(jl) = MAX(zepsec,SQRT(zcfm(jl,klev)*zcdn(jl)*zcons17/zcfnc(jl)))
          zbh(jl) = MAX(zepsec,zch(jl)/zbm(jl)*zcons17)
          zbm(jl) = 1./zbm(jl)
          zbh(jl) = 1./zbh(jl)
       END DO

       IF (lamip2) THEN
          ! Set surface heat capacity and ground heat flux

          DO jl = kidia,kfdia
             zcsurf(jl) = pcapcal(jl,jrow)
             zgrhfx(jl) = pfluxgrd(jl,jrow)
          END DO
       END IF

       ! Initialize surface emission for tracers

       IF (ktrac > 0) THEN
          DO jt = 1, ktrac
             DO jl = kidia, kfdia
                zxtems(jl,jt) = 0.
             END DO
          END DO

          ! Surface emissions and dry deposition

          IF (lemis) THEN
             DO jl = kidia, kfdia
                z1mxtm1(jl) = papm1(jl,klev) &
                            / (ptm1(jl,klev)*rd*(1.+vtmpc1*pqm1(jl,klev)))
             END DO

             CALL xtemiss ( klon,   klev,     irow,    cvdifts,  dtime, &
                            pxtm1,  zxtems,   z1mxtm1,                  &
                            laland, pforestm, psnm1m )
          END IF
       END IF

       !-- 2.3 Equivalent evapotranspiration efficiency coefficient

       DO jl = kidia, kfdia
          IF (lamip2) THEN    
             zeveff(jl) = 1./(1.+zcfh(jl,klev)*zwet(jl))
             zeveff(jl) = MERGE(zeveff(jl),1.,laland(jl))
             ! Relative humidity at the surface
             zrelhum(jl) = zhum(jl)
          END IF
          zwet(jl)  = pcvs(jl) + (1.-pcvs(jl))*(pcvw(jl)+(1.-pcvw(jl))/ &
                                 (1.+zcfh(jl,klev)*zwet(jl)))
          zwet(jl)  = MERGE(zwet(jl),1.,laland(jl))
          lo        = zhum(jl) <= pqm1(jl,klev)/zqs(jl)
          zcsat(jl) = pcvs(jl) + (1.-pcvs(jl))*(pcvw(jl)+(1.-pcvw(jl))* &
                                 MERGE(0.,zhum(jl),lo))
          zcair(jl) = pcvs(jl) + (1.-pcvs(jl))*(pcvw(jl)+(1.-pcvw(jl))* &
                                 MERGE(0.,1.,lo))
          zcsat(jl) = MERGE(zcsat(jl),1.,laland(jl))
          zcair(jl) = MERGE(zcair(jl),1.,laland(jl))
          lo        = pqm1(jl,klev) > zqs(jl)
          zcsat(jl) = MERGE(1.,zcsat(jl),lo)
          zcair(jl) = MERGE(1.,zcair(jl),lo)
          zcsat(jl) = pvgrat(jl)*zwet(jl) + (1.-pvgrat(jl))*zcsat(jl)
          zcair(jl) = pvgrat(jl)*zwet(jl) + (1.-pvgrat(jl))*zcair(jl)
          zcpts(jl) = ptsm1m(jl)*cpd*(1.+vtmpc2*(zcsat(jl)*zqs(jl)+ &
                                 (1.-zcair(jl))*pqm1(jl,klev)))
          ! INCLUDED FOR THE IMPLICIT SOIL SCHEME --------
          IF (lamip2) zcpq(jl) = zcpts(jl)/ptsm1m(jl)
       END DO

       !-- 2.4 Computation of the pbl extension

       DO jl = kidia, kfdia
          zdu2       = MAX(zepdu2,pum1(jl,klev)**2+pvm1(jl,klev)**2)
          zcor       = MAX(ABS(coriol(irow)),zepcor)
          lo         = paz0m(jl) > zepz0o
          zcdn2m     = MERGE((zkap/LOG(1.+pgeom1(jl,klev)/(g*zepz0o)))**2,zcdn(jl),lo)
          zcdnr      = zcdn2m/zcdn(jl)
          IF (lamip2) THEN
             ! ECHAM4 bug correction (found by Erik Van Meijgaard)
             zcfm2m  = MERGE(zcfnc(jl)*zcdnr*(1.-zcons8*zri(jl)/ &
                         (1.+zcons11*zcdn2m*SQRT(ABS(zri(jl))* &
                         (1.+pgeom1(jl,klev)/(g*zepz0o))))),zcfm(jl,klev)*zcdnr, &
                          lo .AND. zri(jl)<0.)
          ELSE
             zcfm2m  = MERGE(zcfnc(jl)*zcdnr*(1.-zcons8*zri(jl))/ &
                         (1.+zcons11*zcdn2m*SQRT(ABS(zri(jl))* &
                         (1.+pgeom1(jl,klev)/(g*zepz0o)))),zcfm(jl,klev)*zcdnr, &
                          lo .AND. zri(jl)<0.)
          END IF
          zustar(jl) = SQRT(zcfm2m*SQRT(zdu2)*ptm1(jl,klev)*(1.+vtmpc1* &
                            pqm1(jl,klev)-pxm1(jl,klev))/(zcons12*paphm1(jl,klevp1)))
          zhdyn(jl)  = MIN(pgeom1(jl,1)/g,zchneu*zustar(jl)/zcor)

          ihpblc(jl) = klev
          ihpbld(jl) = klev
       END DO

       DO jk = klevm1, 1, -1
          DO jl = kidia, kfdia
             zds = zcptgz(jl,jk) - zcptgz(jl,klev)
             zdz = pgeom1(jl,jk)/g - zhdyn(jl)
             ihpblc(jl) = MERGE(jk,ihpblc(jl),ihpblc(jl)==klev .AND. zds>0.)
             ihpbld(jl) = MERGE(jk,ihpbld(jl),ihpbld(jl)==klev .AND. zdz>=0.)
          END DO
       END DO

       ! Convective velocity scale, Monin-Obukhov length and
       ! tke boundary condition (Mailhot/Benoit, 1982)

       DO jl = kidia, kfdia
          ihpbl(jl) = MIN(ihpblc(jl),ihpbld(jl),klev-3)
          zghabl    = MIN(50000.,pgeom1(jl,ihpbl(jl)))
          IF (zwst(jl)>zepsr) THEN
             zconvs = (zwst(jl)*zch(jl)*zghabl)**zcons6
             zmonob = (zustar(jl)**3)/(zkap*g*zwst(jl)*zch(jl))
             zstabf = (pgeom1(jl,klev)/(g*zmonob))**(zcons6*2.)
             zstabf = MIN(zustf*3.,zstabf)
          ELSE
             zconvs = 0.
             zstabf = 0.
          END IF
          ztkevn(jl,klev) = (zustf+zstabf)*(zustar(jl)**2) + zwstf*(zconvs**2)
          ztkevn(jl,klev) = MAX(ztkemin,ztkevn(jl,klev))
       END DO

       IF (nstep==nstart) THEN
          DO jl = kidia, kfdia
             ptkem1m(jl,klev) = ztkevn(jl,klev)
             ptkem(jl,klev)   = ztkevn(jl,klev)
          END DO
       END IF

       !-- 2.5 Vertical loop

       DO jk = ktdia, klevm1

          !-- 2.6 Computation of basic quantities: Wind shear,
          !       buoyancy, Richardson number, mixing lengths.

          ! Modified Richardson number in vertical loop

!DIR$ IVDEP
!OCL NOVREC
          DO jl = kidia, kfdia
             zqtmit = zlwcmit(jl,jk) + zqmit(jl,jk)
             zfux   = zfaxen(jl,jk)/(zcpd*ztmitte(jl,jk))
             zfox   = zfaxen(jl,jk)/(zrd*ztmitte(jl,jk))
             zmult1 = 1. + vtmpc1*zqtmit
             zmult2 = zfux*zmult1 - zrvrd
             zmult3 = zrdrv*zfox*zqssm(jl,jk)/(1.+zrdrv*zfux*zfox*zqssm(jl,jk))
             zmult5 = zmult1 - zmult2*zmult3
             zmult4 = zfux*zmult5 - 1.

             zdus1   = zccover(jl,jk)*zmult5 + (1.-zccover(jl,jk))*zmult1
             zdus2   = zccover(jl,jk)*zmult4 + (1.-zccover(jl,jk))*vtmpc1
             zteldif = (zlteta1(jl,jk)-zlteta1(jl,jk+1))/zhh(jl,jk)*g
             zdqtot  = (pqm1(jl,jk)+pxm1(jl,jk)) - (pqm1(jl,jk+1)+pxm1(jl,jk+1))
             zqddif  = zdqtot/zhh(jl,jk)*g
             zbuoy   = (zteldif*zdus1+ztemit(jl,jk)*zdus2*zqddif)*g/ztvirmit(jl,jk)
             zdivv   = (pum1(jl,jk)-pum1(jl,jk+1))**2
             zdivv1  = (pvm1(jl,jk)-pvm1(jl,jk+1))**2

             IF (lamip2) THEN
                zshear  = (zdivv+zdivv1)*(g/zhh(jl,jk))**2
                zri(jl) = zbuoy/MAX(zshear,zepshr)
             ELSE
                zshear  = MAX(zepdu2,(zdivv+zdivv1))*(g/zhh(jl,jk))**2
                zri(jl) = zbuoy/zshear
             END IF

             ! Asymptotic mixing length for momentum and
             ! heat (zlam) above the pbl as a function of height
             ! according to Holtslag and Boville (1992), j. climate.

             zhexp = EXP(1.-pgeom1(jl,jk)/pgeom1(jl,ihpbl(jl)))
             zlam  = zzzlam + (zcons3-zzzlam)*zhexp
             IF (jk >= ihpbl(jl)) THEN
                zcons23 = zcons25
             ELSE
                zcons23 = zcons2/zlam
             END IF

             ! Mixing length (blackadar) + stability dependent function

             z2geomf = pgeom1(jl,jk) + pgeom1(jl,jk+1)
             zz2geo  = zcons2*z2geomf
             zmix    = zz2geo/(1.+zcons23*z2geomf)

             ! Stability functions (Louis, 1979)

             zalh2 = zmix*zmix
             zucf(jl) = 1./(1.+zcons5*zalh2*SQRT(ABS(zri(jl))*(((pgeom1(jl,jk)/ &
                        pgeom1(jl,jk+1))**zcons6-1.)/(pgeom1(jl,jk)-pgeom1(jl,jk+1)))**3/ &
                        (pgeom1(jl,jk+1))))
             IF (zri(jl)<0.) THEN
                zsh = zshn*(1.-zcons9*zri(jl)*zucf(jl))
                zsm = zsmn*(1.-zcons8*zri(jl)*zucf(jl))
             ELSE
                zsh = zshn/(1.+zcons8*zri(jl)*SQRT(1.+zri(jl)))
                IF (lamip2) THEN
                   zsm = zsmn/(1.+zcons8*zri(jl)/SQRT(1.+zri(jl)))
                ELSE
                   zsm = zsmn/(1.+zcons8*zri(jl)/SQRT(1.+zd*zri(jl)))
                END IF
             END IF

             !-- 2.7 Dimensionless coefficients multiplied by pressure
             !       thicknesses for momentum and heat exchange.

             zzb    = zshear*zmix*zsm - zbuoy*zmix*zsh
             zdisl  = zda1*zmix/ztmst
             zktest = 1. + (zzb*ztmst+SQRT(ptkem1m(jl,jk))*2.)/zdisl
             IF (zktest<=1.) THEN
                ztkevn(jl,jk) = ztkemin
             ELSE
                ztkevn(jl,jk) = MAX(ztkemin,(zdisl*(SQRT(zktest)-1.))**2)
             END IF
             IF (nstep==nstart) THEN
                ptkem1m(jl,jk) = ztkevn(jl,jk)
                ptkem(jl,jk)   = ztkevn(jl,jk)
             END IF
             ztkesq = SQRT(MAX(ztkemin,ptkem1m(jl,jk)))

             ! ------  Diffusion coefficients:  --------
             zzzm = zmix*zsm*ztkesq
             zzzh = zmix*zsh*ztkesq

             IF (lamip2) THEN
                zztvm = (ptvm1(jl,jk)+ptvm1(jl,jk+1))*0.5
                zalf  = paphm1(jl,jk+1)/(zztvm*zhh(jl,jk)*zrd)
             ELSE
                zalf  = paphm1(jl,jk+1)/(ztvirmit(jl,jk)*zhh(jl,jk)*zrd)
             END IF

             zcfm(jl,jk)  = zzzm*zcons18*zalf
             zcfh(jl,jk)  = zzzh*zcons18*zalf
             zcdum(jl,jk) = zcfm(jl,jk)/ztkesq*SQRT(ztkevn(jl,jk))
          END DO
       END DO

       !-- 2.8 Diffusion implicit computations for tke

       DO jk = ktdia, klev
          DO jl = kidia, kfdia
             zedif(jl,jk) = ztpfac2*ztkevn(jl,jk)
          END DO
       END DO

       DO jl = kidia, kfdia
          ztcoe(jl) = (zcdum(jl,itop)+zcdum(jl,itopp1))*0.5
          zqdp  = 1./(papm1(jl,itopp1)-papm1(jl,itop))
          zdisc = 1./(1.+(zcdum(jl,itop)+zcdum(jl,itopp1))*0.5*zqdp)
          zebsm(jl,itop) = zdisc*(zcdum(jl,itop)+zcdum(jl,itopp1))*0.5*zqdp
          zedif(jl,itop) = zdisc*zedif(jl,itop)
       END DO

       DO jk = itopp1, klev - 2
          DO jl = kidia, kfdia
             zqdp = 1./(papm1(jl,jk+1)-papm1(jl,jk))
             zfac = ztcoe(jl)*zqdp
             ztcoe(jl) = (zcdum(jl,jk+1)+zcdum(jl,jk))*0.5
             zdisc = 1./(1.+zfac*(1.-zebsm(jl,jk-1))+(zcdum(jl,jk+1)+ &
                     zcdum(jl,jk))*0.5*zqdp)
             zebsm(jl,jk) = zdisc*(zcdum(jl,jk+1)+zcdum(jl,jk))*0.5*zqdp
             zedif(jl,jk) = zdisc*(zedif(jl,jk)+zfac*zedif(jl,jk-1))
          END DO
       END DO

       DO jl = kidia, kfdia
          zqdp = 1./(papm1(jl,klev)-papm1(jl,klevm1))
          zfac = ztcoe(jl)*zqdp
          ztcoe(jl) = (zcdum(jl,klev)+zcdum(jl,klevm1))*0.5
          zdisc = 1./(1.+zfac*(1.-zebsm(jl,klev-2))+(zcdum(jl,klev)+ &
                  zcdum(jl,klevm1))*0.5*zqdp)
          zedif(jl,klevm1) = zdisc*((zcdum(jl,klev)+zcdum(jl, &
                             klevm1))*0.5*zqdp*zedif(jl,klev)+zedif(jl,klevm1)+ &
                             zfac*zedif(jl,klev-2))
       END DO

       DO jk = klev - 2, itop, -1
          DO jl = kidia, kfdia
             zedif(jl,jk) = zedif(jl,jk) + zebsm(jl,jk)*zedif(jl,jk+1)
          END DO
       END DO

       ! Time integration of turbulent kinetic energy and check

       DO jk = itop, klev
          ztest = 0.
          DO jl = kidia, kfdia
             ptke(jl,jk) = zedif(jl,jk) + ztpfac3*ztkevn(jl,jk)
          END DO
!OCL NOVREC
          DO jl = kidia, kfdia
             ztest = ztest + MERGE(1.,0.,ptke(jl,jk)<0.)
          END DO
          IF (ABS(ztest) > 0.) THEN
             WRITE (nout,*) 'vdiff: tke is negative'
             CALL finish('vdiff','Run terminated.')
          END IF
       END DO

       ! Time filter for turbulent kinetic energy

       IF (nstep/=nstart) THEN
          zeps = eps
       ELSE
          zeps = 0.
       END IF
       DO jk = ktdia, klev
          DO jl = kidia, kfdia
             ptkem1(jl,jk) = ptkem(jl,jk) + zeps*(ptkem1m(jl,jk)- &
                             2.*ptkem(jl,jk)+ptke(jl,jk))
          END DO
       END DO

       !-- 3. Diffusion implicit computations for momentum

       !-- 3.1 Setting of right hand sides

       DO jk = itop, klev
          DO jl = kidia, kfdia
             zudif(jl,jk) = ztpfac2*pum1(jl,jk)
             zvdif(jl,jk) = ztpfac2*pvm1(jl,jk)
          END DO
       END DO

       !-- 3.2 Top layer elimination

       DO jl = kidia, kfdia
          zqdp  = 1./(paphm1(jl,itopp1)-paphm1(jl,itop))
          zdisc = 1./(1.+zcfm(jl,itop)*zqdp)
          zebsm(jl,itop) = zdisc*(zcfm(jl,itop)*zqdp)
          zudif(jl,itop) = zdisc*zudif(jl,itop)
          zvdif(jl,itop) = zdisc*zvdif(jl,itop)
       END DO

       !-- 3.3 Elimination for middle layers

       DO jk = itopp1, klevm1
          DO jl = kidia, kfdia
             zqdp  = 1./(paphm1(jl,jk+1)-paphm1(jl,jk))
             zfac  = zcfm(jl,jk-1)*zqdp
             zdisc = 1./(1.+zfac*(1.-zebsm(jl,jk-1))+zcfm(jl,jk)*zqdp)
             zebsm(jl,jk) = zdisc*(zcfm(jl,jk)*zqdp)
             zudif(jl,jk) = zdisc*(zudif(jl,jk)+zfac*zudif(jl,jk-1))
             zvdif(jl,jk) = zdisc*(zvdif(jl,jk)+zfac*zvdif(jl,jk-1))
          END DO
       END DO

       !-- 3.4 Bottom layer elimination

       DO jl = kidia, kfdia
          zqdp      = 1./(paphm1(jl,klevp1)-paphm1(jl,klev))
          zfac      = zcfm(jl,klevm1)*zqdp
          ztcoe(jl) = zcfm(jl,klev)
          zdisc     = 1./(1.+zfac*(1.-zebsm(jl,klevm1))+zcfm(jl,klev)*zqdp)
          zudif(jl,klev) = zdisc*(zudif(jl,klev)+zfac*zudif(jl,klevm1))
          zvdif(jl,klev) = zdisc*(zvdif(jl,klev)+zfac*zvdif(jl,klevm1))
       END DO

       !-- 3.5 Back-substitution

       DO jk = klevm1, itop, -1
          DO jl = kidia, kfdia
             zudif(jl,jk) = zudif(jl,jk) + zebsm(jl,jk)*zudif(jl,jk+1)
             zvdif(jl,jk) = zvdif(jl,jk) + zebsm(jl,jk)*zvdif(jl,jk+1)
          END DO
       END DO

       !-- 3.6 Incrementation of u and v tendencies and storage of
       !       the dissipation.

       DO jl = kidia, kfdia
          zvidis(jl) = 0.
       END DO

       DO jk = itop, klev
          DO jl = kidia, kfdia
             pvom(jl,jk) = pvom(jl,jk) + (zudif(jl,jk)-ztpfac2*pum1(jl,jk))*zcons13
             pvol(jl,jk) = pvol(jl,jk) + (zvdif(jl,jk)-ztpfac2*pvm1(jl,jk))*zcons13
             zdis(jl,jk) = 0.5*((ztpfac2*pum1(jl,jk)-zudif(jl,jk))      &
                               *(ztpfac4*pum1(jl,jk)+zudif(jl,jk))      &
                               +(ztpfac2*pvm1(jl,jk)-zvdif(jl,jk))      &
                               *(ztpfac4*pvm1(jl,jk)+zvdif(jl,jk)))
             zvidis(jl)  = zvidis(jl) + zdis(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
          END DO
       END DO

       !-- 3.7 Updating of z0 for open sea

       DO jl = kidia, kfdia
          paz0(jl) = MERGE(paz0m(jl),MAX(zcons14*ztcoe(jl)*SQRT(zudif(jl,klev)** &
                     2+zvdif(jl,klev)**2)*ptm1(jl,klev)*(1.+vtmpc1*pqm1(jl,klev)- &
                     pxm1(jl,klev))/paphm1(jl,klevp1),zepzzo),laland(jl))
          lo = ( .NOT. laland(jl)) .AND. (pseaice(jl)>0.5)
          paz0(jl) = MERGE(cz0ice,paz0(jl),lo)
          ztaux = zcons15*ztcoe(jl)*zudif(jl,klev)
          ztauy = zcons15*ztcoe(jl)*zvdif(jl,klev)
          pustr(jl)   = pustrm(jl) + zdiagt*ztaux
          pvstr(jl)   = pvstrm(jl) + zdiagt*ztauy
          ztau        = SQRT(ztaux**2+ztauy**2)
          pustar3(jl) = pustar3m(jl) + zdiagt*SQRT(ztau/zrhos)**3
          pvdis(jl)   = pvdism(jl) + zdiagt*zcons15*zvidis(jl)
!ik ocean coupling
          IF (lcouple) THEN
             ustr_local(jl,jrow) = ustr_local(jl,jrow) + zdiagt*ztaux
             vstr_local(jl,jrow) = vstr_local(jl,jrow) + zdiagt*ztauy
          END IF

       END DO

       zbud = budw(irow)
       dvdisz(irow) = zdiagt*zcons15*zbud*SUM(zvidis(1:klon))

       !-- 4. Diffusion implicit computations for heat (s.+l.)

       DO jk = 1, klev
          DO jl = kidia, kfdia
             ztdif(jl,jk) = 0.
             zqdif(jl,jk) = 0.
             zxdif(jl,jk) = 0.
          END DO
       END DO

       DO jt = 1, ntrac
          DO jk = 1, klev
             DO jl = kidia, kfdia
                zxtdif(jl,jk,jt) = 0.
             END DO
          END DO
       END DO

       !-- 4.1 Setting of right hand sides

       DO jk = itop, klev
          DO jl = kidia, kfdia
             ztdif(jl,jk) = ztpfac2*zcptgz(jl,jk)
             zqdif(jl,jk) = ztpfac2*pqm1(jl,jk)
             zxdif(jl,jk) = ztpfac2*pxm1(jl,jk)
          END DO
       END DO

       DO jt = 1, ntrac
          DO jk = itop, klev
             DO jl = kidia, kfdia
                zxtdif(jl,jk,jt) = ztpfac2*pxtm1(jl,jk,jt)
             END DO
          END DO
       END DO

       !-- 4.2 Top layer elimination

       DO jl = kidia, kfdia
          zqdp  = 1./(paphm1(jl,itopp1)-paphm1(jl,itop))
          zdisc = 1./(1.+zcfh(jl,itop)*zqdp)
          zebsh(jl,itop) = zdisc*(zcfh(jl,itop)*zqdp)
          ztdif(jl,itop) = zdisc*ztdif(jl,itop)
          zqdif(jl,itop) = zdisc*zqdif(jl,itop)
          zxdif(jl,itop) = zdisc*zxdif(jl,itop)
       END DO

       DO jt = 1, ntrac
          DO jl = kidia, kfdia
             zqdp  = 1./(paphm1(jl,itopp1)-paphm1(jl,itop))
             zdisc = 1./(1.+zcfh(jl,itop)*zqdp)
             zxtdif(jl,itop,jt) = zdisc*zxtdif(jl,itop,jt)
          END DO
       END DO

       !-- 4.3 Elimination for middle layers

       DO jk = itopp1, klevm1
          DO jl = kidia, kfdia
             zqdp  = 1./(paphm1(jl,jk+1)-paphm1(jl,jk))
             zfac  = zcfh(jl,jk-1)*zqdp
             zdisc = 1./(1.+zfac*(1.-zebsh(jl,jk-1))+zcfh(jl,jk)*zqdp)
             zebsh(jl,jk) = zdisc*(zcfh(jl,jk)*zqdp)
             ztdif(jl,jk) = zdisc*(ztdif(jl,jk)+zfac*ztdif(jl,jk-1))
             zqdif(jl,jk) = zdisc*(zqdif(jl,jk)+zfac*zqdif(jl,jk-1))
             zxdif(jl,jk) = zdisc*(zxdif(jl,jk)+zfac*zxdif(jl,jk-1))
          END DO
       END DO

       DO jt = 1, ntrac
          DO jk = itopp1, klevm1
             DO jl = kidia, kfdia
                zqdp  = 1./(paphm1(jl,jk+1)-paphm1(jl,jk))
                zfac  = zcfh(jl,jk-1)*zqdp
                zdisc = 1./(1.+zfac*(1.-zebsh(jl,jk-1))+zcfh(jl,jk)*zqdp)
                zxtdif(jl,jk,jt) = zdisc*(zxtdif(jl,jk,jt)+zfac*zxtdif(jl,jk-1,jt))
             END DO
          END DO
       END DO

       !-- 4.4 Bottom layer elimination

       DO jl = kidia, kfdia
          zqdp      = 1./(paphm1(jl,klevp1)-paphm1(jl,klev))
          zfac      = zcfh(jl,klevm1)*zqdp
          ztcoe(jl) = zcfh(jl,klev)
          zdisc = 1./(1.+zfac*(1.-zebsh(jl,klevm1))+zcfh(jl,klev)*zqdp)
          zdisq = 1./(1.+zfac*(1.-zebsh(jl,klevm1))+zcair(jl)*zcfh(jl,klev)*zqdp)

          IF (.NOT. lamip2) THEN
             ztdif(jl,klev) = zdisc*(ztdif(jl,klev)+(zcfh(jl,klev)*zqdp) &
                              *ztpfac2*zcpts(jl)+zfac*ztdif(jl,klevm1))
             zqdif(jl,klev) = zdisq*(zqdif(jl,klev)+(zcsat(jl)*zcfh(jl,klev)*zqdp) &
                              *ztpfac2*zqs(jl)+zfac*zqdif(jl,klevm1))
          ENDIF
          zxdif(jl,klev) = zdisc*(zxdif(jl,klev)+zfac*zxdif(jl,klevm1))

          IF (lamip2) THEN
             ! Coefficients of the Richtmeyer-Morton-scheme
             !   xn = en * xs + fn
             ! with xn = atm. value of s or q, and  xs = surface value

             zetn(jl) = zdisc*zcfh(jl,klev)*zqdp
             zftn(jl) = zdisc*(ztdif(jl,klev)+zfac*ztdif(jl,klevm1))*ztpfac1

             zeqn(jl) = zdisq*zcsat(jl)*zcfh(jl,klev)*zqdp
             zfqn(jl) = zdisq*(zqdif(jl,klev)+zfac*zqdif(jl,klevm1))*ztpfac1
          END IF
       END DO

       DO jt = 1, ntrac
          DO jl = kidia, kfdia
             zqdp   = 1./(paphm1(jl,klevp1)-paphm1(jl,klev))
             zfac   = zcfh(jl,klevm1)*zqdp
             zdisxt = 1./(1.+zfac*(1.-zebsh(jl,klevm1)))
             zxtdif(jl,klev,jt) = zdisxt*(zxtdif(jl,klev,jt)+ztmst*g*zqdp* &
                                  zxtems(jl,jt)*ztpfac2+zfac*zxtdif(jl,klevm1,jt))
          END DO
       END DO

       IF (lamip2) THEN

          ! Copy several variables to fields with length klon which
          ! are used when calling subroutine surftemp.

          DO jl = 1,klon
             zsold(jl)  = zcpts(jl)
             zqsold(jl) = zqs(jl)
             zdqsol(jl) = zdqs(jl)
             znetr(jl)  = psrfl(jl)+pemterm(jl,klevp1)*stbo*(ptsm1m(jl)**4)
             zcdrag(jl) = zcons30*zcfh(jl,klev)  ! rho*ch*|V|
          END DO

          ! Calculation of the ground surface temperature pts

          CALL surftemp(klon, klp2, laland, zcons29, cemiss, stbo, zcpq, zcons16,    &
                        alv, als,                                                    &
                        zftn, zetn, zfqn, zeqn, zsold, zqsold, zdqsol,znetr, zgrhfx, &
                        zcdrag, zcair, zcsat, pcvs, zcsurf,                          &
                        zsnew, zqsnew)

          ! Calculation of sklev and qklev using the new surface values
          ! zsnew and zqsnew which were calculated in subroutine surftemp

          DO jl = kidia,kfdia
             ztdif(jl,klev) = ztpfac2*(zetn(jl)*zsnew(jl)  + zftn(jl))
             zqdif(jl,klev) = ztpfac2*(zeqn(jl)*zqsnew(jl) + zfqn(jl))
          END DO
       END IF

       !-- 4.5 Back-substitution

       DO jk = klevm1, itop, -1
          DO jl = kidia, kfdia
             ztdif(jl,jk) = ztdif(jl,jk) + zebsh(jl,jk)*ztdif(jl,jk+1)
             zqdif(jl,jk) = zqdif(jl,jk) + zebsh(jl,jk)*zqdif(jl,jk+1)
             zxdif(jl,jk) = zxdif(jl,jk) + zebsh(jl,jk)*zxdif(jl,jk+1)
          END DO
       END DO

       DO jt = 1, ntrac
          DO jk = klevm1, itop, -1
             DO jl = kidia, kfdia
                zxtdif(jl,jk,jt) = zxtdif(jl,jk,jt) + zebsh(jl,jk)*zxtdif(jl,jk+1,jt)
             END DO
          END DO
       END DO

       IF (lamip2) THEN

          ! Storage of the tilde-tilde-values of temperature and moisture.
          ! In the implicit scheme in subroutine surftemp the surface
          ! fluxes of energy and moisture were computed using these tilde-
          ! tilde-values. for consistency the fluxes have to be calculated
          ! in the same way later in this subroutine (vdiff) and in evapveg.

          DO jl = 1,klon
             zqstt(jl) = zqsnew(jl)
             zsstt(jl) = zsnew(jl)
             zqatt(jl) = ztpfac1 * zqdif(jl,klev)
             zsatt(jl) = ztpfac1 * ztdif(jl,klev)
          END DO

          ! Conversion from zqsnew to zqs(t+dt) and zsnew to pts(t+dt)

          DO jl = 1,klon
             zqsnew(jl)   = ztpfac2*zqsnew(jl) + ztpfac3*zqs(jl)
             zsnew(jl)    = ztpfac2*zsnew(jl) + ztpfac3*zcpts(jl)
             pts(jl)      = zsnew(jl)/zcpq(jl)
             tstiltil(jl) = zsstt(jl)/zcpq(jl)
          END DO
       END IF

       !-- 4.6 Incrementation of t and q tendencies

       DO jk = itop, klev
          DO jl = kidia, kfdia
             zqdif(jl,jk) = zqdif(jl,jk) + ztpfac3*pqm1(jl,jk)
             zdqdt        = (zqdif(jl,jk)-pqm1(jl,jk))*zcons13
             pqte(jl,jk)  = pqte(jl,jk) + zdqdt
             ztdif(jl,jk) = ztdif(jl,jk) + ztpfac3*zcptgz(jl,jk)
             zdtdt        = ((ztdif(jl,jk)+zdis(jl,jk)-pgeom1(jl,jk)) &
                          /  (cpd*(1.+vtmpc2*zqdif(jl,jk)))-ptm1(jl,jk))*zcons13
             ptte(jl,jk)  = ptte(jl,jk) + zdtdt
             zxdif(jl,jk) = zxdif(jl,jk) + ztpfac3*pxm1(jl,jk)
             zdxmdt       = (zxdif(jl,jk)-pxm1(jl,jk))*zcons13
             pxte(jl,jk)  = pxte(jl,jk) + zdxmdt
          END DO
       END DO
       IF (lxtvdiff) THEN
          DO jt = 1, ntrac
             DO jk = itop, klev
                DO jl = kidia, kfdia
                   zxtdif(jl,jk,jt) = zxtdif(jl,jk,jt) + ztpfac3*pxtm1(jl,jk,jt)
                   zdxtdt = (zxtdif(jl,jk,jt)-pxtm1(jl,jk,jt))*zcons13
                   pxtte(jl,jk,jt) = pxtte(jl,jk,jt) + zdxtdt
                END DO
             END DO
          END DO
       END IF

       !-- 4.7 Storage of the surface heat (s.+l.) and moisture
       !       fluxes and their derivatives
       !       against surface variables

       DO jl = kidia, kfdia
          zcoeff = zcons15*ztcoe(jl)

          lo = zhum(jl) <= pqm1(jl,klev)/zqs(jl)
          zhum(jl) = MERGE(0.,zhum(jl),lo)
          zca      = MERGE(0.,1.,lo)
          zhum(jl) = MERGE(zhum(jl),1.,laland(jl))
          zca      = MERGE(zca,1.,laland(jl))
          lo       = pqm1(jl,klev) > zqs(jl)
          zhum(jl) = MERGE(1.,zhum(jl),lo)
          zca      = MERGE(1.,zca,lo)
          zhum(jl) = (1.-pcvs(jl))*(1.-pcvw(jl))*zhum(jl)
          zca      = (1.-pcvs(jl))*(1.-pcvw(jl))*zca

          zqnlev = zqdif(jl,klev) - ztpfac3*pqm1(jl,klev)

          IF (lamip2) THEN
             zzqs = ztpfac2*zqstt(jl)    ! change to implicit calculation
          ELSE
             zzqs = ztpfac2*zqs(jl)
          END IF

          pqhfl(jl) = zcoeff*(zcair(jl)*zqnlev-zcsat(jl)*zzqs)

          ztnlev = ztdif(jl,klev) - ztpfac3*zcptgz(jl,klev)

          IF (lamip2) THEN
             zzcpts = ztpfac2*zsstt(jl)  ! change to implicit calculation
          ELSE
             zzcpts = ztpfac2*zcpts(jl)
          END IF

          pthfl(jl) = zcoeff*(ztnlev-zzcpts)
          pdhft(jl) = -zcons16*pqhfl(jl)
          pthfl(jl) = pthfl(jl) + ptsm1m(jl)*pdhft(jl)

          zxnlev = zxdif(jl,klev) - ztpfac3*pxm1(jl,klev)
          zxhfl  = zcoeff*zxnlev
          prsfl(jl) = MAX(0.,zxhfl)
          pxhfl(jl) = MIN(zxhfl,0.)

          pahfs(jl) = pahfsm(jl) + zdiagt*pthfl(jl)
          pevap(jl) = pevapm(jl) + zdiagw*(pqhfl(jl)+pxhfl(jl))

!ik ocean coupling
          IF (lcouple) THEN
             ahfs_local(jl,jrow) = ahfs_local(jl,jrow) + zdiagt*pthfl(jl)
          END IF

          pdhfqw(jl) = zcoeff*zwlmxi(jl)*(1.-pcvs(jl))*(zqnlev-zzqs)
          pdhfqs(jl) = zcoeff*(zqnlev-zzqs)
       END DO

       dhfsz(irow)  = zdiagt*zbud*SUM(pthfl(1:klon))
       devapz(irow) = zdiagw*zbud*SUM(pqhfl(1:klon))

       DO jl = kidia, kfdia
          IF (lamip2) THEN
             zqhfl(jl) = pqhfl(jl) - pcvs(jl)*pdhfqs(jl)
             pahfl(jl) = alv*zqhfl(jl) + als*pcvs(jl)*pdhfqs(jl)
          ELSE
             pqhfl(jl) = pqhfl(jl) - pcvs(jl)*pdhfqs(jl)
             pahfl(jl) = alv*pqhfl(jl) + als*pcvs(jl)*pdhfqs(jl)
          END IF
          pthfl(jl) = pthfl(jl) + pahfl(jl)
!ik ocean coupling
          IF (lcouple) THEN
             ahfl_local(jl,jrow) = ahfl_local(jl,jrow) + zdiagt*pahfl(jl)
          END IF
          pahfl(jl) = pahflm(jl) + zdiagt*pahfl(jl)

          lo = laland(jl) .OR. ptsm1m(jl) > tmelt
          zalvs = als*pcvs(jl) + alv*(1.-pcvs(jl))
          zalvs = MERGE(zalvs,als,lo)
          pdhft(jl) = -ztpfac2*zcons15*ztcoe(jl)*(zcpts(jl)/ptsm1m(jl)+ &
                      (zalvs-zcons16*ptsm1m(jl))*zcsat(jl)*zdqs(jl)) + pdhft(jl)
       END DO

       IF (lamip2) THEN

          ! Calculation of the four components of evapotranspiration

          CALL evapotr4(klon, klp2, laland, &
               pcvs, pcvw, pvgrat, zeveff, zrelhum, &
               zqnatm, zqs, zcdrag, zqstt, zqatt, &
               zqsnfl, zqlfl, zqvfl, zqgfl)

       END IF

       !-- 4.8 Compute new t2max and min

       DO jl = kidia, kfdia
          lo1 = zricls(jl) >= 0.
          zrat = zhtq/pgeom1(jl,klev)
          zcbn = LOG(1.+(EXP(zbn(jl))-1.)*zrat)
          zcbs = -(zbn(jl)-zbh(jl))*zrat
          zcbu = -LOG(1.+(EXP(zbn(jl)-zbh(jl))-1.)*zrat)
          zred = (zcbn+MERGE(zcbs,zcbu,lo1))/zbh(jl)
          zh2m = zcpts(jl) + zred*(zcptgz(jl,klev)-zcpts(jl))
          zt2 = (zh2m-zhtq)/(cpd*(1.+vtmpc2*pqm1(jl,klev)))
          IF (lamip2) THEN
             ! ------> no accumulation for 2m temperature
             ptemp2(jl) = zt2
          ELSE
             ptemp2(jl) = ptemp2m(jl) + zdiagt*zt2
          END IF
          pt2max(jl) = MAX(pt2maxm(jl),zt2)
          pt2min(jl) = MIN(pt2minm(jl),zt2)

          !-- 4.9 2m dew point

#ifndef NOLOOKUP
          it = INT(ptm1(jl,klev)*1000.)
          IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
          IF (it<jptlucu1 .OR. it>jptlucu2) THEN
             WRITE(*,*) 'vdiff  : 3 it=',it,jl,klev
          END IF
#endif
          it = MAX(MIN(it,jptlucu2),jptlucu1)
          zqs1 = tlucua(it)/papm1(jl,klev)
#else
          lo = ptm1(jl,klev) > tmelt
          zqs1 = c2es*EXP(MERGE(c3les,c3ies,lo)*(ptm1(jl,klev)-tmelt) &
               /(ptm1(jl,klev)-MERGE(c4les,c4ies,lo)))/papm1(jl,klev)
#endif
          zqs1 = zqs1/(1.-vtmpc1*zqs1)
          zrh2m = MAX(zephum,pqm1(jl,klev)/zqs1)

          lo = zt2 > tmelt
          zcvm3 = MERGE(c3les,c3ies,lo)
          zcvm4 = MERGE(c4les,c4ies,lo)
          zaph2m = paphm1(jl,klevp1)*(1.-zhtq/(rd*zt2*(1.+vtmpc1*pqm1(jl,klev))))
#ifndef NOLOOKUP
          it = INT(zt2*1000.)
          IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
#ifdef DEBUG
          IF (it<jptlucu1 .OR. it>jptlucu2) THEN
             WRITE(*,*) 'vdiff  : 4 it=',it,jl,klev
          END IF
#endif
          it = MAX(MIN(it,jptlucu2),jptlucu1)
          zqs2 = tlucua(it)/zaph2m
#else
          zqs2 = c2es*EXP(zcvm3*(zt2-tmelt) / (zt2-zcvm4))/zaph2m
#endif
          zqs2 = zqs2/(1.-vtmpc1*zqs2)
          zq2m = zrh2m*zqs2
          IF (lamip2) THEN
             ! --------> No accumulation for 2m dewpoint
             !           changed to specific humidity (2m)
             !           name is kept as pdew2
             pdew2(jl) = zq2m
          ELSE
             zfrac = LOG(zaph2m*zq2m/(c2es*(1.+vtmpc1*zq2m)))/zcvm3
             pdew2(jl) = pdew2m(jl) + zdiagt*MIN(zt2,(tmelt-zfrac*zcvm4)/(1.-zfrac))
          END IF

          !-- 4.10 10m wind components, max 10m windspeed

          zrat = zhuv/pgeom1(jl,klev)
          zcbn = LOG(1.+(EXP(zbn(jl))-1.)*zrat)
          zcbs = -(zbn(jl)-zbm(jl))*zrat
          zcbu = -LOG(1.+(EXP(zbn(jl)-zbm(jl))-1.)*zrat)
          zred = (zcbn+MERGE(zcbs,zcbu,lo1))/zbm(jl)
          zu10 = zred*pum1(jl,klev)
          zv10 = zred*pvm1(jl,klev)
          zspeed = SQRT(zu10**2+zv10**2)
          IF (lamip2) THEN
             ! ------------> no accumulation for 10m winds !
             pu10(jl) = zu10
             pv10(jl) = zv10
          ELSE
             pu10(jl) = pu10m(jl) + zdiagt*zu10
             pv10(jl) = pv10m(jl) + zdiagt*zv10
          END IF
          pwimax(jl) = MAX(pwimaxm(jl),zspeed)
          pwind10(jl) = pwind10m(jl) + zdiagt*zspeed

!ik ocean coupling
          IF (lcouple) THEN
             wind10_local(jl,jrow) = wind10_local(jl,jrow) + zdiagt*zspeed
          END IF
!fu++ SST diurnal-cycle 
          IF (ldsst) THEN
             wind10_diurnal(jl,jrow) = wind10_diurnal(jl,jrow) + zdiagt*zspeed
          END IF

       END DO


!ik ocean coupling
       IF (lcouple .AND. (jrow == 1)) THEN
          ! count once on each pe
          accu_time = accu_time + zdiagt
       END IF

#ifndef NOLOOKUP
       IF (lookupoverflow) CALL lookuperror ('vdiff (3)   ')
#endif
    ELSE

       !-- 5. Necessary computations if subroutine is by-passed

       DO jl = kidia, kfdia
          paz0(jl)    = paz0m(jl)
          pvdis(jl)   = pvdism(jl)
          pustr(jl)   = pustrm(jl)
          pvstr(jl)   = pvstrm(jl)
          pahfs(jl)   = pahfsm(jl)
          pahfl(jl)   = pahflm(jl)
          pevap(jl)   = pevapm(jl)
          pthfl(jl)   = 0.
          pdhft(jl)   = 0.
          pqhfl(jl)   = 0.
          pxhfl(jl)   = 0.
          pdhfqw(jl)  = 0.
          pdhfqs(jl)  = 0.
          prsfl(jl)   = 0.
          ptemp2(jl)  = ptemp2m(jl)
          pt2max(jl)  = pt2maxm(jl)
          pt2min(jl)  = pt2minm(jl)
          pdew2(jl)   = pdew2m(jl)
          pu10(jl)    = pu10m(jl)
          pv10(jl)    = pv10m(jl)
          pwind10(jl) = pwind10m(jl)
          pustar3(jl) = pustar3m(jl)
          pwimax(jl)  = pwimaxm(jl)
       END DO
       dvdisz(irow) = 0.
       dhfsz(irow)  = 0.
       devapz(irow) = 0.
    END IF

    RETURN
  END SUBROUTINE vdiff

!END MODULE m_vdiff
