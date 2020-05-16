!+ computes large-scale water phase changes and cloud cover.
!+ $Id: m_cond5.f90,v 1.11 2000/09/14 13:14:13 m214003 Exp $

!MODULE m_cond5
!
!CONTAINS

SUBROUTINE cond5 (kidia,kfdia,klon,klp2,ktdia,klev,klevp1,klab,paclcacm, &
&      paphp1,papm1,papp1,pqm1,ptm1,pxlm1,pxim1,pxtec,laland,ktype, &
&      paclcovm,paprlm,pqvim,pxlvim,pxivim,paclc,paclcac,paclcov,paprl,pqvi, &
&      pxlvi,pxivi,pssfl,pqte,ptte,pxlte,pxite,paprs,prsfl)

  ! Description:
  !
  ! Computes large-scale water phase changes and cloud cover.
  !
  ! Method:
  !
  ! This routine computes the physical tendencies of the three
  ! prognostic variables t, q, and x (cloud water) due to water phase
  ! changes (condensation, evaporation of falling precipitation in
  ! unsaturated layers and melting or freezing of the falling water)
  ! and precipitaion formation (coalescense, sedimentation).
  ! and precipitaion formation (autoconversion and accretion in the
  ! warm phase and aggreagation and accretion in the cold phase).
  !
  ! Rain, snowfall, surface fluxes, cloud water and cloud cover
  ! later to be used for soil processes and radiation are stored.
  !
  ! References:
  !
  !    Lohmann and Roeckner, 1995: MPI-Report no. 179
  !    Levkov et al., 1992: beitr. phys. atmos. 35-58. (ice phase)
  !    Beheng, 1994: atmos. res. 193-206.              (warm phase)
  !    ECHAM3/4-description, 1992: MPI-report no. 93   (cloud cover,
  !                                        condensation,evaporation)
  !
  ! Authors:
  !
  ! U. Lohmann, MPI, in 1995
  ! U. Schulzweida, MPI, Dec 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_constants,         ONLY: als, alv, c4ies, c4les, c5ies, &
                                  c5les, cpd, g, rd, rhoh2o, tmelt,    &
                                  vtmpc1, vtmpc2, api
  USE mo_diagnostics_zonal, ONLY: dlserz, dlsesz, dlsgrz, dlsgsz, dlsmsz
  USE mo_control,           ONLY: nrow, twodt
  USE mo_gaussgrid,         ONLY: budw
  USE mo_start_dataset,     ONLY: nstart, nstep
  USE mo_tmp_buffer,        ONLY: loland
  USE mo_machine,           ONLY: prec
#ifndef NOLOOKUP
  USE mo_convect_tables,    ONLY: tlucua, tlucuaw, &! lookup table
                                  jptlucu1, jptlucu2, &
                                  lookuperror, lookupoverflow
#else
  USE mo_constants,         ONLY: c3les, c3ies, c2es
#endif

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: kfdia, kidia, klev, klevp1, klon, klp2, ktdia

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: paclcacm(klp2,klev), paclcovm(klp2),  &
&        paphp1(klp2,klevp1), papm1(klp2,klev), &
&        papp1(klp2,klev), paprlm(klp2),  &
&        pqm1(klp2,klev),  &
&        pqvim(klp2), ptm1(klp2,klev), pxlvim(klp2), pxivim(klp2), &
&        pxlm1(klp2,klev),  pxim1(klon,klev), pxtec(klp2,klev)
  INTEGER, INTENT (IN) :: klab(klp2,klev), ktype(klp2)
  LOGICAL, INTENT (IN) :: laland(klp2)

  ! klab     : convection flag (0: no convection, 1:  
  ! paclcacm : cloud cover, accumulated (old value)
  ! paphp1   : pressure at half levels (t+dt)
  ! papm1    : pressure at full levels (t-dt)
  ! papp1    : pressure at full levels (t+dt)
  ! pqm1     : specific humidity (t-dt)
  ! ptm1     : temperature (t-dt)
  ! pxlm1    : cloud water (t-dt)
  ! pxim1    : cloud ice (t-dt)
  ! pxtec    : tendency of detrained convective cloud water
  ! ktype    : type of convection
  ! laland   : logical mask for land
  ! paclcovm : total cloud cover, accumulated (old value)
  ! paprlm   : stratiform precipitaion, accumulated (old value)
  ! pqvim    : vertically integrated spec. humidity, accumulated (old value)
  ! pxlvim   : vertically integrated cloud water, accumulated (old value)
  ! pxivim   : vertically integrated cloud ice, accumulated (old value)

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: ptte(klp2,klev), pqte(klp2,klev), pxlte(klp2,klev),  &
&                         pxite(klon,klev), paprs(klp2), prsfl(klp2)

  !  ptte     : tendency of temperature
  !  pqte     : tendency of specific humidity
  !  pxlte    : tendency of cloud water
  !  pxite    : tendency of cloud ice
  !  paprs    : snow fall, accumulated
  !  prsfl    : surface rain flux

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: paclc(klp2,klev), paclcac(klp2,klev),  &
&                       paclcov(klp2), paprl(klp2),  &
&                       pqvi(klp2), pssfl(klp2), pxlvi(klp2), pxivi(klp2)

  !  paclc    : cloud cover
  !  paclcac  : cloud cover, accumulated
  !  paclcov  : total cloud cover
  !  paprl    : stratiform precipitation, accumulated
  !  pqvi     : vertically integrated spec. humidity, accumulated
  !  pssfl    : surface snow flux
  !  pxlvi    : vertically integrated cloud water, accumulated
  !  pxivi    : vertically integrated cloud ice, accumulated

  !  Local scalars: 
  REAL :: rdcpd, zbud, zc1, zf1, zqrho0, zdpg, zxrp1, zxsp1, zraut, zracc, &
&      zdt2, zsaut, zsaci, zlamsm, zcolleffi, zzdrr, zzdrs, zzepr, zzeps,  &
&      zclcm1, zclear, zcloia, zclp1, zcons, zcons1, zcons2, zradl, &
&      zcor, zcpn, zcvm4, zcvm5, zdiagt, zdiagw, zprat, zcdnc, &
&      zdt, zdthdp, zepclc, zxletmst, zxlmvn, zxli, zxlmax, zxlps, &
&      zmeltp2, zqpa, zqr, zxlma, zxima, zxietmst, zximvn, zxii, zximax, &
&      zqsm1, zqtmst, zrcp, zrflo, zrhsc, zrs, zrtc, zrtl, zsatsc, &
&      zsflo, zsnmlt, ztc, ztmst, zxips, zrho0, zrieff, zrih, zri, &
&      zzep, zzepcons, zzqt, zzrh, zss, zsr, zlfdcp, zxlp1, zxip1,  &
&      zsec, zmin, zqg, zthomi, zn0s, zrhoi, zcraut, zcsaut, zcsacl, &
&      zxlold, zdxlcor, zxiold, zdxicor
  REAL :: zxlb(klon), zxib(klon), zxsource(klon), zqtvi(klon), zxltvi(klon), &
&      zxitvi(klon)
  INTEGER :: irow, jb, jk, jl, klev2p1, nexc, nexl
#ifndef NOLOOKUP
  INTEGER :: it
#endif
  LOGICAL :: lo, lo1, lo2

  !  Local arrays: 
  REAL :: zclcp1(klon,klev), zqte(klon,klev), ztte(klon,klev), &
&      zdp(klon), zdthmin(klon), zimlt(klon), zsacl(klon), &
&      zepr(klon), zeps(klon), zrpr(klon), zspr(klon), zevp(klon), zsub(klon), &
&      zlsdcp(klon,klev), zlvdcp(klon,klev), zqcon(klon,klev), &
&      zqp1(klon,klev), zqsp1(klon,klev), zrain(klon), zrfl(klon), zfrl(klon), &
&      zrfln(klon), zrhc(klon,klev), zsat(klon,klev), zsfl(klon), zsfln(klon), &
&      zsnow(klon), ztheta(klon,klev), zxite(klon,klev), zcnd(klon), zdep(klon), &
&      ztp1(klon,klev), zsmlt(klon), zxlte(klon,klev)
  INTEGER :: invb(klon)

  !  Intrinsic functions 
#ifdef ECLIB
  REAL :: EXPHF
!DIR$ VFUNCTION EXPHF
#define EXP(x)  EXPHF(x)
#else
  INTRINSIC EXP
#endif
  INTRINSIC MAX, MERGE, MIN, SQRT, SUM


  !  Executable statements 

#ifndef NOLOOKUP
  lookupoverflow = .FALSE.
#endif

  irow=nrow(1)

! security parameters.

  zsec=1.e-12
  zmin=1.e-30
  zepclc=1.e-2

! computational constants.

!  1. input

  zmeltp2=tmelt+2.
  zqg=1./g
  rdcpd=rd/cpd
  zrtc=0.6
  zrtl=0.6
  nexc=4
  nexl=4
  klev2p1=klev/2+1
  zcloia=1.0e+02
  zrs=0.99
  zrhsc=0.4
  zsatsc=0.8
  zthomi=238.16
  zn0s=3.e6
  zrhoi=500.
  
  zzepcons=0.15
  zcraut=15.
  zcsaut=220.
  zcsacl=0.1
  
  ztmst=twodt
  IF (nstep.EQ.nstart) ztmst=0.5*twodt
  zqtmst=1./ztmst
  zdiagt=0.5*twodt
  zdiagw=zdiagt*0.001
  zcons1=cpd*vtmpc2
  zcons2=1./(ztmst*g)

!  2.  quantities needed for condensation and precipitation calculations.

! 2.1     set to zero all tendencies

  DO jl=kidia,kfdia
    zrain(jl)=0.
    zsnow(jl)=0.
    zrfl(jl)=0.
    zsfl(jl)=0.
    zepr(jl)=0.
    zeps(jl)=0.
    zsmlt(jl)=0.
    zrfln(jl)=0.
    zsfln(jl)=0.
    invb(jl)=1
    zdthmin(jl)=0.0
  END DO
  zsr=0.
  zss=0.

!  2.2      calculate potential temperatures

  DO jk=klev,klev/2,-1
    DO jl=kidia,kfdia
      ztheta(jl,jk)=ptm1(jl,jk)*(1.0e5/papm1(jl,jk))**rdcpd
    END DO
  END DO

! 2.3    check for presence of low-level inversion ( sea points ONLY )

  DO jk=klev,klev2p1,-1
    DO jl=kidia,kfdia
      IF((.NOT.laland(jl)) .AND. ktype(jl).EQ.0) THEN
       zdthdp=(ztheta(jl,jk)-ztheta(jl,jk-1))*zcloia/(papm1(jl,jk)  &
     &         -papm1(jl,jk-1))
       lo=zdthdp.LT.zdthmin(jl)
       zdthmin(jl)=MERGE(zdthdp,zdthmin(jl),lo)
       invb(jl)=MERGE(jk,invb(jl),lo)
      END IF
    END DO
  END DO

  DO jk=ktdia,klev
    DO jl=kidia,kfdia

! 2.4   t, q and qs provisional values at t+dt
!       effectives l for vaporisation and sublimation

      ztp1(jl,jk)=ptm1(jl,jk)+ptte(jl,jk)*ztmst
      zqp1(jl,jk)=pqm1(jl,jk)+pqte(jl,jk)*ztmst
      zqp1(jl,jk)=MAX(zqp1(jl,jk),0.)
      zxlte(jl,jk)=0.
      zxite(jl,jk)=0.

      zrcp=1./(cpd+zcons1*zqp1(jl,jk))
      zlvdcp(jl,jk)=alv*zrcp
      zlsdcp(jl,jk)=als*zrcp

      lo2=(ztp1(jl,jk).LT.zthomi).OR.(ztp1(jl,jk).LT.tmelt  &
     & .AND.pxim1(jl,jk).GT.zsec)
      zcvm4=MERGE(c4ies,c4les,lo2)
      zcvm5=MERGE(c5ies*zlsdcp(jl,jk),c5les*zlvdcp(jl,jk),lo2)
      zzqt=1./(ztp1(jl,jk)-zcvm4)
#ifndef NOLOOKUP
      it=NINT(ztp1(jl,jk)*1000.)
      IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
      it = MAX(MIN(it,jptlucu2),jptlucu1)
      zqsp1(jl,jk)=MERGE(tlucua(it),tlucuaw(it),lo2)/papp1(jl,jk)
#else
      lo = ztp1(jl,jk) > tmelt
      zqsp1(jl,jk)=MERGE(c2es*EXP(MERGE(c3les,c3ies,lo)*(ztp1(jl,jk)-tmelt) &
     &        /    (ztp1(jl,jk)-MERGE(c4les,c4ies,lo))), &
     &              c2es*EXP(c3les*(ztp1(jl,jk)-tmelt)*(1./(ztp1(jl,jk)-c4les))), &
     &              lo2)/papp1(jl,jk)
#endif
      zqsp1(jl,jk)=MIN(zqsp1(jl,jk),0.5)
      zcor=1./(1.-vtmpc1*zqsp1(jl,jk))
      zqsp1(jl,jk)=zqsp1(jl,jk)*zcor
      zqcon(jl,jk)=1./(1.+zcvm5*zqsp1(jl,jk)*zcor*zzqt**2)

! 2.5    threshold relative humidity

      jb=invb(jl)
      lo=(jb.GT.klev-6.AND. jb.LT.klev-1)
      IF(klab(jl,jk).EQ.2) THEN
        zrhc(jl,jk)=zrtc+(zrs-zrtc)*EXP(1.-  &
     &    (paphp1(jl,klevp1)/papp1(jl,jk))**nexc)
        zsat(jl,jk)=1.
      ELSE
        IF(jk.GE.jb+1 .AND. lo) THEN
          zrhc(jl,jk)=zrhsc
          zsat(jl,jk)=zsatsc
        ELSE
          zrhc(jl,jk)=zrtl+(zrs-zrtl)*EXP(1.- &
     &      (paphp1(jl,klevp1)/papp1(jl,jk))**nexl)
          zsat(jl,jk)=1.
        END IF
      ENDIF

! 2.6  cloud cover at t+dt

      zqr=zqsp1(jl,jk)*zsat(jl,jk)*zrhc(jl,jk)
      zclcp1(jl,jk)=(zqp1(jl,jk)-zqr)/(zqsp1(jl,jk)*zsat(jl,jk)-zqr)
      zclcp1(jl,jk)=MAX(zclcp1(jl,jk),0.)
      zclcp1(jl,jk)=MIN(zclcp1(jl,jk),1.)
      zclcp1(jl,jk)=1.-SQRT(1.-zclcp1(jl,jk))

    END DO
  END DO
#ifndef NOLOOKUP
  IF (lookupoverflow) CALL lookuperror ('cond5 (1)   ')
#endif

  DO jk=ktdia,klev

    DO jl=kidia,kfdia

      zlfdcp=zlsdcp(jl,jk)-zlvdcp(jl,jk)
      zcnd(jl)=0.
      zdep(jl)=0.
      zfrl(jl)=0.
      zrpr(jl)=0.
      zspr(jl)=0.
      zimlt(jl)=0.
      zsacl(jl)=0.
      
      zxlb(jl)=pxlm1(jl,jk)
      zxib(jl)=pxim1(jl,jk)
      zxlma=zxlb(jl)
      zxima=zxib(jl)
      
      IF(zclcp1(jl,jk).GT.1.e-10.AND.pxlm1(jl,jk).GT.0.) THEN
        zxlb(jl)=pxlm1(jl,jk)/zclcp1(jl,jk)
        zxlma=0.
      ENDIF
      IF(zclcp1(jl,jk).GT.1.e-10.AND.pxim1(jl,jk).GT.0.) THEN
        zxib(jl)=pxim1(jl,jk)/zclcp1(jl,jk)
        zxima=0.
      ENDIF

! 3.1   melting of incoming snow and
!       instantenous melting of cloud ice and
!       instantenous freezing of cloud water

      zdp(jl)=paphp1(jl,jk+1)-paphp1(jl,jk)
      zcons=zcons2*zdp(jl)/(zlsdcp(jl,jk)-zlvdcp(jl,jk))
      zsnmlt=MIN(zsfl(jl),zcons*MAX(0.,(ztp1(jl,jk)-zmeltp2)))
      zsmlt(jl)=zsmlt(jl)+zsnmlt
      zrfln(jl)=zrfl(jl)+zsnmlt
      zsfln(jl)=zsfl(jl)-zsnmlt
      zdt=-zsnmlt/zcons
      ztp1(jl,jk)=ztp1(jl,jk)+zdt
      ptte(jl,jk)=ptte(jl,jk)+zdt*zqtmst
      IF (ztp1(jl,jk).GE.tmelt) THEN
        zimlt(jl)=zxib(jl)*zclcp1(jl,jk)
        ztp1(jl,jk)=ztp1(jl,jk)-zlfdcp*zimlt(jl)
        zxlb(jl)=zxlb(jl)+zxib(jl)
        zxib(jl)=0.
      ENDIF
      IF (ztp1(jl,jk).LE.zthomi) THEN
        zfrl(jl)=zxlb(jl)*zclcp1(jl,jk)
        ztp1(jl,jk)=ztp1(jl,jk)+zlfdcp*zfrl(jl)
        zxib(jl)=zxib(jl)+zxlb(jl)
        zxlb(jl)=0.
      ENDIF

!  3.2   saturation specific humidity for t-dt.

      lo2=(ztp1(jl,jk).LT.zthomi).OR.(ztp1(jl,jk).LT.tmelt  &
     &       .AND.zxib(jl).GT.zsec)
#ifndef NOLOOKUP
      it=NINT(ptm1(jl,jk)*1000.)
      IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
      it = MAX(MIN(it,jptlucu2),jptlucu1)
      zqsm1=MERGE(tlucua(it),tlucuaw(it),lo2)/papp1(jl,jk)
#else
      lo = ptm1(jl,jk) > tmelt
      zqsm1=MERGE(c2es*EXP(MERGE(c3les,c3ies,lo)*(ptm1(jl,jk)-tmelt) &
     &      /   (ptm1(jl,jk)-MERGE(c4les,c4ies,lo))), &
     &          c2es*EXP(c3les*(ptm1(jl,jk)-tmelt)*(1./(ptm1(jl,jk)-c4les))), &
     &          lo2)/papp1(jl,jk)
#endif
      zqsm1=MIN(zqsm1,0.5)
      zcor=1./(1.-vtmpc1*zqsm1)
      zqsm1=zqsm1*zcor

! 3.3  cloud cover at  t-dt

      zqr=zqsm1*zsat(jl,jk)*zrhc(jl,jk)
      zclcm1=(pqm1(jl,jk)-zqr)/(zqsm1*zsat(jl,jk)-zqr)
      zclcm1=MAX(zclcm1,0.)
      zclcm1=MIN(zclcm1,1.)
      zclcm1=1.-SQRT(1.-zclcm1)

! 3.4  humidity and liquid water in the cloudy and cloud-free
!      part for t-1.

      zqpa=pqm1(jl,jk)
      IF(zclcm1.GT.1.e-10) zqpa=MAX(zqsm1,pqm1(jl,jk))

! 3.5.   condensation/evaporation rate + detrained cloud water

      zcpn=(zqpa+pqte(jl,jk)*ztmst-zqsp1(jl,jk))*zqcon(jl,jk)
      lo2=(ztp1(jl,jk).LT.zthomi).OR.(ztp1(jl,jk).LT.tmelt &
     &       .AND.zxib(jl).GT.zsec)
      IF (.NOT.lo2) THEN
        zxletmst=pxlte(jl,jk)*ztmst
        zxlb(jl)=zxlb(jl)+zxletmst
        zxlmvn=zxlma+zxletmst
        lo=zxlb(jl).GE.0.
        lo1=lo.AND.-zcpn.GT.zxlb(jl)
        zcpn=MERGE(-zxlb(jl),zcpn,lo1)
        zxli=-zxlb(jl)
        zxlmax=MAX(zcpn,0.)
        zxlps=zxlb(jl)+zcpn
        zxlb(jl)=MERGE(zxlps,zxlmax,lo)
        zcpn=MERGE(zcpn,zxlb(jl)+zxli,lo)
        zcnd(jl)=zclcp1(jl,jk)*zcpn+(1.-zclcp1(jl,jk))*(-zxlmvn)
        IF(zcnd(jl).GT.0.) zcnd(jl)=MIN(zcnd(jl),zqp1(jl,jk))
        zxlb(jl)=zxlb(jl)+pxtec(jl,jk)*ztmst
        zxlte(jl,jk)=pxtec(jl,jk)
      ELSE
        zxietmst=pxite(jl,jk)*ztmst
        zxib(jl)=zxib(jl)+zxietmst
        zximvn=zxima+zxietmst
        lo=zxib(jl).GE.0.
        lo1=lo.AND.-zcpn.GT.zxib(jl)
        zcpn=MERGE(-zxib(jl),zcpn,lo1)
        zxii=-zxib(jl)
        zximax=MAX(zcpn,0.)
        zxips=zxib(jl)+zcpn
        zxib(jl)=MERGE(zxips,zximax,lo)
        zcpn=MERGE(zcpn,zxib(jl)+zxii,lo)
        zdep(jl)=zclcp1(jl,jk)*zcpn+(1.-zclcp1(jl,jk))*(-zximvn)
        IF(zdep(jl).GT.0.) zdep(jl)=MIN(zdep(jl),zqp1(jl,jk))
        zxib(jl)=zxib(jl)+pxtec(jl,jk)*ztmst
        zxite(jl,jk)=pxtec(jl,jk)
      ENDIF

! 3.8 ice nuclei after ou and liou

      zrho0=papp1(jl,jk)/(rd*ztp1(jl,jk))
      ztc=ztp1(jl,jk)-tmelt
      ztc=MIN(0.,ztc)
      zrieff=10.
!     zri=0.5e-6*(326.3+12.42*ztc+0.197*ztc*ztc+0.0012*ztc*ztc*ztc)
      IF (zxib(jl).GE.zsec) zrieff=83.8*(zxib(jl)*zrho0*1000.)**0.216
      zrih=-2261.236+SQRT(5113188.044+2809.*zrieff*zrieff*zrieff)
      zri=1.e-6*zrih**(1./3.)
      zri=MAX(1.e-5,zri)
      
      !  calculation of the cloud droplet number concentration
      
      zprat=MIN(2.,80000./papp1(jl,jk))**2
      IF (loland(jl)) THEN
        IF (papp1(jl,jk).LT.80000.) THEN
          zcdnc=1.e6*(50.+(135.-50.)*(EXP(1.-zprat)))
        ELSE
          zcdnc=135.e6
        ENDIF
      ELSE
        IF (papp1(jl,jk).LT.80000.) THEN
          zcdnc=1.e6*(50.+(80.-50.)*(EXP(1.-zprat)))
        ELSE
          zcdnc=80.e6
        ENDIF
      ENDIF

      lo=zxlb(jl)*zclcp1(jl,jk).GE.zmin  &
     &   .AND.ztp1(jl,jk).LT.tmelt.AND.ztp1(jl,jk).GT.zthomi
      IF (lo) THEN
        zfrl(jl)=100.*(EXP(0.66*(tmelt-ztp1(jl,jk)))-1.)* &
     &      zrho0*zxlb(jl)*zxlb(jl)*0.001/zcdnc*ztmst
        zradl=(0.75*zxlb(jl)*zrho0/(api*rhoh2o*zcdnc))**(1./3.)
        zf1=4.*api*zradl*zcdnc*2.e5*(270.16-ztp1(jl,jk))/zrho0
        zf1=MAX(0.,zf1)
        zfrl(jl)=zfrl(jl)+ztmst*1.4e-20*zf1
        zfrl(jl)=MAX(0.,MIN(zfrl(jl),zxlb(jl)))
        zxlb(jl)=zxlb(jl)-zfrl(jl)
        zxib(jl)=zxib(jl)+zfrl(jl)
        zfrl(jl)=zclcp1(jl,jk)*zfrl(jl)
        zxlb(jl)=MAX(zxlb(jl),0.)
      ENDIF

      ztp1(jl,jk)=ztp1(jl,jk)+zcnd(jl)*zlvdcp(jl,jk)+zdep(jl)*zlsdcp(jl,jk)+zlfdcp*zfrl(jl)
      zqp1(jl,jk)=zqp1(jl,jk)-zcnd(jl)-zdep(jl)

      zrho0=papp1(jl,jk)/(rd*ztp1(jl,jk))
      zqrho0=1.3/zrho0

!  4.  cloud physics and precipitation fluxes.

!  4.1 warm cloud physic (after beheng, 1994)

      zdpg=zdp(jl)*zqg*zqtmst
      zxrp1=zrfln(jl)/zdpg
      zxsp1=zsfln(jl)/zdpg
      
      IF (zxlb(jl)*zclcp1(jl,jk).GE.zsec) THEN
        zraut=zcraut*1.2e27/zrho0*(zcdnc*1.e-6)**(-3.3)*(zrho0*zxlb(jl)*1.e-3)**4.7*ztmst
        zracc=6.*zrho0*zxrp1*zxlb(jl)*ztmst
        zrpr(jl)=zclcp1(jl,jk)*MAX(0.,MIN(zracc+zraut,zxlb(jl)))
      ELSE
        zrpr(jl)=zclcp1(jl,jk)*zxlb(jl)
      ENDIF

!  4.2  cold cloud physics
!  conversion from cloud ice to snow after levkov et al. 1992

      IF (zxib(jl)*zclcp1(jl,jk).GT.zsec) THEN
        zc1=17.5*zxib(jl)*zrho0/zrhoi*zqrho0**0.33
        zdt2=0.
        IF (zri.LE.0.99e-4) zdt2=-6./zc1*alog10(zri*1.e4)
        IF (zdt2.GE.1.e-3) THEN
          zsaut=zcsaut*zxib(jl)*ztmst/zdt2
        ELSE
          zsaut=zxib(jl)*ztmst
        ENDIF
        zsaci=0.
        IF (zxsp1.GE.zsec) THEN
          zlamsm=(zrho0*zxsp1/(api*rhoh2o*zn0s))**0.8125
          zcolleffi=EXP(0.025*(ztp1(jl,jk)-tmelt))
          zsaci=api*zn0s*3.078*zxib(jl)*zlamsm*zqrho0**0.5*zcolleffi*ztmst
        ENDIF
        zspr(jl)=zclcp1(jl,jk)*MAX(0.,MIN(zsaci+zsaut,zxib(jl)))
      ELSE
        zspr(jl)=zclcp1(jl,jk)*zxib(jl)
      ENDIF
    
      IF (zxlb(jl)*zclcp1(jl,jk).GT.zsec.AND.zxsp1.GE.zsec) THEN
        zlamsm=(zrho0*zxsp1/(api*rhoh2o*zn0s))**0.8125
        zsacl(jl)=api*zn0s*3.078*zxlb(jl)*zlamsm*zqrho0**0.5*ztmst*zcsacl
        zsacl(jl)=zclcp1(jl,jk)*MAX(0.,MIN(zsacl(jl),(zxlb(jl)-zrpr(jl))))
      ENDIF
    END DO
#ifndef NOLOOKUP
    IF (lookupoverflow) CALL lookuperror ('cond5 (2)   ')
#endif

    DO jl=kidia,kfdia
      zevp(jl)=0.
      zsub(jl)=0.

! 5.1   new precipitation fluxes after consideration of
!       coalescence and sedimentation processes.

      zzdrr=MAX(0.,zcons2*zdp(jl)*zrpr(jl))
      zzdrs=MAX(0.,zcons2*zdp(jl)*(zspr(jl)+zsacl(jl)))
      zrfln(jl)=zrfln(jl)+zzdrr
      zrain(jl)=zrain(jl)+zzdrr
      zsfln(jl)=zsfln(jl)+zzdrs
      zsnow(jl)=zsnow(jl)+zzdrs
    END DO

! 5.2 at top layer or if there is no precipitation avoid evaporation.

    IF (ABS(zsr)<prec .AND. jk>1) THEN
      zsr=SUM(zrfln(1:klon))
    END IF

    IF (jk>1 .AND. zsr>0.) THEN

! 5.3     evaporation of precipitation.

      DO jl=kidia,kfdia
        zclear=1.-zclcp1(jl,jk)
        zzrh=zrhc(jl,jk)*zsat(jl,jk)+zclcp1(jl,jk)*(1.-zrhc(jl,jk)*zsat(jl,jk))
        zzep=MIN(0.,(MIN(zqp1(jl,jk),zqsp1(jl,jk)*zzrh)-zqsp1(jl,jk)))
        zzep=zzepcons*zzep
        zzepr=MAX(-zrfl(jl),zzep*zqcon(jl,jk)*zcons2*zdp(jl))
        zevp(jl)=-zzepr/(zcons2*zdp(jl))*zclear
        zrflo=zrfln(jl)
        zrfln(jl)=zrfln(jl)+zzepr*zclear
        zepr(jl)=zepr(jl)+(zrflo-zrfln(jl))
      END DO
    END IF
    
    IF (ABS(zss)<prec .AND. jk>1) THEN
      zss=SUM(zsfln(1:klon))
    END IF

    IF (jk>1 .AND. zss>0. .AND. ztp1(jl,jk)<=zmeltp2) THEN

! 5.4     sublimation of precipitation.

      DO jl=kidia,kfdia
        zclear=1.-zclcp1(jl,jk)
        zzrh=zrhc(jl,jk)*zsat(jl,jk)+zclcp1(jl,jk)*(1.-zrhc(jl,jk)*zsat(jl,jk))
        zzep=MIN(0.,(MIN(zqp1(jl,jk),zqsp1(jl,jk)*zzrh)-zqsp1(jl,jk)))
        zzep=zzepcons*zzep
        zzeps=MAX(-zsfl(jl),zzep*zqcon(jl,jk)*zcons2*zdp(jl))
        zzeps=MIN(zzeps,0.)
        zsub(jl)=-zzeps/(zcons2*zdp(jl))*zclear
        zsflo=zsfln(jl)
        zsfln(jl)=zsfln(jl)+zzeps*zclear
        zeps(jl)=zeps(jl)+(zsflo-zsfln(jl))
      END DO
    END IF

! 6.    incrementation of t,q and x tendencies and surface fluxes

    DO jl=kidia,kfdia

! 6.1     modification of the t,q and x tendencies,
!         cloud cover for diagnostics and radiation.

      zqte(jl,jk)=(zevp(jl)+zsub(jl)-zcnd(jl)-zdep(jl))*zqtmst
      ztte(jl,jk)=(zlvdcp(jl,jk)*(zcnd(jl)-zevp(jl))+zlsdcp(jl,jk)*(zdep(jl)-zsub(jl)) &
     &         +(zlsdcp(jl,jk)-zlvdcp(jl,jk))*(-zimlt(jl)+zfrl(jl)+zsacl(jl)))*zqtmst
      zxlte(jl,jk)=zxlte(jl,jk)+(-zrpr(jl)+zimlt(jl)-zfrl(jl)-zsacl(jl)+zcnd(jl))*zqtmst
      zxite(jl,jk)=zxite(jl,jk)+(-zspr(jl)-zimlt(jl)+zfrl(jl)+zdep(jl))*zqtmst


      ptte(jl,jk)=ptte(jl,jk)+ztte(jl,jk)
      pqte(jl,jk)=pqte(jl,jk)+zqte(jl,jk)
      pxlte(jl,jk)=pxlte(jl,jk)+zxlte(jl,jk)
      pxite(jl,jk)=pxite(jl,jk)+zxite(jl,jk)
      ztp1(jl,jk)=ptm1(jl,jk)+ptte(jl,jk)*ztmst
      zqp1(jl,jk)=pqm1(jl,jk)+pqte(jl,jk)*ztmst
      zxlp1=pxlm1(jl,jk)+pxlte(jl,jk)*ztmst
      zxip1=pxim1(jl,jk)+pxite(jl,jk)*ztmst

      lo2=(ztp1(jl,jk).LT.zthomi).OR.(ztp1(jl,jk).LT.tmelt.AND.zxip1.GT.zsec)
#ifndef NOLOOKUP
      it=NINT(ztp1(jl,jk)*1000.)
      IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
      it = MAX(MIN(it,jptlucu2),jptlucu1)
      zqsp1(jl,jk)=MERGE(tlucua(it),tlucuaw(it),lo2)/papp1(jl,jk)
#else
      lo = ztp1(jl,jk) > tmelt
      zqsp1(jl,jk)=MERGE(c2es*EXP(MERGE(c3les,c3ies,lo)*(ztp1(jl,jk)-tmelt) &
     &         /  (ztp1(jl,jk)-MERGE(c4les,c4ies,lo))),  &
     &            c2es*EXP(c3les*(ztp1(jl,jk)-tmelt)*(1./(ztp1(jl,jk)-c4les))),  &
     &            lo2)/papp1(jl,jk)
#endif
      zqsp1(jl,jk)=MIN(zqsp1(jl,jk),0.5)
      zcor=1./(1.-vtmpc1*zqsp1(jl,jk))
      zqsp1(jl,jk)=zqsp1(jl,jk)*zcor
      
      zqp1(jl,jk)=MAX(zmin,zqp1(jl,jk))
      zqr=zqsp1(jl,jk)*zsat(jl,jk)*zrhc(jl,jk)
      zclp1=(zqp1(jl,jk)-zqr)/(zqsp1(jl,jk)*zsat(jl,jk)-zqr)
      zclp1=MAX(zclp1,0.)
      zclp1=MIN(zclp1,1.)
      zclp1=1.-SQRT(1.-zclp1)
      
      paclc(jl,jk)=zclp1
      paclcac(jl,jk)=paclcacm(jl,jk)+paclc(jl,jk)*zdiagt
      
      zxlold=zxlp1
      lo=(zxlp1.LT.zmin.OR.zclp1.LT.zepclc)
      zxlp1=MERGE(0.,zxlp1,lo)
      zdxlcor=(zxlp1-zxlold)*zqtmst
      
      zxiold=zxip1
      lo=(zxip1.LT.zmin.OR.zclp1.LT.zepclc)
      zxip1=MERGE(0.,zxip1,lo)
      zdxicor=(zxip1-zxiold)*zqtmst


      pxlte(jl,jk)=pxlte(jl,jk)+zdxlcor
      pxite(jl,jk)=pxite(jl,jk)+zdxicor
      pqte(jl,jk)=pqte(jl,jk)-zdxlcor-zdxicor
      ptte(jl,jk)=ptte(jl,jk)+zlvdcp(jl,jk)*zdxlcor+zlsdcp(jl,jk)*zdxicor
    END DO
#ifndef NOLOOKUP
    IF (lookupoverflow) CALL lookuperror ('cond5 (3)   ')
#endif

! 6.4 swap of the fluxes and END of the vertical loop.

    DO jl=kidia,kfdia
      zrfl(jl)=zrfln(jl)
      zsfl(jl)=zsfln(jl)
    END DO

  END DO

! 7. vertically integrated humidity and cloud water and surface fluxes

  DO jl=kidia,kfdia
    prsfl(jl)=prsfl(jl)+zrfl(jl)
    pssfl(jl)=zsfl(jl)
    paprl(jl)=paprlm(jl)+zdiagw*(prsfl(jl)+pssfl(jl))
    paprs(jl)=paprs(jl)+pssfl(jl)*zdiagw
    
    pqvi(jl)=0.
    pxlvi(jl)=0.
    pxivi(jl)=0.
    
    zqtvi(jl)=0.
    zxltvi(jl)=0.
    zxitvi(jl)=0.
    zxsource(jl)=0.
  END DO

  DO jk=ktdia,klev
    DO jl=kidia,kfdia
      zdpg=(paphp1(jl,jk+1)-paphp1(jl,jk))*zqg
      pqvi(jl)=pqvi(jl)+pqm1(jl,jk)*zdpg
      pxlvi(jl)=pxlvi(jl)+pxlm1(jl,jk)*zdpg
      pxivi(jl)=pxivi(jl)+pxim1(jl,jk)*zdpg
      zqtvi(jl)=zqtvi(jl)+zqte(jl,jk)*zdpg
      zxltvi(jl)=zxltvi(jl)+zxlte(jl,jk)*zdpg
      zxitvi(jl)=zxitvi(jl)+zxite(jl,jk)*zdpg
      zxsource(jl)=zxsource(jl)+pxtec(jl,jk)*zdpg
    END DO
  END DO

  DO jl=kidia,kfdia
    pqvi(jl)=pqvim(jl)+zdiagt*pqvi(jl)
    pxlvi(jl)=pxlvim(jl)+zdiagt*pxlvi(jl)
    pxivi(jl)=pxivim(jl)+zdiagt*pxivi(jl)
  END DO

!  7.1    accumulated total cloudcover

  DO jl=kidia,kfdia
    paclcov(jl)=1.-paclc(jl,1)
  END DO

  DO jk=2,klev
    DO jl=kidia,kfdia
      paclcov(jl)=paclcov(jl)*(1.-MAX(paclc(jl,jk),paclc(jl,jk-1))) &
     &                    /(1.-MIN(paclc(jl,jk-1),1.-zsec))
    END DO
  END DO

  DO jl=kidia,kfdia
    paclcov(jl)=1.-paclcov(jl)
    paclcov(jl)=paclcovm(jl)+zdiagt*paclcov(jl)
  END DO

  zbud=budw(irow)
  dlsgrz(irow)=zdiagw*zbud*SUM(zrain(1:klon))
  dlsgsz(irow)=zdiagw*zbud*SUM(zsnow(1:klon))
  dlsmsz(irow)=zdiagw*zbud*SUM(zsmlt(1:klon))
  dlserz(irow)=zdiagw*zbud*SUM(zepr(1:klon))
  dlsesz(irow)=zdiagw*zbud*SUM(zeps(1:klon))

  RETURN
END SUBROUTINE cond5

!END MODULE m_cond5
