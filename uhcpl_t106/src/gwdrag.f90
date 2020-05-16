SUBROUTINE gwdrag ( klon, ktdia, klev, klevm1, klevp1,  klp2,     &
                    paphm1,      papm1,        pgeom1,            &
                    ptm1,        pum1,         pvm1,              &
                    pewov,       pnsov,        pnwov,   pneov,    &
                    pustrgw,     pvstrgw,      pvdisgw,           &
                    ptte,        pvol,         pvom )

  !
  ! Description:
  !
  ! Does the gravity wave parameterisation.
  !
  ! Method:
  !
  ! This routine computes the physical tendencies of the
  ! prognostic variables u,v  and t due to vertical transports by
  ! subgridscale orographically excited gravity waves.
  !
  ! *gwdrag* is called from *physc*.
  !
  ! The scheme consists of two parts, the calculation of gravity
  ! wave stress and the stress profile in the vertical.
  ! The stress is computed using a low-level wind, static stability
  ! and an orographic variance. Four components of variance are
  ! available, the choice determined by the wind direction.
  ! A wave Richardson number is computed at each level and by
  ! requiring that its value is never less than a critical one
  ! a value of stress is determined at each model level.
  !
  ! The gravity wave stress comprises two parts depending on
  ! the critical froude number of the low level flow, and an
  ! orographic anisotropy function (a measure of two-dim) is
  ! used for the supercritical component.
  !
  ! Authors:
  !
  ! B. Ritter, ECMWF, in 1986, original source
  ! M. Miller, ECMWF, in 1986, changed
  ! M. Miller, ECMWF, in 1989, changed
  ! U. Schlese, DKRZ, March 1994, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, March 2000, use unpacked orographic variances
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_control,           ONLY: nn, nrow, twodt, lamip2
  USE mo_gaussgrid,         ONLY: budw
  USE mo_diagnostics_zonal, ONLY: dgwdisz
  USE mo_param_switches,    ONLY: lgwdrag
  USE mo_constants,         ONLY: api, cpd, g, rd
  USE mo_start_dataset,     ONLY: nstart, nstep

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: klev, klevm1, klevp1, klon, ktdia, klp2

  !  Array arguments with intent(In):
  ! Input 1D
  REAL, TARGET, INTENT (IN) :: pewov(klp2) !  E-W  orographic variance
  REAL, TARGET, INTENT (IN) :: pnsov(klp2) !  N-S  orographic variance
  REAL, TARGET, INTENT (IN) :: pnwov(klp2) ! NW-SE orographic variance
  REAL, TARGET, INTENT (IN) :: pneov(klp2) ! NE-SW orographic variance
  ! Input 2D
  REAL, INTENT (IN) :: paphm1(klp2,klevp1) ! half level pressure (t-dt)
  REAL, INTENT (IN) :: papm1(klp2,klev)    ! full level pressure (t-dt)
  REAL, INTENT (IN) :: pgeom1(klp2,klev)   ! geopotential above surface (t-dt)
  REAL, INTENT (IN) :: ptm1(klp2,klev)     ! temperature (t-dt)
  REAL, INTENT (IN) :: pum1(klp2,klev)     ! zonal wind (t-dt)
  REAL, INTENT (IN) :: pvm1(klp2,klev)     ! meridional wind (t-dt)

  !  Array arguments with intent(InOut):
  ! Input 1D
  REAL, INTENT (INOUT) :: pustrgw(klp2)    ! u-gravity wave stress (accumulated, new value)
  REAL, INTENT (INOUT) :: pvstrgw(klp2)    ! v-gravity wave stress (accumulated, new value)
  REAL, INTENT (INOUT) :: pvdisgw(klp2)    ! dissipation by gravity wave drag
  ! Input 2D
  REAL, INTENT (INOUT) :: ptte(klp2,klev)  ! tendency of temperature
  REAL, INTENT (INOUT) :: pvol(klp2,klev)  ! tendency of meridional wind
  REAL, INTENT (INOUT) :: pvom(klp2,klev)  ! tendency of zonal wind

  !  Local scalars: 
  REAL :: r1, r2, v1, v2, zalfa, zanif, zanisot, zbud, zcons1, zcons2,            &
          zcons4, zcons5, zdamp, zdedt, zdelc, zdelp, zdelv, zdiagt, zdis, zdtdt, &
          zdudt, zdvdt, zdwind, zdz2n, zfr, zfrcrit, zhgeo, zhmax, zkdrag,        &
          zkdragl, zmaxd, zr, zrahilo, zrcrit, zriw, zsig1, zsigcr, zsigt,        &
          zsmooth, zsqr, zssec, zst, ztheta, ztl, ztmst, ztsec, zu, zust, zvcrit, &
          zvelo, zvew, zvne, zvns, zvnw, zvsec, zvst, zzvar
  INTEGER :: ikbcrit, ikcrit, iktop, ilevh, ilevm2, irow, jk, jkcr, &
             jkh, jl, jlcbot, jlcrit, jltop
  LOGICAL :: lo

  !  Local arrays: 
  REAL :: zdp(klon,klev), zdz2(klon,klev), zhcrit(klon,klev),                  &
          zncrit(klon,klev), znorm(klon), zrho(klon,klevp1), zri(klon,klevp1), &
          zstab(klon,klevp1), ztau(klon,klevp1), ztfr(klon), zulow(klon),      &
          zvidis(klon), zvlow(klon), zvpf(klon,klev),                          &
          zvph(klon,klevp1)
  REAL :: zvar(klon,4)

  INTEGER :: icrit(klon), isect(klon), jkcrit(klon), jkcrith(klon), &
             jkhlim(klon), jktop(klon)
  LOGICAL :: lo1(klon,klev)

  !  External functions 
  INTEGER, EXTERNAL :: intmax, intmin

  !  Intrinsic functions 
#ifdef ECLIB
  REAL :: EXPHF, SQRTHF
!DIR$ VFUNCTION EXPHF, SQRTHF
#define EXP(x)  EXPHF(x)
#define SQRT(x) SQRTHF(x)
#else
  INTRINSIC EXP, SQRT
#endif
  INTRINSIC ABS, ATAN, INT, MAX, MERGE, MIN, MOD, SUM


  !  Executable statements 

  ! Physical constants:

  ! *zkdrag * drag constant (proportional to wavenumber)
  ! *zkdragl* drag constant for supercritical stress
  ! *zrahilo* ratio of high to low level stress
  ! *zhmax  * limiting height for supercritical stress
  ! *zsigt  * top of layer for stress computation
  ! *zsigcr * top of layer for low level drag
  ! *zfrcrit* critical froude number
  ! *zrcrit*  critical richardson number for onset of wave turbulence.
  ! *zvcrit*  lower limit for low level wind to be used for gwd.
  ! *zmaxd*   maximum velocity tendency in fraction of velocity
  !           to avoid negative dissipation

  IF (nn==21) THEN
     zkdrag = 5.E-6
  ELSE IF (nn==106) THEN
     zkdrag = 3.E-5
  ELSE IF (nn==30) THEN
     zkdrag = 1.E-5
  ELSE
     IF (lamip2) THEN
        zkdrag = 3.3E-5
     ELSE
        zkdrag = 2.5E-5
     END IF  
  END IF

  zkdragl = 4.*zkdrag
  zrahilo = 0.30
  zhmax   = 10000.
  zsigt   = 0.93
  zsigcr  = 0.80
  zrcrit  = 0.25
  zvcrit  = 0.0
  zdamp   = 4.
  zsmooth = 1./(1.+zdamp)
  zfrcrit = 2.0
  zmaxd   = 0.5

  ! Security parameters.

  ! *zvsec* to secure the projection calculation.
  ! *ztsec* to secure the stress calculation.
  ! *zssec* to secure stability.

  zvsec = 0.10
  ztsec = 1.E-07
  zssec = 1.E-12

  ! Computational constants.

  ilevm2 = klev - 2
  ilevh  = klev/2
  ztmst  = twodt
  IF (nstep == nstart) ztmst = 0.5*twodt
  zdiagt = 0.5*twodt

  zcons1 = 1./rd
  zcons2 = g**2/cpd
  zcons4 = 1./(g*ztmst)
  zcons5 = 1.5*api
  irow   = nrow(1)

  IF (lgwdrag) THEN

     !-- 1. Initialize work array with directional orographic variances

     zvar(:,3) = pewov(1:klon)
     zvar(:,1) = pnsov(1:klon)
     zvar(:,2) = pneov(1:klon)
     zvar(:,4) = pnwov(1:klon)

     !-- 2. Precompute basic state variables

     !-- 2.1 Define low level wind, project winds in plane of
     !       low level wind, determine sector in which to take
     !       the variance and set indicator for critical levels.

     jkcrit(:) = klev
     jktop(:)  = klev
     jkhlim(:) = klevm1

     ! Define top of low level flow

     DO 2101 jk = klev, ilevh, -1
        DO 2102 jl = 1, klon
           lo = (paphm1(jl,jk)/paphm1(jl,klevp1)) >= zsigcr
           jkcrit(jl) = MERGE(jk,jkcrit(jl),lo)
2102    END DO
2101 END DO

     DO 2103 jl = 1, klon
        zstab(jl,klevp1) = 0.0
        zri  (jl,klevp1) = 9999.0
        zrho (jl,klevp1) = 0.0
        zulow(jl)        = 0.0
        zvlow(jl)        = 0.0
        jkcrith(jl)      = klev
        icrit(jl)        = 1
2103 END DO

     ! Define top of subcritical lowlevel drag

     DO 2110 jk = klev, ilevh, -1
        DO 2111 jl = 1, klon
           lo = (papm1(jl,jk)/paphm1(jl,klevp1)) >= zsigt
           jktop(jl) = MERGE(jk,jktop(jl),lo)
2111    END DO
2110 END DO

     jltop = intmin(klon,jktop,1)
     iktop = jktop(jltop)

     DO 2113 jk = iktop, klev
        DO 2114 jl = 1, klon
           IF (jk >= jktop(jl)) THEN
              zulow(jl) = zulow(jl) + pum1(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
              zvlow(jl) = zvlow(jl) + pvm1(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
           END IF
2114    END DO
2113 END DO

     DO 211 jl = 1, klon
        zulow(jl) = zulow(jl)/(paphm1(jl,klevp1)-paphm1(jl,jktop(jl)))
        zvlow(jl) = zvlow(jl)/(paphm1(jl,klevp1)-paphm1(jl,jktop(jl)))
        znorm(jl) = MAX(SQRT(zulow(jl)**2+zvlow(jl)**2),zvsec)
        zvph(jl,klevp1) = znorm(jl)
211  END DO

     DO 2112 jl = 1, klon
        lo = ABS(zulow(jl)) < zvsec
        zu = MERGE(zulow(jl)+2.*zvsec,zulow(jl),lo)
        ztheta = ATAN(zvlow(jl)/zu)
        isect(jl) = MOD(INT((ztheta/api+0.625)*4.),4) + 1
2112 END DO

     DO 213 jk = ktdia, klev
        DO 212 jl = 1, klon
           zvpf(jl,jk)   = (zulow(jl)*pum1(jl,jk)+zvlow(jl)*pvm1(jl,jk))/znorm(jl)
           ztau(jl,jk)   = 0.0
           zhcrit(jl,jk) = 0.0
           zncrit(jl,jk) = 0.0
           lo1(jl,jk)    = .FALSE.
212     END DO
213  END DO

     DO 215 jk = 2, klev
        DO 214 jl = 1, klon
           zdp(jl,jk)  = papm1(jl,jk) - papm1(jl,jk-1)
           zvph(jl,jk) = ((paphm1(jl,jk)-papm1(jl,jk-1))*zvpf(jl,jk) +  &
                          (papm1(jl,jk)-paphm1(jl,jk))*zvpf(jl,jk-1))/zdp(jl,jk)
           lo = zvph(jl,jk) < zvsec
           zvph(jl,jk) = MERGE(zvsec,zvph(jl,jk),lo)
           icrit(jl)   = MERGE(jk,icrit(jl),lo)
214     END DO
215  END DO

     !-- 2.2 Brunt-Vaisala frequency and density at half levels

     DO 223 jk = klev, 2, -1
        DO 222 jl = 1, klon
           zrho(jl,jk)  = 2.*paphm1(jl,jk)*zcons1/(ptm1(jl,jk)+ptm1(jl,jk-1))
           zstab(jl,jk) = 2.*zcons2/(ptm1(jl,jk)+ptm1(jl,jk-1))*     &
                         (1.-cpd*zrho(jl,jk)*(ptm1(jl,jk)-ptm1(jl,jk-1))/zdp(jl,jk))
           zstab(jl,jk) = MAX(zstab(jl,jk),zssec)
222     END DO
223  END DO

     DO 2211 jk = iktop + 1, klevm1
        DO 221 jl = 1, klon
           IF (jk > jktop(jl)) THEN
              zst = zcons2/ptm1(jl,jk)*(1.-cpd*zrho(jl,jk)*(ptm1(jl,jk)-  &
                    ptm1(jl,jk-1))/zdp(jl,jk))
              zstab(jl,klevp1) = zstab(jl,klevp1) + zst*zdp(jl,jk)
              zstab(jl,klevp1) = MAX(zstab(jl,klevp1),zssec)
              zrho(jl,klevp1)  = zrho(jl,klevp1) + paphm1(jl,jk)*2.*zdp(jl,jk)*  &
                                 zcons1/(ptm1(jl,jk)+ptm1(jl,jk-1))
           END IF
221     END DO
2211 END DO

     DO 2212 jl = 1, klon
        zstab(jl,klevp1) = zstab(jl,klevp1)/(papm1(jl,klevm1)-papm1(jl,jktop(jl)))
        zrho(jl,klevp1) =  zrho(jl,klevp1)/(papm1(jl,klevm1)-papm1(jl,jktop(jl)))
2212 END DO

     !-- 2.3 Mean flow Richardson number and critical height for froude layer

     DO 232 jk = 2, klev
        DO 231 jl = 1, klon
           zdwind     = MAX(ABS(zvpf(jl,jk)-zvpf(jl,jk-1)),zvsec)
           zri(jl,jk) = zstab(jl,jk)*(zdp(jl,jk)/(g*zrho(jl,jk)*zdwind))**2
           zri(jl,jk) = MAX(zri(jl,jk),0.25)
231     END DO
232  END DO

     DO 234 jk = klevm1, 2, -1
        DO 233 jl = 1, klon
           zhcrit(jl,jk) = zhcrit(jl,jk+1) + zvph(jl,jk)*zdp(jl,jk)
           zncrit(jl,jk) = zncrit(jl,jk+1) + SQRT(zstab(jl,jk))*zdp(jl,jk)
233     END DO
234  END DO

     DO 236 jk = klevm1, 2, -1
        DO 235 jl = 1, klon
           zhcrit(jl,jk) = zcons5*zhcrit(jl,jk)/zncrit(jl,jk)
           zhgeo = pgeom1(jl,jk)/9.81
           lo1(jl,jk) = (zhgeo>zhcrit(jl,jk))
           jkcrith(jl) = MERGE(jk,jkcrith(jl),(lo1(jl,jk) .NEQV. lo1(jl,jk+1)))
           lo = zhgeo < zhmax
           jkhlim(jl) = MERGE(jk,jkhlim(jl),lo)
235     END DO
236  END DO

     !-- 2.4 Lowest level stress

!CDIR NODEP
     DO 241 jl = 1, klon
        jkcr  = jkcrit(jl)
        zzvar = zvar(jl,isect(jl))
        zfr   = SQRT(zstab(jl,klevp1)*zzvar)/zvph(jl,klevp1)
        zdz2n = (zfrcrit*zvph(jl,klevp1))**2/zstab(jl,klevp1)
        zzvar = MIN(zzvar,zdz2n)
        ztau(jl,klevp1) = zrho(jl,klevp1)*zkdrag*SQRT(zstab(jl,klevp1))*zzvar* &
                          zvph(jl,klevp1)
        lo = (ztau(jl,klevp1)<ztsec) .OR. (icrit(jl)>=jkcr) .OR. &
             (zvph(jl,klevp1)<zvcrit)
        icrit(jl)       = MERGE(klevp1,icrit(jl),lo)
        ztau(jl,klevp1) = MERGE(0.,ztau(jl,klevp1),lo)
        ztau(jl,jkcr)   = zrahilo*ztau(jl,klevp1)
        jkcrith(jl)     = MAX(jkcrith(jl),jkhlim(jl))
        jkcrith(jl)     = MAX(jkcrith(jl),MIN(icrit(jl),jkcr))

        ! Compute anisotropy function for froude contribution to stress

        r1 = 1.0
        r2 = 1.0
        zanisot = 1.0
        zvns = pnsov(jl) + 1.0
        zvne = pneov(jl) + 1.0
        zvew = pewov(jl) + 1.0
        zvnw = pnwov(jl) + 1.0
        v1 = MAX(zvew,zvns)
        IF (v1 > 10000.) r1 = zvew/zvns
        v2 = MAX(zvne,zvnw)
        IF (v2 > 10000.) r2 = zvne/zvnw
        IF (r1 < 1.0)    r1 = 1./r1
        IF (r2 < 1.0)    r2 = 1./r2
        zanisot  = MAX(SQRT(r1),SQRT(r2))
        zanif    = 1.0 - EXP(1.0-zanisot)
        ztl      = zkdragl*zrho(jl,klevp1)*zanif*zvph(jl,klevp1)**3* &
                   (zfr-zfrcrit)**2/SQRT(zstab(jl,klevp1))
        ztfr(jl) = MERGE(ztl,0.0,(zfr>zfrcrit))
241  END DO

     !-- 3. Compute stress profile

     jlcbot  = intmax(klon,jkcrit,1)
     ikbcrit = jkcrit(jlcbot)
     ikbcrit = MIN(ikbcrit,klevm1)

     DO 330 jk = ikbcrit - 1, 2, -1

        !-- 3.1 Wave displacement at next level

        DO 311 jl = 1, klon
           IF (jk < jkcrit(jl)) THEN
              znorm(jl)   = zrho(jl,jk)*zkdrag*SQRT(zstab(jl,jk))*zvph(jl,jk)
              zdz2(jl,jk) = ztau(jl,jk+1)/znorm(jl)*2.
           END IF
311     END DO

        !-- 3.2 Wave Richardson number, new wave displacement and stress.

        DO 321 jl = 1, klon
           IF (jk < jkcrit(jl)) THEN
              zsqr  = SQRT(zri(jl,jk))
              zalfa = SQRT(zstab(jl,jk)*zdz2(jl,jk))/zvph(jl,jk)
              zriw  = zri(jl,jk)*(1.-zalfa)/(1.+zalfa*zsqr)**2
              zr    = 2. + 1./zsqr
              zdz2n = (zvph(jl,jk)*(2.*SQRT(zr)-zr))**2/zstab(jl,jk)
              lo    = zriw < zrcrit
              ztau(jl,jk) = znorm(jl)*MERGE(zdz2n,zdz2(jl,jk),lo)*0.5
              lo = (ztau(jl,jk+1)<ztsec) .OR. (jk<=icrit(jl))
              ztau(jl,jk) = MERGE(0.,ztau(jl,jk),lo)
              ztau(jl,jk) = MIN(ztau(jl,jk),ztau(jl,jk+1))
           END IF
321     END DO
330  END DO

     !-- 3.4 Low level stress profile (subcritical component)

     jlcrit = intmin(klon,jkcrit,1)
     ikcrit = jkcrit(jlcrit)

     DO 342 jk = ikcrit + 1, klevp1
!DIR$ IVDEP
!OCL NOVREC
        DO 341 jl = 1, klon
           jkcr = jkcrit(jl)
           IF (jk > jkcr) THEN
              ztau(jl,jk) = ztau(jl,jkcr)*(paphm1(jl,klevp1)-paphm1(jl,jk)) + &
                            ztau(jl,klevp1)*(paphm1(jl,jk)-paphm1(jl,jkcr))
              ztau(jl,jk) = ztau(jl,jk)/(paphm1(jl,klevp1)-paphm1(jl,jkcr))
           END IF
341     END DO
342  END DO

     !-- 3.5 Addition of low level stress profile(super critical)

     DO 352 jk = 2, klevp1
        DO 351 jl = 1, klon
           jkh = jkcrith(jl)
           zdelp = MERGE(paphm1(jl,jk)-paphm1(jl,jkh),0.0,(jk>=jkh))
           ztau(jl,jk) = ztau(jl,jk)+ztfr(jl)*zdelp/(paphm1(jl,klevp1)-paphm1(jl,jkh))
351     END DO
352  END DO

     ! top layer stress set to zero
     ! and smooth stress in upper layers

     DO 361 jl = 1, klon
        ztau(jl,1)  = 0.0
        zstab(jl,1) = 0.0
        zri(jl,1)   = 0.0
        zvph(jl,1)  = 0.0
        zvidis(jl)  = 0.0
        ztau(jl,2)  = zsmooth*ztau(jl,3)
361  END DO

     !-- 4. Tendencies

     DO 402 jk = ktdia, klev

        DO 401 jl = 1, klon
           zdudt = -g*zulow(jl)*(ztau(jl,jk+1)-ztau(jl,jk)) /      &
                   (zvph(jl,klevp1)*(paphm1(jl,jk+1)-paphm1(jl,jk)))
           zdvdt = -g*zvlow(jl)*(ztau(jl,jk+1)-ztau(jl,jk)) /      &
                   (zvph(jl,klevp1)*(paphm1(jl,jk+1)-paphm1(jl,jk)))
           zvelo = SQRT(pum1(jl,jk)**2+pvm1(jl,jk)**2)
           zdelv = SQRT((ztmst*zdudt)**2+(ztmst*zdvdt)**2)
           zdelv = MAX(zdelv,0.001)
           zdelc = MIN(zmaxd*zvelo/zdelv,1.)
           zdudt = zdelc*zdudt
           zdvdt = zdelc*zdvdt
           pvom(jl,jk) = pvom(jl,jk) + zdudt
           pvol(jl,jk) = pvol(jl,jk) + zdvdt
           zust = pum1(jl,jk) + ztmst*zdudt
           zvst = pvm1(jl,jk) + ztmst*zdvdt
           zdis = 0.5*(pum1(jl,jk)**2+pvm1(jl,jk)**2-zust**2-zvst**2)
           zdedt = zdis/ztmst
           zvidis(jl) = zvidis(jl) + zdis*(paphm1(jl,jk+1)-paphm1(jl,jk))
           zdtdt = zdedt/cpd
           ptte(jl,jk) = ptte(jl,jk) + zdtdt
401     END DO
402  END DO

     !-- 4.1 Update stress components and dissipation

     DO 411 jl = 1, klon
        pustrgw(jl) = pustrgw(jl) + zdiagt*ztau(jl,klevp1)*zulow(jl) /  &
                      zvph(jl,klevp1)
        pvstrgw(jl) = pvstrgw(jl) + zdiagt*ztau(jl,klevp1)*zvlow(jl) /  &
                      zvph(jl,klevp1)
        pvdisgw(jl) = pvdisgw(jl) + zdiagt*zcons4*zvidis(jl)
411  END DO

     zsig1 = SUM(zvidis(1:klon))
     zbud  = budw(irow)
     dgwdisz(irow) = zdiagt*zbud*zsig1*zcons4
  ELSE

     !-- 5. Necessary computations if subroutine is by-passed

     dgwdisz(irow) = 0.
  END IF

  RETURN
END SUBROUTINE gwdrag
