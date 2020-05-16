!+ compute statistics on prognostics variables.
!+ $Id: statd.f90,v 1.9 1999/07/20 14:42:35 m214089 Exp $

SUBROUTINE statd

  ! Description:
  !
  ! Computes some statistics for the prognostic dynamical variables.
  !
  ! Method:
  !
  ! *statd* is called from *gpc*.
  !
  ! Results:
  ! *statd* computes and print with a chosen frequency rms of
  ! vorticity,divergence and means of temperature,surface pressure,
  ! humidity and energy(kinetic,potential and total).
  ! All these quantities are accumulated for each level
  ! in arrays located in *mo_stat*.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, September 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_control,     ONLY: nlev, nlevp1, nlon, nrow
  USE mo_gaussgrid,   ONLY: budw, rcsth
  USE mo_constants,   ONLY: als, alv, cpd, g, tmelt, vtmpc2
  USE mo_stat_global, ONLY: dh, qh, th, voh, xh
  USE mo_stat_zonal,  ONLY: gkez, glqz, gpez, gpsz, &
                            delph, dz, psz, qz, tz, voz, xz
  USE mo_scan_buffer, ONLY: d_scb, t_scb, u_scb, v_scb, vo_scb
  USE mo_memory_gl,   ONLY: q, x
  USE mo_memory_g3a,  ONLY: geospm
  USE mo_memory_g3b,  ONLY: aps

  IMPLICIT NONE

  !  Local scalars: 
  REAL :: zgpe, zgps, zke, zlq, zpe, zqnlon, zqpsz, zw, zwke, zwpe, zwsg, &
&      zzdelp
  INTEGER :: irow, jlev, jlon, krow, jrow

  !  Local arrays: 
  REAL :: zdelp(nlev), zdh(nlev), zqh(nlev), zth(nlev), zvoh(nlev), zxh(nlev)
  REAL :: aph(nlon,nlevp1)

  !  Intrinsic functions 
  INTRINSIC DOT_PRODUCT, MERGE, SQRT, SUM


  !  Executable statements 

!-- 1. Compute statistics

  irow = nrow(1)
  jrow = nrow(2)
  krow = (irow+1)/2

!-- 1.1    Set up weighting constants.

  zw = budw(irow)
  zwsg = zw/g
  zwke = zwsg*rcsth(krow)
  zwpe = zwsg*cpd
  zqnlon = 1./nlon

!-- 1.2    Accumulate horizontal statistics.

  voz(irow) = 0.
  dz(irow) = 0.
  qz(irow) = 0.
  xz(irow) = 0.
  tz(irow) = 0.

  aph(:,nlevp1) = aps(1:nlon,jrow)
  call pres (aph,nlon,aph(:,nlevp1),nlon)

  psz(irow) = SUM(aph(1:nlon,nlevp1))*zqnlon
  DO jlev = 1, nlev

    zdelp(jlev) = 0.
    zvoh(jlev) = 0.
    zdh(jlev) = 0.
    zqh(jlev) = 0.
    zxh(jlev) = 0.
    zth(jlev) = 0.

    DO jlon = 1, nlon
      zzdelp = aph(jlon,jlev+1) - aph(jlon,jlev)
      zdelp(jlev) = zdelp(jlev) + zzdelp
      zvoh (jlev) = zvoh (jlev) + zzdelp*vo_scb(jlon,jlev,jrow)**2
      zdh  (jlev) = zdh  (jlev) + zzdelp*d_scb (jlon,jlev,jrow)**2
      zqh  (jlev) = zqh  (jlev) + zzdelp*q     (jlon,jlev,jrow)
      zxh  (jlev) = zxh  (jlev) + zzdelp*x     (jlon,jlev,jrow)
      zth  (jlev) = zth  (jlev) + zzdelp*t_scb (jlon,jlev,jrow)
    END DO

    voz(irow) = voz(irow) + zvoh(jlev)*zqnlon
    dz(irow) = dz(irow) + zdh(jlev)*zqnlon
    qz(irow) = qz(irow) + zqh(jlev)*zqnlon
    xz(irow) = xz(irow) + zxh(jlev)*zqnlon

    tz(irow) = tz(irow) + zth(jlev)*zqnlon
  END DO

  DO jlev = 1, nlev
    delph(jlev) = delph(jlev) + zw*zdelp(jlev)
    voh(jlev) = voh(jlev) + zw*zvoh(jlev)
    dh(jlev) = dh(jlev) + zw*zdh(jlev)
    qh(jlev) = qh(jlev) + zw*zqh(jlev)
    xh(jlev) = xh(jlev) + zw*zxh(jlev)
    th(jlev) = th(jlev) + zw*zth(jlev)
  END DO

  zqpsz = 1./psz(irow)
  voz(irow) = SQRT(voz(irow)*zqpsz)
  dz(irow) = SQRT(dz(irow)*zqpsz)
  qz(irow) = qz(irow)*zqpsz
  xz(irow) = xz(irow)*zqpsz
  tz(irow) = tz(irow)*zqpsz

  zgps = zw*psz(irow)*nlon
  zgpe = zwsg*DOT_PRODUCT(aph(1:nlon,nlevp1),geospm(1:nlon,jrow))

!-- 1.3    Accumulate global statistics.

  zke = 0.
  zpe = 0.
  zlq = 0.

  DO jlev = 1, nlev

    DO jlon = 1, nlon
      zzdelp = aph(jlon,jlev+1) - aph(jlon,jlev)
      zke = zke + zzdelp*(u_scb(jlon,jlev,jrow)**2+v_scb(jlon,jlev,jrow)**2)
      zpe = zpe + zzdelp*(1.+vtmpc2*q(jlon,jlev,jrow))*t_scb(jlon,jlev,jrow)
      zlq = zlq + zzdelp*q(jlon,jlev,jrow)&
                  *MERGE(alv,als,t_scb(jlon,jlev,jrow)>tmelt)
    END DO
  END DO

  gkez(irow) = zwke*zke
  gpez(irow) = zwpe*zpe + zgpe
  glqz(irow) = zwsg*zlq
  gpsz(irow) = zgps

  RETURN
END SUBROUTINE statd
