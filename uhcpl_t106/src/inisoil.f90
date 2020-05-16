!+ initiate soil temperatures at five levels
!+ $Id: inisoil.f90,v 1.11 1999/07/23 08:21:20 m214030 Exp $

SUBROUTINE inisoil

  ! Description:
  !
  ! Initiate soil temperatures at five levels
  !
  ! Method:
  !
  ! Initialize soil parameterization scheme
  !
  ! *inisoil* is called from *gpc* after *clsst* has been called
  !
  ! Starting from the ts field (input file) temperatures are set
  ! in relation to depth of the soil layer and position of the initial
  ! day in the annual cycle over sea all levels equal sea surface
  ! temperature ts.
  ! ts is at 0.07 m
  ! thickness of layers 0.065, 0.254, 0.913, 2.902, 5.700 m
  !
  ! Authors:
  !
  ! L. Dumenil, MPI, June 1989, original source
  ! U. Schlese, DKRZ, January 1993, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_memory_g3a,    ONLY: slmm, td3m, td3m1m, td4m, td4m1m, td5m, td5m1m, &
                              tdclm, tdclm1m, tdm, tdm1m, tsm, tsm1m, tsnm, &
                              tsnm1m
  USE mo_sst,           ONLY: sst
  USE mo_constants,     ONLY: api, dayl, yearl
  USE mo_rad_switches,  ONLY: nmonth
  USE mo_control,       ONLY: ncbase, nrow
  USE mo_year,          ONLY: cd2dat
  USE mo_decomposition, ONLY: dc=>local_decomposition

  IMPLICIT NONE

  !  Local scalars: 
  REAL :: zd, zd3, zd4, zd5, zday, zdcl, zdmax, zkap, zmax, zmin, zpi, zsqrt
  INTEGER :: id, iday, im, jrow, iy, iyday, jk, jl, jm, jmax, jmin, &
&      jmonth, nmomid, nmonthl
  INTEGER :: nglon

  !  Local arrays: 
  REAL :: zcount(dc%nglon), zmonth(dc%nglon,12), zmth(12), &
          znmea(dc%nglon), zrange(dc%nglon)

  !  External functions 
  INTEGER, EXTERNAL :: ismax, ismin

  !  Intrinsic functions 
  INTRINSIC COS, EXP, REAL, SQRT


  !  Executable Statements 

  nglon = dc% nglon ! local number of longitudes
  jrow  = nrow(2)   ! local continuous latitude index

  ! Computational constants

  zpi = api
  zkap = 7.5E-7
  nmonthl = 30.
  nmomid = 15.
  ! Year day for which interpolation is done
  ! zdmax= day of local annual maximum
  ! layer depths
  zd3 = (-0.07) + 0.5*0.065
  zd4 = (-0.07) + 0.065 + 0.5*0.254
  zd5 = (-0.07) + 0.065 + 0.254 + 0.5*0.913
  zd = (-0.07) + 0.065 + 0.254 + 0.913 + 0.5*2.902
  zdcl = (-0.07) + 0.065 + 0.254 + 0.913 + 2.902 + 0.5*5.7

  ! Get initial year day

  ! Some corrections for 365 version needed !!!! LK

  iday = ncbase
  CALL cd2dat(iday,id,im,iy)
  iyday = (im-1)*nmonthl + id
  IF (nmonth/=0) iyday = (nmonth-1)*nmonthl + nmomid

!-- 1. Copy sst

  jm = 12
  DO jmonth = 1, jm
    DO jl = 1, nglon
      zmonth(jl,jmonth) = sst(jl,jrow,jmonth)
    END DO
  END DO

!-- 2. Calculate annual mean temperature

  DO jl = 1, nglon
    znmea(jl) = 0.
  END DO

  DO jk = 1, 12
    DO jl = 1, nglon
      IF (slmm(jl,jrow)>0.5) THEN
        znmea(jl) = znmea(jl) + zmonth(jl,jk)
      END IF
    END DO
  END DO

  DO jl = 1, nglon
    IF (slmm(jl,jrow)>0.5) THEN
      znmea(jl) = znmea(jl)/12.
    END IF
  END DO

!-- 3. Month of annual maximum/minimum

  jm = 12
  DO jl = 1, nglon
    IF (slmm(jl,jrow)>0.5) THEN
      DO jk = 1, jm
        zmth(jk) = zmonth(jl,jk)
      END DO
      jmax = ismax(jm,zmth,1)
      jmin = ismin(jm,zmth,1)
      zmax = zmth(jmax)
      zmin = zmth(jmin)
      zrange(jl) = zmax - zmin
      zcount(jl) = REAL(jmax)
    END IF
  END DO

!-- 4. Algorithm for temperatures at five levels in the soil

  zday = REAL(iyday)
  zsqrt = SQRT(zkap*yearl*dayl/zpi)

  DO jl = 1, nglon
    IF (slmm(jl,jrow)>0.5) THEN
      zdmax = zcount(jl)*nmonthl - nmomid
      td3m(jl,jrow) = znmea(jl) + 0.5*zrange(jl)*EXP(-zd3/zsqrt)*COS(2.*zpi*(zday- &
&          zdmax)/yearl-zd3/zsqrt)
      td4m(jl,jrow) = znmea(jl) + 0.5*zrange(jl)*EXP(-zd4/zsqrt)*COS(2.*zpi*(zday- &
&          zdmax)/yearl-zd4/zsqrt)
      td5m(jl,jrow) = znmea(jl) + 0.5*zrange(jl)*EXP(-zd5/zsqrt)*COS(2.*zpi*(zday- &
&          zdmax)/yearl-zd5/zsqrt)
      tdm(jl,jrow) = znmea(jl) + 0.5*zrange(jl)*EXP(-zd/zsqrt)*COS(2.*zpi*(zday- &
&          zdmax)/yearl-zd/zsqrt)
      tdclm(jl,jrow) = znmea(jl) + 0.5*zrange(jl)*EXP(-zdcl/zsqrt)*COS(2.*zpi*( &
&          zday-zdmax)/yearl-zdcl/zsqrt)
      tsm(jl,jrow) = td3m(jl,jrow)
      tsnm(jl,jrow) = td3m(jl,jrow)

      tsm1m(jl,jrow) = tsm(jl,jrow)
      tsnm1m(jl,jrow) = tsnm(jl,jrow)
      td3m1m(jl,jrow) = td3m(jl,jrow)
      td4m1m(jl,jrow) = td4m(jl,jrow)
      td5m1m(jl,jrow) = td5m(jl,jrow)
      tdm1m(jl,jrow) = tdm(jl,jrow)
      tdclm1m(jl,jrow) = tdclm(jl,jrow)

    ELSE

      tsm(jl,jrow) = tsm(jl,jrow)
      tsnm(jl,jrow) = tsm(jl,jrow)
      td3m(jl,jrow) = tsm(jl,jrow)
      td4m(jl,jrow) = tsm(jl,jrow)
      td5m(jl,jrow) = tsm(jl,jrow)
      tdm(jl,jrow) = tsm(jl,jrow)
      tdclm(jl,jrow) = tsm(jl,jrow)

      tsm1m(jl,jrow) = tsm(jl,jrow)
      tsnm1m(jl,jrow) = tsm(jl,jrow)
      td3m1m(jl,jrow) = tsm(jl,jrow)
      td4m1m(jl,jrow) = tsm(jl,jrow)
      td5m1m(jl,jrow) = tsm(jl,jrow)
      tdm1m(jl,jrow) = tsm(jl,jrow)
      tdclm1m(jl,jrow) = tsm(jl,jrow)

    END IF

  END DO

  RETURN
END SUBROUTINE inisoil
