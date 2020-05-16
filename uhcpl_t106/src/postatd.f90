!+ complete statistics for dynamics.
!+ $Id: postatd.f90,v 1.6 1999/07/20 14:42:34 m214089 Exp $

SUBROUTINE postatd

  ! Description:
  !
  ! Complete statistics for dynamics.
  !
  ! Method:
  !
  ! *postatd* is called from *nnsc1*.
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
  USE mo_gaussgrid
  USE mo_semi_impl
  USE mo_stat_global
  USE mo_forecast_switches
  USE mo_constants
  USE mo_start_dataset
  USE mo_stat_zonal
  USE mo_doctor

  IMPLICIT NONE

  !  Local scalars: 
  REAL :: gqp, gxp, qhd, xhd, zlat, zqdelp
  INTEGER :: jjrow, jlev

  !  Intrinsic functions 
  INTRINSIC ASIN, SQRT, SUM


  !  Executable statements 

!-- 1. Complete statistics

  ! Complete global sums

  gke = SUM(gkez(1:ngl))
  gpe = SUM(gpez(1:ngl))
  glq = SUM(glqz(1:ngl))
  gps = SUM(gpsz(1:ngl))

  DO jlev = 1, nlev

    gvo = gvo + voh(jlev)
    gd = gd + dh(jlev)
    gq = gq + qh(jlev)
    gx = gx + xh(jlev)
    gt = gt + th(jlev)

    IF (lverdia) THEN
      zqdelp = 1./delph(jlev)
      voh(jlev) = SQRT(voh(jlev)*zqdelp)
      dh(jlev) = SQRT(dh(jlev)*zqdelp)
      qh(jlev) = qh(jlev)*zqdelp
      xh(jlev) = xh(jlev)*zqdelp
      th(jlev) = th(jlev)*zqdelp

    END IF

  END DO

  gvo = SQRT(gvo/gps)
  gd = SQRT(gd/gps)
  gq = gq/gps
  gx = gx/gps
  gt = gt/gps
  gte = gke + gpe
  gtpe = gte + glq

!-- 2. Print statistics

  WRITE (nout,'(a,i6,/,2(2(a,1p,e11.4,/),a,0p,f10.2,/),5(a,1p,e11.4,/))') &
&      ' Step: ', nstep,                               & 
&      '  stdVO [1/s]                   : ', gvo,      &
&      '  stdD  [1/s]                   : ', gd,       &
&      '  T [K]                         : ', gt,       &
&      '  q                             : ', gq,       &
&      '  x                             : ', gx,       &
&      '  ps [hPa]                      : ', gps*0.01, &
&      '  kinetic energy [J/m**2]       : ', gke,      &
&      '  potential energy [J/m**2]     : ', gpe,      &
&      '  total energy [J/m**2]         : ', gte,      &
&      '  latent energy [J/m**2]        : ', glq,      &
&      '  total physics energy [J/m**2] : ', gtpe

  IF (lverdia) THEN
    gqp = gq*gps/g
    gxp = gx*gps/g
    WRITE (nout,'(8x,3(11x),a,6x,2(1p,e11.4))') 'q*ps/g, x*ps/g =', gqp, gxp
    DO jlev = 1, nlev
      qhd = qh(jlev)*delph(jlev)/g
      xhd = xh(jlev)*delph(jlev)/g
      WRITE (nout,'(4x,a,i2,a,2(1p,e11.4),0p,f10.2,1x,4(1p,e11.4))') &
&          '(',jlev,')', voh(jlev), dh(jlev), th(jlev), qh(jlev),    &
&          xh(jlev), qhd, xhd
    END DO
  END IF

  IF (lumax) THEN
    WRITE (nout,'(a,f7.2,a,i3,a,f7.1,a,f7.1)') &
&        ' lat = ',ulat, 'deg  level = ', nulev, &
&        '  max sqrt(u**2+v**2) = ',uvmax, '  ul=', ulm
  END IF
  IF (lzondia) THEN
    WRITE (nout,'(a,i7)') 'Zonal statistics, nstep= ', nstep
    DO jjrow = 1, maxrow - 1, 2
      zlat = 180.*ASIN(0.5*twomu(jjrow))/api
      WRITE (nout,'(2x,f5.0,2(2(1p,e11.4),0p,f10.2,1x))') &
&          zlat, voz(jjrow), dz(jjrow), tz(jjrow), qz(jjrow), &
&          xz(jjrow), psz(jjrow)*0.01
    END DO

    DO jjrow = maxrow, 2, -2
      zlat = 180.*ASIN(0.5*twomu(jjrow))/api
      WRITE (nout,'(2x,f5.0,2(2(1p,e11.4),0p,f10.2,1x))') &
&          zlat, voz(jjrow), dz(jjrow), tz(jjrow), qz(jjrow), &
&          xz(jjrow), psz(jjrow)*0.01
    END DO
  END IF

  RETURN
END SUBROUTINE postatd
