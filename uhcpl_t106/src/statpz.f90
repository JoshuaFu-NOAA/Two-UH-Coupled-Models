!+ complete global sums for physcis diagnostics
!+ $Id: statpz.f90,v 1.4 1998/10/28 12:33:43 m214003 Exp $

SUBROUTINE statpz

  ! Description:
  !
  ! Sums up the zonal sums collected in module *mo_diagnostics_zonal*
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, April 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_parameters
  USE mo_control
  USE mo_diagnostics
  USE mo_diagnostics_zonal

  IMPLICIT NONE

  !  Intrinsic functions 
  INTRINSIC SUM


  !  Executable statements 

!-- 1. Complete global sums for physics diagnostics

  dadcon = dadcon + SUM(dadconz(1:ngl))

  dsrad0 = dsrad0 + SUM(dsrad0z(1:ngl))
  dtrad0 = dtrad0 + SUM(dtrad0z(1:ngl))
  dsrads = dsrads + SUM(dsradsz(1:ngl))
  dtrads = dtrads + SUM(dtradsz(1:ngl))

  dvdis = dvdis + SUM(dvdisz(1:ngl),1)
  dhfs = dhfs + SUM(dhfsz(1:ngl),1)
  devap = devap + SUM(devapz(1:ngl),1)

  dgwdis = dgwdis + SUM(dgwdisz(1:ngl),1)

  dcvgr = dcvgr + SUM(dcvgrz(1:ngl))
  dcver = dcver + SUM(dcverz(1:ngl))
  dcvgs = dcvgs + SUM(dcvgsz(1:ngl))
  dcvfr = dcvfr + SUM(dcvfrz(1:ngl))
  dcvms = dcvms + SUM(dcvmsz(1:ngl))

  dlsgr = dlsgr + SUM(dlsgrz(1:ngl))
  dlsgs = dlsgs + SUM(dlsgsz(1:ngl))
  dlsms = dlsms + SUM(dlsmsz(1:ngl))
  dlser = dlser + SUM(dlserz(1:ngl))
  dlses = dlses + SUM(dlsesz(1:ngl))

  dssrad = dssrad + SUM(dssradz(1:ngl))
  dstrad = dstrad + SUM(dstradz(1:ngl))
  dshfl = dshfl + SUM(dshflz(1:ngl))
  dsdtfl = dsdtfl + SUM(dsdtflz(1:ngl))
  dslsr = dslsr + SUM(dslsrz(1:ngl))
  dslss = dslss + SUM(dslssz(1:ngl))
  dscvr = dscvr + SUM(dscvrz(1:ngl))
  dscvs = dscvs + SUM(dscvsz(1:ngl))
  dsevw = dsevw + SUM(dsevwz(1:ngl))
  dsevi = dsevi + SUM(dseviz(1:ngl))
  dssnmt = dssnmt + SUM(dssnmtz(1:ngl))
  dsros = dsros + SUM(dsrosz(1:ngl))

  RETURN
END SUBROUTINE statpz
