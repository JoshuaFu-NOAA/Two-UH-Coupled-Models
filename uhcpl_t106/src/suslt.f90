!+ initializes module slt
!+ $Id: suslt.f90,v 1.11 2003/05/05 11:27:16 m214003 Exp $

SUBROUTINE suslt

  ! Description:
  !
  ! Initializes module *mo_slt* for slt-scheme
  !
  ! Method:
  !
  ! Implicit arguments:   module mo_grid
  !
  ! Authors:
  !
  ! M. Esch, MPI, March 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! T. Diehl, DKRZ, July 1999, parallel version
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_decomposition, ONLY: lc => local_decomposition
  USE mo_tracer,        ONLY: jps, ntraca
  USE mo_grid,          ONLY: plon, plev, plevp1, plat, plato2, pcnst, plond,     &
                              nxpt, platd, jintmx, plevd, pgls, jfirst, joverlap, &
                              jlast, i2pi, plono2, istart, istop, js, jn, jstart, &
                              jstop, pbpts, plevm1

  IMPLICIT NONE


  !  Executable statements 

  !-- 1. Set parameter

  plon   = lc%nglon
  plev   = lc%nlev
  plevp1 = lc%nlev + 1
  plat   = lc%nglat
  plato2 = lc%nglat/2

  pcnst  = jps + ntraca

  plond  = plon + 1 + 2*nxpt
  ! platd calculated per hemisphere
  platd  = plato2 + 2*nxpt + 2*jintmx
  plevd  = plev*(2+pcnst)

  pgls   = plon*plev
  !  jfirst = nxpt + 1
  !  jlast  = platd - nxpt - 1
  ! assume that the interpolant requires only one extra point
  jfirst = 2
  jlast  = platd - 2

  i2pi   = nxpt + plon + 1
  plono2 = plon/2
  istart = nxpt + 1
  istop  = nxpt + plon
  js     = 1 + nxpt + jintmx
  jn     = plato2 + nxpt + jintmx
  jstart = nxpt + jintmx + 1
  jstop  = jstart - 1 + plato2
  pbpts  = plond*plev*pcnst
  plevm1 = plev - 1

END SUBROUTINE suslt
