MODULE mo_diagnostics

  IMPLICIT NONE

  ! ---------------------------------------------------------------------
  !
  ! module *mo_diagnostics* - quantities related to global physical
  !                           diagnostics.
  !
  ! ---------------------------------------------------------------------

  REAL :: cdiats   !  land top layer energy/ts
  REAL :: cdiatd   !  land deep layer energy/td
  REAL :: cdiawd   !  land deep layer water content/wd
  REAL :: gpe0     !  global potential energy.    (nstep=nstart)
  REAL :: gke0     !  global kinetic energy.            "
  REAL :: gqm0     !  equivalent water content.         "
  REAL :: gts0     !  land top layer energy.            "
  REAL :: gtd0     !  land deep layer energy.           "
  REAL :: gws0     !  land top layer water content.     "
  REAL :: gwd0     !  land deep layer water content.    "
  REAL :: gsn0     !  land snow equivalent depth.       "
  REAL :: dsrad0   !  top solar radiation downward     (acc. j/m**2)
  REAL :: dtrad0   !  top thermal radiation downward         "
  REAL :: dsrads   !  surf. solar radiation downward         "
  REAL :: dtrads   !  surf. thermal radiation downward       "
  REAL :: dvdis    !  kin. energ. dissip. by vert. dif.      "
  REAL :: dhfs     !  surface sensible heat flux             "
  REAL :: devap    !  surface evaporation flux         (accumul. m )
  REAL :: dcvfr    !  kin. energ. gen. by cum. friction(acc. j/m**2)
  REAL :: dcvqac
  REAL :: dcvmoi   !  convec. environmental moistening (accumul. m )
  REAL :: dcvgr    !  convec. generation of rain             "
  REAL :: dcvgs    !  convec. generation of snow             "
  REAL :: dcvms    !  atmosph. melting of convec. snow       "
  REAL :: dcver    !  convec. evaporation of rain            "
  REAL :: dcves    !  convec. evaporation of snow            "
  REAL :: dlsgr    !  large scale generation of rain         "
  REAL :: dlsgs    !  large scale generation of snow         "
  REAL :: dlsms    !  atmosph. melting of l.s. snow          "
  REAL :: dlser    !  large scale evaporation of rain        "
  REAL :: dlses    !  large scale evaporation of snow        "
  REAL :: dssrad   !  soil downward solar radiation    (acc. j/m**2)
  REAL :: dstrad   !  soil downward thermal radiation        "
  REAL :: dshfl    !  soil total heat flux (lat.+sens.)      "
  REAL :: dsdtfl   !  upward thermal conduction in soil      "
  REAL :: dslsr    !  soil large scale rain fall       (accumul. m )
  REAL :: dslss    !  soil large scale snow fall             "
  REAL :: dscvr    !  soil convective rain fall              "
  REAL :: dscvs    !  soil convective snow fall              "
  REAL :: dsevw    ! *1st layer soil evaporation of water  "
  REAL :: dsevi    !  soil evaporation of snow/ice           "
  REAL :: dsdwfl   !  upward water conduction in soil        "
  REAL :: dssnmt   !  snow melt on land soil                 "
  REAL :: ddctfl   !  upward thermal conduction in deep soil "
  REAL :: ddcwfl   !  upward water conduction in deep soil   "
  REAL :: dsros    !  soil top layer run-off                 "
  REAL :: dsrod    !  soil deep layer run-off                "
  REAL :: dadcon   !  adiabatic conversion pot. to kin.(acc. j/m**2)
  REAL :: dgwdis   !  kin. energ dissip. by.gwd             ''
  REAL :: dsevdw   ! *2nd layer soil evaporation of water  "
  REAL :: dsevcw   !  3rd layer soil evaporation of water   "
  REAL :: dstsml   !  latent heat effect on ts
                   !  due to snow melting.
  REAL :: dstdml   !  latent heat effect on td
                   !  due to snow melting.

  LOGICAL :: ldiap !  *true for the physical budgets to be printed.

CONTAINS

!+ preset constants in mo_diagnostics.
!+ $Id: mo_diagnostics.f90,v 1.4 1999/01/25 19:30:31 m214030 Exp $

SUBROUTINE inidia

  ! Description:
  !
  ! Preset constants in mo_diagnostics.
  !
  ! Method:
  !
  ! *inidia* is called from *setphys* and from *initial*.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, February 1983, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A. Rhodin, MPI, Jan 1999, subroutine inidia -> module mo_diagnostics
  ! 
  ! for more details see file AUTHORS
  !

!  USE mo_diagnostics,                                                    &
!      ONLY: dadcon, dcver, dcves, dcvfr, dcvgr, dcvgs, dcvmoi, dcvms,    &
!            dcvqac, ddctfl, ddcwfl, devap, dgwdis, dhfs, dlser, dlses,   &
!            dlsgr, dlsgs, dlsms, dscvr, dscvs, dsdtfl, dsdwfl, dsevcw,   &
!            dsevdw, dsevi, dsevw, dshfl, dslsr, dslss, dsrad0, dsrads,   &
!            dsrod, dsros, dssnmt, dssrad, dstdml, dstrad, dstsml, dtrad0,&
!            dtrads, dvdis

  IMPLICIT NONE

  !  Executable statements 

!-- 1. Preset constants

  dsrad0 = 0.
  dtrad0 = 0.
  dsrads = 0.
  dtrads = 0.
  dvdis = 0.
  dgwdis = 0.
  dhfs = 0.
  devap = 0.
  dcvfr = 0.
  dcvmoi = 0.
  dcvgr = 0.
  dcvgs = 0.
  dcvms = 0.
  dcver = 0.
  dcves = 0.
  dcvqac = 0.
  dlsgr = 0.
  dlsgs = 0.
  dlsms = 0.
  dlser = 0.
  dlses = 0.
  dssrad = 0.
  dstrad = 0.
  dshfl = 0.
  dsdtfl = 0.
  ddctfl = 0.
  dslsr = 0.
  dslss = 0.
  dscvr = 0.
  dscvs = 0.
  dsevw = 0.
  dsevdw = 0.
  dstsml = 0.
  dstdml = 0.
  dsevcw = 0.
  dsevi = 0.
  dsdwfl = 0.
  ddcwfl = 0.
  dssnmt = 0.
  dsros = 0.
  dsrod = 0.
  dadcon = 0.

  RETURN
END SUBROUTINE inidia

END MODULE mo_diagnostics
