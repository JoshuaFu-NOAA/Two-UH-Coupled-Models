!+ drives the atmospheric model and an optional ocean model
!+ $Id: drive.f90,v 1.16 1999/09/24 08:02:54 k202057 Exp $

SUBROUTINE drive

  ! Description:
  !
  ! Drives the atmospheric model and an optional ocean model
  !
  ! Method:
  !
  ! Controls the calls to the atmopheric model and to an optional
  ! ocean model.
  ! It provides memory space for the semi Lagrangian transport
  ! scheme on automatic arrrays and prepares semi Lagrangian
  ! calculations.
  !
  ! *drive* is called from *control*
  !
  ! *stepon*  branch to the atmospheric model
  ! *grdini*  initialises grid for the semi Lagrangian scheme
  !

  USE mo_grid,          ONLY: pcnst,platd,plev,plevp1,plond,plat,&
&                             plon,plato2,i1,j1
  USE mo_start_dataset, ONLY: nstep
  USE mo_hyb,           ONLY: ceta, cetah
  USE mo_constants,     ONLY: g
  USE mo_memory_gl,     ONLY: lammp, phimp, sigmp

  IMPLICIT NONE

  !  Local parameters: 
  INTEGER, PARAMETER :: pmap = 20000

  !  Local scalars: 
  REAL    :: cwava, dlam
  INTEGER :: i, k, lat, ih

  ! Arrays depending on platd get an extra dimension for the hemisphere (since
  ! platd is the local number of latitudes on the extended grid per hemisphere)
  !  Local arrays: 
  REAL ::     coslam(plon), detai(plevp1), detam(plev), dphi(platd,2), &
       &      etaint(plevp1), etamid(plev), gauw(plat), lam(plond), &
       &      lbasdy(4,2,platd,2), lbasdz(4,2,plev), &
       &      lbasiy(4,2,platd,2), lbassd(4,2,plevp1), phi(platd,2), &
       &      sinlam(plon)

  ! Memory space for the semi Lagrangian scheme
  REAL    :: ub(plond,plev,platd,2)
  REAL    :: vb(plond,plev,platd,2)
  REAL    :: fb(plond,plev,pcnst,platd,2)
  INTEGER :: kdpmpf(pmap), kdpmph(pmap), kftype(pcnst)

  !  External subroutines 
  EXTERNAL grdini, stepon, sufix

  !  Intrinsic functions 
  INTRINSIC MIN

  !  Executable statements 

  ub(:,:,:,:)   = 0.
  vb(:,:,:,:)   = 0.
  fb(:,:,:,:,:) = 0.

  !-- 1. Prepare semi lagrangian calculations

  ! Set mass fixer parameters

  CALL sufix(kftype)

  ! Define eta coordinates: Used for calculation etadot vertical velocity
  ! for slt.

  DO k = 1, plev
     etamid(k) = ceta(k)
  END DO
  etaint(1) = 0.0001
  IF (plev>19) etaint(1) = MIN(1.E-5,etamid(1)*0.1)
  DO k = 2, plevp1
     etaint(k) = cetah(k)
  END DO

  CALL grdini(pmap,etamid,etaint,g,dlam,lam,phi,dphi,gauw,sinlam,coslam, &
       &      lbasdy,lbasdz,lbassd,lbasiy,detam,detai,kdpmpf,kdpmph,cwava)

  ! Initial guess for trajectory midpoints in spherical coords.
  ! nstep = 0:  Use arrival points as initial guess for trajectory midpoin
  ! nstep > 0:  Use calculated trajectory midpoints from previous time
  ! step as first guess.
  ! NOTE:  Reduce number of iters necessary for convergence after nstep=1

  IF (nstep==0) THEN
     DO ih = 1,2
        DO lat = 1, plato2
           DO k = 1, plev
              DO i = 1, plon
                 lammp(i,k,(ih-1)*plato2 + lat) = lam(i1 + i   - 1)
                 phimp(i,k,(ih-1)*plato2 + lat) = phi(j1 + lat - 1,ih)
                 sigmp(i,k,(ih-1)*plato2 + lat) = etamid(k)
              END DO
           END DO
        END DO
     END DO
  END IF

  !-- 2. Control a coupled run

  ! optional code

  !-- 3. Branch to the atmospheric model

  CALL stepon(pmap,kdpmpf,gauw,cwava,kdpmph,lam,phi,dphi,sinlam,coslam, &
       &      lbasdy,lbasdz,lbassd,lbasiy,detam,detai,dlam,etamid,etaint, &
       &      ub,vb,fb,kftype)

  !-- 4. Final wrapup of a coupled run

  ! optional code.

END SUBROUTINE drive
