!+ grid point computations.

SUBROUTINE gpc

  ! Description:
  !
  ! Grid point computations.
  !
  ! Method:
  !
  ! This subroutine controls the computations in grid points, that is 
  ! the dynamical computations(*dyn*), the time filtering(*tf1* and *tf2*), 
  ! the physical computations(*phys*) and grid point contributions to the 
  ! semi implicit scheme (*si1*).
  !
  ! *gpc* is called from *nnsc1*.
  ! 
  ! Externals:
  !   *dyn*       dynamical computations.
  !   *tf12*      time filtering.
  !   *si1*       grid point contributions to the semi implicit scheme.
  !   *phys*      physical computations.
  !   *statd*     compute statistics of prognostic variables.
  !   *statp*     compute physical diagnostics.
  !   *writup*    controls write up times.
  !
  ! 1-appendix *b1:organisation of the spectral model.
  ! m.j       5/10/81.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, January 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, November 1998, modify nudging
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_start_dataset, ONLY: nstart, nstep
  USE mo_control,       ONLY: lnudge, nlat, nrow, ngl, lalai, lavgrat, lamip2&
                    ,nlon,lcouple  !fu++
  USE mo_nudging,       ONLY: NudgingSSTnew
  USE mo_sst,           ONLY: clsst, clsst2, cplsst
  USE mo_alai,          ONLY: clalai
  USE mo_avgrat,        ONLY: clavgrat
  !fu++
  USE mo_landsea,       ONLY: bzone, sstatm
  USE mo_memory_g3a,    ONLY: tsm, tsm1m
  !fu++end

  IMPLICIT NONE

  !  Local scalars: 
  INTEGER :: irow,jrow
  INTEGER :: i,j   !fu++

  !  External subroutines 
  EXTERNAL atmice, inisoil, physc, si1

  !  Intrinsic functions 
  INTRINSIC MOD


  !  Executable statements 

!-- 1. Dynamical computations

  ! Compute current "geographical" latitude

  irow = nrow(1)
  jrow = nrow(2)

  IF (MOD(irow,2) /= 0) THEN
    nlat(1) = irow/2 + 1
  ELSE
    nlat(1) = (ngl-irow/2) + 1
  END IF

!-- 2. Distribute climate values

  IF (lnudge) THEN
    CALL NudgingSSTnew
  ELSE IF (lamip2) THEN
    CALL clsst2
  !fu++ELSE IF (lcouple) THEN
  !  CALL cplsst
  ELSE
    CALL clsst
  ENDIF

  if(lcouple) CALL cplsst  !fu++
!======merge ocean model SST into AGCM fu++
!	do i=1,nlon
!	if(bzone(i,jrow).eq.1.and.lcouple) then
!        if(sstatm(i,jrow).lt.12.0) goto 377
!	tsm(i,jrow)=sstatm(i,jrow)+273.16
!	tsm1m(i,jrow)=tsm(i,jrow)
!	end if
! 377	end do
!============fu++ end=========================

  IF (lalai)   CALL clalai
  IF (lavgrat) CALL clavgrat

  ! Initialisation of soil temperatures.

  IF (nstep == nstart) CALL inisoil

  ! Compute seaice

  CALL atmice

!-- 3. Parametrisation of diabatic processes

  CALL physc

  RETURN
END SUBROUTINE gpc
